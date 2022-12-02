/// @addtogroup autodiff
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Parser.cpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <iostream>
#include <algorithm> // for remove_if
#include <cctype>    // for isspace

#include "capd/autodiff/NodeType.h"
#include "capd/map/Parser.h"

namespace capd{
  namespace map{

/* _________________________ parseBrackets ____________________________*/

/**
   return a position of '(' which is paired with the last ')'
   or 0, when there are not any parentheses
*/

size_t parseBrackets(const std::string &e, size_t position)
{
  size_t len = e.find_last_of(')', position);
  if (len==std::string::npos) return 0; //  no parentheses found

  size_t cnt=1;  // contains a difference of the numbers of ')' and '('
  do{
    do{
      if(len==0)
      {
        std::string message = "Parse error in the formula ";
        message += e;
        throw std::runtime_error(message);
      }
      len--;
    }
    while((e.at(len)!='(')&&(e.at(len)!=')'));

    if(e.at(len)==')') cnt++;
    if(e.at(len)=='(')  cnt--;

  }while(cnt);
  return len;
}

/* _________________________ searchForOperator _____________________ */

size_t searchForOperator(const std::string& e, unsigned char op, size_t position)
/*
   checks if operator <op> appears in text <e> outside any parentheses
   and returns its position
   Otherwise 0 is returned.

   the search starts at the end
*/
{
  size_t isOperator = e.find_last_of(op,position);
  if (isOperator == std::string::npos) return 0; // no operator found

  size_t lastBracket = e.find_last_of(')',position);
  if(lastBracket == std::string::npos || isOperator>lastBracket) return isOperator; // no brackets found

  size_t pB = parseBrackets(e,position);
  if(pB) return searchForOperator(e,op,pB-1);
  return 0;
}


/* ________________________ searchForFunction _______________________ */

size_t searchForFunction(const std::string &fun, const std::string &e)
/**
   returns position of arguments of function <fun> , or
   0 if <fun> does not appear or is an argument for some other function
*/
{
  size_t pos = e.find(fun);
  if(pos || (e.at(fun.size())!='(')) return 0;  // if given fun does not appear or if it is not on the 0 position
  return parseBrackets(e);
}
/**
 * Checks if \b text begins with \b prefix
 * @param prefix
 * @param text
 * @return true if \b text begins with \b prefix, false otherwise
 */
bool checkPrefix(const std::string & prefix, const std::string & text){
  std::string::const_iterator prefix_it = prefix.begin(),
                              prefix_end = prefix.end(),
                              text_it = text.begin();
  while(prefix_it != prefix_end){
    if(*prefix_it != *text_it)
      return false;
    prefix_it++; text_it++;
  }
  return true;
}
/**
 * Checks if equation is of the form fun(params)
 * @param[in]  fun       function name to be searched for
 * @param[in]  equation  where to search
 * @param[out] args    if equation has form "fun(params)" then on exit it contains "(args)" otherwise it is not changed
 * @return true if equation has correct form, false otherwise
 */
bool searchForFunction(const std::string &fun, const std::string & equation, std::string & args)
{
  if((!checkPrefix(fun, equation))||(equation.at(fun.size())!='(')) // if given fun does not appear or is not on the 0 position
    return false;
  //TODO : we do not check here that parenthesis are matched properly.
  args=equation.substr(fun.size());
  return true;
}


/* __________________________ removeBrackets ______________________ */

std::string& removeBrackets(std::string &eq)
/**
  This function removes exterior brackets
*/
{
  while ((!parseBrackets(eq))&&(eq.at(eq.size()-1)==')')){
    eq.erase(eq.size()-1);
    eq.erase(0,1);
  }
  return eq;
}

/*---------------------- splitVariables -----------------------------*/

/*
 This function searches for text 'what' in 'where', then splits the following text
 separated by: ':' ',' ';' ' ' beginning from word next to 'what'.
 The last word is indicated by ';' after that word.
 These words are written to the array result.
*/

void splitVariables(const std::string &what, const std::string &where, std::vector<std::string>& result)
{
  size_t start = where.find(what);
  if(start==std::string::npos)
  {
    std::string message = "Cannot find '";
    message += what;
    message += "' in ";
    message += where;
    throw std::runtime_error(message);
  }
  start = where.find_first_of(':',start);
  if(start==std::string::npos)
  {
    std::string message = "Cannot find delimiter ':' in ";
    message += where;
    throw std::runtime_error(message);
  }
  size_t last = where.find_first_of(';',start);
  if(last==std::string::npos)
  {
    std::string message = "Cannot find delimiter ';' in ";
    message += where;
    throw std::runtime_error(message);
  }

  while(start<last)
  {
    size_t p = where.find_first_of(", :;",start+1);
    if(p==std::string::npos) break;
    result.push_back(where.substr(start+1,p-start-1));
    start=p;
  }
}

/*---------------------- isConstant -----------------------------*/

bool isConstant(std::string& s, double& value)
{
  removeBrackets(s);
  size_t pos = s.find_first_not_of("0123456789.+-");
  if(pos!=std::string::npos)
    return false;

  pos = s.find_last_of("+-");
  if(pos!=0 && pos!=std::string::npos)
    return false;

  value = std::strtod(s.c_str(),NULL);
  return true;
}

// --- stringToDouble ----------------------------------------------
/**
 * Converts given text to double.
 *
 * @param[in] text double number in the C format
 * @param[out] result on success value of converted text, on failure result is not changed
 * @return true on success, false if given text do not contain correct number or contains some additional characters.
 */
bool stringToDouble(std::string const & text, double & result) {
    char * endptr;
    char const * strc = text.c_str();

    double d = std::strtod(strc, &endptr);
    if(endptr != strc) {
        while(*endptr && std::isspace(*endptr)) // skip trailing whitespace
           endptr++;
        if(*endptr==0){
          result = d;
          return true;
        }
    }
   return false;
}

// ------------------------------------------------------------------------
/// removes all white spaces from \b text
void removeWhiteSpaces(std::string & text){
  text.erase(std::remove_if(text.begin(), text.end(),isspace), text.end());
}

// ------------------------------------------------------------------------
/// a general function for creating new binary node in DAG representing an expression
/// should not been used by a user
static int createBinaryNode(
    const std::string& expression,
    int hp, capd::autodiff::NodeType op,
    const std::vector<std::string>& vars,
    std::vector<capd::autodiff::Node>& dag,
    std::map<std::string,int>& knownNodes
)
{
  std::string left = std::string(expression,0,hp);
  std::string right = std::string(expression,hp+1,std::string::npos);
  capd::autodiff::Node node(0,0,0,op);
  int leftIndex = parseExpression(left,vars,dag,knownNodes);
  int rightIndex = parseExpression(right,vars,dag,knownNodes);
  node.left = leftIndex;
  node.right = rightIndex;
  node.result = dag.size();
  node.isConst = dag[leftIndex].isConst and dag[rightIndex].isConst;
  node.isTimeDependentOnly = dag[leftIndex].isTimeDependentOnly and dag[rightIndex].isTimeDependentOnly;
  dag.push_back(node);
  knownNodes[expression] = node.result;
  return node.result;
}

// ------------------------------------------------------------------------
/// a general function for creating new unary node in DAG representing an expression
/// should not been used by a user
static int createUnaryNode(
    const std::string& expression,
    std::string& args, capd::autodiff::NodeType op,
    const std::vector<std::string>& vars,
    std::vector<capd::autodiff::Node>& dag,
    std::map<std::string,int>& knownNodes
)
{
  capd::autodiff::Node node(0,capd::autodiff::NODE_NULL,0,op);
  int leftIndex = parseExpression(args,vars,dag,knownNodes);
  node.left = leftIndex;
  node.result = dag.size();
  node.isConst = dag[leftIndex].isConst;
  node.isTimeDependentOnly = dag[leftIndex].isTimeDependentOnly;
  dag.push_back(node);
  knownNodes[expression] = node.result;
  return node.result;
}

// ------------------------------------------------------------------------
/// a general function for creating new unary node in DAG representing an expression
/// should not been used by a user
static int createNode(
    const std::string& expression,
    int leftIndex, capd::autodiff::NodeType op,
    std::vector<capd::autodiff::Node>& dag,
    std::map<std::string,int>& knownNodes
)
{
  capd::autodiff::Node node(leftIndex,capd::autodiff::NODE_NULL,dag.size(),op);
  node.isConst = dag[leftIndex].isConst;
  node.isTimeDependentOnly = dag[leftIndex].isTimeDependentOnly;
  dag.push_back(node);
  knownNodes[expression] = node.result;
  return node.result;
}
// ------------------------------------------------------------------------
/// a general function that parses expression from a given string.
/// It returns:
/// @param [out] vars - vector of variables and parameters
/// @param [out] dag - DAG of the expression represented as an array. Each element of the array encodes operation, and index of left and/or right subexpressions.
/// @param [out] knownNodes - this is a temporary variable for recursive calls of this function. It contains an array of already recognized subexpression. Used for optimization of DAG.
int parseExpression(
      std::string& expression,
      const std::vector<std::string>& vars,
      std::vector<capd::autodiff::Node>& dag,
      std::map<std::string,int>& knownNodes
    )
{
  int hp;
  if (expression.at(0)=='+') expression.erase(0,1); // remove '+' from at the begin of <expression>

  removeBrackets(expression);

  std::map<std::string,int>::iterator it = knownNodes.find(expression);
  if(it!=knownNodes.end())
    return it->second;

  // check for constants
  double value;
  if(stringToDouble(expression, value)){
    capd::autodiff::Node node(capd::autodiff::NODE_NULL,capd::autodiff::NODE_NULL,dag.size(),capd::autodiff::NODE_CONST);
    node.val = value;
    node.isConst = true;
    node.isTimeDependentOnly = true;
    knownNodes[expression] = node.result;
    dag.push_back(node);
    return node.result;
  }

  // if <expression> starts with '-', then we insert '0' in front to turn
  // this into an occurrence of a subtraction
  if (expression.at(0)=='-') expression.insert(0,"0");

  hp = searchForOperator(expression,'+');
  if(hp) return createBinaryNode(expression,hp,capd::autodiff::NODE_ADD,vars,dag,knownNodes);

  hp = searchForOperator(expression,'-');
  if(hp) return createBinaryNode(expression,hp,capd::autodiff::NODE_SUB,vars,dag,knownNodes);

  hp = searchForOperator(expression,'*');
  if(hp) return createBinaryNode(expression,hp,capd::autodiff::NODE_MUL,vars,dag,knownNodes);

  hp = searchForOperator(expression,'/');
  if(hp) return createBinaryNode(expression,hp,capd::autodiff::NODE_DIV,vars,dag,knownNodes);

  hp = searchForOperator(expression,'^');
  if(hp) return createBinaryNode(expression,hp,capd::autodiff::NODE_POW,vars,dag,knownNodes);

  // there was'nt any operator found, hence we look for functions
  std::string params;
  if(searchForFunction("exp",expression,params))
    return createUnaryNode(expression,params,capd::autodiff::NODE_EXP,vars,dag,knownNodes);

  if(searchForFunction("log",expression,params))
    return createUnaryNode(expression,params,capd::autodiff::NODE_LOG,vars,dag,knownNodes);

  if(searchForFunction("sqrt",expression,params))
    return createUnaryNode(expression,params,capd::autodiff::NODE_SQRT,vars,dag,knownNodes);

  if(searchForFunction("sqr",expression,params))
    return createUnaryNode(expression,params,capd::autodiff::NODE_SQR,vars,dag,knownNodes);

  if(searchForFunction("sin",expression,params))
  {
    std::string args = params;
    int c = createUnaryNode("cos"+args,args,capd::autodiff::NODE_COS,vars,dag,knownNodes);
    int s = createUnaryNode(expression,params,capd::autodiff::NODE_SIN,vars,dag,knownNodes);
    dag[c].right = s;
    dag[s].right = c;
    return s;
  }

  if(searchForFunction("cos",expression,params))
  {
    std::string args = params;
    int s = createUnaryNode("sin"+args,args,capd::autodiff::NODE_SIN,vars,dag,knownNodes);
    int c = createUnaryNode(expression,params,capd::autodiff::NODE_COS,vars,dag,knownNodes);
    dag[c].right = s;
    dag[s].right = c;
    return c;
  }

  if(searchForFunction("atan", expression, params)){
    std::string args = params;
    int v = createUnaryNode("sqr"+args,args,capd::autodiff::NODE_SQR,vars,dag,knownNodes);
    int a = createUnaryNode(expression,params,capd::autodiff::NODE_ATAN,vars,dag,knownNodes);
    dag[a].right = v;
    return a;
  }

  if(searchForFunction("asin", expression, params)){
    std::string args1 = "1-sqr"+params;
    std::string args2 = "sqrt("+args1+")";
    int u = parseExpression(params,vars,dag,knownNodes);
    int v = createNode(args1,u,capd::autodiff::NODE_ONE_MINUS_SQR,dag,knownNodes);
    v = createNode(args2,v,capd::autodiff::NODE_SQRT,dag,knownNodes);
    int a = createNode(expression,u,capd::autodiff::NODE_ASIN,dag,knownNodes);
    dag[a].right = v;
    return a;
  }

  if(searchForFunction("acos", expression, params)){
    std::string args1 = "1-sqr"+params;
    std::string args2 = "sqrt("+args1+")";
    int u = parseExpression(params,vars,dag,knownNodes);
    int v = createNode(args1,u,capd::autodiff::NODE_ONE_MINUS_SQR,dag,knownNodes);
    v = createNode(args2,v,capd::autodiff::NODE_SQRT,dag,knownNodes);
    int a = createNode(expression,u,capd::autodiff::NODE_ACOS,dag,knownNodes);
    dag[a].right = v;
    return a;
  }

  /*
  // sqrt1mx2(u) = sqrt(1-u^2)
  if(Parser::searchForFunction("sqrt1mx2", eq, params)){
    std::string auxiliaryFunction =  params+"^2";

    m_nodes[eq] = *N = addBinaryNode<capd::diffAlgebra::Sqrt1mx2Node>(params, auxiliaryFunction);
    return;
  }
  //   sqrp1(u)=u^2+1
  hp = Parser::searchForFunction("sqrp1",eq);

  if(hp)
  {
    std::string left = std::string(eq,hp,std::string::npos);
    *N = new capd::diffAlgebra::Sqrp1Node<ScalarType>(m_order,NULL);
    eqnanal(left,&((*N)->left));
    (*N)->left->m_links++;
    m_nodes[eq]=*N;
    return;
  }
  // arc sin
  if(Parser::searchForFunction("asin", eq, params)){
    std::string help = "sqrt1mx2("+params+")";
    m_nodes[eq]= *N = addBinaryNode<capd::diffAlgebra::AsinNode>(params, help);
    return;
  }
  // arc cos
  if(Parser::searchForFunction("acos", eq, params)){
    std::string help = "sqrt1mx2"+params;
    m_nodes[eq]= *N = addBinaryNode<capd::diffAlgebra::AcosNode>(params, help);
    return;
  }
*/


  std::string re = "undefined symbol '" + expression + "' in the expression!\nIf this is a numerical constant in scientific notation please replace it by a parameter.";

  throw std::runtime_error(re);
}

// ------------------------------------------------------------------------

int parseVariables(std::string expression, std::vector<std::string>& var)
{
  var.clear();
  removeWhiteSpaces(expression);
  splitVariables("var:",expression,var);
  int numberOfVariables = var.size();
  size_t timePos = expression.find("time:");
  if(timePos!=std::string::npos)
    splitVariables("time:",expression,var);
  else
   var.push_back("_CAPD_DEFAULT_TIME_IDENTIFIER_");

  size_t parPos = expression.find("par:");
  if(parPos!=std::string::npos)
    splitVariables("par:",expression,var);

  return numberOfVariables;
}

// --------------------- optimizeMul ------------------------------------

void optimizeMul(std::vector<capd::autodiff::Node>& dag, int left, int right, int i)
{
  using namespace capd::autodiff;
  if(dag[i].isConst)    // most restrictive case const*const
  {
    dag[i].op = NODE_MUL_CONST_BY_CONST;
    return;
  }

  if(dag[i].isTimeDependentOnly)
  {
    // we have few subcases: t*const, t*f(t), f(t)*g(t)
    if(dag[left].isConst)
    {
      if(dag[right].op==NODE_TIME)           // case const*time
        dag[i].op = NODE_MUL_CONST_BY_TIME;
      else                                   // case const*f(time)
        dag[i].op = NODE_MUL_CONST_BY_FUNTIME;
      return;
    }
    // the same but in different order
    if(dag[right].isConst)
    {
      if(dag[left].op==NODE_TIME)           // case time*const
        dag[i].op = NODE_MUL_CONST_BY_TIME;
      else                                   // case f(time)*const
        dag[i].op = NODE_MUL_CONST_BY_FUNTIME;
      dag[i].left = right;
      dag[i].right = left;
      return;
    }
    // none of subexpressions is constant but the result is time dependent only
    // check for time*f(time) or f(time)*time
    if(dag[left].op==NODE_TIME){
      dag[i].op = NODE_MUL_TIME_BY_FUNTIME;
      return;
    }
    if(dag[right].op==NODE_TIME){
      dag[i].op = NODE_MUL_TIME_BY_FUNTIME;
      dag[i].left = right;
      dag[i].right = left;
      return;
    }
    // the last case is f(time)*g(time)
    dag[i].op = NODE_MUL_FUNTIME_BY_FUNTIME;
    return;
  } // end of dag[i].isTimeDependentOnly

  // the expression now depends on some variables
  // maybe one of subexpressions is constant
  if(dag[left].isConst){
    dag[i].op = NODE_MUL_CONST_BY_VAR;
    return;
  }
  if(dag[right].isConst)
  {
    dag[i].op = NODE_MUL_CONST_BY_VAR;
    dag[i].left = right;
    dag[i].right = left;
    return;
  }

  // maybe one of subexpressions is just a time?
  if(dag[left].op == NODE_TIME){
    dag[i].op = NODE_MUL_TIME_BY_VAR;
    return;
  }
  if(dag[right].op == NODE_TIME){
    dag[i].left = right;
    dag[i].right = left;
    dag[i].op = NODE_MUL_TIME_BY_VAR;
    return;
  }

  // maybe one of subexpressions depends on time, only?
  // maybe one of subexpressions is just a time?
  if(dag[left].isTimeDependentOnly){
    dag[i].op = NODE_MUL_FUNTIME_BY_VAR;
    return;
  }
  if(dag[right].isTimeDependentOnly){
    dag[i].op = NODE_MUL_FUNTIME_BY_VAR;
    dag[i].left = right;
    dag[i].right = left;
    return;
  }
  // cannot reduce the expression, this is just multiplication of two full expressions
}

// --------------------- optimizeSum ------------------------------------

void optimizeSum(std::vector<capd::autodiff::Node>& dag, int left, int right, int i)
{
  using namespace capd::autodiff;
  if(dag[i].isConst)    // most restrictive case const+const
  {
    dag[i].op = NODE_CONST_PLUS_CONST;
    return;
  }

  if(dag[i].isTimeDependentOnly)
  {
    // we have few subcases: t+const, t+f(t), f(t)+g(t), f(t)+const
    if(dag[left].isConst)
    {
      if(dag[right].op==NODE_TIME)           // case const+time
        dag[i].op = NODE_CONST_PLUS_TIME;
      else                                   // case const+f(time)
        dag[i].op = NODE_CONST_PLUS_FUNTIME;
      return;
    }
    // the same but in different order
    if(dag[right].isConst)
    {
      if(dag[left].op==NODE_TIME)           // case time+const
        dag[i].op = NODE_CONST_PLUS_TIME;
      else                                   // case f(time)+const
        dag[i].op = NODE_CONST_PLUS_FUNTIME;
      dag[i].left = right;
      dag[i].right = left;
      return;
    }
    // none of subexpressions is constant but the result is time dependent only
    // check for time+f(time) or f(time)+time
    if(dag[left].op==NODE_TIME){
      dag[i].op = NODE_TIME_PLUS_FUNTIME;
      return;
    }
    if(dag[right].op==NODE_TIME){
      dag[i].op = NODE_TIME_PLUS_FUNTIME;
      dag[i].left = right;
      dag[i].right = left;
      return;
    }
    // the last case is f(time)+g(time)
    dag[i].op = NODE_FUNTIME_PLUS_FUNTIME;
    return;
  } // end of dag[i].isTimeDependentOnly

  // the expression now depends on some variables
  // maybe one of subexpressions is constant
  if(dag[left].isConst){
    dag[i].op = NODE_CONST_PLUS_VAR;
    return;
  }
  if(dag[right].isConst)
  {
    dag[i].op = NODE_CONST_PLUS_VAR;
    dag[i].left = right;
    dag[i].right = left;
    return;
  }

  // maybe one of subexpressions is just a time?
  if(dag[left].op == NODE_TIME){
    dag[i].op = NODE_TIME_PLUS_VAR;
    return;
  }
  if(dag[right].op == NODE_TIME){
    dag[i].left = right;
    dag[i].right = left;
    dag[i].op = NODE_TIME_PLUS_VAR;
    return;
  }

  // maybe one of subexpressions depends on time only?
  if(dag[left].isTimeDependentOnly){
    dag[i].op = NODE_FUNTIME_PLUS_VAR;
    return;
  }
  if(dag[right].isTimeDependentOnly){
    dag[i].op = NODE_FUNTIME_PLUS_VAR;
    dag[i].left = right;
    dag[i].right = left;
    return;
  }
  // cannot reduce the expression, this is just sum of two full expressions
}

// --------------------- optimizeSub ------------------------------------

void optimizeSub(std::vector<capd::autodiff::Node>& dag, int left, int right, int i)
{
  using namespace capd::autodiff;
  if(dag[i].isConst)    // most restrictive case const-const
  {
    dag[i].op = NODE_CONST_MINUS_CONST;
    return;
  }

  // check all cases for unary minus
  if(dag[left].op == NODE_CONST and dag[left].val==0.)
  {
    dag[i].left = right;
    dag[i].right = left;
    if(dag[right].op==NODE_TIME)
    {
      dag[i].op = NODE_UNARY_MINUS_TIME;
      return;
    }
    if(dag[right].isTimeDependentOnly)
    {
      dag[i].op = NODE_UNARY_MINUS_FUNTIME;
      return;
    }
    dag[i].op = capd::autodiff::NODE_UNARY_MINUS;
    return;
  }

  if(dag[i].isTimeDependentOnly)
  {
    // we have few subcases: t-const, t-f(t), f(t)-g(t), const-t, const-f(t), etc
    if(dag[left].isConst)
    {
      if(dag[right].op==NODE_TIME)           // case const-time
        dag[i].op = NODE_CONST_MINUS_TIME;
      else                                   // case const-f(time)
        dag[i].op = NODE_CONST_MINUS_FUNTIME;
      return;
    }
    // the same but in different order
    if(dag[right].isConst)
    {
      if(dag[left].op==NODE_TIME)           // case time-const
        dag[i].op = NODE_TIME_MINUS_CONST;
      else                                   // case f(time)-const
        dag[i].op = NODE_FUNTIME_MINUS_CONST;
      return;
    }
    // none of subexpressions is constant but the result is time dependent only
    // check for time-f(time) or f(time)-time
    if(dag[left].op==NODE_TIME){
      dag[i].op = NODE_TIME_MINUS_FUNTIME;
      return;
    }
    if(dag[right].op==NODE_TIME){
      dag[i].op = NODE_FUNTIME_MINUS_TIME;
      return;
    }
    // the last case is f(time)-g(time)
    dag[i].op = NODE_FUNTIME_MINUS_FUNTIME;
    return;
  } // end of dag[i].isTimeDependentOnly

  // the expression now depends on some variables
  // maybe one of subexpressions is constant
  if(dag[left].isConst){
    dag[i].op = NODE_CONST_MINUS_VAR;
    return;
  }

  if(dag[right].isConst)
  {
    dag[i].op = NODE_VAR_MINUS_CONST;
    return;
  }

  // maybe one of subexpressions is just a time?
  if(dag[left].op == NODE_TIME){
    dag[i].op = NODE_TIME_MINUS_VAR;
    return;
  }
  if(dag[right].op == NODE_TIME){
    dag[i].op = NODE_VAR_MINUS_TIME;
    return;
  }

  // maybe one of subexpressions depends on time, only?
  // maybe one of subexpressions is just a time?
  if(dag[left].isTimeDependentOnly){
    dag[i].op = NODE_FUNTIME_MINUS_VAR;
    return;
  }
  if(dag[right].isTimeDependentOnly){
    dag[i].op = NODE_VAR_MINUS_FUNTIME;
    return;
  }
  // cannot reduce the expression, this is just subtraction of two full expressions
}


// --------------------- optimizeDiv------------------------------------

void optimizeDiv(std::vector<capd::autodiff::Node>& dag, int left, int right, int i)
{
  using namespace capd::autodiff;
  if(dag[i].isConst)    // most restrictive case const/const
  {
    dag[i].op = NODE_DIV_CONST_BY_CONST;
    return;
  }

  // case 1. constant in the denominator
  if(dag[right].isConst)
  {
    if(dag[left].op==NODE_TIME)           // case time/const
      dag[i].op = NODE_DIV_TIME_BY_CONST;
    else if(dag[i].isTimeDependentOnly)  // case f(time)/const
      dag[i].op = NODE_DIV_FUNTIME_BY_CONST;
    else
      dag[i].op = NODE_DIV_VAR_BY_CONST;
    return;
  }

  // case 2. time in the denominator
  if(dag[right].op == NODE_TIME)
  {
    if(dag[i].isTimeDependentOnly)          // case f(time)/time
      dag[i].op = NODE_DIV_FUNTIME_BY_TIME;
    else                                    // x/time
      dag[i].op = NODE_DIV_VAR_BY_TIME;
    return;
  }

  // case 3. f(time) in the denominator
  if(dag[right].isTimeDependentOnly)
  {
    if(dag[i].isTimeDependentOnly)          // case f(time)/g(time)
      dag[i].op = NODE_DIV_FUNTIME_BY_FUNTIME;
    else                                    // x/f(time)
      dag[i].op = NODE_DIV_VAR_BY_FUNTIME;
    return;
  }

  // denominator depends on some variables. Cannot optimize division
}

void optimizeUnivariateFunction(
    std::vector<capd::autodiff::Node>& dag,
    int left, int i,
    capd::autodiff::NodeType n,
    capd::autodiff::NodeType nConst,
    capd::autodiff::NodeType nTime,
    capd::autodiff::NodeType nFunTime
  )
{
  if(dag[i].op == n)
  {
    if(dag[left].isConst)
      dag[i].op = nConst;
    else if(dag[left].isTimeDependentOnly)
    {
      if(dag[left].op == capd::autodiff::NODE_TIME)
        dag[i].op = nTime;
      else
        dag[i].op = nFunTime;
    }
  }
}

// --------------------- optimizePow------------------------------------

void optimizePow(std::vector<capd::autodiff::Node>& dag, std::vector<int>& pos, int left, int right, int i){
  using namespace capd::autodiff;
  if(!dag[right].isConst)
    throw std::logic_error("Error in the expression: exponent in x^c must be a constant expression (can depend on constants and parameters only)!");
  optimizeUnivariateFunction(dag,left,i,NODE_POW,NODE_POW_CONST,NODE_POW_TIME,NODE_POW_FUNTIME);

  if(dag[right].op == capd::autodiff::NODE_CONST)
  {
    if(dag[right].val == 0.)
      throw std::logic_error("Map constructor error: an expression of the form x^c, with c=0 is not allowed.");
    else if(dag[right].val == 1.){
      dag[i].op = NODE_NULL;
      for(unsigned j=0;j<dag.size();++j)
      {
        if(dag[j].left==(int)i) dag[j].left = left;
        if(dag[j].right==(int)i) dag[j].right = left;
      }
      for(unsigned j=0;j<pos.size();++j)
        if(pos[j]==(int)i) pos[j] = left;
    }
    else  if(dag[right].val == 2.){
      dag[i].op = NODE_SQR;
      optimizeUnivariateFunction(dag,left,i,NODE_SQR,NODE_SQR_CONST,NODE_SQR_TIME,NODE_SQR_FUNTIME);
    }
    else  if(dag[right].val == 3.){
      dag[i].op = NODE_CUBE;
      optimizeUnivariateFunction(dag,left,i,NODE_CUBE,NODE_CUBE_CONST,NODE_CUBE_TIME,NODE_CUBE_FUNTIME);
    }
    else  if(dag[right].val == 4.){
      dag[i].op = NODE_QUARTIC;
      optimizeUnivariateFunction(dag,left,i,NODE_QUARTIC,NODE_QUARTIC_CONST,NODE_QUARTIC_TIME,NODE_QUARTIC_FUNTIME);
    }
    else if((unsigned)dag[right].val == dag[right].val){
      dag[i].op = NODE_NATURAL_POW;
      optimizeUnivariateFunction(dag,left,i,NODE_NATURAL_POW,NODE_NATURAL_POW_CONST,NODE_NATURAL_POW_TIME,NODE_NATURAL_POW_FUNTIME);
    }
    else if((int)dag[right].val == dag[right].val){
      dag[i].op = NODE_NEG_INT_POW;
      optimizeUnivariateFunction(dag,left,i,NODE_NEG_INT_POW,NODE_NEG_INT_POW_CONST,NODE_NEG_INT_POW_TIME,NODE_NEG_INT_POW_FUNTIME);
    }
    else if((int)(2*dag[right].val) == 2*dag[right].val){
      dag[i].op = NODE_HALF_INT_POW;
      optimizeUnivariateFunction(dag,left,i,NODE_HALF_INT_POW,NODE_HALF_INT_POW_CONST,NODE_HALF_INT_POW_TIME,NODE_HALF_INT_POW_FUNTIME);
    }
  }
}

void optimizeDAG(std::vector<capd::autodiff::Node>& dag, std::vector<int>& pos)
{
  using namespace capd::autodiff;
  // optimize expression
  // here we try to simplify computational DAG,
  // recognize nodes that depend on time only
  // recognize constant expressions
  for(unsigned i=0;i<dag.size();++i)
  {
    if( dag[i].op==NODE_NULL or
        dag[i].op==NODE_CONST or
        dag[i].op==NODE_TIME or
        dag[i].op==NODE_PARAM or
        dag[i].op==NODE_VAR or
        dag[i].op==NODE_COS
    ) continue;

    int left = dag[i].left;
    int right = dag[i].right;

    if(dag[i].op == NODE_MUL) optimizeMul(dag,left,right,i);
    if(dag[i].op == NODE_ADD) optimizeSum(dag,left,right,i);
    if(dag[i].op == NODE_SUB) optimizeSub(dag,left,right,i);
    if(dag[i].op == NODE_DIV) optimizeDiv(dag,left,right,i);
    if(dag[i].op == NODE_POW) optimizePow(dag,pos,left,right,i);

    optimizeUnivariateFunction(dag,left,i,NODE_SIN,NODE_SIN_CONST,NODE_SIN_TIME,NODE_SIN_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_EXP,NODE_EXP_CONST,NODE_EXP_TIME,NODE_EXP_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_LOG,NODE_LOG_CONST,NODE_LOG_TIME,NODE_LOG_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_SQR,NODE_SQR_CONST,NODE_SQR_TIME,NODE_SQR_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_ONE_MINUS_SQR,NODE_ONE_MINUS_SQR_CONST,NODE_ONE_MINUS_SQR_TIME,NODE_ONE_MINUS_SQR_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_SQRT,NODE_SQRT_CONST,NODE_SQRT_TIME,NODE_SQRT_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_ATAN,NODE_ATAN_CONST,NODE_ATAN_TIME,NODE_ATAN_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_ASIN,NODE_ASIN_CONST,NODE_ASIN_TIME,NODE_ASIN_FUNTIME);
    optimizeUnivariateFunction(dag,left,i,NODE_ACOS,NODE_ACOS_CONST,NODE_ACOS_TIME,NODE_ACOS_FUNTIME);
  }
}

int parseMap(unsigned numberOfVariables, std::string expression, const std::vector<std::string>& var, std::vector<capd::autodiff::Node>& dag, std::vector<int>& pos)
{

  removeWhiteSpaces(expression);

  unsigned i;
  std::vector<std::string> fun;
  std::map<std::string,int> knownNodes;
  using namespace capd::autodiff;

  dag.clear();
  pos.clear();

  splitVariables("fun:",expression,fun);
  size_t numberOfFunctions = fun.size();
  i=0;
  for(;i< (size_t)numberOfVariables;++i)
  {
    dag.push_back(Node(NODE_NULL,NODE_NULL,i,NODE_VAR));
    knownNodes[var[i]]=i;
  }

  Node timeVar(NODE_NULL,NODE_NULL,i,NODE_TIME);
  timeVar.isTimeDependentOnly = true;
  dag.push_back(timeVar);
  knownNodes[var[i]]=i;
  ++i;

  for(;i<var.size();++i)
  {
    Node param(NODE_NULL,NODE_NULL,i,NODE_PARAM);
    param.isConst = true;
    param.isTimeDependentOnly = true;
    dag.push_back(param);
    knownNodes[var[i]]=i;
  }
  for(i=0;i<(size_t)numberOfFunctions;++i)
  {
    int current = parseExpression(fun[i],var,dag,knownNodes);
    pos.push_back(current);
  }

  optimizeDAG(dag,pos);

  return numberOfFunctions;
}

/*
// --------------------- isVariable ------------------------------------

bool Parser::isVariable(const std::vector<std::string>& var, std::string &f, int i)
{
  removeBrackets(f);
  for(size_t i=0;i<var.size();++i)
    if(f==var[i])
      return true;
  return false;
}

template<typename Scalar>
int BasicFunction<Scalar>::isParam(std::string &f) const
{
  Parser::removeBrackets(f); // removes external brackets

  int hp=m_indexOfFirstParam;
  while((hp<m_dim)&&(f!=m_var[hp])) // looking for a variable
    hp++;

  return hp;
}
*/
}} // namespace capd::map

/// @}
