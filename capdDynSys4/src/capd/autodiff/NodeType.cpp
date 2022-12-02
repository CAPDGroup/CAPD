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

#include <stdexcept>
#include <cstdlib>

#include "capd/autodiff/NodeType.h"

// ------------------------------------------------------------------------
/// a general function for creating new binary node in DAG representing an expression
/// should not been used by a user
capd::autodiff::Node createBinaryNode(
    const capd::autodiff::Node& left,
    const capd::autodiff::Node& right,
    capd::autodiff::NodeType op
)
{
  std::vector<capd::autodiff::Node>& dag = * capd::autodiff::Node::dag;
  capd::autodiff::Node node(0,0,0,op);
  node.left = left.result;
  node.right = right.result;
  node.result = dag.size();
  node.isConst = dag[node.left].isConst and dag[node.right].isConst;
  node.isTimeDependentOnly = dag[node.left].isTimeDependentOnly and dag[node.right].isTimeDependentOnly;
  dag.push_back(node);
  return node;
}

// ------------------------------------------------------------------------
/// a general function for creating new unary node in DAG representing an expression
/// should not been used by a user
capd::autodiff::Node createUnaryNode(
    const capd::autodiff::Node& left,
    capd::autodiff::NodeType op
)
{
  std::vector<capd::autodiff::Node>& dag = * capd::autodiff::Node::dag;
  capd::autodiff::Node node(0,capd::autodiff::NODE_NULL,0,op);
  node.left = left.result;
  node.result = dag.size();
  node.isConst = dag[node.left].isConst;
  node.isTimeDependentOnly = dag[node.left].isTimeDependentOnly;
  dag.push_back(node);
  return node;
}

std::vector<capd::autodiff::Node>* capd::autodiff::Node::dag = 0;

namespace capd{
  namespace autodiff{


Node operator+(const Node& x, const Node& y){
  return createBinaryNode(x,y,NODE_ADD);
}

Node operator+(const Node& x, double y)
{
  Node c(y);
  return x+c;
}

Node operator+(double x, const Node& y){
  Node c(x);
  return c+y;
}

// ###############################################################

Node operator-(const Node& x, const Node& y){
  return createBinaryNode(x,y,NODE_SUB);
}

Node operator-(const Node& x, double y)
{
  Node c(y);
  return x-c;
}

Node operator-(double x, const Node& y){
  Node c(x);
  return c-y;
}

Node operator-(const Node& x){
  return createUnaryNode(x,NODE_UNARY_MINUS);
}

// ###############################################################

Node operator*(const Node& x, const Node& y){
  return createBinaryNode(x,y,NODE_MUL);
}

Node operator*(const Node& x, double y)
{
  Node c(y);
  return x*c;
}

Node operator*(double x, const Node& y){
  Node c(x);
  return c*y;
}

// ###############################################################

Node operator^(const Node& x, double y){
  if(y==1.)
    return x;
  if(y==0.)
    throw std::logic_error("Map constructor error: an expression of the form x^c, where c=0 is not allowed.");
  Node c(y);
  return createBinaryNode(x,c,NODE_POW);
}

// ###############################################################

Node operator/(const Node& x, const Node& y){
  return createBinaryNode(x,y,NODE_DIV);
}

Node operator/(const Node& x, double y)
{
  Node c(y);
  return x/c;
}

Node operator/(double x, const Node& y){
  Node c(x);
  return c/y;
}

// ###############################################################

Node& operator+=(Node& x, const Node& y){
  Node r=x+y;
  x = r;
  return x;
}

Node& operator+=(Node& x, double y){
  Node c(y);
  return x+=c;
}

// ###############################################################

Node& operator-=(Node& x, const Node& y){
  Node r=x-y;
  x = r;
  return x;
}

Node& operator-=(Node& x, double y){
  Node c(y);
  return x-=c;
}

// ###############################################################

Node& operator*=(Node& x, const Node& y){
  Node r=x*y;
  x = r;
  return x;
}

Node& operator*=(Node& x, double y){
  Node c(y);
  return x*=c;
}

// ###############################################################

Node& operator/=(Node& x, const Node& y){
  Node r=x/y;
  x = r;
  return x;
}

Node& operator/=(Node& x, double y){
  Node c(y);
  return x/=c;
}

}} // namespace capd::map

// ###############################################################

capd::autodiff::Node sqr(const capd::autodiff::Node& x){
  return createUnaryNode(x,capd::autodiff::NODE_SQR);
}

capd::autodiff::Node sqrt(const capd::autodiff::Node& x){
  return createUnaryNode(x,capd::autodiff::NODE_SQRT);
}

capd::autodiff::Node exp(const capd::autodiff::Node& x){
  return createUnaryNode(x,capd::autodiff::NODE_EXP);
}

capd::autodiff::Node log(const capd::autodiff::Node& x){
  return createUnaryNode(x,capd::autodiff::NODE_LOG);
}

capd::autodiff::Node sin(const capd::autodiff::Node& x){
  using namespace capd::autodiff;
  std::vector<Node>& dag = *Node::dag;
  Node c = createUnaryNode(x,NODE_COS);
  Node s = createUnaryNode(x,NODE_SIN);
  dag[c.result].right = s.result;
  dag[s.result].right = c.result;
  return s;
}

capd::autodiff::Node cos(const capd::autodiff::Node& x){
  using namespace capd::autodiff;
  std::vector<Node>& dag = *Node::dag;
  Node s = createUnaryNode(x,NODE_SIN);
  Node c = createUnaryNode(x,NODE_COS);
  dag[c.result].right = s.result;
  dag[s.result].right = c.result;
  return c;
}

capd::autodiff::Node atan(const capd::autodiff::Node& x){
  using namespace capd::autodiff;
  Node u = createUnaryNode(x,NODE_SQR);
  return createBinaryNode(x,u,NODE_ATAN);
}

capd::autodiff::Node asin(const capd::autodiff::Node& x){
  using namespace capd::autodiff;
  Node u = createUnaryNode(x,NODE_ONE_MINUS_SQR);
  Node v = createUnaryNode(u,NODE_SQRT);
  return createBinaryNode(x,v,NODE_ASIN);
}

capd::autodiff::Node acos(const capd::autodiff::Node& x){
  using namespace capd::autodiff;
  Node u = createUnaryNode(x,NODE_ONE_MINUS_SQR);
  Node v = createUnaryNode(u,NODE_SQRT);
  return createBinaryNode(x,v,NODE_ACOS);
}

/*
Node pow(const Node&, int);
Node pow(const Node&, double);
*/


/// @}
