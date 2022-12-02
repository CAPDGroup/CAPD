/////////////////////////////////////////////////////////////////////////////
/// @file NodeType.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_NODETYPE_H_
#define _CAPD_AUTODIFF_NODETYPE_H_

#include <vector>
#include "capd/autodiff/DagIndexer.h"
#include "capd/vectalg/Multiindex.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

template<class T>
struct MaskIterator{
  MaskIterator(T* data, const bool* mask) : data(data), mask(mask){}
  void operator++() { ++data; ++mask; }
  void operator++(int) { ++data; ++mask; }
  const T& operator*() const { return *data; }
  T& operator*() { return *data; }
  const T& operator[](unsigned j) const { return data[j]; }
  T& operator[](unsigned j) { return data[j]; }
  void operator+=(unsigned j) { data+=j; mask+=j; }
  MaskIterator operator+(int j) const { return MaskIterator(data+j,mask+j); }

  T* data;
  const bool* mask;
};

template<class T>
inline bool getMask(MaskIterator<T> i){
  return *(i.mask);
}

template<class T>
inline bool getMask(MaskIterator<T> i, unsigned j){
  return i.mask[j];
}

template<class T>
inline bool getMask(T*) { return true; }

template<class T>
inline bool getMask(T*, unsigned) { return true; }

// ANY change in NodeType MUST be synchronized with an array of functions file eval.hpp
enum NodeType {
// ------------------ ADD -------------------------------
    NODE_ADD                    =0, // f(x,y,..)*g(x,y,..)
    NODE_CONST_PLUS_VAR         =1, // const + f(x,y,..)
    NODE_CONST_PLUS_CONST       =2, // const1 + const2
    NODE_CONST_PLUS_TIME        =3, // const + time
    NODE_CONST_PLUS_FUNTIME     =4, // const + f(time)
    NODE_TIME_PLUS_VAR          =5, // time + f(x,y,..)
    NODE_TIME_PLUS_FUNTIME      =6, // time + f(time)
    NODE_FUNTIME_PLUS_VAR       =7, // f(time) + g(x,y,...)
    NODE_FUNTIME_PLUS_FUNTIME   =8, // f(time) + g(time)
// ------------------ SUB -------------------------------
    NODE_SUB                    =10 ,
    NODE_CONST_MINUS_CONST      =11,
    NODE_CONST_MINUS_VAR        =12,
    NODE_CONST_MINUS_TIME       =13,
    NODE_CONST_MINUS_FUNTIME    =14,
    NODE_TIME_MINUS_CONST       =15,
    NODE_TIME_MINUS_FUNTIME     =16,
    NODE_TIME_MINUS_VAR         =17,
    NODE_FUNTIME_MINUS_CONST    =18,
    NODE_FUNTIME_MINUS_TIME     =19,
    NODE_FUNTIME_MINUS_FUNTIME  =20,
    NODE_FUNTIME_MINUS_VAR      =21,
    NODE_VAR_MINUS_CONST        =22,
    NODE_VAR_MINUS_TIME         =23,
    NODE_VAR_MINUS_FUNTIME      =24,
// -------------- UNARY MINUS ---------------------------
    NODE_UNARY_MINUS            =30,
    NODE_UNARY_MINUS_CONST      =31,
    NODE_UNARY_MINUS_TIME       =32,
    NODE_UNARY_MINUS_FUNTIME    =33,
// ------------------- MUL ------------------------------
    NODE_MUL                    =40, // f(x,y,..)*g(x,y,..)
    NODE_MUL_CONST_BY_VAR       =41, // const*f(x,y,...)
    NODE_MUL_CONST_BY_CONST     =42, // const1*const2
    NODE_MUL_CONST_BY_TIME      =43, // const*t
    NODE_MUL_CONST_BY_FUNTIME   =44, // const*f(t)
    NODE_MUL_TIME_BY_VAR        =45, // t*f(x,y,...)
    NODE_MUL_TIME_BY_FUNTIME    =46, // t*f(t)
    NODE_MUL_FUNTIME_BY_VAR     =47, // f(t)*g(x,y,.....)
    NODE_MUL_FUNTIME_BY_FUNTIME =48, // f(t)*g(t)
// ------------------- DIV ------------------------------
    NODE_DIV                    =50, // Division of type const/var cannot be significantly improved.
    NODE_DIV_VAR_BY_CONST       =51, // Formulae for var/var and const/var differ by only one subtraction.
    NODE_DIV_VAR_BY_TIME        =52, // The formula depends mainly on the result and the denominator.
    NODE_DIV_VAR_BY_FUNTIME     =53, // Thus we define special nodes for these cases when denominator is simpler.
    NODE_DIV_TIME_BY_CONST      =54, // cases TIME_BY_FUNTIME and TIME_BY_VAR are covered by FUNTIME_BY_FUNTIME and FUNTIME_BY_VAR, respectively.
    NODE_DIV_FUNTIME_BY_CONST   =55,
    NODE_DIV_FUNTIME_BY_TIME    =56,
    NODE_DIV_FUNTIME_BY_FUNTIME =57,
    NODE_DIV_CONST_BY_CONST     =58,
// --------------- SQUARE AND SQUARE ROOT--------------------
    NODE_SQR                    =60,
    NODE_SQR_CONST              =61,
    NODE_SQR_TIME               =62,
    NODE_SQR_FUNTIME            =63,
    NODE_SQRT                   =64,
    NODE_SQRT_CONST             =65,
    NODE_SQRT_TIME              =66,
    NODE_SQRT_FUNTIME           =67,
// ------------------- VARIOUS POWERS -----------------------
    NODE_POW                      =70, // general power, exponent can be negative or an interval
    NODE_POW_CONST                =71,
    NODE_POW_TIME                 =72,
    NODE_POW_FUNTIME              =73,
    NODE_NATURAL_POW              =74,  // exponent is a natural number
    NODE_NATURAL_POW_CONST        =75,  // exponent is a natural number
    NODE_NATURAL_POW_TIME         =76,  // exponent is a natural number
    NODE_NATURAL_POW_FUNTIME      =77,  // exponent is a natural number
    NODE_NEG_INT_POW              =78,  // exponent is a negative integer
    NODE_NEG_INT_POW_CONST        =79,  // exponent is a negative integer
    NODE_NEG_INT_POW_TIME         =80,  // exponent is a negative integer
    NODE_NEG_INT_POW_FUNTIME      =81,  // exponent is a negative integer
    NODE_HALF_INT_POW             =82,  // exponent is of the form N/2, where N is an integer (positive or negative)
    NODE_HALF_INT_POW_CONST       =83,  // exponent is of the form N/2
    NODE_HALF_INT_POW_TIME        =84,  // exponent is of the form N/2
    NODE_HALF_INT_POW_FUNTIME     =85,  // exponent is of the form N/2
    NODE_CUBE                     =86,
    NODE_CUBE_CONST               =87,
    NODE_CUBE_TIME                =88,
    NODE_CUBE_FUNTIME             =89,
    NODE_QUARTIC                  =90,
    NODE_QUARTIC_CONST            =91,
    NODE_QUARTIC_TIME             =92,
    NODE_QUARTIC_FUNTIME          =93,
//  -------------------- EXP LOG ----------------------------
    NODE_EXP                    =100,
    NODE_EXP_CONST              =101,
    NODE_EXP_TIME               =102,
    NODE_EXP_FUNTIME            =103,
    NODE_LOG                    =104,
    NODE_LOG_CONST              =105,
    NODE_LOG_TIME               =106,
    NODE_LOG_FUNTIME            =107,
//  ------------------- SIN -----------------------------
    NODE_SIN                    =110,
    NODE_SIN_CONST              =111,
    NODE_SIN_TIME               =112,
    NODE_SIN_FUNTIME            =113,
//  ------------------- ATAN -----------------------------
    NODE_ATAN,
    NODE_ATAN_CONST,
    NODE_ATAN_TIME,
    NODE_ATAN_FUNTIME,
//  ------------------- ASIN -----------------------------
    NODE_ONE_MINUS_SQR,
    NODE_ONE_MINUS_SQR_CONST,
    NODE_ONE_MINUS_SQR_TIME,
    NODE_ONE_MINUS_SQR_FUNTIME,
    NODE_ASIN,
    NODE_ASIN_CONST,
    NODE_ASIN_TIME,
    NODE_ASIN_FUNTIME,
    NODE_ACOS,
    NODE_ACOS_CONST,
    NODE_ACOS_TIME,
    NODE_ACOS_FUNTIME,
// ---------------------- VARS, PARAMS, CONST and COS ----------------
    NODE_NULL,
    NODE_CONST,
    NODE_TIME,
    NODE_PARAM,
    NODE_VAR,
    NODE_COS
  };

struct Int4{
  int left, right, result, op;
  Int4(int _left, int _right, int _result, int _op)
    : left(_left), right(_right), result(_result), op(_op)
  {}
};

struct Node : public Int4{
  double val;
  bool isConst;
  bool isTimeDependentOnly;
  static std::vector<Node>* dag;
  Node(int left, int right, int result, int op)
    : Int4(left,right,result,op),
      val(0.0), isConst(false), isTimeDependentOnly(false)
  {}

  explicit Node(double val)
    : Int4(NODE_NULL,NODE_NULL,dag->size(),NODE_CONST),
      val(val), isConst(true),isTimeDependentOnly(false)
  {
    dag->push_back(*this);
  }

  Node()
    : Int4(NODE_NULL,NODE_NULL,0,NODE_NULL),
      val(0.), isConst(false),isTimeDependentOnly(false)
  {}
};

Node operator+(const Node& x, const Node& y);
Node operator+(const Node& x, double y);
Node operator+(double x, const Node& y);

Node operator-(const Node& x, const Node& y);
Node operator-(const Node& x, double y);
Node operator-(double x, const Node& y);

Node operator-(const Node& x);

Node operator*(const Node& x, const Node& y);
Node operator*(const Node& x, double y);
Node operator*(double x, const Node& y);

Node operator/(const Node& x, const Node& y);
Node operator/(const Node& x, double y);
Node operator/(double x, const Node& y);

Node operator^(const Node& x, double);

Node& operator+=(Node& x, const Node& y);
Node& operator+=(Node& x, double y);

Node& operator-=(Node& x, const Node& y);
Node& operator-=(Node& x, double y);

Node& operator*=(Node& x, const Node& y);
Node& operator*=(Node& x, double y);

Node& operator/=(Node& x, const Node& y);
Node& operator/=(Node& x, double y);

template<class T>
class AbstractNode{
public:
  T* left;
  T* right;
  T* result;
  DagIndexer<T>* dag;
  void setDag(DagIndexer<T>* dag){
    this->dag = dag;
  }
  AbstractNode() : left(0), right(0), result(0){}
  virtual void evalC0(const int coeffNo) = 0;
  virtual void eval(const int degree, const int coeffNo) = 0;
  virtual void eval(const int degree, const int coeffNo, const bool* mask) = 0;
  virtual void evalC0HomogenousPolynomial() = 0;
  virtual void evalHomogenousPolynomial(const int degree, const int coeffNo) = 0;
  virtual void evalHomogenousPolynomial(const int degree, const int coeffNo, const bool* mask) = 0;
  virtual const char* name() const = 0;
  virtual ~AbstractNode() {}
};

#define CAPD_MAKE_DAG_NODE(ClassName)\
  template<class T>\
  class ClassName##Node : public AbstractNode<T>\
  { \
    public:\
    void evalC0(const int coeffNo) {\
      ClassName::evalC0(this->left,this->right,this->result,coeffNo);\
    }\
    void eval(const int degree, const int coeffNo) {\
      ClassName::eval(degree,this->left,this->right,this->result,this->dag,coeffNo);\
    }\
    void eval(const int degree, const int coeffNo, const bool* mask) {\
      ClassName::eval(degree,this->left,this->right,MaskIterator<T>(this->result,mask),this->dag,coeffNo);\
    }\
    void evalC0HomogenousPolynomial() {\
      ClassName::evalC0HomogenousPolynomial(this->left,this->right,this->result);\
    }\
    void evalHomogenousPolynomial(const int degree, const int coeffNo) {\
      ClassName::evalHomogenousPolynomial(degree,this->left,this->right,this->result,this->dag,coeffNo);\
    }\
    void evalHomogenousPolynomial(const int degree, const int coeffNo, const bool* mask) {\
      ClassName::evalHomogenousPolynomial(degree,this->left,this->right,MaskIterator<T>(this->result,mask),this->dag,coeffNo);\
    }\
    const char* name() const {\
      return #ClassName;\
    }\
  }

typedef Int4 MyNode;
/// @}
}} // namespace capd::autodiff

capd::autodiff::Node sqr(const capd::autodiff::Node&);
capd::autodiff::Node sqrt(const capd::autodiff::Node&);
capd::autodiff::Node exp(const capd::autodiff::Node&);
capd::autodiff::Node log(const capd::autodiff::Node&);
capd::autodiff::Node sin(const capd::autodiff::Node&);
capd::autodiff::Node cos(const capd::autodiff::Node&);
capd::autodiff::Node atan(const capd::autodiff::Node&);
capd::autodiff::Node asin(const capd::autodiff::Node&);
capd::autodiff::Node acos(const capd::autodiff::Node&);

#endif
