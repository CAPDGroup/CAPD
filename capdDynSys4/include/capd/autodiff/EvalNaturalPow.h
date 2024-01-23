/////////////////////////////////////////////////////////////////////////////
/// @file EvalNaturalPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_NATURALPOW_H_
#define _CAPD_AUTODIFF_EVAL_NATURALPOW_H_

#include <algorithm>
#include "capd/autodiff/EvalPow.h"
#include "capd/autodiff/EvalMul.h"
#include "capd/autodiff/EvalSqr.h"
#include "capd/autodiff/EvalNegIntPow.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- NaturalPow ------------------------------------
/**
 * natural powers x^c for c=2,3,4 are hand optimized and implemented in
 * EvalSqr.h, EvalCubePow.h and EvalQuarticPow.h, respectively.
 * This file implements automatic differentiation for c>4, provided x\neq 0.
 * Special case of C^0-C^1 jet propagation is implemented also in the case x=0.
 */
 
namespace SingularNaturalPow{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    int c = toInt(leftBound(*right));
    unsigned int jetSize = dag->timeJetSize();
    const T* y = nullptr;
    const T* x = left;
    T* p = const_cast<T*>(right);
    
    while(c>3 or (c>1 and y!=nullptr)){
      if((c&1)==1){
        if(y==nullptr){
          y = x;
        } else {
          p += jetSize;
          Mul::evalC0(x,y,p,coeffNo);
          y = p;
        }
      }
      p += jetSize;
      Sqr::evalC0(x,static_cast<const T*>(nullptr),p,coeffNo);
      x = p;
      c >>= 1;
    }
    if(y==nullptr){
      if(c==2){
        Sqr::evalC0(x,static_cast<const T*>(nullptr),result,coeffNo);
      } else {
        p += jetSize;
        Sqr::evalC0(x,static_cast<const T*>(nullptr),p,coeffNo);
        Mul::evalC0(x,p,result,coeffNo);
      }
    } else {
      Mul::evalC0(x,y,result,coeffNo);
    }    
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result, DagIndexer<T>* dag)
  {
    int c = toInt(leftBound(*right));
    unsigned int jetSize = dag->timeJetSize();
    const T* y = nullptr;
    const T* x = left;
    T* p = const_cast<T*>(right);
    
    while(c>3 or (c>1 and y!=nullptr)){
      if((c&1)==1){
        if(y==nullptr){
          y = x;
        } else {
          p += jetSize;
          Mul::evalC0HomogenousPolynomial(x,y,p);
          y = p;
        }
      }
      p += jetSize;
      Sqr::evalC0HomogenousPolynomial(x,static_cast<const T*>(nullptr),p);
      x = p;
      c >>= 1;
    }
    if(y==nullptr){
      if(c==2){
        Sqr::evalC0HomogenousPolynomial(x,static_cast<const T*>(nullptr),result);
      } else {
        p += jetSize;
        Sqr::evalC0HomogenousPolynomial(x,static_cast<const T*>(nullptr),p);
        Mul::evalC0HomogenousPolynomial(x,p,result);
      }
    } else {
      Mul::evalC0HomogenousPolynomial(x,y,result);
    }    
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    int c = toInt(leftBound(*right));
    unsigned int jetSize = dag->timeJetSize();
    const T* y = nullptr;
    const T* x = left;
    T* p = const_cast<T*>(right);
    
    while(c>3 or (c>1 and y!=nullptr)){
      if((c&1)==1){
        if(y==nullptr){
          y = x;
        } else {
          p += jetSize;
          Mul::eval(degree,x,y,p,dag,coeffNo);
          y = p;
        }
      }
      p += jetSize;
      Sqr::eval(degree,x,static_cast<const T*>(nullptr),p,dag,coeffNo);
      x = p;
      c >>= 1;
    }
    if(y==nullptr){
      if(c==2){
        Sqr::eval(degree,x,static_cast<const T*>(nullptr),result,dag,coeffNo);
      } else {
        p += jetSize;
        Sqr::eval(degree,x,static_cast<const T*>(nullptr),p,dag,coeffNo);
        Mul::eval(degree,x,p,result,dag,coeffNo);
      }
    } else {
      Mul::eval(degree,x,y,result,dag,coeffNo);
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    int c = toInt(leftBound(*right));
    unsigned int jetSize = dag->timeJetSize();
    const T* y = nullptr;
    const T* x = left;
    T* p = const_cast<T*>(right);
    
    while(c>3 or (c>1 and y!=nullptr)){
      if((c&1)==1){
        if(y==nullptr){
          y = x;
        } else {
          p += jetSize;
          Mul::evalHomogenousPolynomial(degree,x,y,p,dag,coeffNo);
          y = p;
        }
      }
      p += jetSize;
      Sqr::evalHomogenousPolynomial(degree,x,static_cast<const T*>(nullptr),p,dag,coeffNo);
      x = p;
      c >>= 1;
    }
    if(y==nullptr){
      if(c==2){
        Sqr::evalHomogenousPolynomial(degree,x,static_cast<const T*>(nullptr),result,dag,coeffNo);
      } else {
        p += jetSize;
        Sqr::evalHomogenousPolynomial(degree,x,static_cast<const T*>(nullptr),p,dag,coeffNo);
        Mul::evalHomogenousPolynomial(degree,x,p,result,dag,coeffNo);
      }
    } else {
      Mul::evalHomogenousPolynomial(degree,x,y,result,dag,coeffNo);
    }     
  }
}
 
namespace NaturalPow
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(!(isSingular(*left)))
      NegIntPow::evalC0IntPow(left,toInt(leftBound(*right)),result,coeffNo);
    else 
      SingularNaturalPow::evalC0(left,right,result,dag,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result, DagIndexer<T>* dag)
  {
    if(!(isSingular(*left)))
      *result = power(*left, toInt(leftBound(*right)));
    else
      SingularNaturalPow::evalC0HomogenousPolynomial(left,right,result,dag);
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(!(isSingular(*left))){
      NegIntPow::evalC0IntPow(left,toInt(leftBound(*right)),result,coeffNo);
      if(degree)
        Pow::evalJetWithoutC0(degree,left,*right,result,dag,coeffNo);
    } else {
      SingularNaturalPow::eval(degree,left,right,result,dag,coeffNo);
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(!(isSingular(*left)))
      Pow::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
    else
      SingularNaturalPow::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- NaturalPowFunTime ------------------------------------

namespace NaturalPowFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,dag,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result, DagIndexer<T>* dag)
  {
    NaturalPow::evalC0HomogenousPolynomial(left,right,result,dag);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,dag,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- NaturalPowTime ------------------------------------

namespace NaturalPowTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(!isSingular(*left)){
      NegIntPowTime::evalC0(left,right,result,coeffNo);
    } else {
      int c = toInt(leftBound(*right));
      if(coeffNo<=c){
        result[coeffNo] = power(*left,c-coeffNo);
        for(int i=0;i<coeffNo;++i){
          result[coeffNo] *= typename capd::TypeTraits<T>::Real(c-i);
          result[coeffNo] /= typename capd::TypeTraits<T>::Real(i+1);
        }
      } else {
        result[coeffNo] = capd::TypeTraits<T>::zero();
      }
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    NaturalPowTime::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- NaturalPowConst ------------------------------------

namespace NaturalPowConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    NaturalPowConst::evalC0(left,right,result,coeffNo);
  }


  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

CAPD_MAKE_DAG_NODE(NaturalPowConst);
CAPD_MAKE_DAG_NODE(NaturalPowTime);
// because of hand-optimized code for x^c, 0\in x we cannot use general macro to generate classes. In method evalC0 an extra argument dag is required
//CAPD_MAKE_DAG_NODE(NaturalPow);
//CAPD_MAKE_DAG_NODE(NaturalPowFunTime);

#define CAPD_MAKE_NAT_POW_NODE(ClassName)\
  template<class T>\
  class ClassName##Node : public AbstractNode<T>\
  { \
    public:\
    void evalC0(const int coeffNo) {\
      ClassName::evalC0(this->left,this->right,this->result,this->dag,coeffNo);\
    }\
    void eval(const int degree, const int coeffNo) {\
      ClassName::eval(degree,this->left,this->right,this->result,this->dag,coeffNo);\
    }\
    void eval(const int degree, const int coeffNo, const bool* mask) {\
      ClassName::eval(degree,this->left,this->right,MaskIterator<T>(this->result,mask),this->dag,coeffNo);\
    }\
    void evalC0HomogenousPolynomial() {\
      ClassName::evalC0HomogenousPolynomial(this->left,this->right,this->result,this->dag);\
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

CAPD_MAKE_NAT_POW_NODE(NaturalPow);
CAPD_MAKE_NAT_POW_NODE(NaturalPowFunTime);

/// @}
}} // namespace capd::autodiff

#endif
