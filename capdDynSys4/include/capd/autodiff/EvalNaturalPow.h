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
namespace NaturalPow
{
  /**
   * Auxiliary function.
   * Computes d^i x^c where c is integer, i>0 and x can be zero
   * @param x - array of coefficients
   */

  template<class T>
  inline
  T* evalC0SingularNaturalPow(const unsigned coeffNo, const T* x, T* temp1, T* temp2, const unsigned c)
  {
    for(unsigned i=1;i<c;++i)
    {
      for(unsigned p=0;p<=coeffNo;++p)
      {
        temp2[p] = TypeTraits<T>::zero();
        for(unsigned j=0;j<=p;++j)
          temp2[p] += x[j]*temp1[p-j];
      }
      std::swap(temp1,temp2);
    }
    return temp1;
  }

  template<class T, class R>
  void evalC0SingularNaturalPow(const T* left, T* r, const int c, R result, const unsigned coeffNo){
    T* t = new T[2*(coeffNo+1)];
    std::copy(left,left+coeffNo+1,t);
    T* p = evalC0SingularNaturalPow(coeffNo,left,t,t+coeffNo+1,c-1);
    std::copy(p,p+coeffNo+1,r);
    delete[]t;

    result[coeffNo] = left[0]*r[coeffNo];
    for(unsigned i=1;i<=coeffNo;++i)
      result[coeffNo] += left[i]*r[coeffNo-i];
  }

  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    if(!(isSingular(*left)) or coeffNo==0)
      NegIntPow::evalC0IntPow(left,c,result,coeffNo);
    else {
      T* r = const_cast<T*>(right);
      evalC0SingularNaturalPow(left,r+1,c,result,coeffNo);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  inline void evalC1SingularNaturalPow(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const T* leftDer = left + order;
    R resultDer = result + order;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      if(getMask(resultDer)){
        resultDer[coeffNo] = capd::TypeTraits<T>::zero();
        for(unsigned j=0;j<=coeffNo;++j)
          resultDer[coeffNo] += right[j+1] * leftDer[coeffNo-j];
        resultDer[coeffNo] *= leftBound(*right);
      }
    }
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    if(!(isSingular(*left))){
      NegIntPow::evalC0IntPow(left,c,result,coeffNo);
      if(degree)
        Pow::evalJetWithoutC0(degree,left,*right,result,dag,coeffNo);
    } else {
        T* r = const_cast<T*>(right);
        evalC0SingularNaturalPow(left,r+1,c,result,coeffNo);

        if(degree<=1)
          evalC1SingularNaturalPow(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        else
         throw std::logic_error("Jet propagation of natural power x^c is not implemented if x and singular or c>4 and requested derivative order is biegger than >1");
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(!(isSingular(*left))){
      Pow::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
    } else {
        if(degree<=1)
          evalC1SingularNaturalPow(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        else
         throw std::logic_error("Jet propagation of natural power x^c is not implemented if x and singular or c>4 and requested derivative order is biegger than >1");
    }
  }
}

// -------------------- NaturalPowFunTime ------------------------------------

namespace NaturalPowFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
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
    NaturalPow::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
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
  void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    NaturalPow::evalC0(left,right,result,coeffNo);
  }


  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

CAPD_MAKE_DAG_NODE(NaturalPow);
CAPD_MAKE_DAG_NODE(NaturalPowConst);
CAPD_MAKE_DAG_NODE(NaturalPowTime);
CAPD_MAKE_DAG_NODE(NaturalPowFunTime);

/// @}
}} // namespace capd::autodiff

#endif
