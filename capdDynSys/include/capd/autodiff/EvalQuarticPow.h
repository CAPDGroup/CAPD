/////////////////////////////////////////////////////////////////////////////
/// @file EvalQuarticPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_QUARTIC_POW_H_
#define _CAPD_AUTODIFF_EVAL_QUARTIC_POW_H_

#include "capd/autodiff/NodeType.h"
#include "capd/autodiff/EvalSqr.h"
#include "capd/autodiff/EvalNegIntPow.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Quartic ------------------------------------

namespace Quartic{

  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
     typedef typename TypeTraits<T>::Real Real;
     if(!(TypeTraits<T>::isSingular(*left)))
       NegIntPow::evalC0IntPow(left,4,result,coeffNo);
     else {
      // compute as x*x^2. Use right node to store data for x^2
      if(coeffNo!=0){
        *right = capd::Math<T>::_sqr(*left);
        right[coeffNo] = Real(2.0) * Sqr::sqrProduct(left,coeffNo);
        result[coeffNo] = Real(2.0) * Sqr::sqrProduct(right,coeffNo);
        *right = T(Real(4.));
      } else
        *result = power(*left,4);
     }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = power(*left,4);
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    const static T four = T(Real(4.));
    if(!(TypeTraits<T>::isSingular(*left))) {
      NegIntPow::evalC0IntPow(left,4,result,coeffNo);
      Pow::evalJetWithoutC0(degree,left,four,result,dag,coeffNo);
    } else {
      *right = capd::Math<T>::_sqr(*left);
      Sqr::eval(degree,left,(const T*)0,right,dag,coeffNo);
      Sqr::eval(degree,right,(const T*)0,result,dag,coeffNo);
      *right = four;
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    typedef typename capd::TypeTraits<T>::Real Real;
    const static T four = T(Real(4.));
    if(!(capd::TypeTraits<T>::isSingular(*left))) {
      switch(degree)
      {
        case 1:
          Pow::evalC1HomogenousPolynomial(left,four,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
          break;
        case 2:
          Pow::evalC2HomogenousPolynomial(left,four,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
          break;
        case 3:
          Pow::evalC3HomogenousPolynomial(left,four,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
          break;
        default:
          for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
            Pow::evalMultiindex(left,four,result,i,dag,coeffNo);
      }
    } else {
      *right = capd::Math<T>::_sqr(*left);
      Sqr::evalHomogenousPolynomial(degree,left,(const T*)0,right,dag,coeffNo);
      Sqr::evalHomogenousPolynomial(degree,right,(const T*)0,result,dag,coeffNo);
      *right = four;
    }
  }
}

// -------------------- QuarticFunTime ------------------------------------

namespace QuarticFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, unsigned coeffNo)
  {
    Quartic::evalC0(left,right,result,coeffNo);
  }

 template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    Quartic::evalC0HomogenousPolynomial(left,right,result);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Quartic::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- QuarticTime ------------------------------------

namespace QuarticTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
   typedef typename TypeTraits<T>::Real Real;
    switch(coeffNo)
    {
      case 0: *result = capd::Math<T>::_pow(*left,4); break;
      case 1: result[1] = Real(4.)*capd::Math<T>::_pow(*left,3); break;
      case 2: result[2] = Real(6.)*Math<T>::_sqr(*left); break;
      case 3: result[3] = Real(4.)*(*left); break;
      case 4: result[4] = T(Real(1.));
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = power(*left,4);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- QuarticConst ------------------------------------

namespace QuarticConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = power(*left,4);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = power(*left,4);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_DAG_NODE(Quartic);
CAPD_MAKE_DAG_NODE(QuarticTime);
CAPD_MAKE_DAG_NODE(QuarticFunTime);
CAPD_MAKE_DAG_NODE(QuarticConst);
/// @}
}} // namespace capd::autodiff

#endif
