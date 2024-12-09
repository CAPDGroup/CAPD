/////////////////////////////////////////////////////////////////////////////
/// @file EvalOneMinusSqr.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_AUTODIFF_EVAL_ONE_MINUS_SQR_H
#define CAPD_AUTODIFF_EVAL_ONE_MINUS_SQR_H

#include "capd/autodiff/EvalSqr.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- OneMinusSqr ------------------------------------

namespace OneMinusSqr{

  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    if(coeffNo)
      result[coeffNo] = Real(-2.)*Sqr::sqrProduct(left,coeffNo);
    else
      *result = TypeTraits<T>::one() - Math<T>::_sqr(*left);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = TypeTraits<T>::one() - Math<T>::_sqr(*left);
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      case 1:
        evalC0(left,right,result,coeffNo);
        Sqr::evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        break;
      case 2:
        evalC0(left,right,result,coeffNo);
        Sqr::evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        Sqr::evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        break;
      case 3:
        evalC0(left,right,result,coeffNo);
        Sqr::evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        Sqr::evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        Sqr::evalC3HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        break;
      default:
        evalC0(left,right,result,coeffNo);
        Sqr::evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        Sqr::evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        Sqr::evalC3HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        for(const MultiindexData* i = dag->getIndexArray().begin(0,4);i<dag->getIndexArray().end(0,degree);++i)
          Sqr::evalMultiindex(left,result,i,coeffNo,-2.);
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        Sqr::evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        break;
      case 2:
        Sqr::evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        break;
      case 3:
        Sqr::evalC3HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,-2.);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          Sqr::evalMultiindex(left,result,i,coeffNo,-2.);
    }
  }
}

// -------------------- OneMinusSqrFunTime ------------------------------------

namespace OneMinusSqrFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, unsigned coeffNo)
  {
    OneMinusSqr::evalC0(left,right,result,coeffNo);
  }

 template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = TypeTraits<T>::one() - Math<T>::_sqr(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    OneMinusSqr::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- OneMinusSqrTime ------------------------------------

namespace OneMinusSqrTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    switch(coeffNo)
    {
      case 2: result[2] = T(Real(-1.)); break;
      case 1: result[1] = Real(-2.)*(*left); break;
      case 0: *result = TypeTraits<T>::one()-Math<T>::_sqr(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = TypeTraits<T>::one()-Math<T>::_sqr(*left);
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

// -------------------- OneMinusSqrConst ------------------------------------

namespace OneMinusSqrConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = TypeTraits<T>::one()-Math<T>::_sqr(*left);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = TypeTraits<T>::one()-Math<T>::_sqr(*left);
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

CAPD_MAKE_DAG_NODE(OneMinusSqr);
CAPD_MAKE_DAG_NODE(OneMinusSqrTime);
CAPD_MAKE_DAG_NODE(OneMinusSqrFunTime);
CAPD_MAKE_DAG_NODE(OneMinusSqrConst);
/// @}
}} // namespace capd::autodiff

#endif // CAPD_AUTODIFF_EVAL_ONE_MINUS_SQR_H
