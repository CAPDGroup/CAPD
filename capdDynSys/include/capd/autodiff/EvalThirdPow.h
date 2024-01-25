
/////////////////////////////////////////////////////////////////////////////
/// @file EvalThirdPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_THIRD_POW_H_
#define _CAPD_AUTODIFF_EVAL_THIRD_POW_H_

#include "capd/autodiff/NodeType.h"
#include "capd/autodiff/EvalMul.h"
#include "capd/autodiff/EvalSqr.h"
#include "capd/autodiff/EvalNegIntPow.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Cube ------------------------------------

namespace Cube{

  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
     if(!(isSingular(*left)))
       NegIntPow::evalC0IntPow(left,3,result,coeffNo);
     else {
      // compute as x*x^2. Use right node to store data for x^2
      if(coeffNo!=0){
        *right = sqr(*left);
        right[coeffNo] = 2.*Sqr::sqrProduct(left,coeffNo);
        Mul::evalC0(left,right,result,coeffNo);
        *right = T(3.);
      } else
        *result = power(*left,3);
     }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = power(*left,3);
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(!(isSingular(*left))) {
      NegIntPow::evalC0IntPow(left,3,result,coeffNo);
      if(degree)
        Pow::evalJetWithoutC0(degree,left,T(3.),result,dag,coeffNo);
    } else {
      *right = sqr(*left);
      Sqr::eval(degree,left,(const T*)0,right,dag,coeffNo);
      Mul::eval(degree,left,right,result,dag,coeffNo);
      *right = T(3.);
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(!(isSingular(*left))) {
      switch(degree)
      {
        case 1:
          Pow::evalC1HomogenousPolynomial(left,T(3.),result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
          break;
        case 2:
          Pow::evalC2HomogenousPolynomial(left,T(3.),result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
          break;
        case 3:
          Pow::evalC3HomogenousPolynomial(left,T(3.),result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
          break;
        default:
          for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
            Pow::evalMultiindex(left,T(3.),result,i,dag,coeffNo);
      }
    } else {
      *right = sqr(*left);
      Sqr::evalHomogenousPolynomial(degree,left,(const T*)0,right,dag,coeffNo);
      Mul::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
      *right = T(3.);
    }
  }
}

// -------------------- CubeFunTime ------------------------------------

namespace CubeFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, unsigned coeffNo)
  {
    Cube::evalC0(left,right,result,coeffNo);
  }

 template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    Cube::evalC0HomogenousPolynomial(left,right,result);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Cube::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- CubeTime ------------------------------------

namespace CubeTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 0: *result = power(*left,3); break;
      case 1: result[1] = 3.*sqr(*left); break;
      case 2: result[2] = 3.*(*left); break;
      case 3: result[3] = 1.;
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = power(*left,3);
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

// -------------------- CubeConst ------------------------------------

namespace CubeConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = power(*left,3);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = power(*left,3);
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

CAPD_MAKE_DAG_NODE(Cube);
CAPD_MAKE_DAG_NODE(CubeTime);
CAPD_MAKE_DAG_NODE(CubeFunTime);
CAPD_MAKE_DAG_NODE(CubeConst);
/// @}
}} // namespace capd::autodiff

#endif
