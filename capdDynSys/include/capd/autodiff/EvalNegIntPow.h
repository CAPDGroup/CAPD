/////////////////////////////////////////////////////////////////////////////
/// @file EvalNegIntPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_NEG_INT_POW_H_
#define _CAPD_AUTODIFF_EVAL_NEG_INT_POW_H_

#include <algorithm>
#include "capd/autodiff/EvalPow.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff 
/// @{
// -------------------- NegIntPow ------------------------------------

namespace NegIntPow
{
  template<class T, class R>
  inline void evalC0IntPow(const T* left, const int c, R result, const unsigned coeffNo)
  {
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += (c*((int)coeffNo-j)-j) * result[j]* left[coeffNo-j];
      result[coeffNo] = temp/((double)coeffNo * (*left));
    }
    else
      *result = power(*left,c);
  }

  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    evalC0IntPow(left,toInt(leftBound(*right)),result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    Pow::evalJetWithoutC0(degree,left,*right,result,dag,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    Pow::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- NegIntPowFunTime ------------------------------------

namespace NegIntPowFunTime
{

  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    NegIntPow::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
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

// -------------------- NegIntPowTime ------------------------------------

namespace NegIntPowTime
{

  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    const int c = toInt(leftBound(*right));
    if(coeffNo)
      result[coeffNo] = (c-((int)coeffNo-1)) * result[coeffNo-1]/((double)coeffNo * (*left));
    else
      *result = power(*left,c);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
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

// -------------------- NegIntPowConst ------------------------------------

namespace NegIntPowConst
{

  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = power(*left,toInt(leftBound(*right)));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = power(*left, toInt(leftBound(*right)));
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

CAPD_MAKE_DAG_NODE(NegIntPow);
CAPD_MAKE_DAG_NODE(NegIntPowConst);
CAPD_MAKE_DAG_NODE(NegIntPowTime);
CAPD_MAKE_DAG_NODE(NegIntPowFunTime);
/// @}
}} // namespace capd::autodiff

#endif
