/////////////////////////////////////////////////////////////////////////////
/// @file EvalHalfIntPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_AUTODIFF_EVAL_HALF_INT_POW_H
#define CAPD_AUTODIFF_EVAL_HALF_INT_POW_H

#include <algorithm>
#include "capd/autodiff/EvalPow.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- HalfIntPow ------------------------------------

namespace HalfIntPow
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    double c = TypeTraits<T>::_double(*right);
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += Real(c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/(Real(coeffNo) * (*left));
    }
    else
      *result = Math<T>::_pow(Math<T>::_sqrt(*left),(int)(2*c));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    typedef typename TypeTraits<T>::Real Real;
    const int c = TypeTraits<T>::_int(Real(2.)*(*right));
    *result = Math<T>::_pow(Math<T>::_sqrt(*left),c);
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
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

// -------------------- HalfIntPowFunTime ------------------------------------


namespace HalfIntPowFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    const double c = TypeTraits<T>::_double(*right);
    
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += Real(c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/(Real(coeffNo) * (*left));
    }
    else
      *result = Math<T>::_pow(Math<T>::_sqrt(*left),static_cast<int>(2*c));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    typedef typename TypeTraits<T>::Real Real;
    const int c = capd::TypeTraits<T>::_int(Real(2.)*(*right));
    *result = Math<T>::_pow(Math<T>::_sqrt(*left),c);
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

// -------------------- HalfIntPowTime ------------------------------------


namespace HalfIntPowTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    const double c = TypeTraits<T>::_double(*right);
    if(coeffNo){
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += Real(c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/(Real(coeffNo) * (*left));
    }
    else
      *result = Math<T>::_pow(Math<T>::_sqrt(*left),static_cast<int>(2*c));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    typedef typename TypeTraits<T>::Real Real;
    const int c = TypeTraits<T>::_int(Real(2.)*(*right));
    *result = Math<T>::_pow(Math<T>::_sqrt(*left),c);
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

// -------------------- HalfIntPowConst ------------------------------------


namespace HalfIntPowConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    if(coeffNo==0){
      const int c = TypeTraits<T>::_int(Real(2.)*(*right));
      *result = Math<T>::_pow(Math<T>::_sqrt(*left),c);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    typedef typename TypeTraits<T>::Real Real;
    const int c = TypeTraits<T>::_int(Real(2.)*(*right));
    *result = Math<T>::_pow(Math<T>::_sqrt(*left),c);
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

CAPD_MAKE_DAG_NODE(HalfIntPow);
CAPD_MAKE_DAG_NODE(HalfIntPowConst);
CAPD_MAKE_DAG_NODE(HalfIntPowTime);
CAPD_MAKE_DAG_NODE(HalfIntPowFunTime);
/// @}
}} // namespace capd::autodiff

#endif // CAPD_AUTODIFF_EVAL_HALF_INT_POW_H
