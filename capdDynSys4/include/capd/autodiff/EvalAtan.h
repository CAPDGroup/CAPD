/////////////////////////////////////////////////////////////////////////////
/// @file EvalAtan.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_ATAN_H_
#define _CAPD_AUTODIFF_EVAL_ATAN_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Atan ------------------------------------

namespace Atan
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned j=1;j<coeffNo;++j)
        temp += double(j) * result[j]*right[coeffNo-j];
      result[coeffNo] = (left[coeffNo] - temp/(double)coeffNo)/(TypeTraits<T>::one() + *right);
    }else{
      *result = atan(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = atan(*left);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo, const T r)
  {
    R resultDer = result + order;
    const T* leftDer = left + order + coeffNo;

    for(unsigned derNo=0;derNo<dim;++derNo,resultDer+=order, leftDer+=order)
    {
      if(getMask(resultDer)){
        T temp = *leftDer;
        for(unsigned j=0;j<coeffNo;++j)
          temp -= resultDer[j] * right[coeffNo-j];
        resultDer[coeffNo] = temp/r;
      }
    }
  }

  template<class T, class R>
  void evalMultiindex(const T* left, const T* right, R result, DagIndexer<T>* dag, const MultiindexData* i, const unsigned coeffNo, const T r)
  {
    if(getMask(result,i->index))
    {
      T t = capd::TypeTraits<T>::zero();
      int h = i->getConvolutionPairsFromEpToK(coeffNo).size();

      for(int q=0;q<h-1;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairsFromEpToK(coeffNo)[q];
        t += dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p]*result[p.first]*right[p.second];
      }
      t= t/(double)i->k[i->p];
      result[i->index+coeffNo] = (left[i->index+coeffNo]-t)/r;
    }
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    const T r = (TypeTraits<T>::one()+*right);
    switch(degree)
    {
      case 0:
        evalC0(left,right,result,coeffNo);
        break;
      case 1:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,r);
        break;
      default:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,r);
        for(const MultiindexData* i = dag->getIndexArray().begin(0,2);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,dag,i,coeffNo,r);
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    const T r = (TypeTraits<T>::one()+*right);
    switch(degree)
    {
      case 1:
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo,r);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,dag,i,coeffNo,r);
    }
  }
}

namespace AtanFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Atan::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = atan(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Atan::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

namespace AtanTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 0:
        *result = atan(*left); break;
      case 1:
        result[1] = left[1]/(TypeTraits<T>::one() + *right); break;
      default:
        const T tmp =
            double(coeffNo-1)*result[coeffNo-1]*right[1] +
            double(coeffNo-2)*result[coeffNo-2]*right[2];

        result[coeffNo] = (left[coeffNo] - tmp/double(coeffNo))/(TypeTraits<T>::one() + *right);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = atan(*left);
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

namespace AtanConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = atan(*left);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = atan(*left);
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

CAPD_MAKE_DAG_NODE(Atan);
CAPD_MAKE_DAG_NODE(AtanConst);
CAPD_MAKE_DAG_NODE(AtanTime);
CAPD_MAKE_DAG_NODE(AtanFunTime);
/// @}
}} // namespace capd::autodiff

#endif
