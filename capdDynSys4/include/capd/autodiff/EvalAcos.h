/////////////////////////////////////////////////////////////////////////////
/// @file EvalAcos.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_ACOS_H_
#define _CAPD_AUTODIFF_EVAL_ACOS_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Acos ------------------------------------

namespace Acos
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T temp = TypeTraits<T>::zero();
      for(unsigned j=1;j<coeffNo;++j)
        temp -= double(j) * result[j]*right[coeffNo-j];
      result[coeffNo] = (temp/(double)coeffNo - left[coeffNo])/(*right);
    }else{
      *result = acos(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = acos(*left);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    R resultDer = result + order;
    const T* leftDer = left + order + coeffNo;

    for(unsigned derNo=0;derNo<dim;++derNo,resultDer+=order, leftDer+=order)
    {
      if(getMask(resultDer)){
        T temp = *leftDer;
        for(unsigned j=0;j<coeffNo;++j)
          temp += resultDer[j] * right[coeffNo-j];
        resultDer[coeffNo] = -temp/(*right);
      }
    }
  }

  template<class T, class R>
  void evalMultiindex(const T* left, const T* right, R result, DagIndexer<T>* dag, const MultiindexData* i, const unsigned coeffNo)
  {
    if(getMask(result,i->index))
    {
      T t = capd::TypeTraits<T>::zero();
      int h = i->getConvolutionPairsFromEpToK(coeffNo).size();

      for(int q=0;q<h-1;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairsFromEpToK(coeffNo)[q];
        t -= dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p]*result[p.first]*right[p.second];
      }
      t = t/(double)i->k[i->p];
      result[i->index+coeffNo] = (t-left[i->index+coeffNo])/(*right);
    }
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
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        for(const MultiindexData* i = dag->getIndexArray().begin(0,2);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,dag,i,coeffNo);
    }
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,dag,i,coeffNo);
    }
  }
}

namespace AcosFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Acos::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = acos(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Acos::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

namespace AcosTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    // TODO: implement optimized code
    Acos::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = acos(*left);
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

namespace AcosConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = acos(*left);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = acos(*left);
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

CAPD_MAKE_DAG_NODE(Acos);
CAPD_MAKE_DAG_NODE(AcosConst);
CAPD_MAKE_DAG_NODE(AcosTime);
CAPD_MAKE_DAG_NODE(AcosFunTime);
/// @}
}} // namespace capd::autodiff

#endif
