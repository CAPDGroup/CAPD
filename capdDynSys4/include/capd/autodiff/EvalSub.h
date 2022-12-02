/////////////////////////////////////////////////////////////////////////////
/// @file EvalSub.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_SUB_H_
#define _CAPD_AUTODIFF_EVAL_SUB_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Sub ------------------------------------

namespace Sub
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    left += shift;
    right += shift;
    result += shift;
    const T* end = left + dataSize*order;
    for(;left!=end; left+=order,right+=order,result+=order)
      if(getMask(result))
        *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalHelper(left,right,result,binomial(dag->domainDimension()+degree,degree),dag->getOrder()+1,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    const unsigned shift = binomial(dim+degree-1,dim)*order;
    const unsigned dataSize = binomial(dim+degree-1,degree);
    evalHelper(left,right,result,dataSize,order,coeffNo+shift);
  }
}

// -------------------- ConstMinusVar  --------------------------

namespace ConstMinusVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = -right[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void evalHelper(const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    right += shift;
    result += shift;
    const T* end = right + dataSize*order;
    for(;right!=end; right+=order,result+=order)
      if(getMask(result))
        *result = -(*right);
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    evalHelper(right,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* /*left*/, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo){
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    const unsigned shift = binomial(dim+degree-1,dim)*order;
    const unsigned dataSize = binomial(dim+degree-1,degree);
    evalHelper(right,result,dataSize,order,coeffNo+shift);
  }
}

// -------------------- ConstMinusFunTime  --------------------------

namespace ConstMinusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = -right[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
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

// -------------------- ConstMinusTime  --------------------------

namespace ConstMinusTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==1)
      result[coeffNo] = -TypeTraits<T>::one();
    else if (coeffNo==0)
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
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

// -------------------- ConstMinusConst  --------------------------

namespace ConstMinusConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
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

// -------------------- TimeMinusConst  --------------------------

namespace TimeMinusConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==1)
      result[coeffNo] = TypeTraits<T>::one();
    else if (coeffNo==0)
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo){
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- TimeMinusFunTime  --------------------------

namespace TimeMinusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo>1)
      result[coeffNo] = - right[coeffNo];
    else
      result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
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

// -------------------- TimeMinusVar  --------------------------

namespace TimeMinusVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    TimeMinusFunTime::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    ConstMinusVar::evalHelper(right,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    ConstMinusVar::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- FunTimeMinusVar  --------------------------

namespace FunTimeMinusVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T,class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    ConstMinusVar::evalHelper(right,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    ConstMinusVar::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- FunTimeMinusFunTime  --------------------------

namespace FunTimeMinusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- FunTimeMinusTime  --------------------------

namespace FunTimeMinusTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo>1)
      result[coeffNo] = left[coeffNo];
    else
      result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
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

// -------------------- FunTimeMinusConst  --------------------------

namespace FunTimeMinusConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = left[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo){
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- VarMinusConst --------------------------

namespace VarMinusConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = left[coeffNo];
    else
      *result = *left - *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void evalHelper(const T* left, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    left += shift;
    result += shift;
    const T* end = left + dataSize*order;
    for(;left!=end; left+=order,result+=order)
      if(getMask(result))
        *result = *left;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    evalHelper(left,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, R result, DagIndexer<T>* dag, const unsigned coeffNo){
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    const unsigned shift = binomial(dim+degree-1,dim)*order;
    const unsigned dataSize = binomial(dim+degree-1,degree);
    evalHelper(left,result,dataSize,order,coeffNo+shift);
  }
}

// -------------------- VarMinusFunTime  --------------------------

namespace VarMinusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] - right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    VarMinusConst::evalHelper(left,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    VarMinusConst::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- VarMinusTime  --------------------------

namespace VarMinusTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    FunTimeMinusTime::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left - *right;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    VarMinusConst::evalHelper(left,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    VarMinusConst::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// ----------------------------------------------------------------------------------

//use macro to define classes/*

CAPD_MAKE_DAG_NODE(Sub);
CAPD_MAKE_DAG_NODE(ConstMinusVar);
CAPD_MAKE_DAG_NODE(ConstMinusFunTime);
CAPD_MAKE_DAG_NODE(ConstMinusTime);
CAPD_MAKE_DAG_NODE(ConstMinusConst);
CAPD_MAKE_DAG_NODE(TimeMinusConst);
CAPD_MAKE_DAG_NODE(TimeMinusFunTime);
CAPD_MAKE_DAG_NODE(TimeMinusVar);
CAPD_MAKE_DAG_NODE(FunTimeMinusConst);
CAPD_MAKE_DAG_NODE(FunTimeMinusTime);
CAPD_MAKE_DAG_NODE(FunTimeMinusFunTime);
CAPD_MAKE_DAG_NODE(FunTimeMinusVar);
CAPD_MAKE_DAG_NODE(VarMinusConst);
CAPD_MAKE_DAG_NODE(VarMinusTime);
CAPD_MAKE_DAG_NODE(VarMinusFunTime);
/// @}
}} // namespace capd::autodiff

#endif
