
/////////////////////////////////////////////////////////////////////////////
/// @file EvalAdd.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_ADD_H_
#define _CAPD_AUTODIFF_EVAL_ADD_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// ---------------------------- Add  ------------------------------------
namespace Add
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] + right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
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
        *result = *left + *right;
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

// -------------------- ConstPlusVar --------------------------
namespace ConstPlusVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = right[coeffNo];
    else
      *result = *left + *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
  }

  template<class T, class R>
  inline void evalHelper(const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    right += shift;
    result += shift;
    const T* end = right + dataSize*order;
    for(;right!=end; right+=order,result+=order)
      if(getMask(result))
        *result = *right;
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

// -------------------- ConstPlusConst --------------------------
namespace ConstPlusConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = *left + *right;
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
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

// -------------------- ConstPlusTime --------------------------
namespace ConstPlusTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 1:
        result[1] = TypeTraits<T>::one();
        break;
      case 0:
        *result = (*left) + (*right);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
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

// -------------------- ConstPlusFunTime --------------------------
namespace ConstPlusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = right[coeffNo];
    else
      *result = (*left) + (*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
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

// -------------------- TimePlusVar --------------------------
namespace TimePlusVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
    case 0:
      *result = (*left) + (*right);
      break;
    case 1:
      result[1] = TypeTraits<T>::one() + right[1];
      break;
    default:
      result[coeffNo] = right[coeffNo];
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    ConstPlusVar::evalHelper(right,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    ConstPlusVar::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- TimePlusFunTime --------------------------
namespace TimePlusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    TimePlusVar::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
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

// -------------------- FunTimePlusVar --------------------------
namespace FunTimePlusVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] + right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(T* left, T* right, R result)
  {
    *result = *left + *right;
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    ConstPlusVar::evalHelper(right,result,binomial(dim+degree,degree)-1,order,coeffNo+order);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    ConstPlusVar::evalHomogenousPolynomial(degree,left,right,result,dag,coeffNo);
  }
}

// -------------------- FunTimePlusFunTime --------------------------

namespace FunTimePlusFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] + right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = *left + *right;
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

//use macro to define classe

CAPD_MAKE_DAG_NODE(Add);
CAPD_MAKE_DAG_NODE(ConstPlusVar);
CAPD_MAKE_DAG_NODE(ConstPlusConst);
CAPD_MAKE_DAG_NODE(ConstPlusTime);
CAPD_MAKE_DAG_NODE(ConstPlusFunTime);
CAPD_MAKE_DAG_NODE(TimePlusVar);
CAPD_MAKE_DAG_NODE(TimePlusFunTime);
CAPD_MAKE_DAG_NODE(FunTimePlusVar);
CAPD_MAKE_DAG_NODE(FunTimePlusFunTime);
/// @}
}} // namespace capd::autodiff

#endif
