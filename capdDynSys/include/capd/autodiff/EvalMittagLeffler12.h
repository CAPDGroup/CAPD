/////////////////////////////////////////////////////////////////////////////
/// @file EvalMittagLeffler12.h
///
/// @author Jonathan Jaquette
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_AUTODIFF_EVAL_MITTAG_LEFFLER_12_H
#define CAPD_AUTODIFF_EVAL_MITTAG_LEFFLER_12_H

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- MittagLeffler12 ------------------------------------
//
// Computes E_{1,2}(u) = (exp(u) - 1) / u,  E_{1,2}(0) = 1
//
// Data layout:
//   left  = input u coefficients
//   right = exp(u) companion coefficients (written by this node)
//   result = E_{1,2}(u) coefficients
//
// Identity: u * E_{1,2}(u) = exp(u) - 1
// Recurrences (PDF equations 4 and 5):
//   e_0 = exp(u_0),  e_k = (1/k) sum_{j=1}^{k} j * u_j * e_{k-j}
//   v_0 = (e_0 - 1) / u_0,
//   v_k = (1/u_0) * (e_k - sum_{j=1}^{k} u_j * v_{k-j})
//
// When u_0 = 0 (degenerate case), evalC0 uses the shifted division:
//   v(t) = (exp(u(t)) - 1) / u(t) rewritten as a formal power series
//   division with numerator and denominator both starting at degree 1.

namespace MittagLeffler12
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    if(coeffNo == 0)
    {
      *right = Math<T>::_exp(*left);
      *result = Math<T>::_mittagLeffler12(*left);
    }
    else
    {
      // Step 1: compute exp companion coefficient (standard exp recurrence)
      T tempExp = TypeTraits<T>::zero();
      for(unsigned j=1;j<=coeffNo;++j)
        tempExp += Real(j) * left[j] * right[coeffNo-j];
      right[coeffNo] = tempExp / Real(coeffNo);

      // Step 2: compute E_{1,2} coefficient
      // From u * v = exp(u) - 1, Cauchy product at index k:
      //   left[0]*result[k] = right[k] - sum_{j=1}^{k} left[j]*result[k-j]
      T temp = right[coeffNo];
      for(unsigned j=1;j<=coeffNo;++j)
        temp -= left[j] * result[coeffNo-j];
      result[coeffNo] = temp / left[0];
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *right = Math<T>::_exp(*left);
    *result = Math<T>::_mittagLeffler12(*left);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    // Pointer setup: leftDer[k] = du/dx_derNo at time k
    // rightDer[k] and resultDer[k] similarly for exp companion and E_{1,2}
    const T* leftDer = left + order;
    T* rightDer = right + order;
    R resultDer = result + order;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      if(getMask(resultDer,coeffNo)){
        // Step 1: exp companion C1 (chain rule: d[exp(u)]/dx_i = exp(u) * du/dx_i)
        T tempExp = TypeTraits<T>::zero();
        for(unsigned j=0;j<=coeffNo;++j)
          tempExp += leftDer[coeffNo-j] * right[j];
        rightDer[coeffNo] = tempExp;

        // Step 2: E_{1,2} C1 from Cauchy product identity u*v = exp(u) - 1
        // Differentiating w.r.t. x_i at time coefficient k:
        //   left[0]*resultDer[k] = rightDer[k] - sum_{j=0}^{k} leftDer[j]*result[k-j]
        //                          - sum_{j=1}^{k} left[j]*resultDer[k-j]
        T temp = tempExp;
        for(unsigned j=0;j<=coeffNo;++j)
          temp -= leftDer[j] * result[coeffNo-j];
        for(unsigned j=1;j<=coeffNo;++j)
          temp -= left[j] * resultDer[coeffNo-j];
        resultDer[coeffNo] = temp / left[0];
      }
    }
  }

  template<class T, class R>
  void evalMultiindex(const T* left, T* right, R result, DagIndexer<T>* dag, const MultiindexData* i, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    if(getMask(result,i->index))
    {
      // Step 1: compute exp companion at this multiindex (Faa di Bruno formula)
      T te = TypeTraits<T>::zero();
      int he = i->getConvolutionPairsFromEpToK(coeffNo).size();
      for(int q=0;q<he;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairsFromEpToK(coeffNo)[q];
        te += Real(dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p])*left[p.first]*right[p.second];
      }
      right[i->index+coeffNo] = te/Real(i->k[i->p]);

      // Step 2: compute E_{1,2} at this multiindex (Cauchy product: left*result = right)
      // left[0]*result_alpha = right_alpha - left_alpha*result_0
      //                        - sum_{0<beta<alpha} (left_beta*result_{alpha-beta} + left_{alpha-beta}*result_beta)
      //                        - (middle term if odd count)
      T t = right[i->index+coeffNo] - left[i->index+coeffNo]*(*result);
      int h = i->getConvolutionPairs(coeffNo).size();
      for(int q=1;q<h/2;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        t -= left[p.first]*result[p.second] + left[p.second]*result[p.first];
      }
      if(h & 1){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[h/2];
        t -= left[p.first]*result[p.second];
      }
      result[i->index+coeffNo] = t / (*left);
    }
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
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
  void evalHomogenousPolynomial(const unsigned degree, const T* left, T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
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

// -------------------- MittagLeffler12Const ------------------------------------

namespace MittagLeffler12Const
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
    {
      *right = Math<T>::_exp(*left);
      *result = Math<T>::_mittagLeffler12(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *right = Math<T>::_exp(*left);
    *result = Math<T>::_mittagLeffler12(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- MittagLeffler12Time ------------------------------------
// Argument is the time variable itself: u = t
// E_{1,2}(t) coefficients are computed via the recurrence with left[0] = t_0 and left[k] = delta(k,1)

namespace MittagLeffler12Time
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    typedef typename TypeTraits<T>::Real Real;
    if(coeffNo == 0)
    {
      *right = Math<T>::_exp(*left);
      *result = Math<T>::_mittagLeffler12(*left);
    }
    else
    {
      // For time variable: left[k] = delta(k,1) when u = t
      // Exp recurrence: right[k] = (1/k) * sum j*left[j]*right[k-j]
      // Only j=1 contributes (left[j]=0 for j!=1 in the time variable case,
      // but we use the general formula since left could be shifted time t+c)
      T tempExp = TypeTraits<T>::zero();
      for(unsigned j=1;j<=coeffNo;++j)
        tempExp += Real(j) * left[j] * right[coeffNo-j];
      right[coeffNo] = tempExp / Real(coeffNo);

      T temp = right[coeffNo];
      for(unsigned j=1;j<=coeffNo;++j)
        temp -= left[j] * result[coeffNo-j];
      result[coeffNo] = temp / left[0];
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *right = Math<T>::_exp(*left);
    *result = Math<T>::_mittagLeffler12(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- MittagLeffler12FunTime ------------------------------------
// Argument depends only on time (not spatial variables)

namespace MittagLeffler12FunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    MittagLeffler12::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *right = Math<T>::_exp(*left);
    *result = Math<T>::_mittagLeffler12(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    MittagLeffler12::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- ExpForML12 ------------------------------------
// Companion node: stores exp(u) for the MittagLeffler12 node.
// This node is NOT evaluated independently — it's a data slot written by the ML12 node.
// It is filtered from the eval path in BasicFunction::createEvalPath().
// We still need a node class so that eval.hpp can create it without throwing.

namespace ExpForML12
{
  template<class T, class R>
  inline void evalC0(const T* /*left*/, T* /*right*/, R /*result*/, const unsigned /*coeffNo*/)
  {}

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* /*left*/, T* /*right*/, R /*result*/)
  {}

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* /*left*/, T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_DAG_NODE(MittagLeffler12);
CAPD_MAKE_DAG_NODE(MittagLeffler12Const);
CAPD_MAKE_DAG_NODE(MittagLeffler12Time);
CAPD_MAKE_DAG_NODE(MittagLeffler12FunTime);
CAPD_MAKE_DAG_NODE(ExpForML12);
/// @}
}} // namespace capd::autodiff

#endif // CAPD_AUTODIFF_EVAL_MITTAG_LEFFLER_12_H
