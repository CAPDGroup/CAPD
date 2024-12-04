/////////////////////////////////////////////////////////////////////////////
/// @file EvalMul.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_MUL_H_
#define _CAPD_AUTODIFF_EVAL_MUL_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Mul ------------------------------------

namespace Mul
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result += coeffNo;
    *result = (*left)*right[coeffNo];
    for(unsigned i=1;i<=coeffNo;++i)
      *result += left[i]* right[coeffNo-i];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    (*result) = (*left)*(*right);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const T* leftDer = left+order;
    const T* rightDer = right+order;
    result += (order+coeffNo);
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,result+=order)
    {
      if(getMask(result)){
        *result = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          *result += (left[i]*rightDer[j] + leftDer[i]*right[j]);
      }
    }
  } // evalC1HomogenousPolynomial

  /// hand optimized code for second order jet propagation of multiplication
  template<class T, class R>
  inline void evalC2HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    // begin of C^1
    const T* leftDer = left+order;
    const T* rightDer = right+order;

    // begin of C^2
    const unsigned s = (dim+1)*order;
    const T* leftHess = left + s;
    const T* rightHess = right + s;
    result += coeffNo + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order)
    {
      // case dx^2
      if(getMask(result)){
        *result  = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          *result += (left[i]*rightHess[j] + leftDer[i]*rightDer[j] + leftHess[i]*right[j]);
      }

      leftHess += order;
      rightHess += order;
      result += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      const T* rightDer2 = rightDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,rightDer2+=order,leftHess+=order,rightHess+=order,result+=order)
      {
        if(getMask(result)){
          *result = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            *result += left[i] * rightHess[j];
            *result += leftDer[i] * rightDer2[j];
            *result += leftDer2[i] * rightDer[j];
            *result += leftHess[i] * right[j];
          }
        }
      }
    }
  }  // evalC2HomogenousPolynomial

  /// hand optimized code for third order jet propagation of multiplication
  template<class T, class R>
  inline void evalC3HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    result += coeffNo;
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      // case dxdxdx
      if(getMask(result,i111)){
        T temp = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          temp += left[i]* right[i111+j];
          temp += left[i1+i] * right[i11+j];
          temp += left[i11+i] * right[i1+j];
          temp += left[i111+i]* right[j];
        }
        result[i111] = temp;
      }

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111+order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;

        if(getMask(result,i112)){
          T temp = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            temp += left[i]*right[i112+j];    // 0,xxy
            temp += left[i1+i]*right[i12+j];  // x,xy
            temp += left[i2+i]*right[i11+j];  // y,xx
            temp += left[i11+i]*right[i2+j];  // xx,y
            temp += left[i12+i]*right[i1+j];  // xy,x
            temp += left[i112+i]*right[j];    // xxy,0
          }
          result[i112] = temp;
        }

        if(getMask(result,i122)){
          T temp = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {

            temp += left[i]*right[i122+j];   // 0,xyy
            temp += left[i1+i]*right[i22+j]; // x,yy
            temp += left[i2+i]*right[i12+j]; // y,xy
            temp += left[i12+i]*right[i2+j]; // xy,y
            temp += left[i22+i]*right[i1+j]; // yy,x
            temp += left[i122+i]*right[j];   // xyy,0
          }
          result[i122] = temp;
        }

        // case dxdydz, assume x<y<z

        unsigned i3 = i2+order;
        unsigned i23 = i22+order;
        unsigned i13 = i12+order;
        unsigned i123 = i122+order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,
                      i3+=order,i123+=order,i23+=order,i13+=order)
        {
          if(getMask(result,i123+coeffNo)){
            T temp = TypeTraits<T>::zero();
            for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
            {
              temp += left[i] * right[i123+j];    // 0,xyz
              temp += left[i1+i] * right[i23+j];  // x,yz
              temp += left[i2+i] * right[i13+j];  // y,xz
              temp += left[i3+i] * right[i12+j];  // z,xy
              temp += left[i12+i] * right[i3+j];  // xy,z
              temp += left[i13+i] * right[i2+j];  // xz,y
              temp += left[i23+i]*right[i1+j];    // yz,x
              temp += left[i123+i]*right[j];      // xyz,0
            }
            result[i123] = temp;
          }
        }
      }
    }
  }  // evalC3HomogenousPolynomial

  template<class T, class R>
  void evalMultiindex(const T* left, const T* right, R result, const MultiindexData* i, const unsigned coeffNo)
  {
    if(getMask(result,i->index))
    {
      T t = capd::TypeTraits<T>::zero();
      int h = i->getConvolutionPairs(coeffNo).size();
      int q=0;
      for(;q<h/2;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        t += left[p.first]*right[p.second] + left[p.second]*right[p.first];
      }
      if(h & 1){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        t += left[p.first]*right[p.second];
      }
      result[i->index+coeffNo] = t;
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
      case 2:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC3HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC3HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        for(const MultiindexData* i = dag->getIndexArray().begin(0,4);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,i,coeffNo);
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
      case 2:
        evalC2HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC3HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,i,coeffNo);
    }
  }
}

// -------------------- MulConstByVar -------------------------------

namespace MulConstByVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = (*left) * right[coeffNo];
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    right += shift;
    result += shift;
    const T* end = right + dataSize*order;
    for(;right!=end; right+=order,result+=order)
      if(getMask(result))
        *result = (*left) * (*right);
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

// -------------------- MulConstByConst -------------------------------

namespace MulConstByConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo){
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- MulConstByFunTime -------------------------------

namespace MulConstByFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = (*left) * right[coeffNo];
  }


  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo){
    evalC0(left,right,result,coeffNo);
  }


  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- MulConstByTime -------------------------------

namespace MulConstByTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 0:
        *result = (*left) * (*right); break;
      case 1:
        result[1] = *left;
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo){
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- MulTimeByVar -------------------------------

namespace MulTimeByVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
    case 0:
      *result = (*left) * (*right);
      break;
    default:
      result[coeffNo] = (*left) * right[coeffNo] + right[coeffNo-1];
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    right += shift;
    result += shift;
    const T* end = right + dataSize*order;
    for(;right!=end; right+=order,result+=order)
      if(getMask(result))
          *result = (*left) * (*right) + *(right-1);
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(coeffNo)
      evalHelper(left,right,result,binomial(dag->domainDimension()+degree,degree),dag->getOrder()+1,coeffNo);
    else
      MulConstByVar::eval(degree,left,right,result,dag,0);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(coeffNo) {
      const unsigned dim = dag->domainDimension();
      const unsigned order = dag->getOrder()+1;
      const unsigned shift = binomial(dim+degree-1,dim)*order;
      const unsigned dataSize = binomial(dim+degree-1,degree);
      evalHelper(left,right,result,dataSize,order,coeffNo+shift);
    }
    else
      MulConstByVar::evalHomogenousPolynomial(degree,left,right,result,dag,0);
  }
}

// -------------------- MulTimeByFunTime -------------------------------

namespace MulTimeByFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    MulTimeByVar::evalC0(left,right,result,coeffNo);
  }


  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    MulTimeByVar::evalC0(left,right,result,coeffNo);
  }


  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}


// -------------------- MulFunTimeByVar -------------------------------

namespace MulFunTimeByVar
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Mul::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift, const unsigned coeffNo)
  {
    right += shift;
    result += shift+coeffNo;
    const T* end = right + dataSize*order;
    for(;right!=end; right+=order,result+=order){
      if(getMask(result)){
        *result = TypeTraits<T>::zero();
        for(unsigned i=0;i<=coeffNo;++i)
          *result += left[i]* right[coeffNo-i];
      }
    }
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    evalHelper(left,right,result,binomial(dim+degree,degree),order,0,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo){
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    const unsigned shift = binomial(dim+degree-1,dim)*order;
    const unsigned dataSize = binomial(dim+degree-1,degree);
    evalHelper(left,right,result,dataSize,order,shift,coeffNo);
  }
}

// -------------------- MulFunTimeByFunTime -------------------------------

namespace MulFunTimeByFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Mul::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) * (*right);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Mul::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_DAG_NODE(Mul);
CAPD_MAKE_DAG_NODE(MulConstByVar);
CAPD_MAKE_DAG_NODE(MulConstByConst);
CAPD_MAKE_DAG_NODE(MulConstByTime);
CAPD_MAKE_DAG_NODE(MulConstByFunTime);
CAPD_MAKE_DAG_NODE(MulTimeByVar);
CAPD_MAKE_DAG_NODE(MulTimeByFunTime);
CAPD_MAKE_DAG_NODE(MulFunTimeByVar);
CAPD_MAKE_DAG_NODE(MulFunTimeByFunTime);
/// @}
}} // namespace capd::autodiff

#endif
