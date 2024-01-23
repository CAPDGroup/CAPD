/////////////////////////////////////////////////////////////////////////////
/// @file EvalSinCos.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_SIN_COS_H_
#define _CAPD_AUTODIFF_EVAL_SIN_COS_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Sin ------------------------------------

namespace Sin
{
  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      T tempSin = capd::TypeTraits<T>::zero();
      T tempCos = capd::TypeTraits<T>::zero();
      for(unsigned j=1;j<=coeffNo;++j)
      {
        tempSin += double(j) * right[coeffNo-j] * left[j];
        tempCos -= double(j) * result[coeffNo-j] * left[j];
      }
      result[coeffNo] = tempSin/(double)coeffNo;
      right[coeffNo] = tempCos/(double)coeffNo;
    }
    else
    {
      (*result) = sin(*left);
      (*right) = cos(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    (*result) = sin(*left);
    (*right) = cos(*left);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const T* leftDer = left+order;
    T* rightDer = right+order+coeffNo;
    R resultDer = result+order+coeffNo;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      if(getMask(resultDer)){
        T tempSin = capd::TypeTraits<T>::zero();
        T tempCos = capd::TypeTraits<T>::zero();
        for(unsigned j=0;j<=coeffNo;++j)
        {
          tempSin += leftDer[coeffNo-j]*right[j];
          tempCos -= leftDer[coeffNo-j]*result[j];
        }
        *resultDer = tempSin;
        *rightDer = tempCos;
      }
    }
  } // evalC1HomogenousPolynomial

  template<class T, class R>
  inline void evalC2HomogenousPolynomial(const T* left, T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const unsigned s = dim*order;
    const T* leftDer = left + order;
    T* rightDer = right + order;
    R resultDer = result + order;
    const T* leftHess = leftDer + s;
    T* rightHess = rightDer + s;
    R resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      // case dx^2
      if(getMask(resultHess)){
        T tempSin1 = TypeTraits<T>::zero();
        T tempSin2 = TypeTraits<T>::zero();
        T tempCos1 = TypeTraits<T>::zero();
        T tempCos2 = TypeTraits<T>::zero();
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          tempSin1 += leftDer[i]*rightDer[j];
          tempSin2 += leftHess[i]*right[j];
          tempCos1 -= leftDer[i]*resultDer[j];
          tempCos2 -= leftHess[i]*result[j];
        }
        resultHess[coeffNo] = 0.5*tempSin1 + tempSin2;
        rightHess[coeffNo]  = 0.5*tempCos1 + tempCos2;
      }

      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      T* rightDer2 = rightDer + order;
      R resultDer2 = resultDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,rightDer2+=order,resultDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
      {
        if(getMask(resultHess)){
          T tempSin1 = TypeTraits<T>::zero();
          T tempCos1 = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            tempSin1 += leftDer[i] * rightDer2[j];
            tempSin1 += leftHess[i] * right[j];
            tempCos1 -= leftDer[i] * resultDer2[j];
            tempCos1 -= leftHess[i] * result[j];
          }
          resultHess[coeffNo] = tempSin1;
          rightHess[coeffNo] = tempCos1;
        }
      }
    }
  }  // evalC2HomogenousPolynomial

  /// hand optimized code for third order jet propagation of sine
  template<class T, class R>
  inline void evalC3HomogenousPolynomial(const T* left, T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      T tempSin1 = TypeTraits<T>::zero(),
        tempSin2 = TypeTraits<T>::zero(),
        tempSin3 = TypeTraits<T>::zero();
      T tempCos1 = TypeTraits<T>::zero(),
        tempCos2 = TypeTraits<T>::zero(),
        tempCos3 = TypeTraits<T>::zero();
      // case dxdxdx
      if(getMask(result,i111)){
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
        {
          tempSin1 += left[i1+i] * right[i11+j];
          tempSin2 += left[i11+i] * right[i1+j];
          tempSin3 += left[i111+i]* right[j];

          tempCos1 -= left[i1+i] * result[i11+j];
          tempCos2 -= left[i11+i] * result[i1+j];
          tempCos3 -= left[i111+i]* result[j];
        }
        result[i111+coeffNo] = (tempSin1+2.*tempSin2)/3. + tempSin3;
        right[i111+coeffNo] = (tempCos1+2.*tempCos2)/3. + tempCos3;
      }

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        if(getMask(result,i122)){
          result[i122+coeffNo] = TypeTraits<T>::zero();
          right[i122+coeffNo] = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            result[i122+coeffNo] += left[i1+i]*right[i22+j];  // x,yy
            result[i122+coeffNo] += left[i12+i]*right[i2+j];  // xy,y
            result[i122+coeffNo] += left[i122+i]*right[j];    // xyy,0

            right[i122+coeffNo] -= left[i1+i]*result[i22+j];  // x,yy
            right[i122+coeffNo] -= left[i12+i]*result[i2+j];  // xy,y
            right[i122+coeffNo] -= left[i122+i]*result[j];    // xxy,0
          }
        }

        if(getMask(result,i112)){
          result[i112+coeffNo] = TypeTraits<T>::zero();
          right[i112+coeffNo] = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            result[i112+coeffNo] += left[i2+i]*right[i11+j]; // y,xx
            result[i112+coeffNo] += left[i12+i]*right[i1+j]; // xy,x
            result[i112+coeffNo] += left[i112+i]*right[j];   // xxy,0

            right[i112+coeffNo] -= left[i2+i]*result[i11+j]; // y,xx
            right[i112+coeffNo] -= left[i12+i]*result[i1+j]; // xy,x
            right[i112+coeffNo] -= left[i112+i]*result[j];   // xxy,0
          }
        }

        // case dxdydz, assume x<y<z
        unsigned i3 = i2 + order;
        unsigned i123 = i122 + order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          if(getMask(result,i123)){
            result[i123+coeffNo] = TypeTraits<T>::zero();
            right[i123+coeffNo] = TypeTraits<T>::zero();
            for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
            {
              result[i123+coeffNo] += left[i1+i]   * right[i23+j];  // x,yz
              result[i123+coeffNo] += left[i12+i]  * right[i3+j];  // xy,z
              result[i123+coeffNo] += left[i13+i]  * right[i2+j];  // xz,y
              result[i123+coeffNo] += left[i123+i] * right[j];      // xyz,0

              right[i123+coeffNo] -= left[i1+i]   * result[i23+j];  // x,yz
              right[i123+coeffNo] -= left[i12+i]  * result[i3+j];  // xy,z
              right[i123+coeffNo] -= left[i13+i]  * result[i2+j];  // xz,y
              right[i123+coeffNo] -= left[i123+i] * result[j];      // xyz,0
            }
          }
        }
      }
    }
  }  // evalC3HomogenousPolynomial

  template<class T, class R>
  void evalMultiindex(const T* left, T* right, R result, DagIndexer<T>* dag, const MultiindexData* i, const unsigned coeffNo)
  {
    if(getMask(result,i->index))
    {
      T ts = capd::TypeTraits<T>::zero();
      T tc = capd::TypeTraits<T>::zero();
      int h = i->getConvolutionPairsFromEpToK(coeffNo).size();

      for(int q=0;q<h;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairsFromEpToK(coeffNo)[q];
        ts += ((double)dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p])*left[p.first]*right[p.second];
        tc -= ((double)dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p])*left[p.first]*result[p.second];
      }
      right[i->index+coeffNo] = tc/(double)i->k[i->p];
      result[i->index+coeffNo] = ts/(double)i->k[i->p];
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
      case 2:
        evalC2HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC3HomogenousPolynomial(left,right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,right,result,dag,i,coeffNo);
    }
  }
}

namespace SinConst{

  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
    {
      *result = sin(*left);
      *right = cos(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *result = sin(*left);
    *right = cos(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

namespace SinTime{

  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {
      result[coeffNo] = right[coeffNo-1]/(double)coeffNo;
      right[coeffNo] = -result[coeffNo-1]/(double)coeffNo;
    }
    else
    {
      (*result) = sin(*left);
      (*right) = cos(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *result = sin(*left);
    *right = cos(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

namespace SinFunTime{

  template<class T, class R>
  inline void evalC0(const T* left, T* right, R result, const unsigned coeffNo)
  {
    Sin::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, T* right, R result)
  {
    *result = sin(*left);
    *right = cos(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// ----------------------------------------------------------------------------------

//use macro to define classes

CAPD_MAKE_DAG_NODE(Sin);
CAPD_MAKE_DAG_NODE(SinConst);
CAPD_MAKE_DAG_NODE(SinTime);
CAPD_MAKE_DAG_NODE(SinFunTime);
/// @}
}} // namespace capd::autodiff

#endif
