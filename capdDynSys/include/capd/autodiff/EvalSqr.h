/////////////////////////////////////////////////////////////////////////////
/// @file EvalSqr.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_SQR_H_
#define _CAPD_AUTODIFF_EVAL_SQR_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Sqr ------------------------------------

namespace Sqr{

  template<class T>
  inline T sqrProduct(const T* x, const unsigned n)
  {
    T temp = TypeTraits<T>::zero();
    const unsigned p = n/2;
    for(unsigned j=0;j<p;++j)
      temp += x[j] * x[n-j];
    return temp + ( (n%2) ? x[p]*x[n-p] : 0.5*sqr(x[p]) );
  }

  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    result[coeffNo] = 2.*sqrProduct(left,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    (*result) = sqr(*left);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, R result, const unsigned dim, const unsigned order, const unsigned coeffNo, typename TypeTraits<T>::Real const c=2.)
  {
    const T* leftDer = left+order;
    result += order+coeffNo;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,result+=order)
    {
      if(getMask(result)){
        T temp = TypeTraits<T>::zero();
        for(unsigned i=0;i<=coeffNo;++i)
          temp += left[i] * leftDer[coeffNo-i];
        *result = c*temp;
      }
    }
  }

  template<class T, class R>
  inline void evalC2HomogenousPolynomial(const T* left, R result, const unsigned dim, const unsigned order, const unsigned coeffNo, typename TypeTraits<T>::Real const c=2.)
  {
    const unsigned s = (dim+1)*order;
    const T* leftDer = left + order;
    const T* leftHess = left + s;
    result += coeffNo + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order)
    {
      // case dx^2
      if(getMask(result)){
        *result = TypeTraits<T>::zero();
        for(unsigned i=0;i<=coeffNo;++i)
          *result += left[i]*leftHess[coeffNo-i];
        *result += sqrProduct(leftDer,coeffNo);
        *result *= c;
      }

      leftHess += order;
      result += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,leftHess+=order,result+=order)
      {
        if(getMask(result)){
          *result = TypeTraits<T>::zero();
          for(unsigned i=0;i<=coeffNo;++i)
            *result += left[i] * leftHess[coeffNo-i] + leftDer[i] * leftDer2[coeffNo-i];
          *result *= c;
        }
      }
    }
  }

  template<class T, class R>
  inline void evalC3HomogenousPolynomial(const T* left, R result, const unsigned dim, const unsigned order, const unsigned coeffNo, typename TypeTraits<T>::Real const c=2.)
  {
    result += coeffNo;
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      if(getMask(result,i111)){
        T temp = TypeTraits<T>::zero();
        // case dxdxdx
        for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          temp += left[i]* left[i111+j] + left[i1+i] * left[i11+j];
        result[i111] = c*temp;
      }

      // cases dxdxdy and dxdydy, assume that x<y
      unsigned i2 = i1+order;
      unsigned i12 = i11+order;
      unsigned i112 = i111 + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,i2+=order,i12+=order,i112+=order)
      {
        const unsigned i22 = capd::autodiff::index(dim,derNo2,derNo2)*order;
        const unsigned i122 = capd::autodiff::index(dim,derNo,derNo2,derNo2)*order;
        if(getMask(result,i112)){
          T temp = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            temp += left[i]*left[i112+j];    // 0,xxy
            temp += left[i1+i]*left[i12+j];  // x,xy
            temp += left[i2+i]*left[i11+j];  // y,xx
          }
          result[i112] = c*temp;
        }

        if(getMask(result,i122)){
          T temp = TypeTraits<T>::zero();
          for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
          {
            temp += left[i]*left[i122+j];   // 0,xyy
            temp += left[i1+i]*left[i22+j]; // x,yy
            temp += left[i2+i]*left[i12+j]; // y,xy
          }
          result[i122] = c*temp;
        }

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          if(getMask(result,i123)){
            T temp = TypeTraits<T>::zero();
            for(unsigned i=0,j=coeffNo;i<=coeffNo;++i,--j)
            {
              temp += left[i] * left[i123+j];    // 0,xyz
              temp += left[i1+i] * left[i23+j];  // x,yz
              temp += left[i2+i] * left[i13+j];  // y,xz
              temp += left[i3+i] * left[i12+j];  // z,xy
            }
            result[i123] = c*temp;
          }
        }
      }
    }
  }

  template<class T, class R>
  void evalMultiindex(const T* left, R result, const MultiindexData* i, const unsigned coeffNo, typename TypeTraits<T>::Real const c=2.)
  {
    if(getMask(result,i->index))
    {
      T t = capd::TypeTraits<T>::zero();
      int h = i->getConvolutionPairs(coeffNo).size();
      int q=0;
      for(;q<h/2;++q)
      {
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        t += left[p.first]*left[p.second];
      }
      if(h & 1)
      {
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        t += 0.5*sqr(left[p.first]);
      }
      result[i->index+coeffNo] = c*t;
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
        evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 2:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC3HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        evalC0(left,right,result,coeffNo);
        evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC3HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        for(const MultiindexData* i = dag->getIndexArray().begin(0,4);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,result,i,coeffNo);
    }
  }


  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* /*right*/, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        evalC1HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 2:
        evalC2HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC3HomogenousPolynomial(left,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,result,i,coeffNo);
    }
  }
}

// -------------------- SqrFunTime ------------------------------------

namespace SqrFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, unsigned coeffNo)
  {
    Sqr::evalC0(left,right,result,coeffNo);
  }

 template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = sqr(*left);
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Sqr::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- SqrTime ------------------------------------

namespace SqrTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    switch(coeffNo)
    {
      case 2: result[2] = 1.; break;
      case 1: result[1] = 2.*(*left); break;
      case 0: *result = sqr(*left);
    }
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = sqr(*left);
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

// -------------------- SqrConst ------------------------------------

namespace SqrConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* /*right*/, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = sqr(*left);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* /*right*/, R result)
  {
    *result = sqr(*left);
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

CAPD_MAKE_DAG_NODE(Sqr);
CAPD_MAKE_DAG_NODE(SqrTime);
CAPD_MAKE_DAG_NODE(SqrFunTime);
CAPD_MAKE_DAG_NODE(SqrConst);
/// @}
}} // namespace capd::autodiff

#endif
