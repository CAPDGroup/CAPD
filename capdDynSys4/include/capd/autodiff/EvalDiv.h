/////////////////////////////////////////////////////////////////////////////
/// @file EvalDiv.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_DIV_H_
#define _CAPD_AUTODIFF_EVAL_DIV_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Div  -------------------------------

namespace Div
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    T temp = left[coeffNo];
    for(unsigned j=0;j<coeffNo;++j)
      temp -= result[j] * right[coeffNo-j];
    result[coeffNo] = temp/(*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left)/(*right);
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const T* leftDer = left + order + coeffNo;
    const T* rightDer = right + order;
    R resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,rightDer+=order,resultDer+=order)
    {
      if(getMask(resultDer)){
        T temp = (*leftDer) - result[coeffNo] * (*rightDer);
        for(unsigned j=0;j<coeffNo;++j)
          temp -= ( result[j]* rightDer[coeffNo-j] + resultDer[j] * right[coeffNo-j] );
        resultDer[coeffNo] = temp/(*right);
      }
    }
  }

  /// hand optimized code for second order jet propagation of division
  template<class T, class R>
  inline void evalC2HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const unsigned s = dim*order;
    const T* rightDer = right+order;
    R resultDer = result + order;
    const T* leftHess = left + order + s + coeffNo;
    const T* rightHess = rightDer + s;
    R resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,rightDer+=order,resultDer+=order)
    {
      // case dxdx
      if(getMask(resultHess)){
        T temp = (*leftHess) - result[coeffNo]*(*rightHess) - resultDer[coeffNo]*(*rightDer);
        for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          temp -= ( result[i]* rightHess[j] + resultDer[i]*rightDer[j] + resultHess[i] * right[j] );
        resultHess[coeffNo] = temp/(*right);
      }
      leftHess += order;
      rightHess += order;
      resultHess += order;

      // case dxdy
      R resultDer2 = resultDer + order;
      const T* rightDer2 = rightDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,resultDer2+=order,rightDer2+=order,leftHess+=order,rightHess+=order,resultHess+=order)
      {
        if(getMask(resultHess)){
          T temp = (*leftHess) - result[coeffNo]*(*rightHess) - resultDer2[coeffNo]*(*rightDer) - resultDer[coeffNo]*(*rightDer2);
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp -= result[i] * rightHess[j];
            temp -= resultDer[i] * rightDer2[j];
            temp -= resultDer2[i] * rightDer[j];
            temp -= resultHess[i] * right[j];
          }
          resultHess[coeffNo] = temp/(*right);
        }
      }
    }
  }

  /// hand optimized code for third order jet propagation of division
  template<class T, class R>
  inline void evalC3HomogenousPolynomial(const T* left, const T* right, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      if(getMask(result,i111)){
        T temp = left[i111+coeffNo] - result[coeffNo]*right[i111] - result[i1+coeffNo]*right[i11] - result[i11+coeffNo]*right[i1];
        // case dxdxdx
        for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
        {
          temp -= result[i]* right[i111+j];
          temp -= result[i1+i] * right[i11+j];
          temp -= result[i11+i] * right[i1+j];
          temp -= result[i111+i]* right[j];
        }
        result[i111+coeffNo] = temp/(*right);
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
          T temp = left[i112+coeffNo]
                 - result[coeffNo]     * right[i112]
                 - result[i1+coeffNo]  * right[i12]
                 - result[i2+coeffNo]  * right[i11]
                 - result[i12+coeffNo] * right[i1]
                 - result[i11+coeffNo] * right[i2];
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp -= result[i]*right[i112+j]    // 0,xxy
                  + result[i1+i]*right[i12+j]  // x,xy
                  + result[i2+i]*right[i11+j]  // y,xx
                  + result[i11+i]*right[i2+j]  // xx,y
                  + result[i12+i]*right[i1+j]  // xy,x
                  + result[i112+i]*right[j];    // xxy,0
          }
          result[i112+coeffNo] = temp/(*right);
        }

        if(getMask(result,i122)){
          T temp2 = left[i122+coeffNo]
                  - result[coeffNo]     * right[i122]
                  - result[i1+coeffNo]  * right[i22]
                  - result[i2+coeffNo]  * right[i12]
                  - result[i12+coeffNo] * right[i2]
                  - result[i22+coeffNo] * right[i1];
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp2 -= result[i]*right[i122+j]   // 0,xyy
                   + result[i1+i]*right[i22+j] // x,yy
                   + result[i2+i]*right[i12+j] // y,xy
                   + result[i12+i]*right[i2+j] // xy,y
                   + result[i22+i]*right[i1+j] // yy,x
                   + result[i122+i]*right[j];  // xyy,0
          }
          result[i122+coeffNo] = temp2/(*right);
        }

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          if(getMask(result,i123)){
            T temp = left[i123+coeffNo]
                   - result[coeffNo]     * right[i123]
                   - result[i1+coeffNo]  * right[i23]
                   - result[i2+coeffNo]  * right[i13]
                   - result[i3+coeffNo]  * right[i12]
                   - result[i12+coeffNo] * right[i3]
                   - result[i13+coeffNo] * right[i2]
                   - result[i23+coeffNo] * right[i1];
            for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
            {
              temp -= result[i] * right[i123+j]    // 0,xyz
                    + result[i1+i] * right[i23+j]  // x,yz
                    + result[i2+i] * right[i13+j]  // y,xz
                    + result[i3+i] * right[i12+j]  // z,xy
                    + result[i12+i] * right[i3+j]  // xy,z
                    + result[i13+i] * right[i2+j]  // xz,y
                    + result[i23+i]*right[i1+j]    // yz,x
                    + result[i123+i]*right[j];     // xyz,0
            }
            result[i123+coeffNo] = temp/(*right);
          }
        }
      }
    }
  }

  template<class T, class R>
  void evalMultiindex(const T* left, const T* right, R result, const MultiindexData* i, const unsigned coeffNo)
  {
    if(getMask(result,i->index))
    {
      T t = left[i->index+coeffNo]-(*result)*right[i->index+coeffNo];
      int h = i->getConvolutionPairs(coeffNo).size();

      for(int q=1;q<h/2;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        t -= result[p.first]*right[p.second] + result[p.second]*right[p.first];
      }
      if(h & 1){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[h/2];
        t -= result[p.first]*right[p.second];
      }
      result[i->index+coeffNo] = t/(*right);
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

// -------------------- DivVarByConst  -------------------------------

namespace DivVarByConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo]/(*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left)/(*right);
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned shift)
  {
    left += shift;
    result += shift;
    const T* end = left + dataSize*order;
    for(;left!=end; left+=order,result+=order)
      if(getMask(result))
        *result = (*left) / (*right);
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

// -------------------- DivVarByFunTime  -------------------------------

namespace DivVarByFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Div::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left)/(*right);
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned coeffNo)
  {
    const T* end = left + dataSize*order;
    for(;left!=end; left+=order,result+=order)
      if(getMask(result))
        Div::evalC0(left,right,result,coeffNo);
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
    evalHelper(left+shift,right,result+shift,dataSize,order,coeffNo);
  }
}

// -------------------- DivVarByTime  -------------------------------

namespace DivVarByTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = (left[coeffNo]-result[coeffNo-1])/(*right);
    else
      *result = (*left)/(*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left)/(*right);
  }

  template<class T, class R>
  inline void evalHelper(const T* left, const T* right, R result, const unsigned dataSize, const unsigned order, const unsigned coeffNo)
  {
    const T* end = left + dataSize*order;
    for(;left!=end; left+=order,result+=order)
      if(getMask(result))
        result[coeffNo] = (left[coeffNo]-result[coeffNo-1])/(*right);
  }

  template<class T, class R>
  inline void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(coeffNo)
      evalHelper(left,right,result,binomial(dag->domainDimension()+degree,degree),dag->getOrder()+1,coeffNo);
    else
      DivVarByConst::eval(degree,left,right,result,dag,0);
  }


  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    const unsigned dim = dag->domainDimension();
    const unsigned order = dag->getOrder()+1;
    const unsigned shift = binomial(dim+degree-1,dim)*order;
    const unsigned dataSize = binomial(dim+degree-1,degree);
    if(coeffNo)
      evalHelper(left+shift,right,result+shift,dataSize,order,coeffNo);
    else
      DivVarByConst::evalHelper(left,right,result,dataSize,order,shift);
  }
}

// -------------------- DivTimeByConst -------------------------------

namespace DivTimeByConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo<=1)
      result[coeffNo] = left[coeffNo] / (*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) / (*right);
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
// -------------------- DivFunTimeByConst -------------------------------

namespace DivFunTimeByConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    result[coeffNo] = left[coeffNo] / (*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) / (*right);
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

// -------------------- DivFunTimeByTime -------------------------------

namespace DivFunTimeByTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    DivVarByTime::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) / (*right);
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

// -------------------- DivFunTimeByFunTime -------------------------------

namespace DivFunTimeByFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Div::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) / (*right);
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

// -------------------- DivConstByConst -------------------------------

namespace DivConstByConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
    {}
    else
      *result = (*left) / (*right);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = (*left) / (*right);
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

CAPD_MAKE_DAG_NODE(Div);
CAPD_MAKE_DAG_NODE(DivVarByConst);
CAPD_MAKE_DAG_NODE(DivVarByTime);
CAPD_MAKE_DAG_NODE(DivVarByFunTime);
CAPD_MAKE_DAG_NODE(DivTimeByConst);
CAPD_MAKE_DAG_NODE(DivFunTimeByConst);
CAPD_MAKE_DAG_NODE(DivFunTimeByTime);
CAPD_MAKE_DAG_NODE(DivFunTimeByFunTime);
CAPD_MAKE_DAG_NODE(DivConstByConst);
/// @}
}} // namespace capd::autodiff

#endif
