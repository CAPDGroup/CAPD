/////////////////////////////////////////////////////////////////////////////
/// @file EvalPow.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_AUTODIFF_EVAL_POW_H_
#define _CAPD_AUTODIFF_EVAL_POW_H_

#include "capd/autodiff/NodeType.h"

namespace capd{
namespace autodiff{
/// @addtogroup autodiff
/// @{

// -------------------- Pow ------------------------------------

namespace Pow
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    const T c = *right;
    if(coeffNo)
    {
      if(!((*left)>capd::TypeTraits<T>::zero()))
      	throw std::runtime_error("Map evaluation error: negative argument in x^c, where c is not an integer.");
      T temp = capd::TypeTraits<T>::zero();
      for(int j=0;j<(int)coeffNo;++j)
        temp += (c*((int)coeffNo-j)-j) * result[j] *left[coeffNo-j];
      result[coeffNo] = temp/((double)coeffNo * (*left));
    }
    else
      *result = exp(c * log(*left));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    if(!((*left)>capd::TypeTraits<T>::zero()))
    	throw std::runtime_error("Map evaluation error: negative argument in x^c, where c is not an integer.");
    *result = exp((*right) * log(*left));
  }

  template<class T, class R>
  inline void evalC1HomogenousPolynomial(const T* left, const T exponent, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const T* leftDer = left + order;
    R resultDer = result + order;
    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      if(getMask(resultDer)){
        T temp1 = result[coeffNo] * (*leftDer);
        T temp2 = capd::TypeTraits<T>::zero();
        for(unsigned j=0;j<coeffNo;++j)
        {
          temp1 += result[j] * leftDer[coeffNo-j];
          temp2 += resultDer[j] * left[coeffNo-j];
        }
        resultDer[coeffNo] = (temp1*exponent-temp2)/(*left);
      }
    }
  } // evalC1HomogenousPolynomial

  template<class T, class R>
  inline void evalC2HomogenousPolynomial(const T* left, const  T exponent, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    const unsigned s = dim*order;
    // begin of C^1
    const T* leftDer = left+order;
    R resultDer = result + order;
    // begin of C^2
    const T* leftHess = leftDer + s;
    R resultHess = resultDer + s;

    for(unsigned derNo=0;derNo<dim;++derNo,leftDer+=order,resultDer+=order)
    {
      // case dx^2
      if(getMask(resultHess)){
        T temp1 = result[coeffNo] * (*leftHess);
        T temp2 = resultDer[coeffNo] * (*leftDer);
        T temp3 = capd::TypeTraits<T>::zero();
        for(unsigned j=0;j<coeffNo;++j)
        {
          temp1 += result[j] * leftHess[coeffNo-j];
          temp2 += resultDer[j] * leftDer[coeffNo-j];
          temp3 += resultHess[j] * left[coeffNo-j];
        }
        resultHess[coeffNo] = (exponent*temp1+0.5*(exponent-1.)*temp2-temp3)/(*left);
      }

      leftHess += order;
      resultHess += order;

      // case dxdy
      const T* leftDer2 = leftDer + order;
      R resultDer2 = resultDer + order;
      for(unsigned derNo2=derNo+1;derNo2<dim;++derNo2,leftDer2+=order,resultDer2+=order,leftHess+=order,resultHess+=order)
      {
        if(getMask(resultHess)){
          T temp1 = result[coeffNo] * (*leftHess) + resultDer[coeffNo] *(*leftDer2);
          T temp2 = resultDer2[coeffNo] *(*leftDer);
          for(unsigned j=0;j<coeffNo;++j)
          {
            temp1 += result[j] * leftHess[coeffNo-j] + resultDer[j] * leftDer2[coeffNo-j];
            temp2 += resultDer2[j] * leftDer[coeffNo-j] + resultHess[j]*left[coeffNo-j];
          }
          resultHess[coeffNo] = (exponent*temp1-temp2)/(*left);
        }
      }
    }
  } // evalC2HomogenousPolynomial

  template<class T, class R>
  inline void evalC3HomogenousPolynomial(const T* left, const T exponent, R result, const unsigned dim, const unsigned order, const unsigned coeffNo)
  {
    unsigned i1 = order;
    for(unsigned derNo=0;derNo<dim;++derNo,i1+=order)
    {
      const unsigned i11 = capd::autodiff::index(dim,derNo,derNo)*order;
      const unsigned i111 = capd::autodiff::index(dim,derNo,derNo,derNo)*order;
      if(getMask(result,i111)){
        T temp1 = result[coeffNo]*left[i111];
        T temp2 = result[i1+coeffNo]*left[i11];
        T temp3 = result[i11+coeffNo]*left[i1];
        T temp4 = capd::TypeTraits<T>::zero();
        // case dxdxdx
        for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
        {
          temp1 += result[i]*left[i111+j];
          temp2 += result[i1+i]*left[i11+j];
          temp3 += result[i11+i]*left[i1+j];
          temp4 += result[i111+i]*left[j];
        }
        result[i111+coeffNo] = (exponent*temp1 + (2.*exponent-1.)*temp2/3. + (exponent-2.)*temp3/3. - temp4)/(*left);
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
          T temp1 = result[coeffNo]* left[i112] + result[i1+coeffNo]*left[i12] + result[i11+coeffNo]*left[i2];
          T temp2 = result[i12+coeffNo]  * left[i1] + result[i2+coeffNo]* left[i11];
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp1 += result[i]* left[i112+j] + result[i1+i]*left[i12+j] + result[i11+i]*left[i2+j];
            temp2 += result[i12+i]  * left[i1+j] + result[i2+i]* left[i11+j] + result[i112+i]*left[j];
          }
          result[i112+coeffNo] = (exponent*temp1-temp2)/(*left);
        }

        if(getMask(result,i122)){
          T temp3 = result[coeffNo]* left[i122] + result[i2+coeffNo]*left[i12] + result[i22+coeffNo]*left[i1];
          T temp4 = result[i12+coeffNo]  * left[i2] + result[i1+coeffNo]* left[i22];
          for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
          {
            temp3 += result[i]* left[i122+j] + result[i2+i]*left[i12+j] + result[i22+i]*left[i1+j];
            temp4 += result[i12+i]  * left[i2+j] + result[i1+i]* left[i22+j] + result[i122+i]*left[j];
          }
          result[i122+coeffNo] = (exponent*temp3-temp4)/(*left);
        }

        // case dxdydz, assume x<y<z
        unsigned i3 = i2+order;
        unsigned i123 = i122+order;
        unsigned i23 = i22 + order;
        unsigned i13 = i12 + order;
        for(unsigned derNo3=derNo2+1;derNo3<dim;++derNo3,i3+=order,i123+=order,i23+=order,i13+=order)
        {
          if(getMask(result,i123)){
            T temp1 = result[i12+coeffNo]*left[i3] + result[i2+coeffNo]*left[i13] + result[i1+coeffNo]*left[i23] + result[coeffNo]*left[i123];
            T temp2 = result[i13+coeffNo]*left[i2] + result[i23+coeffNo]*left[i1] + result[i3+coeffNo]*left[i12];
            for(unsigned i=0,j=coeffNo;i<coeffNo;++i,--j)
            {
              temp1 += result[i12+i]*left[i3+j] + result[i2+i]*left[i13+j] + result[i1+i]*left[i23+j] + result[i]*left[i123+j];
              temp2 += result[i13+i]*left[i2+j] + result[i23+i]*left[i1+j] + result[i3+i]*left[i12+j] + result[i123+i]*left[j];
            }
            result[i123+coeffNo] = (exponent*temp1-temp2)/(*left);
          }
        }
      }
    }
  } // evalC3HomogenousPolynomial

  template<class T, class R>
  void evalMultiindex(const T* left, const T exponent, R result, const MultiindexData* i, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    if(getMask(result,i->index))
    {
      const double k = (double)i->k[i->p];
      T t = capd::TypeTraits<T>::zero();
      int h = i->getConvolutionPairs(coeffNo).size();
      int q=1;
      for(q=1;q<h/2;++q){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        unsigned v1 = dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p];
        unsigned v2 = dag->getIndexArray()[p.second/(dag->getOrder()+1)].k[i->p];
        t += (exponent*(k-v1)-v1)*result[p.first]*left[p.second] + (exponent*(k-v2)-v2)*result[p.second]*left[p.first];
      }
      if (h & 1){
        const MultiindexData::IndexPair p = i->getConvolutionPairs(coeffNo)[q];
        unsigned v = dag->getIndexArray()[p.first/(dag->getOrder()+1)].k[i->p];
        t += (exponent*(k-v)-v)*result[p.first]*left[p.second];
      }
      result[i->index+coeffNo] = (exponent*k*left[i->index+coeffNo]*(*result) + t)/(*left*k);
    }
  }

  template<class T, class R>
  void evalJetWithoutC0(const unsigned degree, const T* left, const T exponent, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        evalC1HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 2:
        evalC1HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC1HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC3HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        evalC1HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC2HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        evalC3HomogenousPolynomial(left,exponent,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        for(const MultiindexData* i = dag->getIndexArray().begin(0,4);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,exponent,result,i,dag, coeffNo);
    }
  }

  template<class T, class R>
  void eval(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    evalC0(left,right,result,coeffNo);
    if(degree)
      evalJetWithoutC0(degree,left,*right,result,dag,coeffNo);
  }

  template<class T, class R>
  void evalHomogenousPolynomial(const unsigned degree, const T* left, const T* right, R result, DagIndexer<T>* dag, const unsigned coeffNo)
  {
    switch(degree)
    {
      case 1:
        evalC1HomogenousPolynomial(left,*right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 2:
        evalC2HomogenousPolynomial(left,*right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      case 3:
        evalC3HomogenousPolynomial(left,*right,result,dag->domainDimension(),dag->getOrder()+1,coeffNo);
        break;
      default:
        for(const MultiindexData* i = dag->getIndexArray().begin(0,degree);i<dag->getIndexArray().end(0,degree);++i)
          evalMultiindex(left,*right,result,i,dag,coeffNo);
    }
  }
}

// -------------------- PowFunTime ------------------------------------

namespace PowFunTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    Pow::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = exp((*right) * log(*left));
  }

  template<class T, class R>
  inline void eval(const unsigned /*degree*/, const T* left, const T* right, R result, DagIndexer<T>* /*dag*/, const unsigned coeffNo)
  {
    Pow::evalC0(left,right,result,coeffNo);
  }

  template<class T, class R>
  inline void evalHomogenousPolynomial(const unsigned /*degree*/, const T* /*left*/, const T* /*right*/, R /*result*/, DagIndexer<T>* /*dag*/, const unsigned /*coeffNo*/)
  {}
}

// -------------------- PowTime ------------------------------------

namespace PowTime
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo)
      result[coeffNo] = ((*right)-(coeffNo-1)) * result[coeffNo-1]/((double)coeffNo * (*left));
    else
      *result = exp((*right) * log((*left)));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = exp((*right) * log(*left));
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

// -------------------- PowConst ------------------------------------

namespace PowConst
{
  template<class T, class R>
  inline void evalC0(const T* left, const T* right, R result, const unsigned coeffNo)
  {
    if(coeffNo==0)
      *result = exp((*right) * log(*left));
  }

  template<class T, class R>
  inline void evalC0HomogenousPolynomial(const T* left, const T* right, R result)
  {
    *result = exp((*right) * log(*left));
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


CAPD_MAKE_DAG_NODE(Pow);
CAPD_MAKE_DAG_NODE(PowConst);
CAPD_MAKE_DAG_NODE(PowTime);
CAPD_MAKE_DAG_NODE(PowFunTime);
/// @}
}} // namespace capd::autodiff

#endif
