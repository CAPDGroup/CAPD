/////////////////////////////////////////////////////////////////////////////
/// @file BasicC2Curve.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_BASICC2CURVE_HPP_
#define _CAPD_DIFFALGEBRA_BASICC2CURVE_HPP_

#include <stdexcept>
#include "capd/diffAlgebra/Hessian.hpp"
#include "capd/diffAlgebra/BasicC2Curve.h"
#include "capd/diffAlgebra/BasicCurve.hpp"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template<class MatrixT>
BasicC2Curve<MatrixT>::BasicC2Curve(size_type dimension, size_type order, size_type degree)
  : BaseCurve(dimension, order, degree)
{
  this->c2Allocate();
}

template<class MatrixT>
BasicC2Curve<MatrixT>::BasicC2Curve(const BasicC2Curve& c)
  : BaseCurve(c)
{
  this->c2Allocate();
  this->copyData(c);
}

template<class MatrixT>
BasicC2Curve<MatrixT>& BasicC2Curve<MatrixT>::operator=(const BasicC2Curve& c)
{
  if(&c == this)
    return *this;
  this->deallocate();
  this->c2Deallocate();
  this->allocate();
  this->c2Allocate();
  this->copyData(c);
  return *this;
}

template<class MatrixT>
void BasicC2Curve<MatrixT>::setOrder(size_type order) {
  if(order==this->m_order)
    return;
  if(order> m_allocatedOrder) {
    this->deallocate();
    this->c2Deallocate();
    this->m_order = order;
    this->allocate();
    this->c2Allocate();
  } else {
    this->m_order = order;
  }
}

template<class MatrixT>
BasicC2Curve<MatrixT>::~BasicC2Curve(){
  this->c2Deallocate();
}

template<class MatrixT>
void BasicC2Curve<MatrixT>::clearCoefficients(){

  BaseCurve::clearCoefficients();
  for(size_type i=0;i<=this->m_order+1;++i){
    this->m_hessianCoefficients[i].clear();
    this->m_hessianRemainderCoefficients[i].clear();
  }
}

template<class MatrixT>
void BasicC2Curve<MatrixT>::c2Allocate(){
  m_hessianCoefficients = HessianType::makeArray(this->m_order+2,this->m_dimension,this->m_dimension);
  m_hessianRemainderCoefficients = HessianType::makeArray(this->m_order+2,this->m_dimension,this->m_dimension);
}

template<class MatrixT>
void BasicC2Curve<MatrixT>::c2Deallocate(){

  delete [] m_hessianCoefficients;
  delete [] m_hessianRemainderCoefficients;
}

template<class MatrixT>
void BasicC2Curve<MatrixT>::copyData(const BasicC2Curve& c){

  BaseCurve::copyData(c);
  for(size_type i=0;i<=this->m_order+1;++i){
    this->m_hessianCoefficients[i] = c.m_hessianCoefficients[i];
    this->m_hessianRemainderCoefficients[i] = c.m_hessianRemainderCoefficients[i];
  }
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICC2CURVE_HPP_
