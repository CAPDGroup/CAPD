/////////////////////////////////////////////////////////////////////////////
/// @file BasicCurve.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_BASICCURVE_HPP_
#define _CAPD_DIFFALGEBRA_BASICCURVE_HPP_

#include <stdexcept>
#include "capd/diffAlgebra/BasicCurve.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template<class MatrixT>
BasicCurve<MatrixT>::BasicCurve(size_type dimension, size_type order, size_type /*degree*/)
  : m_order(order), m_allocatedOrder(0), m_dimension(dimension)
{
  this->allocate();
}

template<class MatrixT>
BasicCurve<MatrixT>::BasicCurve(const BasicCurve& c)
  :m_order(c.m_order), m_allocatedOrder(0), m_dimension(c.m_dimension)
{
  this->allocate();
  this->copyData(c);
}

template<class MatrixT>
BasicCurve<MatrixT>& BasicCurve<MatrixT>::operator=(const BasicCurve& c)
{
  if(&c == this)
    return *this;
  this->deallocate();
  this->m_dimension = c.m_dimension;
  this->m_order = c.m_order;
  this->allocate();
  this->copyData(c);
  return *this;
}

template<class MatrixT>
void BasicCurve<MatrixT>::setOrder(size_type order) {
  if(order==this->m_order)
    return;
  if(order> m_allocatedOrder) {
    this->deallocate();
    this->m_order = order;
    this->allocate();
  } else {
    this->m_order = order;
  }
}

template<class MatrixT>
BasicCurve<MatrixT>::~BasicCurve(){
  this->deallocate();
}

template<class MatrixT>
void BasicCurve<MatrixT>::clearCoefficients(){

  for(size_type i=0;i<=this->m_order+1;++i){
    this->m_coefficientsAtCenter[i].clear();
    this->m_coefficients[i].clear();
    this->m_remainderCoefficients[i].clear();
    this->m_matrixCoefficients[i].clear();
    this->m_matrixRemainderCoefficients[i].clear();
  }
}

template<class MatrixT>
void BasicCurve<MatrixT>::allocate(){

  m_coefficientsAtCenter = VectorType::makeArray(this->m_order+2,this->m_dimension);
  m_coefficients = VectorType::makeArray(this->m_order+2,this->m_dimension);
  m_remainderCoefficients = VectorType::makeArray(this->m_order+2,this->m_dimension);
  m_matrixCoefficients = MatrixType::makeArray(this->m_order+2,this->m_dimension,this->m_dimension);
  m_matrixRemainderCoefficients = MatrixType::makeArray(this->m_order+2,this->m_dimension,this->m_dimension);
  this->m_allocatedOrder = m_order;
}

template<class MatrixT>
void BasicCurve<MatrixT>::deallocate(){

  delete [] m_coefficientsAtCenter;
  delete [] m_coefficients;
  delete [] m_remainderCoefficients;
  delete [] m_matrixCoefficients;
  delete [] m_matrixRemainderCoefficients;
}

template<class MatrixT>
void BasicCurve<MatrixT>::copyData(const BasicCurve& c){

  for(size_type i=0;i<=this->m_order+1;++i){
    this->m_coefficientsAtCenter[i] = c.m_coefficientsAtCenter[i];
    this->m_coefficients[i] = c.m_coefficients[i];
    this->m_remainderCoefficients[i] = c.m_remainderCoefficients[i];
    this->m_matrixCoefficients[i] = c.m_matrixCoefficients[i];
    this->m_matrixRemainderCoefficients[i] = c.m_matrixRemainderCoefficients[i];
  }
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICCURVE_HPP_
