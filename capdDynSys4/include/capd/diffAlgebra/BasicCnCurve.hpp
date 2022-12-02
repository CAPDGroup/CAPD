/////////////////////////////////////////////////////////////////////////////
/// @file BasicCnCurve.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_BASICCNCURVE_HPP_
#define _CAPD_DIFFALGEBRA_BASICCNCURVE_HPP_

#include <stdexcept>
#include "capd/diffAlgebra/BasicCnCurve.h"
#include "capd/diffAlgebra/Hessian.hpp"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template<class MatrixT>
BasicCnCurve<MatrixT>::BasicCnCurve(size_type dimension, size_type order, size_type degree)
  : m_order(order), m_allocatedOrder(0)
{
  this->allocate(dimension,degree);
}

template<class MatrixT>
BasicCnCurve<MatrixT>::BasicCnCurve(const BasicCnCurve& c)
  : m_order(c.m_order), m_allocatedOrder(0)
{
  this->allocate(c.dimension(),c.degree());
  this->copyData(c);
}

template<class MatrixT>
BasicCnCurve<MatrixT>& BasicCnCurve<MatrixT>::operator=(const BasicCnCurve& c)
{
  if(&c == this)
    return *this;
  this->deallocate();
  this->m_order = c.m_order;
  this->allocate(c.dimension(),c.degree());
  this->copyData(c);
  return *this;
}

template<class MatrixT>
void BasicCnCurve<MatrixT>::setDegree(size_type degree) {
  if(degree!=this->degree())
  {
    // save dim
    size_type dim = this->dimension();
    this->deallocate();
    this->allocate(dim,degree);
  }
}

template<class MatrixT>
void BasicCnCurve<MatrixT>::setOrder(size_type order) {
  if(order> m_allocatedOrder) {
    // save dim and degree
    size_type degree = this->degree();
    size_type dim = this->dimension();
    this->deallocate();
    this->m_order = order;
    this->allocate(dim,degree);
  } else {
    this->m_order = order;
  }
}

template<class MatrixT>
BasicCnCurve<MatrixT>::~BasicCnCurve(){
  this->deallocate();
}

template<class MatrixT>
void BasicCnCurve<MatrixT>::clearCoefficients()
{
  for(size_type i=0;i<=this->m_order+1;++i)
  {
    this->m_coefficientsAtCenter[i].clear();
    this->m_coefficients[i].clear();
    this->m_remainderCoefficients[i].clear();
  }
}

template<class MatrixT>
void BasicCnCurve<MatrixT>::allocate(size_type dimension, size_type degree){

  m_coefficientsAtCenter = VectorType::makeArray(this->m_order+2,dimension);
  m_coefficients = new JetType[this->m_order+2];
  m_remainderCoefficients = new JetType[this->m_order+2];

  this->m_allocatedOrder = m_order;
  for(size_type i=0;i<=this->m_order+1;++i)
  {
    m_coefficients[i] = JetType(dimension,dimension,degree);
    m_remainderCoefficients[i] = JetType(dimension,dimension,degree);
  }
}

template<class MatrixT>
void BasicCnCurve<MatrixT>::deallocate(){

  delete [] m_coefficientsAtCenter;
  delete [] m_coefficients;
  delete [] m_remainderCoefficients;
}

template<class MatrixT>
void BasicCnCurve<MatrixT>::copyData(const BasicCnCurve& c)
{
  for(size_type i=0;i<=this->m_order+1;++i)
  {
    this->m_coefficientsAtCenter[i] = c.m_coefficientsAtCenter[i];
    this->m_coefficients[i] = c.m_coefficients[i];
    this->m_remainderCoefficients[i] = c.m_remainderCoefficients[i];
  }
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICCNCURVE_HPP_
