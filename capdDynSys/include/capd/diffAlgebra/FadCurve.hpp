/////////////////////////////////////////////////////////////////////////////
/// @file FadCurve.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_FADCURVE_HPP_
#define _CAPD_DIFFALGEBRA_FADCURVE_HPP_

#include "capd/diffAlgebra/FadCurve.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template<class MatrixT>
FadCurve<MatrixT>::FadCurve(size_type dimension, size_type order, size_type /*degree*/)
  : m_order(order), m_dimension(dimension),
    m_center(dimension,true), m_in(dimension,true),
    m_rem(dimension,true),  m_jacRem(dimension,true)
{}

template<class MatrixT>
void FadCurve<MatrixT>::setOrder(size_type order)
{
  if(order<getAllocatedOrder()-1) // MaxLength is a define from tadiff.h. It is set to 40 by default
    m_order = order;
  else
    throw std::out_of_range("FadCurve<MatrixT>::setOrder - to large order. Change constant MaxLength in fadbad/tadiff.h if necessary.");
}


template<class MatrixT>
void FadCurve<MatrixT>::clearCoefficients()
{
  const size_type dim = this->dimension();
  const size_type order = this->getOrder();

  for(size_type i=0;i<dim;++i)
    for(size_type r=0;r<=order+1;++r)
    {
      this->centerCoefficient(i,r) = TypeTraits<ScalarType>::zero();
      this->remainderCoefficient(i,r) = TypeTraits<ScalarType>::zero();
      this->coefficient(i,r) = TypeTraits<ScalarType>::zero();
      for(size_type j=0;j<dim;++j)
      {
        this->coefficient(i,j,r) = TypeTraits<ScalarType>::zero();
        this->remainderCoefficient(i,j,r) = TypeTraits<ScalarType>::zero();
      }
    }
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_CURVE_HPP_

