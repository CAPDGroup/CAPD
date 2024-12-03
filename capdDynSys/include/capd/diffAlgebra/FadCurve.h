/////////////////////////////////////////////////////////////////////////////
/// @file FadCurve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_FADCURVE_H_
#define _CAPD_DIFFALGEBRA_FADCURVE_H_

#include <stdexcept>
#include "capd/basicalg/TypeTraits.h"
#include "capd/diffAlgebra/CurveInterface.h"
#include "capd/fadbad/fadbad.h"
#include "capd/fadbad/fadiff.h"
#include "capd/fadbad/tadiff.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template < class MatrixT >
class FadCurve : public CurveInterface<MatrixT>{
public:
  // general typedefs required for the interface
  typedef MatrixT MatrixType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  typedef typename TypeTraits<ScalarType>::Real Real;

  // own typedefs
  typedef fadbad::T<ScalarType> TScalar;
  typedef fadbad::F<ScalarType,VectorType::csDim> FScalar;
  typedef fadbad::T<FScalar> TFScalar;
  typedef typename VectorType::template rebind<TFScalar>::other TFVector;
  typedef typename VectorType::template rebind<TScalar>::other TVector;

  FadCurve(size_type dimension, size_type order, size_type degree);

  virtual void setOrder(size_type order); ///< Sets the order of Taylor interpolation
  size_type getOrder() const;             ///< Returns the order of Taylor interpolation
  size_type getAllocatedOrder() const;    ///< Returns maximal allocated order - used to avoid memory reallocation
  size_type dimension() const;            ///< Returns the dimension in which the parametric curve is embedded

  void clearCoefficients();             ///< sets all coefficients to zero

  const ScalarType& centerCoefficient(size_type i, size_type j) const;
  const ScalarType& coefficient(size_type i, size_type j) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j) const;
  const ScalarType& coefficient(size_type i, size_type j, size_type k) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j, size_type k) const;

  ScalarType& centerCoefficient(size_type i, size_type j);
  ScalarType& coefficient(size_type i, size_type j);
  ScalarType& remainderCoefficient(size_type i, size_type j);
  ScalarType& coefficient(size_type i, size_type j, size_type k);
  ScalarType& remainderCoefficient(size_type i, size_type j, size_type k);

protected:

  size_type m_order;
  size_type m_dimension;

  mutable TVector m_center;
  mutable TFVector m_in;
  mutable TVector m_rem;
  mutable TFVector m_jacRem;

};

//###########################################################//

template<class MatrixT>
inline
typename FadCurve<MatrixT>::size_type FadCurve<MatrixT>::getOrder() const
{
  return m_order;
}

template<class MatrixT>
inline
typename FadCurve<MatrixT>::size_type FadCurve<MatrixT>::getAllocatedOrder() const
{
  return MaxLength;
}

template<class MatrixT>
inline
typename FadCurve<MatrixT>::size_type FadCurve<MatrixT>::dimension() const
{
  return m_dimension;
}

// -----------------------------------------------------------------------------

template<class MatrixT>
inline
const typename MatrixT::ScalarType& FadCurve<MatrixT>::centerCoefficient(size_type i, size_type r) const{
  return m_center[i][r];
}

template<class MatrixT>
inline
const typename MatrixT::ScalarType& FadCurve<MatrixT>::coefficient(size_type i, size_type r) const{
  return m_in[i][r].x();
}

template<class MatrixT>
inline
const typename MatrixT::ScalarType& FadCurve<MatrixT>::remainderCoefficient(size_type i, size_type r) const{
  return m_rem[i][r];
}

template<class MatrixT>
inline
const typename MatrixT::ScalarType& FadCurve<MatrixT>::coefficient(size_type i, size_type j, size_type r) const{
  return m_in[i][r].d(j);
}

template<class MatrixT>
inline
const typename MatrixT::ScalarType& FadCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type r) const{
  return m_jacRem[i][r].d(j);
}


// -----------------------------------------------------------------------------

template<class MatrixT>
inline
typename MatrixT::ScalarType& FadCurve<MatrixT>::centerCoefficient(size_type i, size_type r) {
  return m_center[i][r];
}

template<class MatrixT>
inline
typename MatrixT::ScalarType& FadCurve<MatrixT>::coefficient(size_type i, size_type r) {
  return m_in[i][r].x();
}

template<class MatrixT>
inline
typename MatrixT::ScalarType& FadCurve<MatrixT>::remainderCoefficient(size_type i, size_type r) {
  return m_rem[i][r];
}

template<class MatrixT>
inline
typename MatrixT::ScalarType& FadCurve<MatrixT>::coefficient(size_type i, size_type j, size_type r) {
  return m_in[i][r].d(j);
}

template<class MatrixT>
inline
typename MatrixT::ScalarType& FadCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type r) {
  return m_jacRem[i][r].d(j);
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_FADCURVE_H_
