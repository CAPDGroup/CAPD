/////////////////////////////////////////////////////////////////////////////
/// @file BasicCnCurve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_BASICCNCURVE_H_
#define _CAPD_DIFFALGEBRA_BASICCNCURVE_H_

#include <stdexcept>
#include "capd/basicalg/TypeTraits.h"
#include "capd/diffAlgebra/CurveInterface.h"
#include "capd/diffAlgebra/Jet.h"
#include "capd/diffAlgebra/Hessian.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

/**
 * This class is a data structure for storing of a parametric curve
 * together with its partial derivatives with respect to initial point
 * up to desired order.
 *
 * More precisely, a curve c(t,x) is represented as
 * c(t,x_0) + d/dx c(t,x)(x-x_0) + smallRemainder(t,x).
 *
 * Derivatives are represented as d^i c(t,x) + smallRemainder(t,x).
 *
 * This class provides a member functions for accessing of all coefficients.
 */

template <class MatrixT>
class BasicCnCurve : public CurveInterface<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Jet<MatrixT,0> JetType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  BasicCnCurve(size_type dimension, size_type order, size_type degree);
  BasicCnCurve(const BasicCnCurve& c);
  virtual ~BasicCnCurve();

  BasicCnCurve& operator=(const BasicCnCurve& c);

  virtual void setOrder(size_type order); ///< Sets the order of Taylor interpolation
  size_type getOrder() const;             ///< Returns the order of Taylor interpolation
  size_type getAllocatedOrder() const;    ///< Returns maximal allocated order - used to avoid memory reallocation
  size_type dimension() const;            ///< Returns the dimension in which the parametric curve is embedded
  size_type degree() const;               ///< Returns maximal degree of the jet
  void setDegree(size_type degree);  ///< Sets new maximal degree of the jet and reallocates memory

  void clearCoefficients();                 ///< sets all coefficients to zero

  const ScalarType& centerCoefficient(size_type i, size_type k) const;
  const ScalarType& coefficient(size_type i, size_type k) const;
  const ScalarType& coefficient(size_type i, size_type j, size_type k) const;
  const ScalarType& coefficient(size_type i, size_type j, size_type c, size_type k) const;
  const ScalarType& coefficient(size_type i, size_type j, size_type c, size_type s, size_type k) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j, size_type k) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j, size_type c, size_type k) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j, size_type c, size_type s, size_type k) const;

  ScalarType& centerCoefficient(size_type i, size_type k);
  ScalarType& coefficient(size_type i, size_type k);
  ScalarType& coefficient(size_type i, size_type j, size_type k);
  ScalarType& coefficient(size_type i, size_type j, size_type c, size_type k);
  ScalarType& coefficient(size_type i, size_type j, size_type c, size_type s, size_type k);
  ScalarType& remainderCoefficient(size_type i, size_type j);
  ScalarType& remainderCoefficient(size_type i, size_type j, size_type k);
  ScalarType& remainderCoefficient(size_type i, size_type j, size_type c, size_type k);
  ScalarType& remainderCoefficient(size_type i, size_type j, size_type c, size_type s, size_type k);


//protected:
  void allocate(size_type dimension, size_type degree);
  void deallocate();
  void copyData(const BasicCnCurve& c);

  const VectorType* getCoefficientsAtCenter() const;
  const JetType* getCoefficients() const;
  const JetType* getRemainderCoefficients() const;

  VectorType* getCoefficientsAtCenter();
  JetType* getCoefficients();
  JetType* getRemainderCoefficients();

  VectorType *m_coefficientsAtCenter;
  JetType *m_coefficients;
  JetType *m_remainderCoefficients;

  size_type m_order;
  size_type m_allocatedOrder;
};

// ----------------- inline definitions ------------------

template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::size_type BasicCnCurve<MatrixT>::getOrder() const{
  return this->m_order;
}


template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::size_type BasicCnCurve<MatrixT>::getAllocatedOrder() const{
  return this->m_allocatedOrder;
}


template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::size_type BasicCnCurve<MatrixT>::dimension() const{
  return this->m_coefficientsAtCenter[0].dimension();
}

template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::size_type BasicCnCurve<MatrixT>::degree() const{
  return this->m_coefficients[0].degree();
}

template <class MatrixT>
inline const typename BasicCnCurve<MatrixT>::VectorType* BasicCnCurve<MatrixT>::getCoefficientsAtCenter() const {
  return this->m_coefficientsAtCenter;
}

template <class MatrixT>
inline typename BasicCnCurve<MatrixT>::VectorType* BasicCnCurve<MatrixT>::getCoefficientsAtCenter() {
  return this->m_coefficientsAtCenter;
}

template <class MatrixT>
inline const typename BasicCnCurve<MatrixT>::JetType* BasicCnCurve<MatrixT>::getCoefficients() const {
  return this->m_coefficients;
}

template <class MatrixT>
inline typename BasicCnCurve<MatrixT>::JetType* BasicCnCurve<MatrixT>::getCoefficients(){
  return this->m_coefficients;
}

template <class MatrixT>
inline const typename BasicCnCurve<MatrixT>::JetType* BasicCnCurve<MatrixT>::getRemainderCoefficients() const {
  return this->m_remainderCoefficients;
}

template <class MatrixT>
inline typename BasicCnCurve<MatrixT>::JetType* BasicCnCurve<MatrixT>::getRemainderCoefficients(){
  return this->m_remainderCoefficients;
}

// -----------------------------------------------------------------------------

template <class MatrixT>
inline const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::centerCoefficient(size_type i, size_type r) const {
  return this->m_coefficientsAtCenter[r][i];
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type r) const {
  return this->m_coefficients[r](i);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type r) const {
  return this->m_remainderCoefficients[r](i);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type j, size_type r) const {
  return this->m_coefficients[r](i,j);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type r) const {
  return this->m_remainderCoefficients[r](i,j);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type j, size_type c, size_type r) const {
  return this->m_coefficients[r](i,j,c);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type c, size_type r) const {
  return this->m_remainderCoefficients[r](i,j,c);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type j, size_type c, size_type s, size_type r) const {
  return this->m_coefficients[r](i,j,c,s);
}

template <class MatrixT>
inline
const typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type c, size_type s, size_type r) const {
  return this->m_remainderCoefficients[r](i,j,c,s);
}

// -----------------------------------------------------------------------------

template <class MatrixT>
inline  typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::centerCoefficient(size_type i, size_type r)  {
  return this->m_coefficientsAtCenter[r][i];
}

template <class MatrixT>
inline
 typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type r)  {
  return this->m_coefficients[r](i);
}

template <class MatrixT>
inline
 typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type r)  {
  return this->m_remainderCoefficients[r](i);
}

template <class MatrixT>
inline
 typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type j, size_type r)  {
  return this->m_coefficients[r](i,j);
}

template <class MatrixT>
inline
 typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type r)  {
  return this->m_remainderCoefficients[r](i,j);
}

template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type j, size_type c, size_type r)  {
  return this->m_coefficients[r](i,j,c);
}

template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type c, size_type r) {
  return this->m_remainderCoefficients[r](i,j,c);
}

template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::coefficient(size_type i, size_type j, size_type c, size_type s, size_type r)  {
  return this->m_coefficients[r](i,j,c,s);
}

template <class MatrixT>
inline
typename BasicCnCurve<MatrixT>::ScalarType& BasicCnCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type c, size_type s, size_type r) {
  return this->m_remainderCoefficients[r](i,j,c,s);
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICCNCURVE_H_
