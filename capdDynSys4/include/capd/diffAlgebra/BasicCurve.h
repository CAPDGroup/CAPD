/////////////////////////////////////////////////////////////////////////////
/// @file BasicCurve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_BASICCURVE_H_
#define _CAPD_DIFFALGEBRA_BASICCURVE_H_

#include <stdexcept>
#include "capd/diffAlgebra/CurveInterface.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

/**
 * This class is a data structure for storing of a parametric curve
 * together with first order derivatives with respect to initial point.
 *
 * More precisely, a curve c(t,x) is represented as
 * c(t,x_0) + d/dx c(t,x)(x-x_0) + smallRemainder(t,x)
 *
 * This class provides a member functions for accessing of all coefficients.
 *
 * This class is a basic class for more general Curve.
 */

template <class MatrixT >
class BasicCurve : public CurveInterface<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef __size_type size_type;
  typedef __difference_type difference_type;
  typedef typename VectorType::iterator iterator;
  typedef typename VectorType::const_iterator const_iterator;

  BasicCurve(size_type dimension, size_type order, size_type degree);
  BasicCurve(const BasicCurve& c);
  virtual ~BasicCurve();

  BasicCurve& operator=(const BasicCurve& c);

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

  const VectorType* getCoefficientsAtCenter() const;
  const VectorType* getCoefficients() const;
  const VectorType* getRemainderCoefficients() const;

  const MatrixType* getMatrixCoefficients() const;
  const MatrixType* getMatrixRemainderCoefficients() const;

  VectorType* getCoefficientsAtCenter();
  VectorType* getCoefficients();
  VectorType* getRemainderCoefficients();

  MatrixType* getMatrixCoefficients();
  MatrixType* getMatrixRemainderCoefficients();

  // iterators
  iterator beginCenterCoefficients(size_type j)           { return m_coefficientsAtCenter[j].begin(); }
  iterator beginCoefficients(size_type j)                 { return m_coefficients[j].begin(); }
  iterator beginRemainderCoefficients(size_type j)        { return m_remainderCoefficients[j].begin(); }
  iterator beginMatrixCoefficients(size_type j)           { return m_matrixCoefficients[j].begin(); }
  iterator beginMatrixRemainderCoefficients(size_type j)  { return m_matrixRemainderCoefficients[j].begin(); }
  iterator endCenterCoefficients(size_type j)             { return m_coefficientsAtCenter[j].end(); }
  iterator endCoefficients(size_type j)                   { return m_coefficients[j].end(); }
  iterator endRemainderCoefficients(size_type j)          { return m_remainderCoefficients[j].end(); }
  iterator endMatrixCoefficients(size_type j)             { return m_matrixCoefficients[j].end(); }
  iterator endMatrixRemainderCoefficients(size_type j)    { return m_matrixRemainderCoefficients[j].end(); }

  const_iterator beginCenterCoefficients(size_type j) const          { return m_coefficientsAtCenter[j].begin(); }
  const_iterator beginCoefficients(size_type j) const                { return m_coefficients[j].begin(); }
  const_iterator beginRemainderCoefficients(size_type j) const       { return m_remainderCoefficients[j].begin(); }
  const_iterator beginMatrixCoefficients(size_type j) const          { return m_matrixCoefficients[j].begin(); }
  const_iterator beginMatrixRemainderCoefficients(size_type j) const { return m_matrixRemainderCoefficients[j].begin(); }
  const_iterator endCenterCoefficients(size_type j) const            { return m_coefficientsAtCenter[j].end(); }
  const_iterator endCoefficients(size_type j) const                  { return m_coefficients[j].end(); }
  const_iterator endRemainderCoefficients(size_type j) const         { return m_remainderCoefficients[j].end(); }
  const_iterator endMatrixCoefficients(size_type j) const            { return m_matrixCoefficients[j].end(); }
  const_iterator endMatrixRemainderCoefficients(size_type j) const   { return m_matrixRemainderCoefficients[j].end(); }

  protected:
  void allocate();
  void deallocate();
  void copyData(const BasicCurve& c);

  VectorType *m_coefficientsAtCenter;
  VectorType *m_coefficients;
  VectorType *m_remainderCoefficients;

  MatrixType *m_matrixCoefficients;
  MatrixType *m_matrixRemainderCoefficients;

  size_type m_order;
  size_type m_allocatedOrder;
  size_type m_dimension;
};

// ----------------- inline definitions ------------------

template<class MatrixT>
inline
typename BasicCurve<MatrixT>::size_type BasicCurve<MatrixT>::getOrder() const{
  return this->m_order;
}


template<class MatrixT>
inline
typename BasicCurve<MatrixT>::size_type BasicCurve<MatrixT>::getAllocatedOrder() const{
  return this->m_allocatedOrder;
}


template<class MatrixT>
inline
typename BasicCurve<MatrixT>::size_type BasicCurve<MatrixT>::dimension() const{
  return this->m_dimension;
}

// -----------------------------------------------------------------------------

template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::VectorType* BasicCurve<MatrixT>::getCoefficientsAtCenter() const{
  return this->m_coefficientsAtCenter;
}


template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::VectorType* BasicCurve<MatrixT>::getCoefficients() const{
  return this->m_coefficients;
}


template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::VectorType* BasicCurve<MatrixT>::getRemainderCoefficients() const{
  return this->m_remainderCoefficients;
}


template<class MatrixT>
inline
const MatrixT* BasicCurve<MatrixT>::getMatrixCoefficients() const{
  return this->m_matrixCoefficients;
}

template<class MatrixT>
inline const MatrixT* BasicCurve<MatrixT>::getMatrixRemainderCoefficients() const{
  return this->m_matrixRemainderCoefficients;
}

// -----------------------------------------------------------------------------


template<class MatrixT>
inline
typename BasicCurve<MatrixT>::VectorType* BasicCurve<MatrixT>::getCoefficientsAtCenter(){
  return this->m_coefficientsAtCenter;
}


template<class MatrixT>
inline
typename BasicCurve<MatrixT>::VectorType* BasicCurve<MatrixT>::getCoefficients(){
  return this->m_coefficients;
}


template<class MatrixT>
inline
typename BasicCurve<MatrixT>::VectorType* BasicCurve<MatrixT>::getRemainderCoefficients(){
  return this->m_remainderCoefficients;
}


template<class MatrixT>
inline MatrixT* BasicCurve<MatrixT>::getMatrixCoefficients(){
  return this->m_matrixCoefficients;
}


template<class MatrixT>
inline MatrixT* BasicCurve<MatrixT>::getMatrixRemainderCoefficients(){
  return this->m_matrixRemainderCoefficients;
}

// -----------------------------------------------------------------------------

template<class MatrixT>
inline const typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::centerCoefficient(size_type i, size_type r) const {
  return this->m_coefficientsAtCenter[r][i];
}

template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::coefficient(size_type i, size_type r) const {
  return this->m_coefficients[r][i];
}

template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::remainderCoefficient(size_type i, size_type r) const {
  return this->m_remainderCoefficients[r][i];
}

template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::coefficient(size_type i, size_type j, size_type r) const {
  return this->m_matrixCoefficients[r](++i,++j);
}

template<class MatrixT>
inline
const typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type r) const {
  return this->m_matrixRemainderCoefficients[r](++i,++j);
}

// -----------------------------------------------------------------------------

template<class MatrixT>
inline typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::centerCoefficient(size_type i, size_type r){
  return this->m_coefficientsAtCenter[r][i];
}

template<class MatrixT>
inline
typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::coefficient(size_type i, size_type r) {
  return this->m_coefficients[r][i];
}

template<class MatrixT>
inline
typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::remainderCoefficient(size_type i, size_type r) {
  return this->m_remainderCoefficients[r][i];
}

template<class MatrixT>
inline
typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::coefficient(size_type i, size_type j, size_type r) {
  return this->m_matrixCoefficients[r](++i,++j);
}

template<class MatrixT>
inline
typename BasicCurve<MatrixT>::ScalarType& BasicCurve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type r) {
  return this->m_matrixRemainderCoefficients[r](++i,++j);
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICCURVE_H_
