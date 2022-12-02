/////////////////////////////////////////////////////////////////////////////
/// @file BasicC2Curve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_BASICC2CURVE_H_
#define _CAPD_DIFFALGEBRA_BASICC2CURVE_H_

#include "capd/diffAlgebra/BasicCurve.h"
#include "capd/diffAlgebra/Hessian.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

/**
 * This class is a data structure for storing of a parametric curve
 * together with first and second order derivatives with respect to initial point.
 *
 * More precisely, a curve c(t,x) is represented as
 * c(t,x_0) + d/dx c(t,x)(x-x_0) + smallRemainder(t,x)
 *
 * This class provides a member functions for accessing of all coefficients.
 *
 * This class is a basic class for more general C2Curve.
 */

template <class MatrixT >
class BasicC2Curve : public BasicCurve<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef BasicCurve<MatrixT> BaseCurve;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef typename HessianType::iterator iterator;
  typedef typename HessianType::const_iterator const_iterator;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  BasicC2Curve(size_type dimension, size_type order, size_type degree);
  BasicC2Curve(const BasicC2Curve& c);
  ~BasicC2Curve();

  BasicC2Curve& operator=(const BasicC2Curve& c);

  virtual void setOrder(size_type order); ///< Sets the order of Taylor interpolation
  using BaseCurve::getOrder;
  using BaseCurve::getAllocatedOrder;
  using BaseCurve::dimension;

  void clearCoefficients();

  using BaseCurve::centerCoefficient;
  using BaseCurve::coefficient;
  using BaseCurve::remainderCoefficient;

  ScalarType& coefficient(size_type i, size_type j, size_type c, size_type k);
  ScalarType& remainderCoefficient(size_type i, size_type j, size_type c, size_type k);

  const ScalarType& coefficient(size_type i, size_type j, size_type c, size_type k) const;
  const ScalarType& remainderCoefficient(size_type i, size_type j, size_type c, size_type k) const;

  using BaseCurve::getCoefficientsAtCenter;
  using BaseCurve::getCoefficients;
  using BaseCurve::getRemainderCoefficients;

  using BaseCurve::getMatrixCoefficients;
  using BaseCurve::getMatrixRemainderCoefficients;

  HessianType* getHessianCoefficients();
  HessianType* getHessianRemainderCoefficients();

  const HessianType* getHessianCoefficients() const;
  const HessianType* getHessianRemainderCoefficients() const;

  HessianType& getHessianCoefficients(size_type p);
  HessianType& getHessianRemainderCoefficients(size_type p);

  const HessianType& getHessianCoefficients(size_type p) const;
  const HessianType& getHessianRemainderCoefficients(size_type p) const;

  // iterators
  iterator beginHessianCoefficients(size_type j)           { return m_hessianCoefficients[j].begin(); }
  iterator beginHessianRemainderCoefficients(size_type j)  { return m_hessianRemainderCoefficients[j].begin(); }
  iterator endHessianCoefficients(size_type j)             { return m_hessianCoefficients[j].end(); }
  iterator endHessianRemainderCoefficients(size_type j)    { return m_hessianRemainderCoefficients[j].end(); }
  const_iterator beginHessianCoefficients(size_type j) const           { return m_hessianCoefficients[j].begin(); }
  const_iterator beginHessianRemainderCoefficients(size_type j) const  { return m_hessianRemainderCoefficients[j].begin(); }
  const_iterator endHessianCoefficients(size_type j) const             { return m_hessianCoefficients[j].end(); }
  const_iterator endHessianRemainderCoefficients(size_type j) const    { return m_hessianRemainderCoefficients[j].end(); }

  protected:

  HessianType* m_hessianCoefficients;
  HessianType* m_hessianRemainderCoefficients;


  void c2Allocate();
  void c2Deallocate();
  void copyData(const BasicC2Curve& c);

  using BaseCurve::m_order;
  using BaseCurve::m_allocatedOrder;
  using BaseCurve::m_dimension;
};

// ----------------- inline definitions ------------------

template<class MatrixT>
inline typename BasicC2Curve<MatrixT>::HessianType*
BasicC2Curve<MatrixT>::getHessianCoefficients() {
  return this->m_hessianCoefficients;
}

template<class MatrixT>
inline const typename BasicC2Curve<MatrixT>::HessianType*
BasicC2Curve<MatrixT>::getHessianCoefficients() const{
  return this->m_hessianCoefficients;
}

template<class MatrixT>
inline typename BasicC2Curve<MatrixT>::HessianType*
BasicC2Curve<MatrixT>::getHessianRemainderCoefficients() {
  return this->m_hessianRemainderCoefficients;
}

template<class MatrixT>
inline const typename BasicC2Curve<MatrixT>::HessianType*
BasicC2Curve<MatrixT>::getHessianRemainderCoefficients() const{
  return this->m_hessianRemainderCoefficients;
}

// ----------------------------------------------------------

template<class MatrixT>
inline typename BasicC2Curve<MatrixT>::HessianType&
BasicC2Curve<MatrixT>::getHessianCoefficients(size_type p) {
  return this->m_hessianCoefficients[p];
}

template<class MatrixT>
inline const typename BasicC2Curve<MatrixT>::HessianType&
BasicC2Curve<MatrixT>::getHessianCoefficients(size_type p) const{
  return this->m_hessianCoefficients[p];
}

template<class MatrixT>
inline typename BasicC2Curve<MatrixT>::HessianType&
BasicC2Curve<MatrixT>::getHessianRemainderCoefficients(size_type p) {
  return this->m_hessianRemainderCoefficients[p];
}

template<class MatrixT>
inline const typename BasicC2Curve<MatrixT>::HessianType&
BasicC2Curve<MatrixT>::getHessianRemainderCoefficients(size_type p) const{
  return this->m_hessianRemainderCoefficients[p];
}

// ----------------------------------------------------------

template<class MatrixT>
inline
const typename BasicC2Curve<MatrixT>::ScalarType& BasicC2Curve<MatrixT>::coefficient(size_type i, size_type j, size_type c, size_type r) const {
  return this->m_hessianCoefficients[r](i,j,c);
}

template<class MatrixT>
inline
const typename BasicC2Curve<MatrixT>::ScalarType& BasicC2Curve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type c, size_type r) const {
  return this->m_hessianRemainderCoefficients[r](i,j,c);
}

template<class MatrixT>
inline
typename BasicC2Curve<MatrixT>::ScalarType& BasicC2Curve<MatrixT>::coefficient(size_type i, size_type j, size_type c, size_type r) {
  return this->m_hessianCoefficients[r](i,j,c);
}

template<class MatrixT>
inline
typename BasicC2Curve<MatrixT>::ScalarType& BasicC2Curve<MatrixT>::remainderCoefficient(size_type i, size_type j, size_type c, size_type r) {
  return this->m_hessianRemainderCoefficients[r](i,j,c);
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICCURVE_H_
