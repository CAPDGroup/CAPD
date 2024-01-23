/////////////////////////////////////////////////////////////////////////////
/// @file ParametricCurve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_PARAMETRICCURVE_H_
#define _CAPD_DIFFALGEBRA_PARAMETRICCURVE_H_

#include <stdexcept>
#include <vector>
#include "capd/basicalg/TypeTraits.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/diffAlgebra/Jet.h"


namespace capd{
namespace diffAlgebra{

/**
 * This file defines an abstract class that represents parametric curve in \f$R^n\f$.
 * We do not specify how it is represented in the memory.
 *
 * The curve can be evaluated at given time 't' as well as differentiated with respect to 't'.
 *
 * The curve can store rigorous enclosures of curves.
 *
 * The main template parameter is an abstract class that represent a piece of curve as one polynomial.
 */
template<class MatrixT, class VectorT = typename MatrixT::RowVectorType>
class ParametricCurve{
public:
  typedef MatrixT MatrixType;
  typedef VectorT VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;

  ParametricCurve(Real left, Real right);

  virtual ~ParametricCurve(){}

  virtual VectorType operator()(const ScalarType& h) const = 0;
  virtual MatrixType derivative(const ScalarType& h) const = 0;
  virtual MatrixType operator[](const ScalarType& h) const = 0;
  virtual HessianType hessian(const ScalarType& /*h*/) const {
    throw std::logic_error("ParametricCurve::hessian is not implemented in abstract class. Use C2Curve or CnCurve instead.");
  }
  virtual JetType jet(const ScalarType& /*h*/) const{
    throw std::logic_error("ParametricCurve::jet is not implemented in abstract class. Use CnCurve instead.");
  }

  virtual void eval(ScalarType /*h*/, JetType& /*v*/) const{
    throw std::logic_error("ParametricCurve::eval is not implemented in abstract class. Use CnCurve instead.");
  }
  virtual void setDomain(Real left, Real right);
  virtual Real getLeftDomain() const;
  virtual Real getRightDomain() const;

protected:
  Real m_left, m_right; ///< domain
};

// ----------------- inline definitions ------------------

template<class MatrixT, class VectorT>
ParametricCurve<MatrixT,VectorT>::ParametricCurve(Real left, Real right)
  : m_left(left), m_right(right)
{}

template<class MatrixT, class VectorT>
inline
void ParametricCurve<MatrixT,VectorT>::setDomain(Real left, Real right)
{
  this->m_left = left;
  this->m_right = right;
}

template<class MatrixT, class VectorT>
inline
typename ParametricCurve<MatrixT,VectorT>::Real ParametricCurve<MatrixT,VectorT>::getLeftDomain() const
{
  return this->m_left;
}

template<class MatrixT, class VectorT>
inline
typename ParametricCurve<MatrixT,VectorT>::Real ParametricCurve<MatrixT,VectorT>::getRightDomain() const
{
  return this->m_right;
}

}} // namespace capd::diffAlgebra

#endif
