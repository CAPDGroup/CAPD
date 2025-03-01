/////////////////////////////////////////////////////////////////////////////
/// @file Curve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_DIFFALGEBRA_CURVE_H
#define CAPD_DIFFALGEBRA_CURVE_H
#include <stdexcept>
#include "capd/basicalg/TypeTraits.h"
#include "capd/diffAlgebra/ParametricCurve.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{
/**
 * This class provides methods for evaluation of the parametric curve
 * for a given parameter value.
 *
 * Template parameter BaseCurveT must provide storage and access for coefficients.
 */

template <class BaseCurveT, bool isInterval = capd::TypeTraits<typename BaseCurveT::ScalarType>::isInterval >
class Curve : public BaseCurveT, public ParametricCurve<typename BaseCurveT::MatrixType>{
public:
  typedef typename BaseCurveT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  Curve(Real left, Real right, size_type dimension, size_type order, size_type degree);
  VectorType timeDerivative(const ScalarType& h) const;

  VectorType operator()(const ScalarType& h) const;
  MatrixType derivative(const ScalarType& h) const;
  MatrixType operator[](const ScalarType& h) const { return this->derivative(h); }
};

/// Specialization for interval types
template <class BaseCurveT>
class Curve<BaseCurveT,true> : public BaseCurveT, public ParametricCurve<typename BaseCurveT::MatrixType>{
public:
  typedef typename BaseCurveT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  Curve(Real left, Real right, size_type dimension, size_type order, size_type degree);
  VectorType timeDerivative(const ScalarType& h) const;

  VectorType operator()(const ScalarType& h) const;
  MatrixType derivative(const ScalarType& h) const { return this->oneStepDerivative(h)*initMatrix; }
  MatrixType operator[](const ScalarType& h) const { return this->oneStepDerivative(h)*initMatrix; }

  void setInitMatrix(const MatrixType& M){
    initMatrix = M;
  }

  MatrixType oneStepDerivative(const ScalarType& h) const;
  MatrixType oneStepDerivativeOfNumericalMethod(const ScalarType& h) const;
  VectorType valueAtCenter(const ScalarType& h) const;
  VectorType remainder(const ScalarType& h) const;
  VectorType getCenter() const;
protected:
  MatrixType initMatrix;
};

///@}
}} // namespace capd::diffAlgebra

#endif // CAPD_DIFFALGEBRA_CURVE_H
