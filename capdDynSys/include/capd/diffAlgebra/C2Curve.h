/////////////////////////////////////////////////////////////////////////////
/// @file C2Curve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_C2CURVE_H_
#define _CAPD_DIFFALGEBRA_C2CURVE_H_

#include <stdexcept>
#include "capd/diffAlgebra/Curve.h"

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

template <class BaseCurveT , bool isInterval = capd::TypeTraits<typename BaseCurveT::ScalarType>::isInterval >
class C2Curve : public Curve<BaseCurveT, isInterval>{
public:
  typedef typename BaseCurveT::HessianType HessianType;
  typedef typename BaseCurveT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  C2Curve(Real left, Real right, size_type dimension, size_type order, size_type degree)
    : Curve<BaseCurveT>(left,right,dimension,order,degree)
  {}
  HessianType hessian(const ScalarType& h) const;

};

template <class BaseCurveT>
class C2Curve<BaseCurveT,true> : public Curve<BaseCurveT, true>{
public:
  typedef typename BaseCurveT::HessianType HessianType;
  typedef typename BaseCurveT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  C2Curve(Real left, Real right, size_type dimension, size_type order, size_type degree)
    : Curve<BaseCurveT>(left,right,dimension,order,degree),
      initHessian(dimension)
  {}

  HessianType hessian(const ScalarType& h) const{
    return this->oneStepHessian(h)*this->initMatrix + this->oneStepDerivative(h)*this->initHessian;
  }

  void setInitHessian(const HessianType& H){
    initHessian = H;
  }
protected:
  HessianType oneStepHessian(const ScalarType& h) const;
  HessianType initHessian;
};

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_C2CURVE_H_
