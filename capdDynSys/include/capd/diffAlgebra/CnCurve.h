/////////////////////////////////////////////////////////////////////////////
/// @file CnCurve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_CNCURVE_H_
#define _CAPD_DIFFALGEBRA_CNCURVE_H_

#include <stdexcept>
#include "capd/diffAlgebra/Curve.h"
#include "capd/diffAlgebra/Jet.h"

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
class CnCurve : public Curve<BaseCurveT, isInterval>{
public:
  typedef typename BaseCurveT::HessianType HessianType;
  typedef typename BaseCurveT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef capd::diffAlgebra::Jet<MatrixType,0> JetType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  CnCurve(Real left, Real right, size_type dimension, size_type order, size_type degree)
    : Curve<BaseCurveT>(left,right,dimension,order,degree)
  {}

  HessianType hessian(const ScalarType& h) const;
  JetType jet(const ScalarType& h) const;
  void eval(ScalarType h, JetType& v) const;
};

template <class BaseCurveT>
class CnCurve<BaseCurveT,true> : public Curve<BaseCurveT, true>{
public:
  typedef typename BaseCurveT::HessianType HessianType;
  typedef typename BaseCurveT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef capd::diffAlgebra::Jet<MatrixType,0> JetType;
  typedef __size_type size_type;
  typedef __difference_type difference_type;

  CnCurve(Real left, Real right, size_type dimension, size_type order, size_type degree)
    : Curve<BaseCurveT>(left,right,dimension,order,degree),
      initHessian(dimension),
      initJet(dimension,degree)
  {
    initJet.setMatrix(this->initMatrix);
  }

  HessianType hessian(const ScalarType& h) const;

  void setInitHessian(const HessianType& H){
    initHessian = H;
  }

  void setInitJet(const JetType& jet){
    initJet = jet;
  }

  void eval(ScalarType h, JetType& v) const;
  JetType jet(const ScalarType& h) const;
  
  HessianType initHessian;
  JetType initJet;
};

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_CNCURVE_H_
