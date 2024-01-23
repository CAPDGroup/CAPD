/////////////////////////////////////////////////////////////////////////////
/// @file CurveInterface.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_CURVEINTERFACE_H_
#define _CAPD_DIFFALGEBRA_CURVEINTERFACE_H_

#include <stdexcept>
#include <sstream>
#include "capd/basicalg/TypeTraits.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/diffAlgebra/Jet.h"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

/**
 * This class provides common interface for all types of curves.
 * Some of functions are provided for compilation purposes only and their call at runtime means internal error in the library.
 * The library user should not use directly this class in own code.
 */
template<class MatrixT>
class CurveInterface{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;

  void setInitMatrix(const MatrixType&){
    throw std::logic_error("CurveInterface::setInitMatrix is not implemented. Use Curve, C2Curve or CnCurve instead.");
  }
  void setInitHessian(const HessianType&){
    throw std::logic_error("CurveInterface::setInitHessian is not implemented. Use C2Curve or CnCurve instead.");
  }
  void setInitJet(const JetType&){
    throw std::logic_error("CurveInterface::setInitJet is not implemented. Use CnCurve instead.");
  }

  std::runtime_error domainErrorMessage(std::string msg, ScalarType h, Real left, Real right) const{
    std::ostringstream out;
    out << "capd::diffAlgebra::" << msg << " error:\nargument "
        << h << " is out of domain=[" << left << "," << right << "]\n";
    return std::runtime_error(out.str());
  }
};

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_BASICCURVE_H_
