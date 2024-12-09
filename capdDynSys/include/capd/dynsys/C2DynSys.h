

/////////////////////////////////////////////////////////////////////////////
/// @file C2DynSys.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2015 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C2DYNSYS_H_
#define _CAPD_DYNSYS_C2DYNSYS_H_

#include <string>
#include "capd/dynsys/C1DynSys.h"
#include "capd/diffAlgebra/Hessian.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<typename MatrixT>
class C2DynSys : public virtual capd::dynsys::C1DynSys<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef capd::diffAlgebra::Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;

  virtual void encloseC2Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc,
      HessianType& o_hessianPhi,
      HessianType& o_hessianRem,
      HessianType& o_hessianEnc
  ) = 0;
};
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_C2DYNSYS_H_


