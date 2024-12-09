

/////////////////////////////////////////////////////////////////////////////
/// @file CnDynSys.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_CNDYNSYS_H_
#define _CAPD_DYNSYS_CNDYNSYS_H_

#include "capd/dynsys/C2DynSys.h"
#include "capd/diffAlgebra/Jet.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
template<typename MatrixT>
class CnDynSys : public virtual capd::dynsys::C2DynSys<MatrixT>{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;

  /**
   *
   * @param t - time in ODE
   * @param x - a centre of current set
   * @param xx - the set
   * @param[out] out_phi - should contain coefficients of solution jet at (t,xx)
   * @param[out] out_rem - should contain remainder coefficients computed at enclosure
   * @param[out] out_enc - should contain enclosure for the main part and variational equations
   * @return approximate solution (without remainder) at (t,x)
   */
  template<class JetT>
  VectorType encloseCnMap(
        const ScalarType& /*t*/,
        const VectorType& /*x*/,
        const VectorType& /*xx*/,
        JetT& /*out_phi*/,
        JetT& /*out_rem*/,
        JetT& /*out_enc*/
      ){
    throw std::logic_error("template method CnDynSys::encloseMap cannot be implemented! Use overridden version from an inherited class.\n");
  }

// from DynSys
  using C1DynSys<MatrixType>::Phi;
  using C1DynSys<MatrixType>::JacPhi;
  using C1DynSys<MatrixType>::Lipschitz;
  using C1DynSys<MatrixType>::enclosure;
};
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_CNDYNSYS_H_


