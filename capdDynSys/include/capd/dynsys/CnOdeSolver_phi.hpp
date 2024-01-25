

/////////////////////////////////////////////////////////////////////////////
/// @file CnOdeSolver_phi.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_CNODESOLVER_PHI_HPP_
#define _CAPD_DYNSYS_CNODESOLVER_PHI_HPP_

#include "capd/dynsys/CnOdeSolver.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
template<typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
typename CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::VectorType
CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::Phi(const ScalarType&t, const VectorType &iv)
{
  this->setCurrentTime(t);
  VectorType* coeff = this->getCoefficientsAtCenter();
  coeff[0] = iv;
  this->m_vField->computeODECoefficients(coeff,this->getOrder());

  // summation of the Taylor series
  VectorType v = coeff[this->getOrder()];
  for(int r = this->getOrder() - 1; r >= 0; --r)
    capd::vectalg::multiplyAssignObjectScalarAddObject(v,this->m_step,coeff[r]);
  return v;
}

//###########################################################//

template<typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
typename CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::MatrixType
CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::JacPhi(const ScalarType& t, const VectorType &iv)
{
  MatrixType result(dimension(),dimension());
  this->setCurrentTime(t);
  BaseTaylor::operator()(iv,result);
  return result;
}

/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_CNODESOLVER_PHI_HPP_

