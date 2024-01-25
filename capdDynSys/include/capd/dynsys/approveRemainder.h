

/////////////////////////////////////////////////////////////////////////////
/// @file approveRemainder.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2015 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#include "capd/poincare/SaveStepControl.h"

#ifndef _CAPD_DYNSYS_APPROVEREMAINDER_H_
#define _CAPD_DYNSYS_APPROVEREMAINDER_H_

namespace capd{
namespace dynsys {
/// @addtogroup dynsys
/// @{

template<class Solver, class RemainderType>
void computeAndApproveRemainder(
      Solver& solver,
      const typename Solver::ScalarType& t,   //< @param[in] current time of ODE
      const typename Solver::VectorType& xx,  //< @param[in] set to be moved along the trajectories of ODE
      RemainderType& o_rem,                   //< @param[out] bound for the error of numerical method over the time step
      RemainderType& o_enc                    //< @param[out] bound for the enclosure over the time step
  )
{
  typedef typename Solver::ScalarType ScalarType;
  typedef typename Solver::VectorType VectorType;
  typedef typename ScalarType::BoundType Real;
  
  capd::poincare::SaveStepControl<Solver> saveStepControl(solver);

  // IMPORTANT: this MUST be computed BEFORE call to computeReamainder!
  Real eps = Real(2.)*Solver::getEffectiveTolerance(solver,xx);

  // DW: do not decrease the factor below to 2-3. 
  // My simulation shows that the restriction in many cases worse the estimations.
  // The probability of significant reduction of predicted step is rather rare.
  if(!isSingular(solver.getStep()))
    solver.setMaxStep(capd::min(solver.getMaxStep(),Real(2)*solver.getStep()));
  solver.computeTimeStep(t,xx);
  solver.computeRemainder(t,xx,o_enc,o_rem);

  // if changing step is allowed, validate remainder and adjust time step if necessary
  if(solver.isStepChangeAllowed() and !isSingular(solver.getStep())) {
    Real remSize = rightBound(maxWidth((VectorType)o_rem));
    while( remSize > eps and solver.getMaxStep() > solver.getStepControl().getMinStepAllowed() )  {
      solver.adjustTimeStep((solver.getStep()*0.95).leftBound());
      solver.computeRemainder(t,xx,o_enc,o_rem);
      remSize = rightBound(maxWidth((VectorType)o_rem));
    }
  }
}

/// @}
}} // namespace capd::dynsys

#endif

