/////////////////////////////////////////////////////////////////////////////
/// @file BasicPoincareMap_operator.hpp
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_BASIC_POINCARE_MAP_OPERATOR_HPP_
#define _CAPD_POINCARE_BASIC_POINCARE_MAP_OPERATOR_HPP_

#include "capd/poincare/BasicPoincareMap.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/diffAlgebra/CnTimeJet.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(const VectorType& v)
{
  ScalarType returnTime = 0.;
  return this->operator()(v,returnTime);
}

// ------------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(VectorType v, ScalarType& in_out_time)
{
  SaveStepControl<SolverT> ssc(this->m_solver);

  // constructor sets currentTime=0 by default
  capd::diffAlgebra::C0TimeJet<VectorType> coeff(&v);
  coeff.setCurrentTime(in_out_time);
  this->integrateUntilSectionCrossing(coeff, 1);
  ScalarType timeJustBeforeSection = coeff.getCurrentTime();
  const CurveType& c = this->m_solver.getCurve();
  ScalarType t0 = findRelativeCrossingTime(timeJustBeforeSection,c);
  in_out_time = timeJustBeforeSection + t0;

  return c(t0);
}

// ------------------------------------------------------------------------------------------


template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(VectorType v, VectorType& pointAfterSection, ScalarType& in_out_time)
{
  SaveStepControl<SolverT> ssc(this->m_solver);

  // constructor sets currentTime=0 by default
  capd::diffAlgebra::C0TimeJet<VectorType> coeff(&v);
  coeff.setCurrentTime(in_out_time);

  pointAfterSection = this->integrateUntilSectionCrossing(coeff, 1);
  ScalarType timeJustBeforeSection = coeff.getCurrentTime();
  const CurveType& c = this->m_solver.getCurve();
  ScalarType t0 = findRelativeCrossingTime(timeJustBeforeSection,c);
  in_out_time = timeJustBeforeSection + t0;

  return c(t0);
}

// ------------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(const VectorType& v, MatrixType& dF)
{
  ScalarType returnTime = 0.;
  return this->operator()(v,dF,returnTime);
}

// ------------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(VectorType v, MatrixType& dF, ScalarType& in_out_time)
{
  SaveStepControl<SolverT> ssc(this->m_solver);
  dF.setToIdentity();
  capd::diffAlgebra::C1TimeJet<MatrixType> coeffs(&v,&dF);
  coeffs.setCurrentTime(in_out_time);

  integrateUntilSectionCrossing(coeffs, 1);
  ScalarType timeJustBeforeSection = coeffs.getCurrentTime();
  const CurveType& c = this->m_solver.getCurve();

  ScalarType t0 = findRelativeCrossingTime(timeJustBeforeSection,c);
  in_out_time = timeJustBeforeSection + t0;

  dF = c[t0];
  return c(t0);
}

// ------------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(const VectorType& v, MatrixType& dF, HessianType& h)
{
  ScalarType returnTime = 0.;
  return this->operator()(v,dF,h,returnTime);
}

// ------------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
typename BasicPoincareMap<SolverT, FunT>::VectorType
BasicPoincareMap<SolverT, FunT>::operator()(VectorType v, MatrixType& dF, HessianType& h, ScalarType& in_out_time)
{
  SaveStepControl<SolverT> ssc(this->m_solver);
  dF.setToIdentity();
  h.clear();
  capd::diffAlgebra::C2TimeJet<MatrixType> coeffs(&v,&dF,&h);
  coeffs.setCurrentTime(in_out_time);

  integrateUntilSectionCrossing(coeffs, 1);
  ScalarType timeJustBeforeSection = coeffs.getCurrentTime();

  const CurveType& c = this->m_solver.getCurve();
  ScalarType t0 = findRelativeCrossingTime(timeJustBeforeSection,c);
  dF = c[t0];
  h = c.hessian(t0);
  in_out_time = timeJustBeforeSection + t0;
  return c(t0);
}

// -----------------------------------------------------------------------------------------

/**
 * Nonrigorous Poincare map
 *
 * A point just after the section on the nonrigorous trajectory is returned
 * The result contains also approximate values of the partial derivatives
 * of the flow
 */

template <typename SolverT, typename SectionT>
typename BasicPoincareMap<SolverT, SectionT>::VectorType BasicPoincareMap<SolverT, SectionT>::operator()(JetType& x)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(x,returnTime);
}

// -----------------------------------------------------------------------------------------

/**
 * Nonrigorous Poincare map
 *
 * A point just after the section on the nonrigorous trajectory is returned
 * The result contains also approximate values of the partial derivatives
 * of the flow
 */

template <typename SolverT, typename SectionT>
typename BasicPoincareMap<SolverT, SectionT>::VectorType BasicPoincareMap<SolverT, SectionT>::operator()(JetType& x, ScalarType& in_out_time)
{
  SaveStepControl<SolverT> ssc(this->m_solver);
  capd::diffAlgebra::CnTimeJet<MatrixType,0> coeffs(&x);
  coeffs.setCurrentTime(in_out_time);
  this->integrateUntilSectionCrossing(coeffs, 1);
  ScalarType t0 = findRelativeCrossingTime(coeffs.getCurrentTime(),this->m_solver.getCurve());
  this->m_solver.getCurve().eval(t0,x);
  return x();
}

/// @}
}} // namespace capd::poincare

#endif // #define _CAPD_POINCARE_BASIC_POINCARE_MAP_OPERATOR_HPP_

