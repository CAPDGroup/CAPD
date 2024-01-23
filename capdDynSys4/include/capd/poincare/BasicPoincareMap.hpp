/////////////////////////////////////////////////////////////////////////////
/// @file BasicPoincareMap.hpp
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.


#ifndef _CAPD_POINCARE_BASIC_POINCARE_MAP_HPP_
#define _CAPD_POINCARE_BASIC_POINCARE_MAP_HPP_

#include <stdexcept>
#include "capd/poincare/BasicPoincareMap.h"
#include "capd/poincare/BasicPoincareMap_operator.hpp"
#include "capd/basicalg/factrial.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
// ------------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
BasicPoincareMap<SolverT, SectionT>::BasicPoincareMap(
   Solver& solver,
   SectionType & section,
   CrossingDirection direction,
   const RealType & errorTolerance
 ) : m_solver(solver),
     m_section(&section),
     m_crossingDirection(direction),
     m_sectionFactor(errorTolerance),
     m_blowUpMaxNorm(1e+10),
     m_maxReturnTime(1000.)
{}

// ------------------------------------------------------------------------------------------

/// Crosses the section and returns the value of Poincare Map.
template <typename SolverT, typename SectionT>
typename BasicPoincareMap<SolverT, SectionT>::ScalarType
BasicPoincareMap<SolverT, SectionT>::findRelativeCrossingTime(ScalarType timeJustBeforeSection, const CurveType & c)
{
  ScalarType t0 = (c.getLeftDomain()+c.getRightDomain())/2;
  int maxNumberOfNewtonIterations= 15;

  // now we resolve for the return time
  for(int i=0; i<maxNumberOfNewtonIterations; ++i) {
    VectorType v = c(t0);
    ScalarType s_t = (*m_section)(v);
    if(abs(s_t) < this->m_sectionFactor)
      break;

    // note: c'(t0) = f(t0,c(t0))
    ScalarType Ds_t = this->m_section->gradient(v) * this->getVectorField()(timeJustBeforeSection+t0,v);
    ScalarType N = t0 - s_t/Ds_t;

    // we know that zero is in the domain
    if(N < c.getLeftDomain())
      N = c.getLeftDomain();
    if(N > c.getRightDomain())
      N = c.getRightDomain();
    t0 = N;
  }
  return t0;
}

// -------------------------------------------------------------------

template <typename SolverT, typename SectionT>
void BasicPoincareMap<SolverT, SectionT>::computeDP(
      const VectorType& Px,
      const MatrixType& derivativeOfFlow,
      const HessianType& hessianOfFlow,
      MatrixType& DP,
      HessianType& D2P,
      VectorType& dT,
      MatrixType& d2T,
      ScalarType returnTime
    )
{
  size_type i;
  const size_type dim = this->m_solver.getVectorField().dimension();
  this->m_solver.computeCoefficientsAtCenter(returnTime,Px,2);
  VectorType fieldOnPx(dim);
  VectorType d2Phidt2(dim);
  for(i=0;i<dim;++i){
    fieldOnPx[i] = this->m_solver.centerCoefficient(i,1);
    d2Phidt2[i] = this->m_solver.centerCoefficient(i,2);
  }
  MatrixType derOfVectorFieldOnPx = this->m_solver.getVectorField().derivative(returnTime,Px);
  this->m_section->computeDP(Px,derivativeOfFlow,hessianOfFlow,fieldOnPx,d2Phidt2,derOfVectorFieldOnPx,DP,D2P,dT,d2T);
}

/// @}
}} // namespace capd::poincare

#endif // #define _CAPD_POINCARE_BASIC_POINCARE_MAP_HPP_

