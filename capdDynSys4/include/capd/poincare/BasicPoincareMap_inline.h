/////////////////////////////////////////////////////////////////////////////
/// @file BasicPoincareMap_inline.h
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_BASIC_POINCARE_MAP_INLINE_H_
#define _CAPD_POINCARE_BASIC_POINCARE_MAP_INLINE_H_

#include "capd/poincare/BasicPoincareMap.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline const typename BasicPoincareMap<SolverT, SectionT>::Solver&
BasicPoincareMap<SolverT, SectionT>::getSolver() const
{
  return m_solver;
}
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::Solver&
BasicPoincareMap<SolverT, SectionT>::getSolver()
{
  return m_solver;
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline const typename BasicPoincareMap<SolverT, SectionT>::Solver&
BasicPoincareMap<SolverT, SectionT>::getDynamicalSystem() const
{
  return m_solver;
}
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::Solver&
BasicPoincareMap<SolverT, SectionT>::getDynamicalSystem()
{
  return m_solver;
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::VectorFieldType const&
BasicPoincareMap<SolverT, SectionT>::getVectorField() const
{
  return m_solver.getVectorField();
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::VectorFieldType&
BasicPoincareMap<SolverT, SectionT>::getVectorField()
{
  return m_solver.getVectorField();
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline const typename BasicPoincareMap<SolverT, SectionT>::SectionType&
BasicPoincareMap<SolverT, SectionT>::getSection() const
{
  return *m_section;
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::size_type
BasicPoincareMap<SolverT, SectionT>::getOrder() const
{
  return m_solver.getOrder();
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::ScalarType
BasicPoincareMap<SolverT, SectionT>::getStep() const
{
  return m_solver.getStep();
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline void BasicPoincareMap<SolverT, SectionT>::setOrder(size_type newOrder)
{
  m_solver.setOrder(newOrder);
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline void BasicPoincareMap<SolverT, SectionT>::setStep(const ScalarType& newStep)
{
  m_solver.setStep(newStep);
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline void BasicPoincareMap<SolverT, SectionT>::setFactor(double newFactor)
{
  m_sectionFactor = newFactor;
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline void BasicPoincareMap<SolverT, SectionT>::setSection(const SectionType& newSection)
{
  m_section = &newSection;
}

// -------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::MatrixType
BasicPoincareMap<SolverT, SectionT>::computeDP(
      const VectorType& Px,
      const MatrixType& derivativeOfFlow,
      VectorType& dT,
      ScalarType returnTime
   )
{
  VectorType fieldOnPx = this->m_solver.getVectorField()(returnTime,Px);
  return this->m_section->computeDP(Px,derivativeOfFlow,fieldOnPx,dT);
}

// -------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline typename BasicPoincareMap<SolverT, SectionT>::MatrixType
BasicPoincareMap<SolverT, SectionT>::computeDP(
      const VectorType& Px,
      const MatrixType& derivativeOfFlow,
      ScalarType returnTime
   )
{
  VectorType dT(Px.dimension(),false);
  return this->computeDP(Px,derivativeOfFlow,dT,returnTime);
}

// -------------------------------------------------------------------

template <typename SolverT, typename SectionT>
inline void BasicPoincareMap<SolverT, SectionT>::computeDP(
      const VectorType& Px,
      const MatrixType& derivativeOfFlow,
      const HessianType& hessianOfFlow,
      MatrixType& DP,
      HessianType& D2P,
      ScalarType returnTime
    )
{
  const size_type dim = this->m_solver.getVectorField().dimension();
  VectorType dT(dim);
  MatrixType d2T(dim,dim);
  this->computeDP(Px,derivativeOfFlow,hessianOfFlow,DP,D2P,dT,d2T,returnTime);
}

// -----------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
inline void BasicPoincareMap<SolverT, FunT>::turnOnStepControl()
{
  m_solver.turnOnStepControl();
}
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
inline void BasicPoincareMap<SolverT, FunT>::turnOffStepControl()
{
  m_solver.turnOffStepControl();
}
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
inline void BasicPoincareMap<SolverT, FunT>::onOffStepControl(bool sc)
{
  m_solver.onOffStepControl(sc);
}
 
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
inline void BasicPoincareMap<SolverT, FunT>::setMaxReturnTime(double maxReturnTime)
{
  this->m_maxReturnTime = maxReturnTime;
}
 
// -----------------------------------------------------------------------------------------

template <typename SolverT, typename FunT>
inline void BasicPoincareMap<SolverT, FunT>::setBlowUpMaxNorm(double blowUpMaxNorm)
{
  this->m_blowUpMaxNorm = blowUpMaxNorm;
}

}}

#endif
