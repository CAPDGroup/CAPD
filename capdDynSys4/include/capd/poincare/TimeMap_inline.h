/////////////////////////////////////////////////////////////////////////////
/// @file TimeMap.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file is part of the CAPD library.  This library is free software;
// you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this software; see the file "license.txt".  If not, write to the
// Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA.

#ifndef _CAPD_POINCARE_TIME_MAP_INLINE_H_
#define _CAPD_POINCARE_TIME_MAP_INLINE_H_

#include "capd/poincare/TimeMap.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
template<typename SolverT>
inline const SolverT & TimeMap<SolverT>::getDynamicalSystem() const{
  return m_solver;
}

template<typename SolverT>
inline SolverT & TimeMap<SolverT>::getDynamicalSystem(){
   return m_solver;
}

template<typename SolverT>
inline const SolverT & TimeMap<SolverT>::getSolver() const{
  return m_solver;
}

template<typename SolverT>
inline SolverT & TimeMap<SolverT>::getSolver(){
   return m_solver;
}

template<typename SolverT>
inline const typename TimeMap<SolverT>::VectorFieldType& TimeMap<SolverT>::getVectorField() const
{
  return this->m_solver.getVectorField();
}

template<typename SolverT>
inline typename TimeMap<SolverT>::VectorFieldType& TimeMap<SolverT>::getVectorField()
{
  return this->m_solver.getVectorField();
}

template<typename SolverT>
inline typename TimeMap<SolverT>::size_type TimeMap<SolverT>::getOrder() const
{
  return this->m_solver.getOrder();
}

template<typename SolverT>
inline void TimeMap<SolverT>::setOrder(size_type newOrder)
{
  this->m_solver.setOrder(newOrder);
}

template<typename SolverT>
inline typename TimeMap<SolverT>::ScalarType TimeMap<SolverT>::getStep() const
{
  return this->m_solver.getStep();
}

template<typename SolverT>
inline void TimeMap<SolverT>::setStep(const ScalarType& newStep)
{
  this->m_solver.setStep(newStep);
}

template<typename SolverT>
inline void TimeMap<SolverT>::turnOnStepControl()
{
  this->getSolver().turnOnStepControl();
}

template<typename SolverT>
inline void TimeMap<SolverT>::turnOffStepControl()
{
  this->getSolver().turnOffStepControl();
}

template<typename SolverT>
inline void TimeMap<SolverT>::onOffStepControl(bool sc)
{
  this->getSolver().onOffStepControl(sc);
}

template<typename SolverT>
inline void TimeMap<SolverT>::stopAfterStep(bool b)
{
  this->m_oneStep = b;
}

template<typename SolverT>
inline bool TimeMap<SolverT>::completed() const
{
  return this->m_completed;
}

template<typename SolverT>
inline const typename SolverT::ScalarType& TimeMap<SolverT>::getCurrentTime() const
{
  return this->m_currentTime;
}

/// @}
}} // namespace capd::poincare

#endif
