/////////////////////////////////////////////////////////////////////////////
/// @file TimeMap_template.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_TIME_MAP_TEMPLATE_H_
#define _CAPD_POINCARE_TIME_MAP_TEMPLATE_H_

#include "capd/poincare/TimeMap.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
template <typename SolverT>
template<typename SetType>
void TimeMap<SolverT>::initComputations(SetType &originalSet, void* current)
{
  m_currentTime = originalSet.getCurrentTime();
  m_currentSet = current;
  m_completed = false;
  this->m_solver.getStepControl().init(this->m_solver,this->m_currentTime,originalSet);
}

// ----------------------------------------------------------------------------------------------

template <typename SolverT>
template<typename SetType>
void TimeMap<SolverT>::moveSet(ScalarType Time, SetType &originalSet, void* current)
{
  SaveStepControl<SolverT> ssc(m_solver);

  ScalarType maxStep;
  if(this->m_currentSet != current)
    initComputations(originalSet,current);

  if(!(this->m_currentTime < Time)){
    this->m_currentSet = 0;
    this->m_solver.adjustTimeStep(TypeTraits<ScalarType>::zero());
    this->m_completed = true;
    return;
  }

  const bool loop = !this->m_oneStep;
  do
  {
    maxStep = Time-this->m_currentTime;

    while(true)
    {
      try{
        this->m_solver.setMaxStep(maxStep);
        this->m_solver(originalSet);
        if(m_solution)
          m_solution->add(this->m_solver.getCurve());
        break;
      } catch(std::exception& e)
      {
        maxStep /= typename TypeTraits<ScalarType>::Real(1.5);
        this->m_solver.clearCoefficients();
        if(!this->m_solver.isStepChangeAllowed() or maxStep<this->m_solver.getStepControl().getMinStepAllowed()){
          this->m_currentSet = 0;
          this->m_completed = true;
          throw;
        }
      }
    }

    this->m_currentTime += this->m_solver.getStep();
    if(!(this->m_currentTime < Time)){
      this->m_currentSet = 0;
      this->m_completed = true;
      return;
    }
  }while(loop);
}

// ----------------------------------------------------------------------------------------------

template <typename SolverT>
template<typename SetType>
typename SolverT::VectorType TimeMap<SolverT>::operator()(ScalarType Time, SetType &originalSet)
{
  m_solution = 0;
  this->moveSet(Time,originalSet,&originalSet);
  return VectorType(originalSet);
}

// ----------------------------------------------------------------------------------------------

template <typename SolverT>
template<typename SetType>
typename SolverT::VectorType TimeMap<SolverT>::operator()(ScalarType Time, SetType &originalSet, SolutionCurve& solution)
{
  m_solution = &solution;
  this->moveSet(Time,originalSet,&originalSet);
  return VectorType(originalSet);
}

// ---------------------------------------------------------

template <typename SolverT>
template<typename SetType>
typename SolverT::VectorType TimeMap<SolverT>::operator()(ScalarType Time, SetType &originalSet, MatrixType& derivative)
{
  m_solution = 0;
  this->moveSet(Time,originalSet,&originalSet);
  derivative = (MatrixType) originalSet;
  return VectorType(originalSet);
}

// ---------------------------------------------------------

template <typename SolverT>
template <typename CnCoeffType>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, CnCoeffType& result)
// the point on trajectory just after time 'time'
{
  if(this->m_currentSet==&result)
  {
    this->moveSet(Time, result,&result);
    return VectorType(result);
  }

  result.clear();
  result() = v;
  result.setMatrix(MatrixType::Identity(v.dimension()));
  result.setCurrentTime(TypeTraits<ScalarType>::zero());
  this->moveSet(Time, result,&result);
  return VectorType(result);
}

/// @}
}} // namespace capd::poincare


#endif

