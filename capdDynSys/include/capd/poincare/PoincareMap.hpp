/////////////////////////////////////////////////////////////////////////////
/// @file  PoincareMap.hpp
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_POINCARE_MAP_HPP_
#define _CAPD_POINCARE_POINCARE_MAP_HPP_

#include <cassert>
#include <exception>
#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/BasicPoincareMap.hpp"
#include "capd/dynset/lib.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
/**
 *   Constructs PoincareMap for given dynamical system and section
 *
 *  @param ds        dynamical system
 *  @param section   Poincare section
 *  @param crossing  section crossing direction
 *  @param factor    time step correction factor during section crossing (should be in interval (0, 1))
**/
template <typename SolverT, typename SectionT>
PoincareMap<SolverT,SectionT>::PoincareMap(
       Solver  & ds,
       SectionType & section,
       typename BasicPoincareMap<SolverT,SectionT>::CrossingDirection crossing,
       const BoundType & errorTolerance
  ) : BasicPoincareMap<SolverT,SectionT>(ds, section, crossing, errorTolerance),
  timeStepCorrectionFactor(0.9),
  maxCorrectionAttempts(10),
  minCrossingTimeStep( capd::TypeTraits< ScalarType >::epsilon()*4),
  sectionCoordinateSystem(MatrixType::Identity(ds.dimension()))
{}

template <typename SolverT, typename SectionT>
typename PoincareMap<SolverT,SectionT>::VectorType 
PoincareMap<SolverT,SectionT>::operator()(VectorType x, int n)
{ 
  DefaultC0Set<MatrixType> s(x); 
  return (*this)(s,n); 
}

template <typename SolverT, typename SectionT>
typename PoincareMap<SolverT,SectionT>::VectorType 
PoincareMap<SolverT,SectionT>::operator()(const VectorType& x, const VectorType& x0, const MatrixType& A, ScalarType& out_returnTime, int n)
{ 
  DefaultC0Set<MatrixType> s(x); 
  return (*this)(s,x0,A,out_returnTime,n); 
}

template <typename SolverT, typename SectionT>
typename PoincareMap<SolverT,SectionT>::VectorType 
PoincareMap<SolverT,SectionT>::operator()(const VectorType& x, ScalarType& out_returnTime, int n)
{ 
  DefaultC0Set<MatrixType> s(x); 
  return (*this)(s,out_returnTime,n); 
}

template <typename SolverT, typename SectionT>
typename PoincareMap<SolverT,SectionT>::VectorType 
PoincareMap<SolverT,SectionT>::operator()(const VectorType& x, MatrixType& D, int n)
{ 
  DefaultC1Set<MatrixType> s(x); 
  return (*this)(s,D,n); 
}

template <typename SolverT, typename SectionT>
typename PoincareMap<SolverT,SectionT>::VectorType 
PoincareMap<SolverT,SectionT>::operator()(const VectorType& x, MatrixType& D, ScalarType& out_returnTime, int n)
{ 
  DefaultC1Set<MatrixType> s(x); 
  return (*this)(s,D,out_returnTime,n); 
}

/// @}
}} // namespace capd::poincare

#endif // _CAPD_POINCARE_POINCARE_MAP_HPP_

