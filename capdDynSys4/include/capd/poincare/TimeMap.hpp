/////////////////////////////////////////////////////////////////////////////
/// @file TimeMap.hpp
///
/// @author Tomasz Kapela, Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_TIME_MAP_HPP_
#define _CAPD_POINCARE_TIME_MAP_HPP_

#include "capd/poincare/TimeMap.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/diffAlgebra/CnTimeJet.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
template <typename SolverT>
TimeMap<SolverT>::TimeMap(SolverT &solver)
  : m_solver(solver),
    m_currentSet(0),
    m_oneStep(false),
    m_currentTime(0.),
    m_completed(true),
    m_solution(0)
{}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v)
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,t);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, ScalarType& in_out_time)
{
  m_solution = 0;
  capd::diffAlgebra::C0TimeJet<VectorType> c(&v);
  c.setCurrentTime(in_out_time);
  this->moveSet(Time, c,(void*)&v);
  in_out_time = c.getCurrentTime();
  return v;
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, SolutionCurve& solution)
{
  capd::diffAlgebra::C0TimeJet<VectorType> c(&v);
  c.setCurrentTime(solution.getLeftDomain());
  m_solution = &solution;
  this->moveSet(Time, c,(void*)&v);
  return v;
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, MatrixType &der)
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,MatrixType::Identity(v.dimension()),der,t);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, MatrixType &der, ScalarType& in_out_time)
{
  return (*this)(Time,v,MatrixType::Identity(v.dimension()),der,in_out_time);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, MatrixType &der, SolutionCurve& solution)
{
  return (*this)(Time,v,MatrixType::Identity(v.dimension()),der,solution);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, const MatrixType& initMatrix, MatrixType &der)
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,initMatrix,der,t);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, const MatrixType& initMatrix, MatrixType& der, ScalarType& in_out_time)
{
  m_solution = 0;
  der = initMatrix;
  capd::diffAlgebra::C1TimeJet<MatrixType> c(&v,&der);
  c.setCurrentTime(in_out_time);

  this->moveSet(Time, c,(void*)&v);
  in_out_time = c.getCurrentTime();
  return v;
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, const MatrixType& initMatrix, MatrixType& der, SolutionCurve& solution)
{
  m_solution = &solution;
  der = initMatrix;
  capd::diffAlgebra::C1TimeJet<MatrixType> c(&v,&der);
  c.setCurrentTime(solution.getLeftDomain());
  this->moveSet(Time, c,(void*)&v);
  return v;
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, MatrixType& der, HessianType& hessian)
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,MatrixType::Identity(v.dimension()),HessianType(v.dimension()), der, hessian, t);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, MatrixType& der, HessianType& hessian, ScalarType& in_out_time)
{
  return (*this)(Time,v,MatrixType::Identity(v.dimension()),HessianType(v.dimension()),der,hessian,in_out_time);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, MatrixType& der, HessianType& hessian, SolutionCurve& solution)
{
  return (*this)(Time,v,MatrixType::Identity(v.dimension()),HessianType(v.dimension()),der,hessian,solution);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(
    ScalarType Time,
    VectorType& v, const MatrixType& initMatrix, const HessianType& initHessian,
    MatrixType& der, HessianType& hessian
    )
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,initMatrix,initHessian,der,hessian,t);
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(
    ScalarType Time,
    VectorType& v, const MatrixType& initMatrix, const HessianType& initHessian,
    MatrixType &der, HessianType& hessian,
    ScalarType& in_out_time
    )
{
  m_solution = 0;
  der = initMatrix;
  hessian = initHessian;
  capd::diffAlgebra::C2TimeJet<MatrixType> c(&v,&der,&hessian);
  c.setCurrentTime(in_out_time);
  this->moveSet(Time, c,(void*)&v);
  in_out_time = c.getCurrentTime();
  return v;
}

// ---------------------------------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time,
    VectorType& v, const MatrixType& initMatrix, const HessianType& initHessian,
    MatrixType& der, HessianType& hessian,
    SolutionCurve& solution)
{
  m_solution = &solution;
  der = initMatrix;
  hessian = initHessian;
  capd::diffAlgebra::C2TimeJet<MatrixType> c(&v,&der,&hessian);
  c.setCurrentTime(solution.getLeftDomain());
  this->moveSet(Time, c,(void*)&v);
  return v;
}

// ---------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, JetType& result)
// the point on trajectory just after time 'time'
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,result,t);
}

// ---------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, JetType& result, ScalarType& in_out_time)
// the point on trajectory just after time 'time'
{
  m_solution = 0;

  result.clear();
  result() = v;
  result.setMatrix(MatrixType::Identity(v.dimension()));
  capd::diffAlgebra::CnTimeJet<MatrixType,0> c(&result);
  c.setCurrentTime(in_out_time);
  this->moveSet(Time, c,(void*)&v);
  in_out_time = c.getCurrentTime();
  return (VectorType)result;
}

// ---------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, VectorType& v, JetType& result, SolutionCurve& solution)
// the point on trajectory just after time 'time'
{
  m_solution = &solution;
  result.clear();
  result() = v;
  result.setMatrix(MatrixType::Identity(v.dimension()));
  capd::diffAlgebra::CnTimeJet<MatrixType,0> c(&result);
  c.setCurrentTime(solution.getLeftDomain());
  this->moveSet(Time, c,(void*)&v);
  return (VectorType)result;
}

// ---------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, const JetType& v, JetType& result)
// the point on trajectory just after time 'time'
{
  ScalarType t = TypeTraits<ScalarType>::zero();
  return (*this)(Time,v,result,t);
}

// ---------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, const JetType& v, JetType& result, ScalarType& in_out_time)
// the point on trajectory just after time 'time'
{
  m_solution = 0;
  result = v;
  capd::diffAlgebra::CnTimeJet<MatrixType,0> c(&result);
  c.setCurrentTime(in_out_time);
  this->moveSet(Time, c,(void*)&v);
  in_out_time = c.getCurrentTime();
  return (VectorType)result;
}

// ---------------------------------------------------------

template <typename SolverT>
typename TimeMap<SolverT>::VectorType
TimeMap<SolverT>::operator()(ScalarType Time, const JetType& v, JetType& result, SolutionCurve& solution)
// the point on trajectory just after time 'time'
{
  m_solution = &solution;
  result = v;
  capd::diffAlgebra::CnTimeJet<MatrixType,0> c(&result);
  c.setCurrentTime(solution.getLeftDomain());
  this->moveSet(Time, c,(void*)&v);
  return (VectorType)result;
}

/// @}
}} // namespace capd::poincare


#endif // #define _CAPD_POINCARE_TIME_MAP_HPP_

