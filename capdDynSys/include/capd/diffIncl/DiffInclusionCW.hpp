/// @addtogroup diffIncl
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusionCW.hpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Tomasz Kapela, 2007 */

#ifndef _CAPD_DIFFINCL_DIFFINCLUSIONCW_HPP_
#define _CAPD_DIFFINCL_DIFFINCLUSIONCW_HPP_

#include <sstream>
#include <string>
#include <stdexcept>

#include "capd/diffIncl/DiffInclusionCW.h"
#include "capd/diffIncl/DiffInclusion.hpp"
#include "capd/vectalg/iobject.hpp"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
namespace diffIncl{


template <typename MapT, typename DynSysT>
DiffInclusionCW<MapT, DynSysT>::DiffInclusionCW(
           MultiMapType& A_diffIncl,
           size_type A_order,
           NormType const & A_norm,
           ScalarType const & expErrorTolerance
)
  : BaseClass(A_diffIncl, A_order, A_norm), m_errorTolerance(expErrorTolerance) {
}

/**
 * Bounds for perturbation of solution of selected ODE after one time step
 *
 * We use Component Wise algorithm
 *
 * @param x  set before the the step
 * @return   rigorous bounds for perturbations
 */
template <typename MapT, typename DynSysT>
typename DiffInclusionCW<MapT, DynSysT>::VectorType DiffInclusionCW<MapT, DynSysT>::perturbations(const ScalarType & time, const VectorType& x){

  VectorType W_1 = dynamicalSystemEnclosure(time, x);
  VectorType W_2 = diffInclusionEnclosure(time, x);

  MatrixType J = m_diffIncl.derivative(time, W_2);
  VectorType delta = m_diffIncl.perturbations(time, W_1);
  
  if(!capd::vectalg::containsZero(delta) )
    throw std::runtime_error("DiffInclusionCW::perturbations error: perturbations of the vector field do not contain 0 (selection not contained in the differential inclusion).");

  VectorType C = rightVector(abs(delta));

  size_type i, j;
  for( i=0; i < J.numberOfRows(); ++i)
    for( j = 0; j < J.numberOfColumns(); ++j)
      if(i != j)
        J[i][j] = (capd::abs(J[i][j])).right();
      else
        J[i][j] = (J[i][j]).right();
  MatrixType At = getStep() * J;
  MatrixType A = MatrixType::Identity(J.numberOfRows());
  MatrixType Sum = A;
  ScalarType AtNorm = right((*m_norm)(At)),
      AnNorm = right((*m_norm)(A));

  ScalarType n = typename capd::TypeTraits<ScalarType>::Real(2.0);  // n = i + 2
  ScalarType q = AtNorm/n;

  while(true){
    A = A * At / n;
    Sum += A;
    AnNorm *= q;
    n += ScalarType(1.0);
    q = AtNorm / n;
    if(q < 1){
      // remainder = |A| * |At/(N+2)| / (1 - |At/(N+2)|)   (the sum of geometric series from N to infinity)
      ScalarType remainder = right(AnNorm * AtNorm / (n - AtNorm ));
      if(remainder < m_errorTolerance)
        break;
    }
  }
  //std::cout << "\n N = " << n-2 << "\n remainder approx : " << remainder;
  // we recompute remainder because norm of A can be smaller than approximation : AnNorm
  ScalarType remainder = right((*m_norm)(A) * AtNorm / (n  - AtNorm ));
  //std::cout << "\n remainder real   : " << remainder;

  for(i=0; i < J.numberOfRows(); ++i)
    for(j = 0; j < J.numberOfColumns(); ++j)
      Sum[i][j] += remainder * ScalarType(-1.0, 1.0);
  VectorType D = getStep() * (Sum * C);
  VectorType result(D.dimension());

  for(i=0; i< D.dimension(); ++i)
    result[i] = ScalarType(-D[i].rightBound(), D[i].rightBound());
  return result;
}

}} //namespace capd::diffIncl

#endif // _CAPD_DIFFINCL_DIFFINCLUSIONCW_HPP_

/// @}
