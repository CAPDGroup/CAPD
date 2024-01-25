/// @addtogroup diffIncl
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusionLN.hpp
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Tomasz Kapela, 2007 */

#ifndef _CAPD_DIFFINCL_DIFFINCLUSIONLN_HPP_ 
#define _CAPD_DIFFINCL_DIFFINCLUSIONLN_HPP_ 

#include "capd/diffIncl/DiffInclusionLN.h"
#include "capd/diffIncl/DiffInclusion.hpp"

namespace capd{
namespace diffIncl{


template <typename MapT, typename DynSysT>
DiffInclusionLN<MapT, DynSysT>::DiffInclusionLN(
    MultiMapType& A_diffIncl,
    size_type A_order,
    NormType const & A_norm
) 
: BaseClass(A_diffIncl, A_order, A_norm) {
}

template <typename MapT, typename DynSysT>
typename DiffInclusionLN<MapT, DynSysT>::VectorType DiffInclusionLN<MapT, DynSysT>::perturbations(const ScalarType & time, const VectorType& x){

  VectorType W_1 = dynamicalSystemEnclosure(time, x);
  VectorType W_2 = diffInclusionEnclosure(time, x);

  MatrixType J = m_diffIncl.getVectorField().derivative(time, W_2);
  VectorType delta = m_diffIncl.perturbations(time, W_2); 

  if(!capd::vectalg::containsZero(delta) )
    throw std::runtime_error("DiffInclusionCW::perturbations error: perturbations of the vector field do not contain 0 (selection not contained in the differential inclusion).");

  ScalarType C = right((*m_norm)(delta));
  ScalarType l = right((*m_norm)(J));

  ScalarType D = (l.contains(0.0))? C*getStep() : (C*(exp(l*getStep())-1))/l;
  VectorType result(x.dimension());

  for(size_type i=0; i< result.dimension(); ++i)
    result[i] = ScalarType(-D.rightBound(), D.rightBound());

  return result;
}

}} //namespace capd::diffIncl

#endif // _CAPD_DIFFINCL_DIFFINCLUSIONLN_HPP_ 

/// @}
