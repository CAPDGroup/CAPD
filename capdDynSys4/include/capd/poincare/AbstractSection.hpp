

/////////////////////////////////////////////////////////////////////////////
/// @file AbstractSection.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_ABSTRACT_SECTION_HPP_
#define _CAPD_POINCARE_ABSTRACT_SECTION_HPP_

#include "capd/poincare/AbstractSection.h"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
// -----------------------------------------------------------------------------------------

template <class MatrixT>
void AbstractSection<MatrixT>::computeDT(
      const MatrixType& derivativeOfFlow,
      const VectorType& gradientOnPx,
      const ScalarType& denominator,
      VectorType& result
   ) const
{
  for (size_type j = 0; j < gradientOnPx.dimension(); ++j)
  {
    ScalarType val = TypeTraits<ScalarType>::zero();
    for (size_type k = 0; k < gradientOnPx.dimension(); ++k)
      val -= gradientOnPx[k] * derivativeOfFlow(k + 1, j + 1);
    result[j] = val / denominator;
  }
}

// -------------------------------------------------------------------

template <class MatrixT>
typename AbstractSection<MatrixT>::MatrixType
AbstractSection<MatrixT>::computeDP(
      const VectorType& Px,
      const MatrixType& derivativeOfFlow,
      const VectorType& fieldOnPx,
      VectorType& dT
   ) const
{
  const size_type dim = Px.dimension();
  MatrixType result(dim,dim);
  VectorType gradientOnPx = this->gradient(Px);
  ScalarType denominator = fieldOnPx*gradientOnPx;

  this->computeDT(derivativeOfFlow,gradientOnPx,denominator,dT);

  for(size_type i=1;i<=dim;++i)
    for(size_type j=1;j<=dim;++j)
      result(i,j) = fieldOnPx[i-1]*dT[j-1] + derivativeOfFlow(i,j);
   return result;
}

// -------------------------------------------------------------------

template <class MatrixT>
void AbstractSection<MatrixT>::computeDP(
      const VectorType& Px,
      const MatrixType& derivativeOfFlow,
      const HessianType& hessianOfFlow,
      const VectorType& fieldOnPx,
      const VectorType& d2Phidt2,
      const MatrixType& derOfVectorFieldOnPx,
      MatrixType& DP,
      HessianType& D2P,
      VectorType& dT,
      MatrixType& d2T
    ) const
{
  size_type i,j,k;
  const size_type dim = Px.dimension();
  VectorType gradientOnPx = this->gradient(Px);
  ScalarType denominator = fieldOnPx*gradientOnPx;

  // here we compute first order derivatives of return time.
  this->computeDT(derivativeOfFlow,gradientOnPx,denominator,dT);

  // first derivative of Poincare map
  for(j=0;j<dim;++j)
    for(i=0;i<dim;++i)
      DP(j+1,i+1) = fieldOnPx[j]*dT[i] + derivativeOfFlow(j+1,i+1);

  // second derivatives of return time
  MatrixType dPhiVectorFieldOnPPx = derOfVectorFieldOnPx*derivativeOfFlow;

  for(i=0;i<dim;++i)
  {
    ScalarType d2Tii = TypeTraits<ScalarType>::zero();
    for(k=0;k<dim;++k)
    {
      D2P(k,i,i) = hessianOfFlow(k,i,i) + dT[i]*( dPhiVectorFieldOnPPx(k+1,i+1) + d2Phidt2[k] * dT[i]) ;
      d2Tii -= gradientOnPx[k]* D2P(k,i,i);
    }
    d2T(i+1,i+1) = d2Tii/denominator;

    for(j=i+1;j<dim;++j)
    {
      ScalarType d2Tij = TypeTraits<ScalarType>::zero();
      for(k=0;k<dim;++k)
      {
        D2P(k,i,j) = hessianOfFlow(k,i,j)
                   + dPhiVectorFieldOnPPx(k+1,i+1) * dT[j]
                   + dPhiVectorFieldOnPPx(k+1,j+1) * dT[i]
                   + 2.*d2Phidt2[k] * dT[i]*dT[j]
                                               ;
        d2Tij -= gradientOnPx[k]*D2P(k,i,j);
      }
      d2T(j+1,i+1) = d2T(i+1,j+1) = d2Tij/denominator;
    }
  }

  // second derivatives of Poincare map
  for(k=0;k<dim;++k)
    for(i=0;i<dim;++i)
      for(j=i;j<dim;++j)
        D2P(k,i,j) += fieldOnPx[k] * d2T(i+1,j+1);

}
/// @}
}} // namespace capd::poincare

#endif  /* _CAPD_POINCARE_ABSTRACT_SECTION_HPP_ */


