

/////////////////////////////////////////////////////////////////////////////
/// @file CnOdeSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_CNODESOLVER_HPP_
#define _CAPD_DYNSYS_CNODESOLVER_HPP_

#include "capd/dynsys/DynSys.hpp"
#include "capd/dynsys/BasicCnOdeSolver.hpp"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/dynsys/CnOdeSolver.h"
#include "capd/dynsys/CnOdeSolver_phi.hpp"
#include "capd/dynsys/CnOdeSolver_remainder.hpp"
#include "capd/dynsys/approveRemainder.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::CnOdeSolver(MapT& vectorField, size_type order, const StepControlPolicy& stepControl)
  : BasicCnOdeSolver<MapT,StepControlPolicy,CurveT>(vectorField,order,stepControl){}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::setInitialCondition(ScalarType t, const VectorType& x, const VectorType& xx)
{
  this->setCurrentTime(t);

  for(size_type i=0;i<this->dimension();++i){
    this->getCoefficientsAtCenter()[0][i] = x[i];
    this->coefficient(i,0) = xx[i];
    for(size_type j=0;j<this->dimension();++j) {
      this->coefficient(i,j,0) =
          (i==j) ? TypeTraits<ScalarType>::one() : TypeTraits<ScalarType>::zero();
    }
  }
}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::encloseC2Map(
    const ScalarType& t,
    const VectorType& x,
    const VectorType& xx,
    VectorType& o_phi,
    VectorType& o_rem,
    VectorType& o_enc,
    MatrixType& o_jacPhi,
    MatrixType& o_jacRem,
    MatrixType& o_jacEnc,
    HessianType& o_hessianPhi,
    HessianType& o_hessianRem,
    HessianType& o_hessianEnc
  )
{
  this->setInitialCondition(t,x,xx);
  for(size_type i=0;i<this->dimension();++i)
    for(size_type j=0;j<this->dimension();++j)
      for(size_type c=j;c<this->dimension();++c)
        this->coefficient(i,j,c,0) = TypeTraits<ScalarType>::zero();

  this->m_vField->computeODECoefficients(this->getCoefficientsAtCenter(),this->getOrder());
  this->m_vField->computeODECoefficients(this->getCoefficients(),2,this->getOrder());

  capd::diffAlgebra::C2TimeJet<MatrixType> rem(&o_rem,&o_jacRem,&o_hessianRem);
  capd::diffAlgebra::C2TimeJet<MatrixType> enc(&o_enc,&o_jacEnc,&o_hessianEnc);

  capd::dynsys::computeAndApproveRemainder(*this,t,xx,rem,enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,o_hessianPhi);
}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::encloseC1Map(
    const ScalarType& t,
    const VectorType& x,
    const VectorType& xx,
    VectorType& o_phi,
    VectorType& o_rem,
    VectorType& o_enc,
    MatrixType& o_jacPhi,
    MatrixType& o_jacRem,
    MatrixType& o_jacEnc
  )
{
  this->setInitialCondition(t,x,xx);
  this->m_vField->computeODECoefficients(this->getCoefficientsAtCenter(),this->getOrder());
  this->m_vField->computeODECoefficients(this->getCoefficients(),1,this->getOrder());

  capd::diffAlgebra::C1TimeJet<MatrixType> rem(&o_rem,&o_jacRem);
  capd::diffAlgebra::C1TimeJet<MatrixType> enc(&o_enc,&o_jacEnc);

  capd::dynsys::computeAndApproveRemainder(*this,t,xx,rem,enc);
  this->sumTaylorSeries(o_phi,o_jacPhi);
}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::encloseC0Map(
    const ScalarType& t,
    const VectorType& x,
    const VectorType& xx,
    VectorType& o_phi,
    VectorType& o_rem,
    VectorType& o_enc,
    MatrixType& o_jacPhi
  )
{
  this->setInitialCondition(t,x,xx);
  this->m_vField->computeODECoefficients(this->getCoefficientsAtCenter(),this->getOrder());
  this->m_vField->computeODECoefficients(this->getCoefficients(),1,this->getOrder());

  capd::dynsys::computeAndApproveRemainder(*this,t,xx,o_rem,o_enc);

  this->sumTaylorSeries(o_phi,o_jacPhi);
}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::sumTaylorSeries(VectorType& o_phi, MatrixType& o_jacPhi)
{
  size_type i,j;
  const int order = this->getOrder();
  for(i=0;i<this->dimension();++i){
    o_phi[i] = this->getCoefficientsAtCenter()[order][i];
    for(j=0;j<this->dimension();++j)
      o_jacPhi(i+1,j+1) = this->coefficient(i,j,order);
  }

  for(int r=order-1;r>=0;--r){
    for(i=0;i<this->dimension();++i){
      o_phi[i] = o_phi[i]*this->m_step + this->getCoefficientsAtCenter()[r][i];
      for(j=0;j<this->dimension();++j)
        o_jacPhi(i+1,j+1) = o_jacPhi(i+1,j+1)*this->m_step + this->coefficient(i,j,r);
    }
  }
}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::sumTaylorSeries(VectorType& o_phi, MatrixType& o_jacPhi, HessianType& o_hessianPhi)
{
  size_type i,j,c;
  const int order = this->getOrder();

  this->sumTaylorSeries(o_phi,o_jacPhi);

  for(i=0;i<this->dimension();++i)
    for(j=0;j<this->dimension();++j)
      for(c=j;c<this->dimension();++c)
        o_hessianPhi(i,j,c) = this->coefficient(i,j,c,order);

  for(int r=order-1;r>=0;--r)
    for(i=0;i<this->dimension();++i)
      for(j=0;j<this->dimension();++j)
        for(c=j;c<this->dimension();++c)
          o_hessianPhi(i,j,c) = o_hessianPhi(i,j,c)*this->m_step + this->coefficient(i,j,c,r);
}

// ####################################################################

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
typename CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::ScalarType
CnOdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT>::getCoeffNorm(size_type r, size_type degree) const
{
  typename TypeTraits<ScalarType>::Real result = 0;
  size_type i;
  for(i=0;i<this->dimension();++i){
    result = capd::max(result,rightBound(this->remainderCoefficient(i,r)) - leftBound(this->remainderCoefficient(i,r)));
    result = capd::max(result,rightBound(this->coefficient(i,r)) - leftBound(this->coefficient(i,r)));
  }

  if(degree)
  {
    const typename CurveT::JetType& c = this->getCoefficients()[r];
    const typename CurveT::JetType& rem = this->getRemainderCoefficients()[r];

    for(i=0;i<this->dimension();++i)
    {
      typename CurveT::JetType::const_iterator b = c.begin(i,1), e=c.end(i,degree);
      typename CurveT::JetType::const_iterator b1 = rem.begin(i,1);
      while(b!=e)
      {
        result = capd::max(result,rightBound(*b)-leftBound(*b));
        result = capd::max(result,rightBound(*b1)-leftBound(*b1));
        ++b;
        ++b1;
      }
    }
  }

  return ScalarType(result);
}
/// @}
}}

#endif // _CAPD_DYNSYS_CNODESOLVER_HPP_


