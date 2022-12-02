

/////////////////////////////////////////////////////////////////////////////
/// @file C2OdeSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C2ODESOLVER_HPP_
#define _CAPD_DYNSYS_C2ODESOLVER_HPP_

#include <string>
#include <stdexcept>

#include "capd/dynsys/C2OdeSolver.h"
#include "capd/dynsys/OdeSolver.hpp"
#include "capd/dynsys/BasicC2OdeSolver.hpp"
#include "capd/vectalg/Norm.hpp"
#include "approveRemainder.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
//###########################################################//

template<typename MapType,typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
C2OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveT>::C2OdeSolver(VectorFieldType& vectorField,size_type order)
  : BasicOdeSolver<VectorFieldType,StepControlPolicy,CurveT>(vectorField,order),
    OdeSolver<VectorFieldType,StepControlPolicy,EnclosurePolicy,CurveT>(vectorField,order),
    BasicC2OdeSolver<VectorFieldType,StepControlPolicy,CurveT>(vectorField,order)
{}

// ####################################################################

template<typename MapType,typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void C2OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveT>::computeRemainderCoefficients(const VectorType& x, const MatrixType& M, const HessianType& H)
{
  VectorType* remCoeff = this->getRemainderCoefficients();
  MatrixType* matrixRemCoeff = this->getMatrixRemainderCoefficients();
  HessianType* hessianRemCoeff = this->getHessianRemainderCoefficients();
  remCoeff[0] = x;
  matrixRemCoeff[0] = M;
  hessianRemCoeff[0] = H;
  this->m_vField->computeODECoefficients(remCoeff,matrixRemCoeff,hessianRemCoeff,this->getOrder()+1);
}

// ####################################################################

template<typename MapType,typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void C2OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveT>::encloseC2Map(
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
  VectorType* coeff = this->getCoefficients();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  HessianType* hessianCoeff = this->getHessianCoefficients();

  const int order=this->getOrder();
  coeff[0] = xx;
  matrixCoeff[0].setToIdentity();
  hessianCoeff[0].clear();
  this->m_vField->computeODECoefficients(coeff,matrixCoeff,hessianCoeff,order);
  this->computeCoefficientsAtCenter(x,order);

  capd::diffAlgebra::C2TimeJet<MatrixType> rem(&o_rem,&o_jacRem,&o_hessianRem);
  capd::diffAlgebra::C2TimeJet<MatrixType> enc(&o_enc,&o_jacEnc,&o_hessianEnc);

  capd::dynsys::computeAndApproveRemainder(*this,t,xx,rem,enc);
  this->sumTaylorSeries(o_phi,o_jacPhi,o_hessianPhi,this->getCoefficientsAtCenter(),matrixCoeff,hessianCoeff,order);
}

// ####################################################################

template<typename MapType,typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
void C2OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveT>::computeRemainder(ScalarType t, const VectorType& xx, C2TimeJetType& o_enc, C2TimeJetType& o_rem)
{
  EnclosurePolicy::computeEnclosureAndRemainder(*this,t,xx,o_enc,o_rem);
}
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_C2ODESOLVER_HPP_


