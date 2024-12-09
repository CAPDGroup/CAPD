

/////////////////////////////////////////////////////////////////////////////
/// @file CnOdeSolver_remainder.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_CNODESOLVER_REMAINDER_HPP_
#define _CAPD_DYNSYS_CNODESOLVER_REMAINDER_HPP_

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<typename MapType, typename StepControlPolicy, typename EnclosurePolicy,typename CurveType>
typename CnOdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::VectorType
CnOdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::Remainder(const ScalarType& t, const VectorType &iv, VectorType& o_enc)
{
  o_enc = this->enclosure(t,iv);
  const static ScalarType I(TypeTraits<ScalarType>::zero().leftBound(),TypeTraits<ScalarType>::one().rightBound());
  this->computeRemainderCoefficients(t + I*this->m_step,o_enc);
  size_type r=this->getOrder()+1;

  ScalarType fac = power(this->m_step,r);
  VectorType result(this->dimension(),true);
  for(size_type i=0;i<this->dimension();++i)
    result[i] = this->remainderCoefficient(i,r)*fac;

  return result;
}

//###########################################################//

template<typename MapType, typename StepControlPolicy, typename EnclosurePolicy,typename CurveType>
void CnOdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::JacRemainder(
      const ScalarType& t,
      const VectorType& vecEnclosure,
      const MatrixType& jacEnclosure,
      VectorType& Remainder,
      MatrixType& jacRemainder
    )
{
  const static ScalarType I(TypeTraits<ScalarType>::zero().leftBound(),TypeTraits<ScalarType>::one().rightBound());
  size_type r=this->getOrder()+1;
  this->computeRemainderCoefficients(t + I*this->m_step,vecEnclosure,jacEnclosure);
  ScalarType fac = power(this->m_step,r);

  for(size_type i=0;i<this->dimension();++i)
  {
    Remainder[i] = fac*this->remainderCoefficient(i,r);
    for(size_type j=0;j<this->dimension();++j)
      jacRemainder(i+1,j+1) = fac*this->remainderCoefficient(i,j,r);
  }
}

// ####################################################################

template <typename MapType, typename StepControlPolicy, typename EnclosurePolicy, typename CurveType>
void CnOdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::c2Remainder(
      const VectorType& Enc,
      const MatrixType& jacEnc,
      const HessianType& hessianEnc,
      VectorType& o_Rem,
      MatrixType& o_jacRem,
      HessianType& o_hessianRem
  )
{
  size_type i,j,c;
  // set initial condition
  for(i=0;i<this->dimension();++i)
  {
    this->remainderCoefficient(i,0) = Enc[i];
    for(j=0;j<this->dimension();++j)
    {
      this->remainderCoefficient(i,j,0) = jacEnc(i+1,j+1);
      for(c=j;c<this->dimension();++c)
        this->remainderCoefficient(i,j,c,0) = hessianEnc(i,j,c);
    }
  }

  this->m_vField->computeODECoefficients(this->getRemainderCoefficients(),2,this->getOrder()+1);

  ScalarType fac = power(this->m_step,this->getOrder()+1);
  for(i=0;i<this->dimension();++i)
  {
    o_Rem[i] = fac*this->remainderCoefficient(i,this->getOrder()+1);
    for(j=0;j<this->dimension();++j)
    {
      o_jacRem(i+1,j+1) = fac*this->remainderCoefficient(i,j,this->getOrder()+1);
      for(c=j;c<this->dimension();++c)
        o_hessianRem(i,j,c) = fac*this->remainderCoefficient(i,j,c,this->getOrder()+1);
    }
  }
}

// ---------------------------------------------------------------------------------------

template<typename MapType, typename StepControlPolicy, typename EnclosurePolicy,typename CurveType>
void CnOdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainder(ScalarType t, const VectorType& xx, VectorType& o_enc, VectorType& o_rem)
{
  o_rem = this->Remainder(t,xx,o_enc);
}

// ---------------------------------------------------------------------------------------

template<typename MapType, typename StepControlPolicy, typename EnclosurePolicy,typename CurveType>
void CnOdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainder(ScalarType t, const VectorType& xx, C1TimeJetType& o_enc, C1TimeJetType& o_rem)
{
  o_enc.vector() = this->enclosure(t,xx);
  o_enc.matrix() = EnclosurePolicy::jacEnclosure(this->getVectorField(),t,this->getStep(),o_enc.vector(),vectalg::EuclLNorm<VectorType,MatrixType>());
  this->JacRemainder(t,o_enc.vector(),o_enc.matrix(),o_rem.vector(),o_rem.matrix());
}

// ---------------------------------------------------------------------------------------

template<typename MapType, typename StepControlPolicy, typename EnclosurePolicy,typename CurveType>
void CnOdeSolver<MapType, StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainder(ScalarType t, const VectorType& xx, C2TimeJetType& o_enc, C2TimeJetType& o_rem)
{
  o_enc.vector() = this->enclosure(t,xx);
  EnclosurePolicy::c2Enclosure(this->getVectorField(),this->getStep(),o_enc.vector(),o_enc.matrix(),o_enc.hessian());
  this->c2Remainder(o_enc.vector(),o_enc.matrix(),o_enc.hessian(),o_rem.vector(),o_rem.matrix(),o_rem.hessian());
}

/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_CNODESOLVER_REMAINDER_HPP_


