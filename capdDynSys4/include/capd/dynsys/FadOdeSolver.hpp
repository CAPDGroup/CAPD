

/////////////////////////////////////////////////////////////////////////////
/// @file FadOdeSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2017 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_FADODESOLVER_HPP_
#define _CAPD_DYNSYS_FADODESOLVER_HPP_

#include "capd/basicalg/power.h"
#include "capd/dynsys/FadOdeSolver.h"
#include "capd/dynsys/BasicFadOdeSolver.hpp"
#include "capd/vectalg/Norm.hpp"
#include "capd/dynsys/FirstOrderEnclosure.hpp"
#include "capd/dynsys/DynSys.hpp"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<class FadMapT, typename StepControlT>
FadOdeSolver<FadMapT,StepControlT>::FadOdeSolver(VectorFieldType& f, size_type _order, const StepControlT& _stepControl)
  : BaseTaylor(f,_order,_stepControl)
{}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename FadOdeSolver<FadMapT,StepControlT>::VectorType
FadOdeSolver<FadMapT,StepControlT>::Phi(const ScalarType& t, const VectorType& v)
{
  VectorType u(v.dimension());
  this->setCurrentTime(t);
  this->setInitialCondition(v,this->m_center);
  this->computeCoeff(this->m_center,this->m_centerOut,this->getOrder());
  this->sumTaylorSeries(u,this->m_center,this->getOrder());
  return u;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename FadOdeSolver<FadMapT,StepControlT>::VectorType
FadOdeSolver<FadMapT,StepControlT>::enclosure(const ScalarType& t, const VectorType& x)
{
  return FirstOrderEnclosure::enclosure(this->getVectorField(),t,x,this->getStep());
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename FadOdeSolver<FadMapT,StepControlT>::MatrixType
FadOdeSolver<FadMapT,StepControlT>::JacPhi(const ScalarType& t, const VectorType& u)
{
  MatrixType result(this->dimension(),this->dimension());
  this->setCurrentTime(t);
  this->setInitialCondition(u,this->m_in);
  this->computeCoeff(this->m_in,this->m_out,this->getOrder());
  this->sumTaylorSeries(result,this->m_in,this->getOrder());
  return result;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void FadOdeSolver<FadMapT,StepControlT>::computeTaylorCoefficients(ScalarType t, const VectorType& x, const VectorType& xx){
  this->setCurrentTime(t);
  this->setInitialCondition(x,this->m_center);
  this->setInitialCondition(xx,this->m_in);
  this->computeCoeff(this->m_center,this->m_centerOut,this->getOrder());
  this->computeCoeff(this->m_in,this->m_out,this->getOrder());
}

//###########################################################//

template<typename MapType, typename StepControlType>
void FadOdeSolver<MapType,StepControlType>::encloseC0Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi
  )
{
  // here we compute all the coefficients for phi(t) and DPhi(t)
  this->computeTaylorCoefficients(t,x,xx);

  this->computeTimeStep(t,xx);

  // in the following function the time step can be adjusted
  o_rem = this->Remainder(t,xx,o_enc);

  this->sumTaylorSeries(o_phi,this->m_center,this->getOrder());
  this->sumTaylorSeries(o_jacPhi,this->m_in,this->getOrder());
}

//###########################################################//

template<typename MapType, typename StepControlType>
void FadOdeSolver<MapType,StepControlType>::encloseC1Map(
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
  // here we compute all the coefficients for phi(t) and DPhi(t)
  this->computeTaylorCoefficients(t,x,xx);

  this->computeTimeStep(t,xx);

  // in the following function the time step can be adjusted
  o_enc = this->enclosure(t,xx);
  o_jacEnc = this->jacEnclosure(t,o_enc);
  this->JacRemainder(t,o_enc,o_jacEnc,o_rem,o_jacRem);

  this->sumTaylorSeries(o_phi,this->m_center,this->getOrder());
  this->sumTaylorSeries(o_jacPhi,this->m_in,this->getOrder());
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename FadOdeSolver<FadMapT,StepControlT>::VectorType
FadOdeSolver<FadMapT,StepControlT>::Remainder(const ScalarType& t, const VectorType& v, VectorType& o_enc)
{
  typedef typename ScalarType::BoundType Real;
  const static ScalarType I(TypeTraits<Real>::zero(),TypeTraits<Real>::one());

  VectorType result(v.dimension());
  o_enc = this->enclosure(t,v);
  this->setCurrentTime(t+this->getStep()*I);
  this->computeRemainderCoefficients(o_enc);
  ScalarType fac = power(this->m_step,this->getOrder()+1);
  for(size_type i=0;i<v.dimension();++i)
    result[i] = fac*this->remainderCoefficient(i,this->getOrder()+1);
  return result;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
void FadOdeSolver<FadMapT,StepControlT>::JacRemainder(
        const ScalarType& t,
        const VectorType& vecEnclosure,
        const MatrixType& jacEnclosure,
        VectorType& Remainder,
        MatrixType& jRemainder
){
  typedef typename ScalarType::BoundType Real;
  const static ScalarType I(TypeTraits<Real>::zero(),TypeTraits<Real>::one());

  this->setCurrentTime(t+this->getStep()*I);
  this->computeRemainderCoefficients(vecEnclosure,jacEnclosure);
  ScalarType fac = power(this->m_step,this->getOrder()+1);
  for(size_type i=0;i<vecEnclosure.dimension();++i)
  {
    Remainder[i] = fac*this->remainderCoefficient(i,this->getOrder()+1);
    for(size_type j=0;j<vecEnclosure.dimension();++j)
      jRemainder(i+1,j+1) = fac*this->remainderCoefficient(i,j,this->getOrder()+1);
  }
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename FadOdeSolver<FadMapT,StepControlT>::MatrixType
FadOdeSolver<FadMapT,StepControlT>::jacEnclosure(const ScalarType& t, const VectorType& enc)
{
  return FirstOrderEnclosure::jacEnclosure(this->m_vectorField,t,this->m_step,enc,capd::vectalg::EuclLNorm<VectorType,MatrixType>());
}

//###########################################################//

template<class FadMapT, typename StepControlT>
typename FadOdeSolver<FadMapT,StepControlT>::ScalarType
FadOdeSolver<FadMapT,StepControlT>::getCoeffNorm(size_type r, size_type degree) const
{
  typename TypeTraits<ScalarType>::Real result = 0;
  for(size_type i=0;i<this->dimension();++i){
    result = capd::max(result,rightBound(this->remainderCoefficient(i,r)) - leftBound(this->remainderCoefficient(i,r)));
    result = capd::max(result,rightBound(this->coefficient(i,r)) - leftBound(this->coefficient(i,r)));
  }
  if(degree)
  {
    for(size_type i=0;i<this->dimension();++i){
      for(size_type j=0;j<this->dimension();++j){
        result = capd::max(result,rightBound(this->remainderCoefficient(i,j,r)) - leftBound(this->remainderCoefficient(i,j,r)));
        result = capd::max(result,rightBound(this->coefficient(i,j,r)) - leftBound(this->coefficient(i,j,r)));
      }
    }
  }
  return ScalarType(result);
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_FADODESOLVER_HPP_


