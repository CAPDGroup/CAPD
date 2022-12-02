

/////////////////////////////////////////////////////////////////////////////
/// @file BasicCnOdeSolver.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICCNODESOLVER_HPP_
#define _CAPD_DYNSYS_BASICCNODESOLVER_HPP_

#include "capd/basicalg/power.h"
#include "capd/dynsys/BasicCnOdeSolver.h"
#include "capd/dynsys/BasicOdeSolver.hpp"
#include "capd/diffAlgebra/BasicCnCurve.hpp"
#include "capd/diffAlgebra/CnCurve.hpp"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
template <typename MapT, typename StepControlT,typename CurveT>
BasicCnOdeSolver<MapT,StepControlT,CurveT>::~BasicCnOdeSolver(){}

// ---------------------------- CONSTRUCTORS ---------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
BasicCnOdeSolver<MapT,StepControlT,CurveT>::BasicCnOdeSolver(
      VectorFieldType& vectorField,
      size_type order,
      const StepControlType& stepControl
    )   : capd::dynsys::StepControlInterface<StepControlType,ScalarType>(stepControl),
          SolutionCurve(0.0,0.0,vectorField.dimension(),order,vectorField.degree()),
          m_vField(&vectorField),
          m_defaultC1InitialCondition(MatrixType::Identity(vectorField.dimension())),
          m_defaultC2InitialCondition(vectorField.dimension(),vectorField.dimension()),
          m_fixedTimeStep(TypeTraits<ScalarType>::zero()),
          m_step(TypeTraits<ScalarType>::zero())
{
  this->m_vField->setOrder(order + 1);
  this->m_vField->differentiateTime();
  Multiindex::generateList(this->dimension(),this->degree(),this->m_listIndices);
  if(order<=0) throw std::logic_error("BasicCnTaylor constructor: order must be a positive integer");
}

// -----------------------------------------------------------------------------


template <typename MapType, typename StepControlType, typename CurveT>
typename BasicCnOdeSolver<MapType, StepControlType, CurveT>::VectorType
BasicCnOdeSolver<MapType, StepControlType, CurveT>::operator()(VectorType v)
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  coeff[0] = v;
  this->evalAndSum(v);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnOdeSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT,StepControlT,CurveT>::operator()(VectorType v,MatrixType& der)
{
  this->setInitialCondition(v,m_defaultC1InitialCondition);
  this->evalAndSum(v,der);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnOdeSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT,StepControlT,CurveT>::operator()(
      VectorType v, const MatrixType& D, MatrixType& out_der
    )
{
  this->setInitialCondition(v,D);
  this->evalAndSum(v,out_der);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnOdeSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT,StepControlT,CurveT>::operator()(VectorType v,MatrixType& der, HessianType& hessian)
{
  this->setInitialCondition(v,m_defaultC1InitialCondition,m_defaultC2InitialCondition);
  this->evalAndSum(v,der,hessian);
  return v;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnOdeSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT,StepControlT,CurveT>::operator()(
      VectorType v, const MatrixType& D, const HessianType& H,
      MatrixType& out_der, HessianType& out_hessian
    )
{
  this->setInitialCondition(v,D,H);
  this->evalAndSum(v,out_der,out_hessian);
  return v;
}


// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
typename BasicCnOdeSolver<MapT,StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT,StepControlT,CurveT>::operator()(ScalarType& t, JetType& coeff)
{
  this->setCurrentTime(t);
  this->setInitialCondition(coeff);
  this->evalAndSum(coeff);
  t += this->m_step;
  return VectorType(coeff);
}


// -----------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
void BasicCnOdeSolver<MapType, StepControlType, CurveT>::setOrder(size_type order)
{
  if(order > this->getAllocatedOrder())
  {
    this->m_vField->setOrder(order + 1);
    this->m_vField->differentiateTime();
  }
  SolutionCurve::setOrder(order);
}

// ---------------------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
typename BasicCnOdeSolver<MapType, StepControlType, CurveT>::ScalarType
BasicCnOdeSolver<MapType, StepControlType, CurveT>::getCoeffNorm(size_type r, size_type degree) const
{
  ScalarType result = TypeTraits<ScalarType>::zero();
  size_type i;
  for(i=0;i<this->dimension();++i)
  {
    result = capd::max(result,capd::abs(this->centerCoefficient(i,r)));
    result = capd::max(result,capd::abs(this->coefficient(i,r)));
  }
  if(degree)
  {
    const typename CurveT::JetType& c = this->getCoefficients()[r];

    for(i=0;i<this->dimension();++i)
    {
      typename CurveT::JetType::const_iterator b = c.begin(i,1), e=c.end(i,degree);
      while(b!=e)
        result = capd::max(result,capd::abs(*b++));
    }
  }

  return result;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT,StepControlT,CurveT>::evalAndSum(VectorType& v,MatrixType& der, HessianType& hessian)
{
  size_type i,j,c;

  int r=this->getOrder();
  this->m_vField->computeODECoefficients(this->getCoefficients(),2,r);
  this->computeTimeStep(v);
  this->sum(v);

  hessian.clear();
  const bool maskDissabled = (this->m_vField->getMask()==0);
  for(j=0;j<this->dimension();++j){
    if(maskDissabled or this->getMask(j)){
      for(i=0;i<this->dimension();++i){
        der(i+1,j+1) = this->coefficient(i,j,this->getOrder());
        for(int r=this->getOrder()-1;r>=0;--r)
          der(i+1,j+1) = der(i+1,j+1)*this->m_step + this->coefficient(i,j,r);
      }
      for(c=j;c<this->dimension();++c){
        if(maskDissabled or this->getMask(j,c)){
          for(i=0;i<this->dimension();++i){
            hessian(i,j,c) = this->coefficient(i,j,c,this->getOrder());
            for(int r=this->getOrder()-1;r>=0;--r)
              hessian(i,j,c) = hessian(i,j,c)*this->m_step + this->coefficient(i,j,c,r);
          }
        }
      }
    }
    else
      der.column(j).clear();
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT,StepControlT,CurveT>::evalAndSum(VectorType& v,MatrixType& der)
{
  size_type i,j;

  this->m_vField->computeODECoefficients(this->getCoefficients(),1,this->getOrder());
  this->computeTimeStep(v);
  this->sum(v);

  const bool maskDissabled = (this->m_vField->getMask()==0);
  for(j=0;j<this->dimension();++j){
    if(maskDissabled or this->getMask(j)){
      for(i=0;i<this->dimension();++i){
        der(i+1,j+1) = this->coefficient(i,j,this->getOrder());
        for(int r=this->getOrder()-1;r>=0;--r)
          der(i+1,j+1) = der(i+1,j+1)*this->m_step + this->coefficient(i,j,r);
      }
    }
    else
      der.column(j).clear();
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT,StepControlT,CurveT>::sum(VectorType& v)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    int r = this->getOrder();
    v[i] = this->centerCoefficient(i,r) = this->coefficient(i,r);
    for(--r;r>=0;--r)
    {
      this->centerCoefficient(i,r) = this->coefficient(i,r);
      v[i] = v[i]*this->m_step + this->coefficient(i,r);
    }
  }
}
// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT,StepControlT,CurveT>::evalAndSum(VectorType& v)
{
  this->m_vField->computeODECoefficients(this->getCoefficientsAtCenter(),this->getOrder());
  this->computeTimeStep(v);
  capd::vectalg::evalPolynomial(this->getCoefficientsAtCenter(),this->m_step,v,this->getOrder());
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT,StepControlT,CurveT>::evalAndSum(JetType& v)
{
  this->m_vField->computeODECoefficients(this->getCoefficients(),v.degree(),this->getOrder());
  this->computeTimeStep(v());

  const size_type vs = binomial(v.dimension()+v.degree(),v.degree());
  const size_type s = binomial(this->dimension()+this->degree(),this->degree());
  const bool* m = this->getMask();

  if(m==0){
    for(size_type j=0;j<vs;++j){
      for(size_type i=0;i<this->dimension();++i){
        v[i*vs+j] = this->getCoefficients()[this->getOrder()][i*s+j];
        for(int r=this->getOrder()-1;r>=0;--r)
          v[i*vs+j] = v[i*vs+j]*this->m_step + this->getCoefficients()[r][i*s+j];
      }
    }
  } else {
    for(size_type j=0;j<vs;++j,m+=this->getOrder()+2){
      if(*m){
        for(size_type i=0;i<this->dimension();++i){
          v[i*vs+j] = this->getCoefficients()[this->getOrder()][i*s+j];
          for(int r=this->getOrder()-1;r>=0;--r)
            v[i*vs+j] = v[i*vs+j]*this->m_step + this->getCoefficients()[r][i*s+j];

        }
      } else {
        for(size_type i=0;i<this->dimension();++i)
          v[i*vs+j] = TypeTraits<ScalarType>::zero();
      }
    }
  }
  for(size_type r=0;r<=this->getOrder();++r)
    this->getCoefficientsAtCenter()[r] = this->getCoefficients()[r]();
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT,StepControlT,CurveT>::setInitialCondition(const JetType& v)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    typename JetType::const_iterator b = v.begin(i), e=v.end(i);
    ScalarType* p = &this->coefficient(i,0);
    while(b!=e)
    {
      *p = *b;
      b++;
      p++;
    }
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
const typename BasicCnOdeSolver<MapType, StepControlType, CurveT>::SolutionCurve&
BasicCnOdeSolver<MapType, StepControlType, CurveT>::getCurve()
{
  this->setDomain(0.,rightBound(this->m_step));
  return *this;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT, StepControlT,CurveT>::setInitialCondition(const VectorType& x, const MatrixType& M)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    this->coefficient(i,0) = x[i];
    for(size_type j=0;j<this->dimension();++j)
      this->coefficient(i,j,0) = M(i+1,j+1);
  }
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
void BasicCnOdeSolver<MapT, StepControlT,CurveT>::setInitialCondition(const VectorType& x, const MatrixType& M, const HessianType& H)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    this->coefficient(i,0) = x[i];
    for(size_type j=0;j<this->dimension();++j){
      this->coefficient(i,j,0) = M(i+1,j+1);
      for(size_type c=j;c<this->dimension();++c)
        this->coefficient(i,j,c,0) = H(i,j,c);
    }
  }
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICCNODESOLVER_HPP_


