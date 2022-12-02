/////////////////////////////////////////////////////////////////////////////
/// @file CnCurve.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_CNCURVE_HPP_
#define _CAPD_DIFFALGEBRA_CNCURVE_HPP_

#include <stdexcept>
#include "capd/diffAlgebra/CnCurve.h"
#include "capd/diffAlgebra/Curve.hpp"
#include "capd/diffAlgebra/BasicCnCurve.hpp"
#include "capd/diffAlgebra/Jet.hpp"

namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra
/// @{

template<class CurveT, class HessianT>
void oneStepHessian(typename CurveT::ScalarType h, const CurveT& curve, HessianT& result)
{
  typedef typename CurveT::size_type size_type;
  const size_type dim=curve.dimension();
  for(size_type i=0;i<dim;++i)
    for(size_type j=0;j<dim;++j)
      for(size_type c=j;c<dim;++c){
        int p = curve.getOrder();
        result(i,j,c) = curve.coefficient(i,j,c,p);
        for(--p;p>=0;--p)
          result(i,j,c) = result(i,j,c)*h + curve.coefficient(i,j,c,p);
      }
}

template <class BaseCurveT, bool isInterval >
typename CnCurve<BaseCurveT,isInterval>::HessianType CnCurve<BaseCurveT,isInterval>::hessian(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    HessianType result(this->dimension());
    oneStepHessian(h,*this,result);
    return result;
  }
  throw this->domainErrorMessage("CnCurve::hessian(h)",h,this->m_left,this->m_right);
}

template <class BaseCurveT>
typename CnCurve<BaseCurveT,true>::HessianType CnCurve<BaseCurveT,true>::hessian(const ScalarType& h) const
{
  const size_type dim=this->dimension();
  if((h>=this->m_left) and (h<=this->m_right))
  {
    HessianType H(dim);
    oneStepHessian(h,*this,H);
    ScalarType r = power(h,this->getOrder()+1);
    for(size_type i=0;i<dim;++i)
      for(size_type j=0;j<dim;++j)
        for(size_type c=j;c<dim;++c)
          H(i,j,c) += r*this->remainderCoefficient(i,j,c,this->getOrder()+1);

    return H*this->initMatrix + this->oneStepDerivative(h)*this->initHessian;
  }
  throw this->domainErrorMessage("CnCurve::hessian(h)",h,this->m_left,this->m_right);
}

template<class Curve, class Jet>
void oneStepJet(typename Curve::ScalarType h, Curve& curve, Jet& v)
{
  for(typename Curve::size_type i=0;i<v.dimension();++i)
  {
    const typename Curve::ScalarType* p = &curve.coefficient(i,curve.getOrder());
    typename Jet::iterator b = v.begin(i), e = v.end(i);
    for(;b!=e;++b,++p)
      *b = *p;
    for(int r = curve.getOrder()-1;r>=0;--r)
    {
      const typename Curve::ScalarType* p = &curve.coefficient(i,r);
      typename Jet::iterator b = v.begin(i), e = v.end(i);
      for(;b!=e;++b,++p)
        *b = (*b)*h + (*p);
    }
  }
}

template <class BaseCurveT, bool isInterval >
void CnCurve<BaseCurveT,isInterval>::eval(ScalarType h, JetType& v) const{
  if((h>=this->m_left) and (h<=this->m_right)){
    oneStepJet(h,*this,v);
    return;
  }
  throw this->domainErrorMessage("CnCurve::eval(h,jet)",h,this->m_left,this->m_right);
}

template <class BaseCurveT, bool isInterval >
typename CnCurve<BaseCurveT,isInterval>::JetType CnCurve<BaseCurveT,isInterval>::jet(const ScalarType& h) const{
  JetType jet(this->dimension(),this->degree());
  this->eval(h,jet);
  return jet;
}

template <class BaseCurveT>
void CnCurve<BaseCurveT,true>::eval(ScalarType h, JetType& v) const{
  if((h>=this->m_left) and (h<=this->m_right)){
    const size_type dim = this->dimension();
    v() = this->operator()(h);

    JetType jet(dim,v.degree());
    oneStepJet(h,*this,jet);

    // add bound for Lagrange remainder
    ScalarType c = power(h,this->getOrder()+1);

    // computation of the Taylor coefficients for homogenous polynomials of degree >1
    for(size_type i=0;i<dim;++i)
    {
      const ScalarType* p = &this->remainderCoefficient(i,this->getOrder()+1);
      typename JetType::iterator b = jet.begin(i), e = jet.end(i);
      for(;b!=e;++b,++p)
        *b += c*(*p);
    }

    if(v.degree()>0)
      substitutionPowerSeries(jet,this->initJet,v,false);
  } else
    throw this->domainErrorMessage("CnCurve::eval(h,jet)",h,this->m_left,this->m_right);
}

template <class BaseCurveT>
typename CnCurve<BaseCurveT,true>::JetType CnCurve<BaseCurveT,true>::jet(const ScalarType& h) const{
  JetType jet(this->dimension(),this->degree());
  this->eval(h,jet);
  return jet;
}

///@}
}} // namespace capd::diffAlgebra

#endif // _CAPD_DIFFALGEBRA_CNCURVE_H_
