/////////////////////////////////////////////////////////////////////////////
/// @file Curve.hpp
///
/// @author Daniel Wilczak, Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef CAPD_DIFFALGEBRA_CURVE_HPP
#define CAPD_DIFFALGEBRA_CURVE_HPP
#include <sstream>

#include "capd/diffAlgebra/Curve.h"
#include "capd/vectalg/iobject.hpp"

namespace capd{
namespace diffAlgebra{

template <class BaseCurveT, bool isInterval >
Curve<BaseCurveT,isInterval>::Curve(Real left, Real right, size_type dimension, size_type order,size_type degree)
  : BaseCurveT(dimension,order,degree),
    ParametricCurve<typename BaseCurveT::MatrixType>(left,right)
{}

template <class BaseCurveT>
Curve<BaseCurveT,true>::Curve(Real left, Real right, size_type dimension, size_type order, size_type degree)
  : BaseCurveT(dimension,order,degree),
    ParametricCurve<typename BaseCurveT::MatrixType>(left,right),
    initMatrix(MatrixType::Identity(dimension))
{}

// -------------------------------------------------------------------------

template <class BaseCurveT , bool isInterval>
typename Curve<BaseCurveT,isInterval>::VectorType
Curve<BaseCurveT,isInterval>::operator()(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    VectorType phi(this->dimension());
    for(size_type d=0;d<this->dimension();++d){
      ScalarType x = this->centerCoefficient(d,this->getOrder());
      for(int i=this->getOrder()-1;i>=0;--i)
          x = x*h + this->centerCoefficient(d,i);
      phi[d] = x;
    }
    return phi;
  }
  throw this->domainErrorMessage("Curve::operator(h)",h,this->m_left,this->m_right);
}  // operator()

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::VectorType Curve<BaseCurveT,true>::operator()(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    VectorType deltaX(this->dimension());
    VectorType phi(this->dimension());
    VectorType xx(this->dimension());
    VectorType rem(this->dimension());
    MatrixType jacPhi(this->dimension(),this->dimension());

    ScalarType c = power(h,this->getOrder()+1);
    for(size_type d=0;d<this->dimension();++d)
    {
      phi[d] = this->centerCoefficient(d,this->getOrder());
      xx[d] = this->coefficient(d,this->getOrder());
      rem[d] = c * this->remainderCoefficient(d,this->getOrder()+1);
      deltaX[d] = this->coefficient(d,0) - this->centerCoefficient(d,0);
      for(size_type k=0;k<this->dimension();++k){
        jacPhi(d+1,k+1) = this->coefficient(d,k,this->getOrder());
      }
    }

    // sum Taylor series
    for(int i=this->getOrder()-1;i>=0;--i)
    {
      for(size_type d=0;d<this->dimension();++d)
      {
        xx[d] = xx[d]*h + this->coefficient(d,i);
        phi[d] = phi[d]*h + this->centerCoefficient(d,i);
        for(size_type k=0;k<this->dimension();++k){
          jacPhi(d+1,k+1) = jacPhi(d+1,k+1)*h + this->coefficient(d,k,i);
        }
      }
    }
    if(!intersection(xx,phi + jacPhi*deltaX,xx))
      throw std::runtime_error("capd::diffAlgebra::Curve::operator() error - empty intersection. Report this error to CAPD developers!");
    return xx + rem;
  }

  throw this->domainErrorMessage("Curve::operator(h)",h,this->m_left,this->m_right);
}  // operator()

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::VectorType Curve<BaseCurveT,true>::remainder(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    VectorType rem(this->dimension());

    ScalarType c = power(h,this->getOrder()+1);
    for(size_type d=0;d<this->dimension();++d)
      rem[d] = c * this->remainderCoefficient(d,this->getOrder()+1);
    return rem;
  }

  throw this->domainErrorMessage("Curve::remainder(h)",h,this->m_left,this->m_right);
}  // remainder

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::VectorType Curve<BaseCurveT,true>::valueAtCenter(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    VectorType phi(this->dimension());
    for(size_type d=0;d<this->dimension();++d)
      phi[d] = this->centerCoefficient(d,this->getOrder());

    // sum Taylor series
    for(int i=this->getOrder()-1;i>=0;--i)
      for(size_type d=0;d<this->dimension();++d)
        phi[d] = phi[d]*h + this->centerCoefficient(d,i);
    return phi;
  }

  throw this->domainErrorMessage("Curve::valueAtCenter(h)",h,this->m_left,this->m_right);
}  // valueAtCenter

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::VectorType Curve<BaseCurveT,true>::getCenter() const
{
  VectorType center(this->dimension());
  for(size_type d=0;d<this->dimension();++d)
    center[d] = this->centerCoefficient(d,0);

  return center;
}  // getCenter

// -------------------------------------------------------------------------

template<class BaseCurveT, bool isInterval>
typename Curve<BaseCurveT,isInterval>::VectorType
Curve<BaseCurveT,isInterval>::timeDerivative(const ScalarType& h) const
{
  size_type i = this->getOrder();
  if((h>=this->m_left) and (h<=this->m_right))
  {
    VectorType phi(this->dimension());
    for(size_type d=0;d<this->dimension();++d)
    {
      ScalarType x = Real(i)*this->centerCoefficient(d,i);
      for(size_type i=this->getOrder()-1;i>=1;--i)
        x = x*h + Real(i)*this->centerCoefficient(d,i);
      phi[d] = x;
    }
    return phi;
  }

  throw this->domainErrorMessage("Curve::timeDerivative(h)",h,this->m_left,this->m_right);
}  // operator()

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::VectorType
Curve<BaseCurveT,true>::timeDerivative(const ScalarType& h) const
{
  size_type i = this->getOrder();
  if((h>=this->m_left) and (h<=this->m_right))
  {
    VectorType deltaX(this->dimension());
    VectorType phi(this->dimension());
    VectorType rem(this->dimension());
    MatrixType jacPhi(this->dimension(),this->dimension());

    ScalarType c = power(h,i);
    for(size_type d=0;d<this->dimension();++d)
    {
      phi[d] = Real(i)*this->centerCoefficient(d,i);
      rem[d] = Real(i+1) * c * this->remainderCoefficient(d,i+1);
      deltaX[d] = this->coefficient(d,1) - this->centerCoefficient(d,1);
      for(size_type k=0;k<this->dimension();++k)
        jacPhi(d+1,k+1) = Real(i)*this->coefficient(d,k,i);
    }

    // sum Taylor series
    for(size_type i=this->getOrder()-1;i>=1;--i)
    {
      for(size_type d=0;d<this->dimension();++d)
      {
        phi[d] = phi[d]*h + Real(i)*this->centerCoefficient(d,i);
        for(size_type k=0;k<this->dimension();++k)
          jacPhi(d+1,k+1) = jacPhi(d+1,k+1)*h + Real(i)*this->coefficient(d,k,i);
      }
    }

    return (phi + jacPhi*deltaX + rem);
  }

  throw this->domainErrorMessage("Curve::timeDerivative(h)",h,this->m_left,this->m_right);
}  // operator()

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::MatrixType Curve<BaseCurveT,true>::oneStepDerivative(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    const size_type dim = this->dimension();
    MatrixType jacPhi(dim,dim);
    ScalarType c = power(h,this->getOrder()+1);

    for(size_type i=0;i<dim;++i)
      for(size_type j=0;j<dim;++j)
      {
        ScalarType x = this->coefficient(i,j,this->getOrder());
        for(int r=this->getOrder()-1;r>=0;--r)
          x = x*h + this->coefficient(i,j,r);
        jacPhi(i+1,j+1) = x + c*this->remainderCoefficient(i,j,this->getOrder()+1);
      }
    return jacPhi;
  }

  throw this->domainErrorMessage("Curve::oneStepDerivative(h)",h,this->m_left,this->m_right);
}

// -------------------------------------------------------------------------

template<class BaseCurveT>
typename Curve<BaseCurveT,true>::MatrixType Curve<BaseCurveT,true>::oneStepDerivativeOfNumericalMethod(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    const size_type dim = this->dimension();
    MatrixType jacPhi(dim,dim);

    for(size_type i=0;i<dim;++i)
      for(size_type j=0;j<dim;++j)
      {
        ScalarType x = this->coefficient(i,j,this->getOrder());
        for(int r=this->getOrder()-1;r>=0;--r)
          x = x*h + this->coefficient(i,j,r);
        jacPhi(i+1,j+1) = x;
      }
    return jacPhi;
  }

  throw this->domainErrorMessage("Curve::oneStepDerivativeOfNumerivalMethod(h)",h,this->m_left,this->m_right);
}

// -------------------------------------------------------------------------

template<class BaseCurveT, bool isInterval>
typename Curve<BaseCurveT,isInterval>::MatrixType Curve<BaseCurveT,isInterval>::derivative(const ScalarType& h) const
{
  if((h>=this->m_left) and (h<=this->m_right))
  {
    const size_type dim = this->dimension();
    MatrixType jacPhi(dim,dim);

    for(size_type i=0;i<dim;++i)
      for(size_type j=0;j<dim;++j)
      {
        ScalarType x = this->coefficient(i,j,this->getOrder());
        for(int r=this->getOrder()-1;r>=0;--r)
          x = x*h + this->coefficient(i,j,r);
        jacPhi(i+1,j+1) = x;
      }
    return jacPhi;
  }

  throw this->domainErrorMessage("Curve::derivative(h)",h,this->m_left,this->m_right);
}

}} // namespace capd::diffAlgebra



#endif // CAPD_DIFFALGEBRA_CURVE_HPP

