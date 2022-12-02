/////////////////////////////////////////////////////////////////////////////
/// @file  PoincareMap_templateOperator.h
///
/// @author Daniel Wilczak
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_OPERATOR_H_
#define _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_OPERATOR_H_

#include <cassert>
#include "capd/poincare/PoincareMap.h"
#include "capd/poincare/BasicPoincareMap.hpp"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::computePoincareMap(T& theSet, int n)
{
  T setAfterSection = theSet;
  this->integrateUntilSectionCrossing(theSet,setAfterSection,n);
  VectorType bound = (VectorType)theSet;
  ScalarType oneStepReturnTime;

  // try to cross section in one step and use Newton method
  if(this->crossSectionInOneStep(theSet,setAfterSection,oneStepReturnTime,bound)){
    this->sectionDerivativesEnclosure.computeOneStepSectionEnclosure(theSet,this->m_solver,bound,oneStepReturnTime);
    theSet = setAfterSection;
    return bound;
  }
  // if not possible, try to cross section by the old method.
  return this->crossSection(theSet, setAfterSection);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& theSet, const VectorType& c, const MatrixType& A, ScalarType& out_returnTime, int n)
{
  // The originalSet after this function contains a set just after section.
  this->sectionDerivativesEnclosure.init(&out_returnTime,nullptr,nullptr,nullptr);

  const size_type D = A.numberOfRows();
  T setAfterSection = theSet;
  this->integrateUntilSectionCrossing(theSet,setAfterSection,n);
  VectorType bound = (VectorType)setAfterSection;
  ScalarType oneStepReturnTime;
  // try to cross section in one step and use Newton method to resolve for return time
  if(this->crossSectionInOneStep(theSet,setAfterSection,oneStepReturnTime,bound)){
    SaveStepControl<Solver> ssc(this->m_solver);
/*  Here we use the following estimates
    t - bound for return time
    A - const matrix
    c - const vector
    coordinates on section are given by A*(x(t)-c)

    t0 := mid(t)
    dt := t-t0

    x0 := mid(x(t0))
    dx := x(t0) - x0

    Using Taylor expansion of order 2 wrt. time and then mean value form we get:

    A*(x(t)-c) = A*(x(t0)-c)  +  A*f(x(t0))*dt  +  0.5*A*Df(x(t))*f(x(t))*dt^2
               = A*(x(t0)-c)  +  A*f(x0)*dt  +  A*Df(x(t0))*dt*deltaX  +  0.5*A*Df(x(t))*f(x(t))*dt^2

   VERY important:
   - first multiply "thin objects", like [ A*f(x0) ] * dt
   - evaluation of A*(x(t0)-c) must take into account representation of x(t0).
     This is hidden in the class that represents the set - usually this is doubleton or tripleton.
*/
    MatrixType M(D,D);
    *(this->sectionDerivativesEnclosure.returnTime) = theSet.getCurrentTime() + oneStepReturnTime;
    ScalarType t0, dt;
    oneStepReturnTime.split(t0,dt);

    // estimate [ A*Df(x(t)) ] * [ f(x(t)) * (0.5*dt^2) ]
    VectorType fx = this->getVectorField()(*(this->sectionDerivativesEnclosure.returnTime),bound,M);
    VectorType result = (A*M)*(fx*(typename capd::TypeTraits<ScalarType>::Real(0.5)*sqr(dt)));
    this->m_solver.setStep(t0);
    this->m_solver(theSet);
    VectorType X = (VectorType)(theSet);
    //VectorType x0(D), dx(D);
    VectorType x0 = X, dx=X;
    split(X,x0,dx);

    // add [ A*Df(x(t0))] *  [ dt*dx ]
    M = A*this->getVectorField().derivative(theSet.getCurrentTime(),X);
    VectorType second = M*dx;

    // add [ A*f(x0) ] * dt. VERY important: first multiply "thin objects", i.e. A and f(x0).
    VectorType first = A*this->getVectorField()(theSet.getCurrentTime(),x0);

    // intersect bounds and add to result
    VectorType Y = intersection(A*this->getVectorField()(theSet.getCurrentTime(),X),first+second);
    result += Y*dt;

    // eventually add A*(x(t0)-c) taking into account representation of the set. This is hidden in the class that represents the set.
    //VectorType Y0 = originalSet.affineTransformation(A,c);
    result += theSet.affineTransformation(A,c);
    // reset the set
    theSet = setAfterSection;

    // intersect with naive estimation
    result = intersection(A*(bound-c),result);
    return result;
  }

  // if not possible, try to cross section by the old method.
  return A*(this->crossSection(theSet, setAfterSection)-c);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, ScalarType& out_returnTime, int n)
{
  // The originalSet after this function contains a set just after section.
  this->sectionDerivativesEnclosure.init(&out_returnTime,nullptr,nullptr,nullptr);
  return this->computePoincareMap(originalSet,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,returnTime,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, ScalarType& out_returnTime, int n)
{
  this->sectionDerivativesEnclosure.init(&out_returnTime,&der,nullptr,nullptr);
  return this->computePoincareMap(originalSet,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,der,returnTime,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, HessianType& hessian, ScalarType& out_returnTime, int n)
{
  this->sectionDerivativesEnclosure.init(&out_returnTime,&der,&hessian,nullptr);
  return this->computePoincareMap(originalSet,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, MatrixType& der, HessianType& hessian, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,der,hessian,returnTime,n);
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, typename T::JetType& jet, ScalarType& out_returnTime, int n)
{
  this->sectionDerivativesEnclosure.init(&out_returnTime,nullptr,nullptr,&jet);
  VectorType r = this->computePoincareMap(originalSet,n);
  jet() = r;
  return r;
}

/*__________________________________________________________________________*/

template <typename DS, typename FunT>
template<typename T>
typename PoincareMap<DS,FunT>::VectorType
PoincareMap<DS,FunT>::operator()(T& originalSet, typename T::JetType& result, int n)
{
  ScalarType returnTime = TypeTraits<ScalarType>::zero();
  return (*this)(originalSet,result,returnTime,n);
}

/// @}
}} // namespace capd::poincare

#endif // _CAPD_POINCARE_POINCARE_MAP_TEMPLATE_OPERATOR_H_

