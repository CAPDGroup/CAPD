/////////////////////////////////////////////////////////////////////////////
/// @file SolutionCurve.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DIFFALGEBRA_SOLUTIONCURVE_H_
#define _CAPD_DIFFALGEBRA_SOLUTIONCURVE_H_

#include <stdexcept>
#include <vector>
#include "capd/basicalg/TypeTraits.h"
#include "capd/diffAlgebra/ParametricCurve.h"
#include "capd/vectalg/Vector_Interval.hpp"


namespace capd{
namespace diffAlgebra{
/// @addtogroup diffAlgebra 
/// @{
/**
 * This file defines class that represents parametric curve in \f$R^n\f$.
 * Internally it is represented by piecewise polynomial approximation.
 *
 * In can also store higher order partial derivatives with respect to space variables, i.e.
 *
 * c(t,x) can be a polynomial both with respect to main parameterization variable 't' and space variable 'x'.
 *
 * The curve can be evaluated at given time 't' as well as differentiated with respect to 't'.
 *
 * The curve can store rigorous enclosures of curves. They are represented using mean value form to reduce some wrapping effect.
 *
 * The main template parameter is an abstract class that represent a piece of curve as one polynomial.
 */
template<class CurveT>
class BaseSolutionCurve : public ParametricCurve<typename CurveT::MatrixType>
{
  public:
  typedef CurveT BaseCurve;
  typedef typename BaseCurve::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;

  typedef std::vector<BaseCurve*> CurveContainer;
  typedef std::vector<ScalarType>Domain;

  BaseSolutionCurve(ScalarType leftDomain)
    : ParametricCurve<typename CurveT::MatrixType>(leftBound(leftDomain),rightBound(leftDomain))
  {
    domain.push_back(leftDomain);
  }

  BaseSolutionCurve(const BaseSolutionCurve<CurveT>& solution_curve)
    : ParametricCurve<typename CurveT::MatrixType>(solution_curve)
  {
    this->domain = solution_curve.domain;
    cloneCurveContainer(solution_curve.pieces);
  }

  BaseSolutionCurve& operator= (const BaseSolutionCurve<CurveT>& solution_curve)
  {
    clearCurveContainer();
    this->domain = solution_curve.domain;
    cloneCurveContainer(solution_curve.pieces);
    return *this;
  }

  ~BaseSolutionCurve()
  {
    clearCurveContainer();  
  }

  Real getLeftDomain() const
  {
    return rightBound(domain[0]);
  }

  Real getRightDomain() const{
    if(pieces.size()>0){
      return leftBound(domain[domain.size()-1]);
    }else
      throw std::domain_error("ParametricCurve::getLeftDomain() - domain is empty.");
  }

  virtual void add(const BaseCurve& c){
    domain.push_back(domain[domain.size()-1]+c.getRightDomain());
    pieces.push_back(new BaseCurve(c));
  }

  typename CurveContainer::size_type getNumberOfPieces() const { return pieces.size(); }
protected:

  Domain domain;
  CurveContainer pieces;

private:
  void clearCurveContainer()
  {
    for(typename CurveContainer::size_type i = 0; i < pieces.size(); ++i)
    {
      delete pieces[i];
      pieces[i] = 0;
    }
  }

  void cloneCurveContainer(const CurveContainer& pieces)
  {
    this->pieces.reserve(pieces.size());
    for(typename CurveContainer::size_type i = 0; i < pieces.size(); ++i)
    {
      this->pieces.push_back( new BaseCurve( *(pieces[i]) ) );
    }
  }
};

// ############################################################################

// this is a template version for all non-interval cases

template<class CurveT, bool isInterval = TypeTraits<typename CurveT::ScalarType>::isInterval >
class SolutionCurve : public BaseSolutionCurve<CurveT>
{
public:
  typedef CurveT BaseCurve;
  typedef typename BaseCurve::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef Jet<MatrixType,0> JetType;

  typedef std::vector<BaseCurve*> CurveContainer;
  typedef std::vector<ScalarType>Domain;

  SolutionCurve(ScalarType leftDomain)
    : BaseSolutionCurve<CurveT>(leftDomain)
  {}

  unsigned checkAndFind(const ScalarType& h, ScalarType& t) const{
    if(!(h>=this->getLeftDomain() and h<=this->getRightDomain() and this->pieces.size()>0))
        throw std::domain_error("SolutionCurve::operator()(ScalarType) - argument is out of the domain.");

    // binary search of subdomain
    int i=0, j = static_cast<int>(this->domain.size())-1;
    int k;
    do{
      k=(i+j)/2;
      if(h>=this->domain[k])
        i=k;
      if(h<=this->domain[k])
        j=k;
    }while(i<j-1);
    t = capd::max(TypeTraits<ScalarType>::zero(),h-this->domain[i]);
    t = capd::min(t,this->pieces[i]->getRightDomain());
    return i;
  }

  VectorType operator()(const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->operator()(t);
  }

  VectorType timeDerivative(const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->timeDerivative(t);
  }

  MatrixType operator[](const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->operator[](t);
  }

  MatrixType derivative(const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->derivative(t);
  }

  HessianType hessian(const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->hessian(t);
  }
  JetType jet(const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->jet(t);
  }
  void eval(ScalarType h, JetType& v) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    this->pieces[i]->eval(t,v);
  }
  /*
  VectorType valueAtCenter(const ScalarType& h) const{
    ScalarType t;
    unsigned i = checkAndFind(h,t);
    return this->pieces[i]->valueAtCenter(t);
  }
*/
};

// ############## specialization for intervals ########################

// this is a template version for all interval cases

template<class CurveT>
class SolutionCurve<CurveT,true> : public BaseSolutionCurve<CurveT>
{
public:
  typedef CurveT BaseCurve;
  typedef typename BaseCurve::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename TypeTraits<ScalarType>::Real Real;
  typedef Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;

  typedef std::vector<BaseCurve*> CurveContainer;
  typedef std::vector<ScalarType>Domain;

  SolutionCurve(ScalarType leftDomain)
    : BaseSolutionCurve<CurveT>(leftDomain)
  {}

  void checkAndFind(const ScalarType& h, unsigned& leftI, unsigned& rightI) const{
    Real left = leftBound(h);
    Real right = rightBound(h);
    if(!(left>=this->getLeftDomain() and right<=this->getRightDomain() and this->pieces.size()>0))
        throw std::domain_error("SolutionCurve::operator()(ScalarType) - argument is out of the domain.");

    // TODO
    // try to implement binary search
    leftI=0;
    while(this->domain[leftI+1]<left) leftI++;
    rightI = leftI;
    while(this->domain[rightI]<right) rightI++;
  }

  template<class ResultType, class MethodPointer>
  ResultType eval(const ScalarType& h,MethodPointer p) const {
    unsigned leftI,rightI;
    checkAndFind(h,leftI,rightI);

    Real left = leftBound(h);
    Real right = rightBound(h);

    Real l = capd::max(TypeTraits<Real>::zero(),leftBound(left-this->domain[leftI]));
    Real r = capd::min(rightBound(right-this->domain[leftI]),this->pieces[leftI]->getRightDomain());

    // evaluate at first subdomain
    ScalarType t(l,r);
    ResultType result = (this->pieces[leftI]->*p)(t);

    while(++leftI<rightI){
      r = capd::min(rightBound(right-this->domain[leftI]),this->pieces[leftI]->getRightDomain());
      ScalarType t(TypeTraits<Real>::zero(),r);
      intervalHull(result,(this->pieces[leftI]->*p)(t),result);
    }
    return result;
  }

  VectorType timeDerivative(const ScalarType& h) const{
    return this->eval<VectorType>(h,&BaseCurve::timeDerivative);
  }

  VectorType operator()(const ScalarType& h) const{
    return this->eval<VectorType>(h,&BaseCurve::operator());
  }

  /// @deprecated
  MatrixType operator[](const ScalarType& h) const{
    return this->eval<MatrixType>(h,&BaseCurve::operator[]);
  }

  MatrixType derivative(const ScalarType& h) const{
    return this->eval<MatrixType>(h,&BaseCurve::derivative);
  }

  HessianType hessian(const ScalarType& h) const{
    return this->eval<HessianType>(h,&BaseCurve::hessian);
  }
};

/// @}
}} // namespace capd::diffAlgebra

#endif
