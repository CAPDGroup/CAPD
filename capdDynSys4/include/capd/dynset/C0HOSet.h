/////////////////////////////////////////////////////////////////////////////
///
/// @file C0HOSet.h
///
/// @author Daniel Wilczak
/// Created on: Dec 28, 2011
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2011 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0HOSET_H_
#define _CAPD_DYNSET_C0HOSET_H_

#include "capd/basicalg/factrial.h"
#include "capd/dynset/AbstractSet.h"
#include "capd/dynset/HOData.h"
#include "capd/dynset/C0Set.h"

namespace capd{
namespace dynset{
// @addtogroup capd
/// @{

/**
 * This class uses representation of subset of R^n inherited from template parameter.
 *
 * The evaluation of the set by an ODE is realized by intersection of two methods:
 * the Taylor method and the Hermite-Obreshkov method.
 *
 * IMPORTANT: present implementation is valid for orders of the Taylor method less or equal 64, only.
 * This is due to capacity of integer type used to store binomial coefficients.
 * Minimal order should be at least 3.
 */
template <class BaseSetT>
class C0HOSet : public C0Set<typename BaseSetT::MatrixType>, public HOData<typename BaseSetT::Data>, public BaseSetT::Policy{

public :
  typedef BaseSetT BaseSet;
  typedef typename BaseSet::BaseSet C0BaseSet;
  typedef typename BaseSet::MatrixType MatrixType;
  typedef typename BaseSet::VectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef typename BaseSet::DynSysType DynSysType;
  typedef typename BaseSet::SetType SetType;
  typedef HOData<typename BaseSetT::Data> Data;

  C0HOSet(const BaseSet & set);
  explicit C0HOSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOSet(const VectorType& x, const MatrixType& C,
      const VectorType& r0, const MatrixType& B,
      const VectorType& r,
      ScalarType t = TypeTraits<ScalarType>::zero()
    );

  /// This method computes value of functor f at interval vector represented by this set.
  /// It computes the value as an intersection of evalAt for two representations of this set: predictor and corrector.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    VectorType gradient = f.gradient(this->m_currentSet);
    ScalarType r = f(this->m_currentSet);
    if(subset(this->predictor.get_x(),this->m_currentSet)){
      intersection( r, f(this->predictor.get_x())+predictor.evalAffineFunctional(gradient,this->predictor.get_x()), r );
    } else{
      VectorType x0 = midVector(this->m_currentSet);
      intersection( r, f(x0) + this->predictor.evalAffineFunctional(gradient,x0), r );
    }
    if(subset(this->corrector.get_x(),this->m_currentSet)){
      intersection( r, f(this->corrector.get_x())+corrector.evalAffineFunctional(gradient,this->corrector.get_x()), r );
    } else{
      VectorType x0 = midVector(this->m_currentSet);
      intersection( r, f(x0) + this->corrector.evalAffineFunctional(gradient,x0), r );
    }
    return r;
  }

  VectorType affineTransformation(const MatrixType& M, const VectorType& x) const {
    return intersection(predictor.affineTransformation(M,x),corrector.affineTransformation(M,x));
  }

  virtual ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& u) const{
    ScalarType r;
    if(!intersection(predictor.evalAffineFunctional(gradient,u),corrector.evalAffineFunctional(gradient,u),r)){
      throw std::logic_error("C0HOSet::evalAffineFunctional - empty intersection. Report this error to CAPD developers!");
    }
    return r;
  }

  template<class Solver>
  void move(Solver& solver){
    this->move(solver,*this);
  }

  template<class Solver>
  void move(Solver& solver, C0HOSet& result) const;
  std::string show() const{
    return BaseSet(predictor).show() + BaseSet(corrector).show();
  }

  void move(DynSysType& /*solver*/){
    throw std::logic_error("C0HOSet: cannot move HO-type sets via abstract C0DynSys.");
  }

  std::string name() const { return "C0HOSet"; }
  C0BaseSet predictor, corrector;
};

// -------------------------------------------


template <class BaseSetT>
template<class Solver>
void C0HOSet<BaseSetT>::move(Solver& solver, C0HOSet& result) const
{
  size_type order = solver.getOrder();
  size_type q = order/2;
  size_type p = order - q;

  // predictor step - Taylor method
  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  if(subset(this->predictor.get_x(),this->m_currentSet)){
      result.x = this->predictor.get_x();
      result.deltaX = this->m_currentSet-this->predictor.get_x();
  }else
    split(this->m_currentSet, result.x, result.deltaX);

  solver.encloseC0Map( this->getCurrentTime(), result.x, this->m_currentSet, result.y, result.rem, result.enc, result.jacPhi );
  BaseSet::move( this->predictor, result.predictor, result.pBound, result );

  // corrector step
  capd::dynset::computePsi(solver.getCurve(),p,q,solver.getStep(),result.psiPlus,result.JPlus);

  // we have to recompute coefficients at the image
  split(result.pBound,result.y,result.deltaY);
  ScalarType nextTime =  this->getCurrentTime() + solver.getStep();
  solver.computeImplicitCoefficients(nextTime,result.y,result.pBound,q);
  capd::dynset::computePsi(solver.getImplicitCurve(),q,p,-solver.getStep(),result.psiMinus,result.JMinus);

  // solve implicit equations
  result.computeC0HORemainder(p,q);
  result.computeC0HOCoefficients();

  result.rem =  result.midJMinusInverse*result.rem;
  subtractAssignMatrixByVector(result.rem,result.T,result.deltaY);
  BaseSet::move( this->corrector, result.corrector, result.cBound, result );

  if(!intersection(result.pBound,result.cBound,result.m_currentSet)){
    throw std::logic_error("C0HOSet: intersection of predictor and corrector is empty! Please report this error to CAPD developers.");
  }

  if(subsetInterior(result.cBound,result.pBound)){
    result.predictor=result.corrector;
  }
  else if(subsetInterior(result.pBound,result.cBound)){
    result.corrector=result.predictor;
  }
  this->reorganizeIfNeeded(result.corrector);
  this->reorganizeIfNeeded(result.predictor);
  result.setLastEnclosure(result.enc);
  result.setCurrentTime(nextTime);
}

/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0HOSET_H_

