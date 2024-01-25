/////////////////////////////////////////////////////////////////////////////
///
/// @file C0HODoubletonSetGeometricTail.h
///
/// @author Daniel Wilczak
/// Created on: Dec 28, 2011
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2011 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0HODOUBLETONSETGEOMETRICTAIL_H_
#define _CAPD_DYNSET_C0HODOUBLETONSETGEOMETRICTAIL_H_

#include "capd/basicalg/factrial.h"
#include "capd/dynset/AbstractSet.h"
#include "capd/dynset/C0HOSet.h"

namespace capd{
namespace pdes{
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
class C0HODoubletonSetGeometricTail : public capd::dynset::C0HOSet<BaseSetT>{
public:
  typedef capd::pdes::GeometricBound<capd::interval> VectorType;
  typedef typename BaseSetT::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType FiniteVectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;
  typedef capd::dynset::C0HOSet<BaseSetT> BaseSet;
  typedef typename capd::dynset::C0Set<MatrixType>::SetType SetType;
  typedef capd::pdes::PdeSolver<VectorType> DynSysType;
  typedef typename BaseSet::Data  Data;

  C0HODoubletonSetGeometricTail(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero())
    : BaseSet(x.getExplicitCoefficients(),t),
      m_currentSeries(x),
      m_enc(x)
  {}
  C0HODoubletonSetGeometricTail(const VectorType& x, const FiniteVectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero())
    : BaseSet(x.projection(r0.dimension()),r0,t),
      m_currentSeries(x),
      m_enc(x)
  {
    for(size_type i=0;i<r0.dimension();++i)
      m_currentSeries[i] += r0[i];
  }

  C0HODoubletonSetGeometricTail(const VectorType& x, const MatrixType& C, const FiniteVectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero())
  : BaseSet(x.projection(r0.dimension()),C,r0,t),
    m_currentSeries(x),
    m_enc(x)
  {
    m_currentSeries.projection(r0.dimension()) += C*r0;
  }

  /// computes image of the set after one step/iterate of the dynamical system
  void move(DynSysType & dynsys) { move(dynsys,*this); }

  /// computes image of the set after one step/iterate of the dynamical system and stores it in result
  void move(DynSysType & dynsys, C0HODoubletonSetGeometricTail& result) const;

  /// This method computes value of functor f at interval vector represented by this set.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    return BaseSet::evalAt(f.getProjection(BaseSet::dimension()));
  }

  //using SetType::operator VectorType;
  operator VectorType() const { return m_currentSeries; }

  std::string name() const { return "C0HODoubletonSetGeometricTail"; }

  const VectorType& getCurrentSeries() const { return m_currentSeries; }
  VectorType m_currentSeries, m_enc;
};

// -----------------------------------------------------------------------------------

template<typename Policies>
void C0HODoubletonSetGeometricTail<Policies>::move(DynSysType & solver, C0HODoubletonSetGeometricTail& result) const
{
  size_type order = solver.getOrder();
  size_type q = order/2;
  size_type p = order - q;

  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  if(subset(this->predictor.get_x(),this->m_currentSet)){
      result.x = this->predictor.get_x();
      result.deltaX = this->m_currentSet - this->predictor.get_x();
  }else
    split(this->m_currentSet, result.x, result.deltaX);
  if(&result!=this)
    result.m_currentSet = this->m_currentSet;
  solver.encloseC0Map(result.x,result.m_currentSeries,result.y, result.rem, result.m_enc, result.jacPhi);
  BaseSet::BaseSet::move(this->predictor,result.predictor,result.pBound,result);

  for(size_type i=0;i<result.m_currentSet.dimension();++i){
    ScalarType t = result.m_currentSeries.getCoefficient(i+1);
    intersection(t,result.pBound[i],result.pBound[i]);
    result.m_currentSeries.setCoefficient(i+1,result.pBound[i]);
  }
  // corrector step
  capd::dynset::computePsi(solver.getCurve(),p,q,solver.getStep(),result.psiPlus,result.JPlus);

  // we have to recompute coefficients at the image
  split(result.pBound,result.y,result.deltaY);
  ScalarType nextTime =  this->getCurrentTime() + solver.getStep();
  solver.computeImplicitCoefficients(result.y,result.m_currentSeries,q);
  capd::dynset::computePsi(solver.getImplicitCurve(),q,p,-solver.getStep(),result.psiMinus,result.JMinus);

  // solve implicit equations
  result.computeC0HORemainder(p,q);
  result.computeC0HOCoefficients();

  result.rem =  result.midJMinusInverse*result.rem;
  subtractAssignMatrixByVector(result.rem,result.T,result.deltaY);
  BaseSet::BaseSet::move( this->corrector, result.corrector, result.cBound, result );

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

  for(size_type i=0;i<result.m_currentSet.dimension();++i){
    result.m_currentSeries.setCoefficient(i+1,result.m_currentSet[i]);
  }

}

/// @}
}} // namespace capd::pdes

#endif // _CAPD_DYNSET_C0HODOUBLETONSETGEOMETRICTAIL_H_


