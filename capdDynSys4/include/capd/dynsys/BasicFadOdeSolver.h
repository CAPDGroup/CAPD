

/////////////////////////////////////////////////////////////////////////////
/// @file BasicFadOdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2017 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICFADODESOLVER_H_
#define _CAPD_DYNSYS_BASICFADODESOLVER_H_

#include "capd/fadbad/fadiff.h"
#include "capd/fadbad/differentiate.h"
#include "capd/basicalg/power.h"
#include "capd/dynsys/StepControl.h"
#include "capd/diffAlgebra/FadCurve.h"
#include "capd/diffAlgebra/Curve.h"
#include "capd/dynsys/Move.h"



namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
// a type FadMapT should specify interface defined in class FadMap - see header file
template<class FadMapT, typename StepControlT = capd::dynsys::DLastTermsStepControl>
class BasicFadOdeSolver :
  public capd::dynsys::StepControlInterface<StepControlT,typename FadMapT::ScalarType>,
  public capd::diffAlgebra::Curve< capd::diffAlgebra::FadCurve<typename FadMapT::MatrixType> >
{
public:
  typedef FadMapT VectorFieldType;
  typedef StepControlT StepControlType;
  typedef typename FadMapT::ScalarType ScalarType;
  typedef typename FadMapT::MatrixType MatrixType;
  typedef typename FadMapT::VectorType VectorType;
  typedef typename FadMapT::FunctionType FunctionType;
  typedef typename MatrixType::size_type size_type;

  typedef fadbad::T<ScalarType> TScalar;
  typedef fadbad::F<ScalarType,FadMapT::VectorType::csDim> FScalar;
  typedef fadbad::T<FScalar> TFScalar;
  typedef typename VectorType::template rebind<TFScalar>::other TFVector;
  typedef typename VectorType::template rebind<TScalar>::other TVector;
  typedef typename VectorType::template rebind<FScalar>::other FVector;

  typedef capd::diffAlgebra::Curve< capd::diffAlgebra::FadCurve<typename FadMapT::MatrixType> > SolutionCurve;

  BasicFadOdeSolver(VectorFieldType& f, size_type _order, const StepControlT& stepControl=StepControlT());
  virtual ~BasicFadOdeSolver() {}

  VectorType operator()(VectorType); ///< Computes image of vector v after one time step.
  VectorType operator()(ScalarType& t, VectorType);  ///< Computes image of vector v after one time step. The argument t is updated in this procedure.

  VectorType operator()(VectorType v, MatrixType& o_resultDerivative) ;  ///< Computes image of vector v and derivatives of the flow with respect to init condition (v,identity). Version for autonomous systems.
  VectorType operator()(ScalarType& t, VectorType v, MatrixType& o_resultDerivative) ;  ///< Computes image of vector v and derivatives of the flow with respect to init condition (v,identity). Version for nonautonomous systems. The argument t is updated in this procedure.

  VectorType operator()(VectorType v, const MatrixType& derivative, MatrixType& o_resultDerivative) ; ///< Computes image of vector v and derivatives of a flow with respect to init condition (v, derivative)
  VectorType operator()(ScalarType& t, VectorType v, const MatrixType& derivative, MatrixType& o_resultDerivative) ; ///< Computes image of vector v and derivatives of a flow with respect to init condition (v, derivative). The argument t is updated in this procedure.

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C1JetMove is used meaning that this solver can integrate C^0 and C^1 jets only.
  template <typename JetT>
  void operator()(JetT& jet){
	  C1JetMove<BasicFadOdeSolver,JetT>::move(jet,*this);
  }

  const VectorFieldType& getVectorField() const; ///< Returns vector field
  VectorFieldType& getVectorField();

  ScalarType getStep() const; ///< Returns current time step
  void setStep(const ScalarType& newStep); ///< Sets time step for the next step of integration and turns off step control.

  virtual ScalarType getCoeffNorm(size_type i,size_type degree) const;
  const SolutionCurve& getCurve();

  // template members

  template <typename Description>
  void setParameter(Description name, const ScalarType& value){
    m_vectorField.setParameter(name,value);
    this->recordDags();
  }


  // the following methods provide an interface for generic algorithms based on an abstract solver
  void computeCoefficientsAtCenter(const VectorType& x, size_type order);
  void computeCoefficientsAtCenter(ScalarType t, const VectorType& x, size_type order);
  void computeCoefficients(const VectorType& x, size_type order);
  void computeCoefficients(ScalarType t, const VectorType& x, size_type order);
  void computeCoefficients(const VectorType& x, const MatrixType& M, size_type order);
  void computeCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, size_type order);

  VectorType enclosure(const ScalarType& /*t*/, const VectorType& /*x*/){
    throw std::logic_error("BasicFadTaylor::enclosure - cannot compute enclosure, this is a nonrigorous solver. Implementation only for satisfying an required interface of StepControl");
  }
  void adjustTimeStep(const ScalarType& newStep); ///< sets time step but does not change step control settings (compare setStep)
protected:

  void setCurrentTime(const ScalarType& a_time) const {
    m_time[0] = a_time;
    m_ftime[0].x() = a_time;
  }
  const ScalarType& getCurrentTime() const {
    return m_time[0];
  }

  void computeTimeStep(VectorType& v);
  void recordDags();
  void reset();

  template<class AVector>
  void computeCoeff(AVector& in, AVector& out, size_type order);

  void setInitialCondition(const VectorType& u, TVector& in);
  void setInitialCondition(const VectorType& u, TFVector& in);
  void setInitialCondition(const VectorType& u, const MatrixType& M, TFVector& in);

  void sumTaylorSeries(VectorType& u, TVector& in, size_type order);
  void sumTaylorSeries(MatrixType& M, TFVector& in, size_type order);
  void sumTaylorSeries(VectorType& u, MatrixType& M, TFVector& in, size_type order);

  VectorFieldType& m_vectorField;
  mutable TVector m_centerOut;
  mutable TFVector m_out;
  mutable TVector m_remOut;
  mutable TFVector m_jacRemOut;
  mutable TScalar m_time;
  mutable TFScalar m_ftime;

  ScalarType m_fixedTimeStep;
  ScalarType m_step;
}; // the end of class FadTaylor

  //###########################################################//

template<class FadMapT, typename StepControlT>
template<class AVector>
void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoeff(AVector& in, AVector& out, size_type order)
{
  for(size_type i=0;i<this->dimension();++i)
    out[i].reset();

  for(size_type r=0;r<order;++r)
  {
    for(size_type i=0;i<this->dimension();++i)
    {
      out[i].eval(r);
      in[i][r+1] = out[i][r]/(r+1);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename FadMapT, typename StepControlType>
inline const FadMapT& BasicFadOdeSolver<FadMapT, StepControlType>::getVectorField() const {
  return m_vectorField;
}

template <typename FadMapT, typename StepControlType>
inline FadMapT& BasicFadOdeSolver<FadMapT, StepControlType>::getVectorField() {
  return m_vectorField;
}

template <typename FadMapT, typename StepControlType>
inline typename BasicFadOdeSolver<FadMapT, StepControlType>::ScalarType BasicFadOdeSolver<FadMapT, StepControlType>::getStep() const {
  return m_step;
}

template <typename FadMapT, typename StepControlType>
inline void BasicFadOdeSolver<FadMapT, StepControlType>::setStep(const ScalarType& newStep) {
  m_fixedTimeStep = newStep;
  this->turnOffStepControl();
}

template <typename FadMapT, typename StepControlType>
inline void BasicFadOdeSolver<FadMapT, StepControlType>::adjustTimeStep(const ScalarType& newStep) {
  m_step = newStep;
}

template <typename FadMapT, typename StepControlType>
inline void BasicFadOdeSolver<FadMapT, StepControlType>::computeTimeStep(VectorType& v) {
  m_step = this->isStepChangeAllowed()
      ? this->getStepControl().computeNextTimeStep(*this,this->getCurrentTime(),v)
      : capd::min(this->m_fixedTimeStep,this->getMaxStep());
}

// --------------------------- operators for autonomous systems --------------------------


template<class FadMapT, typename StepControlT>
inline typename BasicFadOdeSolver<FadMapT,StepControlT>::VectorType
BasicFadOdeSolver<FadMapT,StepControlT>::operator()(ScalarType& t, VectorType u)
{
  this->setCurrentTime(t);
  u = (*this)(u);
  t += this->getStep();
  return u;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline typename BasicFadOdeSolver<FadMapT,StepControlT>::VectorType
BasicFadOdeSolver<FadMapT,StepControlT>::operator()(ScalarType& t, VectorType u, const MatrixType& derivative, MatrixType& o_resultDerivative)
{
  this->setCurrentTime(t);
  u = (*this)(u,derivative,o_resultDerivative);
  t += this->getStep();
  return u;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline typename BasicFadOdeSolver<FadMapT,StepControlT>::VectorType
BasicFadOdeSolver<FadMapT,StepControlT>::operator()(ScalarType& t, VectorType u, MatrixType& o_resultDerivative)
{
  this->setCurrentTime(t);
  u = (*this)(u,o_resultDerivative);
  t += this->getStep();
  return u;
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoefficientsAtCenter(const VectorType& x, size_type order)
{
  this->setInitialCondition(x,this->m_center);
  this->computeCoeff(this->m_center,this->m_centerOut,order);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoefficientsAtCenter(ScalarType t, const VectorType& x, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficientsAtCenter(x,order);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoefficients(const VectorType& x, size_type order)
{
  this->setInitialCondition(x,this->m_in);
  this->computeCoeff(this->m_in,this->m_out,order);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoefficients(ScalarType t, const VectorType& x, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficients(x,order);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoefficients(const VectorType& x, const MatrixType& M, size_type order)
{
  this->setInitialCondition(x,M,this->m_in);
  this->computeCoeff(this->m_in,this->m_out,order);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void BasicFadOdeSolver<FadMapT,StepControlT>::computeCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficients(x,M,order);
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICFADODESOLVER_H_


