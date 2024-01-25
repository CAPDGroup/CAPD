

/////////////////////////////////////////////////////////////////////////////
/// @file BasicOdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICSOLVER_H_
#define _CAPD_DYNSYS_BASICSOLVER_H_

#include <string>
#include "capd/dynsys/StepControl.h"
#include "capd/diffAlgebra/BasicCurve.h"
#include "capd/diffAlgebra/Curve.h"
#include "capd/dynsys/Move.h"

namespace capd {
namespace dynsys {
/// @addtogroup dynsys
/// @{

/** MapT constraints:
 *  type definitions:
 *   - FunctionType
 *   - MatrixType
 *   - MatrixType::RowVectorType
 *   - MatrixType::ScalarType
 *
 *  methods:
 *   - dimension            - dimension of the space
 *   - getOrder, setOrder   - order of Taylor expansion
 *   - operator()(a, i)     - computes coeffs a_i from a = a_{i-1}
 *   - variational(x_coeff, m_F, m_order);
 *
 */

template <
  typename MapT,
  typename StepControlT = capd::dynsys::DLastTermsStepControl,
  typename CurveT = capd::diffAlgebra::Curve< capd::diffAlgebra::BasicCurve<typename MapT::MatrixType> >
>
class BasicOdeSolver :  public capd::dynsys::StepControlInterface<StepControlT,typename MapT::ScalarType>, public CurveT {
public:
  typedef MapT VectorFieldType;
  typedef StepControlT StepControlType;
  typedef typename VectorFieldType::FunctionType FunctionType;
  typedef typename VectorFieldType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef CurveT SolutionCurve;

  typedef typename VectorFieldType::JetType JetType;
  typedef typename JetType::Multipointer Multipointer;
  typedef typename JetType::Multiindex Multiindex;

  BasicOdeSolver(VectorFieldType& field, size_type order, const StepControlT& stepControl = StepControlT());

  VectorType operator()(VectorType);                                              ///< Computes image of vector v after one time step.
  VectorType operator()(ScalarType& t, const VectorType& u);                      ///< Computes image of vector v after one time step. The argument t is updated in this procedure.

  VectorType operator()(VectorType, MatrixType& o_resultDerivative) ;  ///< Computes image of vector v and derivatives of the flow with respect to init condition (v,identity). Version for autonomous systems.
  VectorType operator()(ScalarType& t, const VectorType&, MatrixType& o_resultDerivative) ;  ///< Computes image of vector v and derivatives of the flow with respect to init condition (v,identity). Version for nonautonomous systems. The argument t is updated in this procedure.

  VectorType operator()(VectorType, const MatrixType& derivative, MatrixType& o_resultDerivative) ; ///< Computes image of vector v and derivatives of a flow with respect to init condition (v, derivative)
  VectorType operator()(ScalarType& t, const VectorType& v, const MatrixType& derivative, MatrixType& o_resultDerivative) ; ///< Computes image of vector v and derivatives of a flow with respect to init condition (v, derivative). The argument t is updated in this procedure.

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C1JetMove is used meaning that this solver can integrate C^0 and C^1 jets only.
  template <typename JetT>
  void operator()(JetT& jet){
	  C1JetMove<BasicOdeSolver,JetT>::move(jet,*this);
  }

  const VectorFieldType& getVectorField() const; ///< Returns vector field
  VectorFieldType& getVectorField();

  void setOrder(size_type order); ///< Sets the order of the Taylor method

  ScalarType getStep() const; ///< Returns the time step made in the last call to this solver
  void setStep(const ScalarType& newStep); ///< Sets fixed time step and turns off step control

  // the following methods provide an interface for generic algorithms based on an abstract solver
  void computeCoefficientsAtCenter(const VectorType& x, size_type order);
  void computeCoefficientsAtCenter(ScalarType t, const VectorType& x, size_type order);
  void computeCoefficients(const VectorType& x, size_type order);
  void computeCoefficients(ScalarType t, const VectorType& x, size_type order);
  void computeCoefficients(const VectorType& x, const MatrixType& M, size_type order);
  void computeCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, size_type order);

  const SolutionCurve& getCurve();
  virtual ScalarType getCoeffNorm(size_type i, size_type degree) const;
  size_type dimension() const {
    return m_vField->dimension();
  }

  VectorType enclosure(const ScalarType& /*t*/, const VectorType& /*x*/){
    throw std::logic_error("OdeSolver::enclosure - cannot compute enclosure, this is a nonrigorous solver. Implementation provided to satisfy required interface of StepControl algorithms");
  }
  void adjustTimeStep(const ScalarType& newStep); ///< sets time step but does not change step control settings (compare setStep)

  /**
   * The iterator range [b,e) should contain a range of Multiinideces the user requires to compute.
   * The method automatically adds all the depending partial derivatives to this collection and defines
   * a mask for computation of partial derivtives.
   * @param [b,e) iterator range which contains collection of multiindices
   * @warning The method causes undefined behavior if a multiindex in the collection exceeds limits of the vector field (like dimension, maximal allowed degree).
   */
  template<class Iterator>
  void setMask(Iterator b, Iterator e){
    this->m_vField->setMask(b,e);
  }
  /**
   * Adds new multiindex (along with dependencies) to the \b existing mask.
   * @param mi multiindex to be added to the mask
   * @warning causes undefined behavior if the mask has not been set before call to this method.
   */
  void addMultiindexToMask(const Multiindex& mi){
    this->m_vField->addMultiindexToMask(mi);
  }
  /**
   * Resets the mask of derivatives. In consequence, full jet of derivatives will be computed after call to any method that computes derivative, hessian or jet.
   */
  void resetMask(){
    this->m_vField->resetMask();
  }

 protected:
  void setCurrentTime(const ScalarType& a_time) const { m_vField->setCurrentTime(a_time); }
  const ScalarType& getCurrentTime() const { return m_vField->getCurrentTime(); }
  const bool* getMask() const { return this->m_vField->getMask(); }
  bool getMask(size_type j) const { return this->m_vField->getMask(j); }
  void computeTimeStep(VectorType& v);
  void sumTaylorSeries(
      VectorType& v,
      VectorType* coeff,
      size_type order
    );

  void sumTaylorSeries(
      VectorType& v,
      MatrixType& der,
      VectorType* coeff,
      MatrixType* matrixCoeff,
      size_type order
    );
  void operator=(const BasicOdeSolver&) {
  } /// we do not allow copying of objects
  BasicOdeSolver(const BasicOdeSolver & solver)
    : capd::dynsys::StepControlInterface<StepControlT,ScalarType>(solver.getStepControl()),
      SolutionCurve(0.,0.,1,1,1)
  { } /// we do not allow copying of objects

  VectorFieldType* m_vField;
  ScalarType m_fixedTimeStep;
  ScalarType m_step;
};

// -----------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
inline const MapType& BasicOdeSolver<MapType, StepControlType, CurveT>::getVectorField() const {
  return *m_vField;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline MapType& BasicOdeSolver<MapType, StepControlType, CurveT>::getVectorField() {
  return *m_vField;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline typename BasicOdeSolver<MapType, StepControlType, CurveT>::ScalarType
BasicOdeSolver<MapType, StepControlType, CurveT>::getStep() const {
  return m_step;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline void BasicOdeSolver<MapType, StepControlType, CurveT>::setStep(const ScalarType& newStep) {
  m_fixedTimeStep = newStep;
  this->turnOffStepControl();
}

template <typename MapType, typename StepControlType, typename CurveT>
inline void BasicOdeSolver<MapType, StepControlType, CurveT>::adjustTimeStep(const ScalarType& newStep) {
  m_step = newStep;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline void BasicOdeSolver<MapType, StepControlType, CurveT>::computeTimeStep(VectorType& v) {
  m_step = this->isStepChangeAllowed()
      ? this->getStepControl().computeNextTimeStep(*this,this->getCurrentTime(),v)
      : capd::min(this->m_fixedTimeStep,this->getMaxStep());
}

//###########################################################//


template <typename MapType, typename StepControlType, typename CurveT>
inline typename BasicOdeSolver<MapType, StepControlType, CurveT>::VectorType
BasicOdeSolver<MapType, StepControlType, CurveT>::operator()(ScalarType& t, const VectorType &iv)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(iv);
  t += this->getStep();
  return r;
}

//###########################################################//


template <typename MapType, typename StepControlType, typename CurveT>
inline typename BasicOdeSolver<MapType, StepControlType, CurveT>::VectorType
BasicOdeSolver<MapType, StepControlType, CurveT>::operator()(ScalarType& t, const VectorType &v, MatrixType& der)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(v,der);
  t += this->getStep();
  return r;
}

//###########################################################//

template <typename MapType, typename StepControlType, typename CurveT>
inline typename BasicOdeSolver<MapType, StepControlType, CurveT>::VectorType
BasicOdeSolver<MapType, StepControlType, CurveT>::operator()(ScalarType& t,  const VectorType &v, const MatrixType& derivative, MatrixType & resultDerivative)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(v,derivative,resultDerivative);
  t += this->getStep();
  return r;
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
inline void BasicOdeSolver<MapType,StepControlType,CurveType>::computeCoefficientsAtCenter(const VectorType& x, size_type order)
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  coeff[0] = x;
  this->m_vField->computeODECoefficients(coeff,order);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
inline void BasicOdeSolver<MapType,StepControlType,CurveType>::computeCoefficientsAtCenter(ScalarType t, const VectorType& x, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficientsAtCenter(x,order);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
inline void BasicOdeSolver<MapType,StepControlType,CurveType>::computeCoefficients(const VectorType& x, size_type order)
{
  VectorType* coeff = this->getCoefficients();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  coeff[0] = x;
  matrixCoeff[0].setToIdentity();
  this->m_vField->computeODECoefficients(coeff,matrixCoeff,order);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
inline void BasicOdeSolver<MapType,StepControlType,CurveType>::computeCoefficients(ScalarType t, const VectorType& x, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficients(x,order);
}

//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
inline void BasicOdeSolver<MapType,StepControlType,CurveType>::computeCoefficients(const VectorType& x, const MatrixType& M, size_type order)
{
  VectorType* coeff = this->getCoefficients();
  MatrixType* matrixCoeff = this->getMatrixCoefficients();
  coeff[0] = x;
  matrixCoeff[0] = M;
  this->m_vField->computeODECoefficients(coeff,matrixCoeff,order);
}


//###########################################################//

template<typename MapType, typename StepControlType,typename CurveType>
inline void BasicOdeSolver<MapType,StepControlType,CurveType>::computeCoefficients(ScalarType t,const VectorType& x, const MatrixType& M, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficients(x,M,order);
}

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
inline void BasicOdeSolver<MapType,StepControlType,CurveT>::sumTaylorSeries(
      VectorType& v,
      VectorType* coeff,
      size_type order
  )
{
  capd::vectalg::evalPolynomial(coeff,this->m_step,v,order);
}
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICODESOLVER_H_

