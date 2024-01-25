

/////////////////////////////////////////////////////////////////////////////
/// @file BasicCnOdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICCNODESOLVER_H_
#define _CAPD_DYNSYS_BASICCNODESOLVER_H_

#include <utility>
#include <algorithm>
#include "capd/dynsys/StepControl.h"
#include "capd/diffAlgebra/CnCurve.h"
#include "capd/diffAlgebra/BasicCnCurve.h"
#include "capd/dynsys/Move.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template <
  typename MapT,
  typename StepControlT = capd::dynsys::DLastTermsStepControl,
  typename CurveT = capd::diffAlgebra::CnCurve< capd::diffAlgebra::BasicCnCurve< typename MapT::MatrixType > >
  >
class BasicCnOdeSolver :  public capd::dynsys::StepControlInterface<StepControlT,typename MapT::ScalarType>, public CurveT {
public:

  typedef MapT VectorFieldType;
  typedef StepControlT StepControlType;
  typedef typename MapT::FunctionType FunctionType;
  typedef typename VectorFieldType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorFieldType::HessianType HessianType;
  typedef typename VectorFieldType::JetType JetType;
  typedef CurveT SolutionCurve;
  typedef typename JetType::Multipointer Multipointer;
  typedef typename JetType::Multiindex Multiindex;
  typedef typename MatrixType::size_type size_type;

  BasicCnOdeSolver(VectorFieldType& a_vectorField, size_type a_order, const StepControlT& stepControl = StepControlT()); // degree will be read form vectorField object
  virtual ~BasicCnOdeSolver();

  VectorType operator()(VectorType);                        ///< Computes image of vector v after one time step.
  VectorType operator()(ScalarType& t, const VectorType&);  ///< Computes image of vector v after one time step. The argument t is updated in this procedure.

  VectorType operator()(VectorType, MatrixType& o_resultDerivative) ;  ///< Computes image of vector v and derivatives of the flow with respect to init condition (v,identity). Version for autonomous systems.
  VectorType operator()(ScalarType& t, const VectorType&, MatrixType& o_resultDerivative) ;  ///< Computes image of vector v and derivatives of the flow with respect to init condition (v,identity). Version for nonautonomous systems. The argument t is updated in this procedure.

  VectorType operator()(VectorType, const MatrixType& derivative, MatrixType& o_resultDerivative) ; ///< Computes image of vector v and derivatives of a flow with respect to init condition (v, derivative)
  VectorType operator()(ScalarType& t, const VectorType& v, const MatrixType& derivative, MatrixType& o_resultDerivative) ; ///< Computes image of vector v and derivatives of a flow with respect to init condition (v, derivative). The argument t is updated in this procedure.

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization CnJetMove is used meaning that this solver can integrate any type of jets.
  template <typename JetT>
  void operator()(JetT& jet){
	  CnJetMove<BasicCnOdeSolver,JetT>::move(jet,*this);
  }

  /// Computes next point on the trajectory,
  /// first and second order derivatives with respect to initial conditions.
  /// Initial conditions for variational equations are Id and zero, respectively.
  VectorType operator()(VectorType, MatrixType&, HessianType&);

  /// Computes next point on the trajectory of a nonautonomous system,
  /// first and second order derivatives with respect to initial conditions.
  /// Initial conditions for variational equations are Id and zero, respectively.
  VectorType operator()(ScalarType& t, const VectorType&, MatrixType&, HessianType&);

  /// The routine computes next point, derivatives and second order derivatives of a flow.
  /// Initial conditions for variational equations are V and H, respectively.
  VectorType operator()(VectorType, const MatrixType& V, const HessianType& H, MatrixType&, HessianType&);

  /// The routine computes next point, derivatives and second order derivatives of a nonautonomous flow.
  /// Initial conditions for variational equations are V and H, respectively.
  VectorType operator()(ScalarType& t, const VectorType& x, const MatrixType& V, const HessianType& H, MatrixType&, HessianType&);

  /// Initial conditions for the trajectory and for variational equations up to given degree are given in argument in_out_coeffs.
  /// The full result is stored in in_out_coeffs and also the new point on the trajectory is returned.
  /// Note: CnCoeffType is a data structure that stores current time.
  VectorType operator()(ScalarType& t, JetType&);

  // the following methods provide an interface for generic algorithms based on an abstract solver
  void computeCoefficientsAtCenter(const VectorType& x, size_type order);
  void computeCoefficientsAtCenter(ScalarType t, const VectorType& x, size_type order);
  void computeCoefficients(const VectorType& x, size_type order);
  void computeCoefficients(ScalarType t, const VectorType& x, size_type order);
  void computeCoefficients(const VectorType& x, const MatrixType& M, size_type order);
  void computeCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, size_type order);

  const VectorFieldType& getVectorField() const; ///< Returns vector field
  VectorFieldType& getVectorField();

  void setOrder(size_type order); ///< Sets the order of the Taylor method

  ScalarType getStep() const; ///< Returns current time step
  void setStep(const ScalarType& newStep); ///< Sets time step and turns off step control

  const SolutionCurve& getCurve();

// inline definitions
  inline size_type degree() const;
  ScalarType getCoeffNorm(size_type i, size_type degree) const;

  /// TODO
//  void clearCoefficients();

  VectorType enclosure(const ScalarType& /*t*/, const VectorType& /*x*/){
    throw std::logic_error("CnOdeTaylor::enclosure - cannot compute enclosure, this is a nonrigorous solver. Implementation provided to satisfy required interface of StepControl algorithms");
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
    this->clearCoefficients();
  }
  /**
   * Adds new multiindex (along with dependencies) to the \b existing mask.
   * @param mi multiindex to be added to the mask
   * @warning causes undefined behavior if the mask has not been set before call to this method.
   */
  void addMultiindexToMask(const capd::vectalg::Multiindex& mi){
    this->m_vField->addMultiindexToMask(mi);
  }
  /**
   * Resets the mask of derivatives. In consequence, full jet of derivatives will be computed after call to any method that computes derivative, hessian or jet.
   */
  void resetMask(){
    this->m_vField->resetMask();
    this->clearCoefficients();
  }

protected:
  void setCurrentTime(const ScalarType& a_time) const { m_vField->setCurrentTime(a_time); }
  const ScalarType& getCurrentTime() const { return m_vField->getCurrentTime(); }
  const bool* getMask() const { return this->m_vField->getMask(); }
  bool getMask(size_type j) const { return this->m_vField->getMask(j); }
  bool getMask(size_type j, size_type c) const { return this->m_vField->getMask(j,c); }

  VectorFieldType* m_vField;
  MatrixType m_defaultC1InitialCondition;
  HessianType m_defaultC2InitialCondition;
  ScalarType m_fixedTimeStep;
  ScalarType m_step;

  typename Multiindex::IndicesSet m_listIndices;

  void operator=(const BasicCnOdeSolver&){}
  BasicCnOdeSolver(const BasicCnOdeSolver& s)
    : capd::dynsys::StepControlInterface<StepControlT,typename MapT::ScalarType>(s.getStepControl()),
      CurveT(0,0,1,1,1)
  {}

  void setInitialCondition(const JetType& coeff);
  void setInitialCondition(const VectorType& v, const MatrixType& D);
  void setInitialCondition(const VectorType& v, const MatrixType& D, const HessianType& H);

  void computeTimeStep(const VectorType& v);

  void evalAndSum(VectorType& v);
  void evalAndSum(VectorType& v, MatrixType& der);
  void evalAndSum(VectorType& v, MatrixType& der, HessianType& hessian);
  void evalAndSum(JetType& v);

  void sum(VectorType& v);
}; // the end of class BasicCnTaylor

// -----------------------------------------------------------------------------

template <typename MapType, typename StepControlType, typename CurveT>
inline const MapType& BasicCnOdeSolver<MapType, StepControlType, CurveT>::getVectorField() const {
  return *m_vField;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline MapType& BasicCnOdeSolver<MapType, StepControlType, CurveT>::getVectorField() {
  return *m_vField;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline typename BasicCnOdeSolver<MapType, StepControlType, CurveT>::ScalarType
BasicCnOdeSolver<MapType, StepControlType, CurveT>::getStep() const {
  return m_step;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline void BasicCnOdeSolver<MapType, StepControlType, CurveT>::setStep(const ScalarType& newStep) {
  m_fixedTimeStep = newStep;
  this->turnOffStepControl();
}

template <typename MapType, typename StepControlType, typename CurveT>
inline void BasicCnOdeSolver<MapType, StepControlType, CurveT>::adjustTimeStep(const ScalarType& newStep) {
  m_step = newStep;
}

template <typename MapType, typename StepControlType, typename CurveT>
inline void BasicCnOdeSolver<MapType, StepControlType, CurveT>::computeTimeStep(const VectorType& v) {
  m_step = this->isStepChangeAllowed()
      ? this->getStepControl().computeNextTimeStep(*this,this->getCurrentTime(),v)
      : capd::min(this->m_fixedTimeStep,this->getMaxStep());
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline typename BasicCnOdeSolver<MapT, StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT, StepControlT,CurveT>::operator()(ScalarType& t, const VectorType &v)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(v);
  t += this->getStep();
  return r;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline typename BasicCnOdeSolver<MapT, StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT, StepControlT,CurveT>::operator()(ScalarType& t, const VectorType &v, MatrixType& der)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(v,der);
  t += this->getStep();
  return r;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline typename BasicCnOdeSolver<MapT, StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT, StepControlT,CurveT>::operator()(ScalarType& t, const VectorType &v, MatrixType& der, HessianType& coeff)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(v,der,coeff);
  t += this->getStep();
  return r;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline typename BasicCnOdeSolver<MapT, StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT, StepControlT,CurveT>::operator()(
      ScalarType& t, const VectorType& x, const MatrixType& D,
      MatrixType& out_der
    )
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(x,D,out_der);
  t += this->getStep();
  return r;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline typename BasicCnOdeSolver<MapT, StepControlT,CurveT>::VectorType
BasicCnOdeSolver<MapT, StepControlT,CurveT>::operator()(
      ScalarType& t, const VectorType& x, const MatrixType& D, const HessianType& H,
      MatrixType& out_der, HessianType& out_hessian
    )
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(x,D,H,out_der,out_hessian);
  t += this->getStep();
  return r;
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline typename BasicCnOdeSolver<MapT, StepControlT,CurveT>::size_type
BasicCnOdeSolver<MapT, StepControlT,CurveT>::degree() const
{
  return this->m_vField->degree();
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline void BasicCnOdeSolver<MapT, StepControlT,CurveT>::computeCoefficientsAtCenter(const VectorType& x, size_type order)
{
  VectorType* coeff = this->getCoefficientsAtCenter();
  coeff[0] = x;
  this->m_vField->computeODECoefficients(coeff,order);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline void BasicCnOdeSolver<MapT, StepControlT,CurveT>::computeCoefficientsAtCenter(ScalarType t, const VectorType& x, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficientsAtCenter(x,order);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline void BasicCnOdeSolver<MapT, StepControlT,CurveT>::computeCoefficients(const VectorType& x, size_type order)
{
  this->setInitialCondition(x,m_defaultC1InitialCondition);
  this->m_vField->computeODECoefficients(this->getCoefficients(),1,order);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline void BasicCnOdeSolver<MapT, StepControlT,CurveT>::computeCoefficients(ScalarType t, const VectorType& x, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficients(x,order);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline void BasicCnOdeSolver<MapT, StepControlT,CurveT>::computeCoefficients(const VectorType& x, const MatrixType& M, size_type order)
{
  this->setInitialCondition(x,M);
  this->m_vField->computeODECoefficients(this->getCoefficients(),1,order);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlT,typename CurveT>
inline void BasicCnOdeSolver<MapT, StepControlT,CurveT>::computeCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, size_type order)
{
  this->setCurrentTime(t);
  this->computeCoefficients(x,M,order);
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICCNODESOLVER_H_


