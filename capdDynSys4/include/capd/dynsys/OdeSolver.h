  /// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file OdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_ODESOLVER_H_
#define _CAPD_DYNSYS_ODESOLVER_H_

#include <string>
#include <stdexcept>
#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/dynset/C0Set.h"
#include "capd/dynset/C1Set.h"
#include "capd/dynsys/C1DynSys.h"
#include "capd/dynsys/BasicOdeSolver.h"
#include "capd/dynsys/FirstOrderEnclosure.h"
#include "capd/dynsys/HighOrderEnclosure.h"
#include "capd/poincare/SaveStepControl.h"

namespace capd{
namespace dynsys{

template <
  typename MapT,
  typename StepControlPolicyT = capd::dynsys::ILastTermsStepControl,
  typename EnclosurePolicy = HighOrderEnclosure,
  typename CurveT = capd::diffAlgebra::Curve< capd::diffAlgebra::BasicCurve<typename MapT::MatrixType> >
>
class OdeSolver: public virtual C1DynSys<typename MapT::MatrixType>, public virtual BasicOdeSolver<MapT,StepControlPolicyT,CurveT>
{
public:
  typedef MapT VectorFieldType;
  typedef StepControlPolicyT StepControlPolicy;
  typedef typename MapT::FunctionType FunctionType;
  typedef typename VectorFieldType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef BasicOdeSolver<MapT,StepControlPolicy,CurveT> BaseTaylor;
  typedef typename MatrixType::size_type size_type;
  typedef diffAlgebra::C1TimeJet<MatrixType> C1TimeJetType;

  OdeSolver(VectorFieldType& vField, size_type order, const StepControlPolicy& stepControl = StepControlPolicy());
  void setOrder(size_type order); ///< Sets the order of the Taylor method


  // implementation of DynSys interface
  VectorType Phi(const ScalarType& t,const VectorType& iv);
  MatrixType JacPhi(const ScalarType& t,const VectorType& iv);

  VectorType enclosure(const ScalarType& t, const VectorType& x);
  VectorType Remainder(const ScalarType& t, const VectorType& iv, VectorType& o_enc);

  void encloseC0Map(
      const ScalarType& t,  //< @param[in] current time of ODE
      const VectorType& x0, //< @param[in] an internal point of x, usually center of x
      const VectorType& x,  //< @param[in] set to be moved along the trajectories of ODE
      VectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
      VectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
      VectorType& o_enc,    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
      MatrixType& o_jacPhi  //< @param[out] bound for derivative Dphi(x)
  );


  void encloseC1Map(
      const ScalarType& t,  //< @param[in] current time of ODE
      const VectorType& x0, //< @param[in] an internal point of x, usually center of x
      const VectorType& x,  //< @param[in] set to be moved along the trajectories of ODE
      VectorType& o_phi,    //< @param[out] bound for phi(x0), where phi is a numerical method
      VectorType& o_rem,    //< @param[out] bound for the error of numerical method over the time step
      VectorType& o_enc,    //< @param[out] enclosure of all trajectories starting from x over the time interval (time step of numerical method)
      MatrixType& o_jacPhi, //< @param[out] bound for derivative Dphi(x)
      MatrixType& o_jacRem, //< @param[out] bound for the error of numerical method over the time step for variational equation
      MatrixType& o_jacEnc  //< @param[out] enclosure of all trajectories of variational equations with initial condition set to Identity over the time interval (time step of numerical method)
  );

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C1SetMove is used meaning that this solver can integrate C^0 and C^1 sets only.
  /// Moreover, it cannot integrate nonrigorous jets (for user safety).
  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
    C1SetMove<OdeSolver,SetType>::move(set,*this);
  }

  /// Computes image of the set (in set's representation) and stores it in the result set.
   /// @param[in]  set       C^0 or C^1  set representing initial conditions
   /// @param[out] result    on return contains image of the set
  template <typename SetType>
  void operator()(SetType& set, SetType& result){
    this->saveCurrentSet(set);
    C1SetMove<OdeSolver,SetType>::move(set, result, *this);
  }

  // the following methods provide an interface for generic algorithms based on an abstract solver
  void computeRemainderCoefficients(const VectorType& x);
  void computeRemainderCoefficients(const VectorType& x, const MatrixType& M);
  void computeRemainderCoefficients(ScalarType t, const VectorType& x);
  void computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M);
  void initRemainderCoefficients(ScalarType t, const VectorType& x, unsigned degree){
    if(degree)
      computeRemainderCoefficients(t,x,MatrixType::Identity(x.dimension()));
    else
      computeRemainderCoefficients(t,x);
  }
  virtual void computeRemainder(ScalarType t, const VectorType& xx, VectorType& o_enc, VectorType& o_rem);
  virtual void computeRemainder(ScalarType t, const VectorType& xx, C1TimeJetType& o_enc, C1TimeJetType& o_rem);

  void computeTaylorCoefficients(ScalarType t, const VectorType& x, const VectorType& xx);
  void computeImplicitCoefficients(ScalarType t, const VectorType& x, const VectorType& xx, size_type order);

  ScalarType getCoeffNorm(size_type, size_type degree) const;
  ScalarType getStep() const{
    return BaseTaylor::getStep();
  }

  const CurveT& getImplicitCurve() const{
    return implicitCurve;
  }

  // @override
  void computeTimeStep(const ScalarType& t, const VectorType& x){
    this->m_step = this->isStepChangeAllowed()
        ? this->getStepControl().computeNextTimeStep(*this,t,x)
        : capd::min(this->m_fixedTimeStep,this->getMaxStep());
  }

protected:

  // TimeRange is base for all types of sets and nonrigorous CxCoeff
  void saveCurrentSet(const capd::diffAlgebra::TimeRange<ScalarType>& /*set*/){
  }

  void saveCurrentSet(const capd::dynset::C1Set<MatrixType>& set){
    this->setInitMatrix((MatrixType)set);
  }


  void operator=(const OdeSolver& ){}
  OdeSolver(const OdeSolver& t) : BaseTaylor(t), implicitCurve(0.,0.,1,1,1){}

  CurveT implicitCurve; ///< an extra storage for Taylor coefficients used in implicit HO-method
};

// --------------- inline definitions -----------------

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainderCoefficients(const VectorType& x)
{
  VectorType* coeff = this->getRemainderCoefficients();
  coeff[0] = x;
  this->m_vField->computeODECoefficients(coeff,this->getOrder()+1);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainderCoefficients(ScalarType t, const VectorType& x)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainderCoefficients(const VectorType& x, const MatrixType& M)
{
  VectorType* coeff = this->getRemainderCoefficients();
  MatrixType* matrixCoeff = this->getMatrixRemainderCoefficients();
  coeff[0] = x;
  matrixCoeff[0] = M;
  this->m_vField->computeODECoefficients(coeff,matrixCoeff,this->getOrder()+1);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x,M);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainder(
      ScalarType t,
      const VectorType& xx,
      VectorType& o_enc,
      VectorType& o_rem
  )
{
  EnclosurePolicy::computeEnclosureAndRemainder(*this,t,xx,o_enc,o_rem);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeRemainder(
      ScalarType t,
      const VectorType& xx,
      C1TimeJetType& o_enc,
      C1TimeJetType& o_rem
  )
{
  EnclosurePolicy::computeEnclosureAndRemainder(*this,t,xx,o_enc,o_rem);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeImplicitCoefficients(
    ScalarType t, const VectorType& x, const VectorType& xx, size_type order
  )
{
  this->setCurrentTime(t);
  VectorType* coeff = this->implicitCurve.getCoefficientsAtCenter();
  coeff[0] = x;
  this->m_vField->computeODECoefficients(coeff,order);

  coeff = this->implicitCurve.getCoefficients();
  MatrixType* matrixCoeff = this->implicitCurve.getMatrixCoefficients();
  coeff[0] = xx;
  matrixCoeff[0].setToIdentity();
  this->m_vField->computeODECoefficients(coeff,matrixCoeff,order);
}

//###########################################################//

template<typename MapType, typename StepControlPolicy,typename EnclosurePolicy, typename CurveType>
inline void OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveType>::computeTaylorCoefficients(
    ScalarType t, const VectorType& x, const VectorType& xx
  )
{
  this->setCurrentTime(t);
  this->computeCoefficientsAtCenter(x,this->getOrder());
  this->computeCoefficients(xx,this->getOrder());
}

}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_ODESOLVER_H_

/// @}
