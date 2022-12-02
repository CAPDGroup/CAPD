

/////////////////////////////////////////////////////////////////////////////
/// @file CnOdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_CNODESOLVER_H_
#define _CAPD_DYNSYS_CNODESOLVER_H_

#include <string>
#include <vector>
#include "capd/dynset/C1Set.h"
#include "capd/dynset/C2Set.h"
#include "capd/dynset/CnSet.h"
#include "capd/vectalg/Norm.h"
#include "capd/dynsys/CnDynSys.h"
#include "capd/dynsys/FirstOrderEnclosure.h"
#include "capd/dynsys/BasicCnOdeSolver.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/diffAlgebra/CnTimeJet.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template <
  typename MapT,
  typename StepControlT = capd::dynsys::IEncFoundStepControl,
  typename EnclosurePolicyT = FirstOrderEnclosure,
  typename CurveT = capd::diffAlgebra::CnCurve< capd::diffAlgebra::BasicCnCurve< typename MapT::MatrixType > >
  >
class CnOdeSolver : public CnDynSys<typename MapT::MatrixType>, public BasicCnOdeSolver<MapT,StepControlT,CurveT>
{
public:
  typedef MapT VectorFieldType;
  typedef StepControlT StepControlPolicy;
  typedef EnclosurePolicyT EnclosurePolicy;
  typedef BasicCnOdeSolver<MapT,StepControlPolicy,CurveT> BaseTaylor;
  typedef typename VectorFieldType::ScalarType ScalarType;
  typedef typename VectorFieldType::VectorType VectorType;
  typedef typename VectorFieldType::MatrixType MatrixType;
  typedef typename VectorFieldType::HessianType HessianType;
  typedef typename VectorFieldType::FunctionType FunctionType;
  typedef typename BaseTaylor::JetType JetType;
  typedef typename MatrixType::size_type size_type;
  typedef diffAlgebra::C1TimeJet<MatrixType> C1TimeJetType;
  typedef diffAlgebra::C2TimeJet<MatrixType> C2TimeJetType;
  typedef diffAlgebra::CnTimeJet<MatrixType,0> CnTimeJetType;

//  typedef typename SeriesContainerType::Multipointer Multipointer;
//  typedef typename SeriesContainerType::Multiindex Multiindex;
  typedef CurveT SolutionCurve;

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization CnSetMove is used meaning that this solver can integrate C^0, C^1 and C^2 sets only.
  /// Moreover, it cannot integrate nonrigorous jets (for user safety).
  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
    CnSetMove<CnOdeSolver,SetType>::move(set,*this);
  }

  template <typename SetType>
  void operator()(SetType& set, SetType & result){
    this->saveCurrentSet(set);
    CnSetMove<CnOdeSolver,SetType>::move(set, result, *this);
  }

  CnOdeSolver(VectorFieldType& vectorField, size_type order, const StepControlT& stepControl=StepControlT());

  // methods for C^0 computations
  VectorType Phi(const ScalarType& t,const VectorType& iv);
  MatrixType JacPhi(const ScalarType& t,const VectorType& iv);
  VectorType enclosure(const ScalarType& t, const VectorType& x) { return EnclosurePolicy::enclosure(*this,t,x); }
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

  // methods for C^1 computations
  void JacRemainder(
         const ScalarType& t,
         const VectorType& vecEnclosure,
         const MatrixType& jacEnclosure,
         VectorType& Remainder,
         MatrixType& jacRemainder
      ) ;
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

  // methods for C^2
  void encloseC2Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi,
      MatrixType& o_jacRem,
      MatrixType& o_jacEnc,
      HessianType& o_hessianPhi,
      HessianType& o_hessianRem,
      HessianType& o_hessianEnc
  );

  void c2Remainder(
        const VectorType& Enc,
        const MatrixType& jacEnc,
        const HessianType& hessianEnc,
        VectorType& o_Rem,
        MatrixType& o_jacRem,
        HessianType& o_hessianRem
    );

  template<class JetT>
  VectorType encloseCnMap(
        const ScalarType& t,
        const VectorType& x,
        const VectorType& xx,
        JetT& phi,
        JetT& rem,
        JetT& enc
      );

  template<class JetT>
  VectorType cnEnclosure(const ScalarType& t, const VectorType& x, JetT& result);

  template<class JetT>
  void cnRemainder(const JetT& enc, JetT& result) ;

  ScalarType getStep() const{
    return BaseTaylor::getStep();
  }

  ScalarType getCoeffNorm(size_type i, size_type degree) const;

  // the following methods provide an interface for generic algorithms based on an abstract solver
  void computeRemainderCoefficients(const VectorType& x);
  void computeRemainderCoefficients(const VectorType& x, const MatrixType& M);
  void computeRemainderCoefficients(ScalarType t, const VectorType& x);
  void computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M);

  virtual void computeRemainder(ScalarType t, const VectorType& xx, VectorType& o_enc, VectorType& o_rem);
  virtual void computeRemainder(ScalarType t, const VectorType& xx, C1TimeJetType& o_enc, C1TimeJetType& o_rem);
  virtual void computeRemainder(ScalarType t, const VectorType& xx, C2TimeJetType& o_enc, C2TimeJetType& o_rem);
  template<class JetT>
  void computeRemainder(ScalarType t, const VectorType& xx, JetT& o_enc, JetT& o_rem);

  void sumTaylorSeries(VectorType& o_phi, MatrixType& o_jacPhi);
  void sumTaylorSeries(VectorType& o_phi, MatrixType& o_jacPhi, HessianType& o_hessianPhi);
  template<class JetT>
  void sumTaylorSeries(JetT& o_phi);

  using BaseTaylor::getVectorField;
  using BaseTaylor::setOrder;
  using BaseTaylor::getOrder;
  using BaseTaylor::setStep;
  using BaseTaylor::dimension;
  using BaseTaylor::setCurrentTime;
  using BaseTaylor::getCurrentTime;
  using BaseTaylor::getCurve;

  // @override
  void computeTimeStep(const ScalarType& t, const VectorType& x){
    this->m_step = this->isStepChangeAllowed()
        ? this->getStepControl().computeNextTimeStep(*this,t,x)
        : capd::min(this->m_fixedTimeStep,this->getMaxStep());
  }

protected:
  void setInitialCondition(ScalarType t, const VectorType& x, const VectorType& xx);
  // TimeRange is base for all types of sets and nonrigorous CxCoeff
  void saveCurrentSet(const capd::diffAlgebra::TimeRange<ScalarType>& /*set*/){
  }

  void saveCurrentSet(capd::dynset::C1Set<MatrixType>& set){
    this->setInitMatrix((MatrixType)set);
  }
  void saveCurrentSet(capd::dynset::C2Set<MatrixType>& set){
    this->setInitMatrix((MatrixType)set);
    this->setInitHessian((HessianType)set);
  }
  void saveCurrentSet(capd::dynset::CnSet<MatrixType>& set){
    this->setInitJet(set.currentSet());
  }

  void operator=(const CnOdeSolver&){}
  CnOdeSolver(const CnOdeSolver& t) : BaseTaylor(t){}

  using BaseTaylor::m_step;
  using BaseTaylor::m_vField;

}; // the end of class CnSolver

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
inline void CnOdeSolver<MapT, StepControlPolicy,EnclosurePolicy, CurveT>::computeRemainderCoefficients(const VectorType& x)
{
  for(size_type i=0;i<this->dimension();++i)
    this->remainderCoefficient(i,0) = x[i];
  this->m_vField->computeODECoefficients(this->getRemainderCoefficients(),0,this->getOrder()+1);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
inline void CnOdeSolver<MapT, StepControlPolicy,EnclosurePolicy, CurveT>::computeRemainderCoefficients(ScalarType t, const VectorType& x)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
inline void CnOdeSolver<MapT, StepControlPolicy,EnclosurePolicy, CurveT>::computeRemainderCoefficients(const VectorType& x, const MatrixType& M)
{
  for(size_type i=0;i<this->dimension();++i)
  {
    this->remainderCoefficient(i,0) = x[i];
    for(size_type j=0;j<this->dimension();++j)
      this->remainderCoefficient(i,j,0) = M(i+1,j+1);
  }
  this->m_vField->computeODECoefficients(this->getRemainderCoefficients(),1,this->getOrder()+1);
}

// ---------------------------------------------------------------------------------------

template <typename MapT, typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
inline void CnOdeSolver<MapT, StepControlPolicy,EnclosurePolicy, CurveT>::computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x,M);
}
/// @}
}} // the end of the namespace capd::dynsys

#include "capd/dynsys/CnOdeSolver_templateMembers.h"

#endif // _CAPD_DYNSYS_CNODESOLVER_H_


