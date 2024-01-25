

/////////////////////////////////////////////////////////////////////////////
/// @file FadOdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2017 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_FADODESOLVER_H_
#define _CAPD_DYNSYS_FADODESOLVER_H_

#include "capd/vectalg/Norm.h"
#include "capd/dynsys/C1DynSys.h"
#include "capd/dynset/C0Set.h"
#include "capd/dynset/C1Set.h"
#include "capd/dynsys/BasicFadOdeSolver.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<class FadMapT, typename StepControlT = capd::dynsys::IEncFoundStepControl >
class FadOdeSolver : public capd::dynsys::C1DynSys<typename FadMapT::MatrixType>, public BasicFadOdeSolver<FadMapT,StepControlT>{
public:
  typedef FadMapT VectorFieldType;
  typedef BasicFadOdeSolver<FadMapT,StepControlT> BaseTaylor;
  typedef StepControlT StepControlType;
  typedef typename BaseTaylor::ScalarType ScalarType;
  typedef typename BaseTaylor::FScalar FScalar;
  typedef typename BaseTaylor::TFScalar TFScalar;
  typedef typename BaseTaylor::MatrixType MatrixType;
  typedef typename BaseTaylor::VectorType VectorType;
  typedef typename BaseTaylor::FVector FVector;
  typedef typename BaseTaylor::TFVector TFVector;
  typedef typename BaseTaylor::FunctionType FunctionType;
  typedef typename MatrixType::size_type size_type;

  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;

  FadOdeSolver(VectorFieldType& f, size_type _order, const StepControlT& _stepControl=StepControlT());

  // methods for C^0 computations
  VectorType Phi(const ScalarType& t,const VectorType& iv);
  MatrixType JacPhi(const ScalarType& t,const VectorType& iv);
  VectorType enclosure(const ScalarType& t, const VectorType& x);
  VectorType Remainder(const ScalarType& t, const VectorType& iv, VectorType& o_enc);

  void computeTaylorCoefficients(ScalarType t, const VectorType& x, const VectorType& xx);

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
  MatrixType jacEnclosure(const ScalarType& t, const VectorType& enc);
  void JacRemainder(
         const ScalarType& t,
         const VectorType &vecEnclosure,
         const MatrixType &jacEnclosure,
         VectorType &Remainder,
         MatrixType &jRemainder
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

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C1SetMove is used meaning that this solver can integrate C^0 and C^1 sets only.
  /// Moreover, it cannot integrate nonrigorous jets (for user safety).
  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
	  C1SetMove<FadOdeSolver,SetType>::move(set,*this);
  }
  /// Computes image of the set (in set's representation) and stores it in the result set.
   /// @param[in]  set       C^0, C^1  set representing initial conditions
   /// @param[out] result    on return contains image of the set
  template <typename SetType>
  void operator()(SetType& set, SetType& result){
    this->saveCurrentSet(set);
	  C1SetMove<FadOdeSolver,SetType>::move(set, result, *this);
  }

  ScalarType getCoeffNorm(size_type i, size_type degree) const;
  ScalarType getStep() const{
    return BaseTaylor::getStep();
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
protected:
  // TimeRange is base for all types of sets and nonrigorous CxCoeff
  void saveCurrentSet(const capd::diffAlgebra::TimeRange<ScalarType>&/* set*/){
  }

  void saveCurrentSet(const capd::dynset::C1Set<MatrixType>& set){
    this->setInitMatrix((MatrixType)set);
  }

  // @override
  void computeTimeStep(const ScalarType& t, const VectorType& x){
    this->m_step = this->isStepChangeAllowed()
        ? this->getStepControl().computeNextTimeStep(*this,t,x)
        : capd::min(this->m_fixedTimeStep,this->getMaxStep());
  }

}; // the end of class FadTaylor

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void FadOdeSolver<FadMapT,StepControlT>::computeRemainderCoefficients(const VectorType& x)
{
  this->setInitialCondition(x,this->m_rem);
  this->computeCoeff(this->m_rem,this->m_remOut,this->getOrder()+1);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void FadOdeSolver<FadMapT,StepControlT>::computeRemainderCoefficients(ScalarType t, const VectorType& x)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x);
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void FadOdeSolver<FadMapT,StepControlT>::computeRemainderCoefficients(const VectorType& x, const MatrixType& M)
{
  this->setInitialCondition(x,M,this->m_jacRem);
  this->computeCoeff(this->m_jacRem,this->m_jacRemOut,this->getOrder()+1);

  for(size_type i=0;i<this->dimension();++i)
    for(size_type j=0;j<=this->getOrder()+1;++j)
      this->remainderCoefficient(i,j) = this->m_jacRem[i][j].x();
}

//###########################################################//

template<class FadMapT, typename StepControlT>
inline void FadOdeSolver<FadMapT,StepControlT>::computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x,M);
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_FADODESOLVER_H_


