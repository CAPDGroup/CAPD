

/////////////////////////////////////////////////////////////////////////////
/// @file C2OdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_C2SOLVER_H_
#define _CAPD_DYNSYS_C2SOLVER_H_

#include <string>
#include <vector>
#include "capd/vectalg/Norm.h"
#include "capd/dynset/C2Set.h"
#include "capd/dynsys/OdeSolver.h"
#include "capd/dynsys/C2DynSys.h"
#include "capd/dynsys/BasicC2OdeSolver.h"
#include "capd/diffAlgebra/C2TimeJet.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
template <
  typename MapT,
  typename StepControlPolicyT = capd::dynsys::ILastTermsStepControl,
  typename EnclosurePolicyT = HighOrderEnclosure,
  typename CurveT = capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<typename MapT::MatrixType> >
>
class C2OdeSolver : public OdeSolver<MapT,StepControlPolicyT,EnclosurePolicyT,CurveT>,
                 public BasicC2OdeSolver<MapT,StepControlPolicyT,CurveT>,
                 public C2DynSys<typename MapT::MatrixType>
{
public:
  typedef EnclosurePolicyT EnclosurePolicy;
  typedef MapT VectorFieldType;
  typedef StepControlPolicyT StepControlPolicy;
  typedef typename MapT::FunctionType FunctionType;
  typedef typename VectorFieldType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorFieldType::HessianType HessianType;
  typedef BasicC2OdeSolver<MapT,StepControlPolicy,CurveT> BaseC2Taylor;
  typedef OdeSolver<MapT,StepControlPolicy,EnclosurePolicy,CurveT> BaseTaylor;
  typedef CurveT SolutionCurve;
  typedef typename MatrixType::size_type size_type;
  typedef diffAlgebra::C2TimeJet<MatrixType> C2TimeJetType;

  C2OdeSolver(VectorFieldType& vectorField, size_type order);

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

  using BaseTaylor::computeRemainder;
  virtual void computeRemainder(ScalarType t, const VectorType& xx, C2TimeJetType& o_enc, C2TimeJetType& o_rem);
  using BaseTaylor::computeRemainderCoefficients;

  void initRemainderCoefficients(ScalarType t, const VectorType& x, unsigned degree){
    switch(degree){
      case 2:
        computeRemainderCoefficients(t,x,MatrixType::Identity(x.dimension()),HessianType(x.dimension(),x.dimension())); break;
      case 1:
        computeRemainderCoefficients(t,x,MatrixType::Identity(x.dimension())); break;
      default:
        computeRemainderCoefficients(t,x);
    }
  }

  void computeRemainderCoefficients(const VectorType& x, const MatrixType& M, const HessianType& H);
  void computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, const HessianType& H);
  ScalarType getStep() const{
    return BaseC2Taylor::getStep();
  }

  using BaseC2Taylor::getVectorField;
  using BaseC2Taylor::setOrder;
  using BaseC2Taylor::getOrder;
  using BaseC2Taylor::setStep;
  using BaseC2Taylor::dimension;

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C2SetMove is used meaning that this solver can integrate C^0, C^1 and C^2 sets only.
  /// Moreover, it cannot integrate nonrigorous jets (for user safety).
  template <typename SetType>
  void operator()(SetType& set){
    this->saveCurrentSet(set);
	  C2SetMove<C2OdeSolver,SetType>::move(set,*this);
  }


   /// Computes image of the set (in set's representation) and stores it in the result set.
   /// @param[in]  set       C^0, C^1  or C^2 set representing initial conditions
   /// @param[out] result    on return contains image of the set
  template <typename SetType>
  void operator()(SetType& set, SetType& result){
    this->saveCurrentSet(set);
	  C2SetMove<C2OdeSolver,SetType>::move(set, result, *this);
  }

protected:
  using BaseTaylor::saveCurrentSet;
  void saveCurrentSet(capd::dynset::C2Set<MatrixType>& set){
    this->setInitMatrix((MatrixType)set);
    this->setInitHessian((HessianType)set);
  }

  void operator=(const C2OdeSolver& ) {}
  C2OdeSolver(const C2OdeSolver& t) : BasicOdeSolver<MapT,StepControlPolicy,CurveT>(t), BaseTaylor(t), BaseC2Taylor(t) {}

  using BaseTaylor::m_step;
  using BaseTaylor::m_vField;
}; // the end of class C2Solver

// ####################################################################

template<typename MapType,typename StepControlPolicy, typename EnclosurePolicy, typename CurveT>
inline void C2OdeSolver<MapType,StepControlPolicy,EnclosurePolicy,CurveT>::computeRemainderCoefficients(ScalarType t, const VectorType& x, const MatrixType& M, const HessianType& H)
{
  this->setCurrentTime(t);
  this->computeRemainderCoefficients(x,M,H);
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_C2ODESOLVER_H_


