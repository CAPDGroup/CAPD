

/////////////////////////////////////////////////////////////////////////////
/// @file BasicC2OdeSolver.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_BASICC2ODESOLVER_H_
#define _CAPD_DYNSYS_BASICC2ODESOLVER_H_

#include <string>
#include "capd/diffAlgebra/BasicC2Curve.h"
#include "capd/diffAlgebra/C2Curve.h"
#include "capd/dynsys/BasicOdeSolver.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{
template <
  typename MapT,
  typename StepControlT = capd::dynsys::DLastTermsStepControl,
  typename CurveT = capd::diffAlgebra::C2Curve< capd::diffAlgebra::BasicC2Curve<typename MapT::MatrixType> >
>
class BasicC2OdeSolver : public virtual BasicOdeSolver<MapT, StepControlT, CurveT>{
public:
  typedef MapT VectorFieldType;
  typedef StepControlT StepControlType;
  typedef typename MapT::FunctionType FunctionType;
  typedef typename VectorFieldType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorFieldType::HessianType HessianType;
  typedef BasicOdeSolver<MapT, StepControlT, CurveT> BaseTaylor;
  typedef CurveT SolutionCurve;
  typedef typename MatrixType::size_type size_type;

  BasicC2OdeSolver(VectorFieldType& vectorField, size_type order, const StepControlT& stepControl = StepControlT());

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

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  /// In this case a specialization C2JetMove is used meaning that this solver can integrate C^0, C^1 and C^2 jets only.
  template <typename JetT>
  void operator()(JetT& jet){
    C2JetMove<BasicC2OdeSolver,JetT>::move(jet,*this);
  }

  using BaseTaylor::getVectorField;
  using BaseTaylor::dimension;
  using BaseTaylor::getCurve;
  using BaseTaylor::setStep;
  using BaseTaylor::getStep;
  using BaseTaylor::getOrder;
  using BaseTaylor::setOrder;
  using BaseTaylor::setCurrentTime;
  using BaseTaylor::getCurrentTime;
  using BaseTaylor::operator();

protected:
  using BaseTaylor::getMask;
  bool getMask(size_type j, size_type c) const {
    return this->m_vField->getMask(j,c);
  }
  BasicC2OdeSolver(const BasicC2OdeSolver& t) : BaseTaylor(t){}
  void operator=(const BasicC2OdeSolver&){}

  using BaseTaylor::m_vField;
  using BaseTaylor::m_step;

  void sumTaylorSeries(
      VectorType& v,
      MatrixType& der,
      HessianType& hessian,
      VectorType* coeff,
      MatrixType* matrixCoeff,
      HessianType* hessianCoeff,
      size_type order
    );
}; // the end of class BasicC2Taylor


// ------------------ inline definitions -------------

template <typename MapType, typename StepControlType,typename CurveT>
inline typename BasicC2OdeSolver<MapType, StepControlType,CurveT>::VectorType
BasicC2OdeSolver<MapType, StepControlType,CurveT>::operator()(ScalarType& t, const VectorType &v, MatrixType& der, HessianType& coeff)
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(v,der,coeff);
  t += this->getStep();
  return r;
}

//###########################################################//

template <typename MapType, typename StepControlType,typename CurveT>
inline typename BasicC2OdeSolver<MapType,StepControlType,CurveT>::VectorType
BasicC2OdeSolver<MapType,StepControlType,CurveT>::operator()(
      ScalarType& t, const VectorType& x, const MatrixType& D, const HessianType& H,
      MatrixType& out_der, HessianType& out_hessian
    )
{
  this->setCurrentTime(t);
  VectorType r = this->operator()(x,D,H,out_der,out_hessian);
  t += this->getStep();
  return r;
}
/// @}
}} // the end of the namespace capd::dynsys

#endif // _CAPD_DYNSYS_BASICC2ODESOLVER_H_


