/// @addtogroup dynsys
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file MpStepControl.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2021- by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_MP_STEP_CONTROL_H_
#define _CAPD_DYNSYS_MP_STEP_CONTROL_H_



namespace capd {
namespace dynsys {
#include "capd/multiPrec/mplib.h"
#include "capd/dynsys/StepControl.h"



///////////////////////////////////////////////////////////////////////////////////
///
/// Computes next time step using already computed Taylor series of the solution.
///
/// We choose time step so that expected error of the next step
/// of the ODE integration is less than epsilon.
///
/// @remark For non rigorous computations it is good to take more that one term
///         to prevent big time steps when the last term happens to be close to zero
///
///////////////////////////////////////////////////////////////////////////////////
template<class Solver>
typename Solver::ScalarType
computeNextStep(
		Solver& solver,                               ///< ODE solver containing Taylor series of the solution for current step
		int numberOfTerms,                            ///< number of terms in Taylor series to use in predictions
		const typename TypeTraits<typename Solver::ScalarType>::Real& epsilon,     ///< expected error for one step of integration
		const typename TypeTraits<typename Solver::ScalarType>::Real& minTimeStep, ///< minimal time step allowed
		int degree,
		bool clearMantissaFlag                  ///< decides if less important mantissa bits are set to zero (improves computations)
) {
	  typedef typename Solver::ScalarType ScalarType;
	  typedef typename TypeTraits<ScalarType>::Real Real;
	  int order = solver.getOrder();
	  Real optStep = 1.5*rightBound(solver.getMaxStep());

	  for(int i = order; (i > order - numberOfTerms) && (i > 0); i--) {
	    Real coeffNorm = rightBound(solver.getCoeffNorm(i,degree));
	    capd::rounding::DoubleRounding::roundNearest();
	    if(isSingular(coeffNorm)) continue;
	    Real errorTerm = epsilon/ coeffNorm;
	    Real step = exp(log(errorTerm)/i);
	    optStep = capd::min(optStep, step);
	  }

	  optStep = capd::max(optStep, minTimeStep);
	  if(clearMantissaFlag)
	    optStep = clearMantissaBits(optStep);

	  // note: solver.getMaxStep() might be a desired nondegenerated interval! Do not change it to leftBound.
	  return optStep < solver.getMaxStep() ? ScalarType(optStep) :  solver.getMaxStep();
}

// ---------------------------------------------------------------------

class MpLastTermsStepControl {
 public:
  typedef capd::MpFloat Real;

  explicit MpLastTermsStepControl(
  	int numberOfTerms = 1,
  	Real minimalTimeStep = 1. / 1048576.
  		) : m_numberOfTerms(numberOfTerms),
  		m_minTimeStep(minimalTimeStep) {
  }

  template <class Solver, class SetType>
  typename Solver::ScalarType computeNextTimeStep(Solver& solver, typename Solver::ScalarType const& /*t*/, const SetType& s) const {
    auto time = computeNextStep(
    	solver,
    	m_numberOfTerms,
    	Solver::getEffectiveTolerance(solver,s),
    	m_minTimeStep,
    	s.degree(),
    	false
    	);
  //  std::cout << "\n Computed time : " << time << std::endl;
    // bitWrite(std::cout, time);
    // hexWrite(std::cout, time);

    return time;
  }

  template <class Solver, class SetType>
  	void init(Solver& solver, typename Solver::ScalarType const& /*t*/, const SetType& s) const {
    typedef typename Solver::ScalarType ScalarType;
    typedef typename Solver::VectorType VectorType;
    typedef typename Solver::MatrixType MatrixType;
    typedef typename TypeTraits<ScalarType>::Real Float;

    VectorType x(s);

    ScalarType currentTime = s.getCurrentTime();
    MatrixType df = solver.getVectorField().derivative(currentTime,x);
    // cannot use EuclNorm because it is not implemented for VectorType=GeometricSeries but this is just an approximate bound
    Float lipConstant = sqrt(matrixAlgorithms::spectralRadiusOfSymMatrix(Transpose(df) * df, TypeTraits<ScalarType>::epsilon()*1.0e-4).rightBound());
    lipConstant = isSingular(lipConstant) ? Float(1.) : Float(1.)/lipConstant;
    Float h = capd::min(lipConstant,Float(1.));
    h = capd::min(h,solver.getMaxStep().leftBound());
    ScalarType I(Float(0.),h);
    x = x+I*solver.getVectorField()(currentTime+I,x);
    solver.initRemainderCoefficients(currentTime+I,x,s.degree());
    solver.adjustTimeStep(h);
  }

  Real getMinStepAllowed() const
  {
    return m_minTimeStep;
  }

  int m_numberOfTerms;
  Real m_minTimeStep;
};


// ---------------------------------------------------------------------
/// This class is a common interface for StepControl
/// used in PoincareMap and TimeMap.
/// Both classes inherit this interface
template<typename Scalar>
class StepControlInterface<MpLastTermsStepControl, Scalar> {
  typedef capd::MpFloat TolScalarType;
 public:
  typedef MpLastTermsStepControl StepControlType;

  StepControlInterface()
	  : m_onOffStepControl(true),
		m_absoluteTolerance(power(TolScalarType(10), -TypeTraits<Scalar>::numberOfDigits() - 3)),
		m_relativeTolerance(power(TolScalarType(10), -TypeTraits<Scalar>::numberOfDigits() - 3)),
		m_maxStep(1.e+100) {}

  StepControlInterface(const StepControlType &stepControl)
	  : m_stepControl(stepControl),
		m_onOffStepControl(true),
		m_absoluteTolerance(power(TolScalarType(10), -TypeTraits<Scalar>::numberOfDigits() - 3)),
		m_relativeTolerance(power(TolScalarType(10), -TypeTraits<Scalar>::numberOfDigits() - 3)),
		m_maxStep(1.e+100) {}

  void turnOnStepControl() {
	m_onOffStepControl = true;
  }

  void turnOffStepControl() {
	m_onOffStepControl = false;
  }

  void onOffStepControl(bool _onOffStepControl) {
	this->m_onOffStepControl = _onOffStepControl;
  }

  const StepControlType &getStepControl() const {
	return m_stepControl;
  }

  void setStepControl(const StepControlType &stepControl) {
	m_stepControl = stepControl;
  }

  bool isStepChangeAllowed() const {
	return m_onOffStepControl;
  }

  void setAbsoluteTolerance(TolScalarType tol) {
	m_absoluteTolerance = abs(tol);
  }

  void setRelativeTolerance(TolScalarType tol) {
	m_relativeTolerance = abs(tol);
  }

  TolScalarType getAbsoluteTolerance() const {
	return m_absoluteTolerance;
  }

  TolScalarType getRelativeTolerance() const {
	return m_relativeTolerance;
  }

  Scalar getMaxStep() const {
	return m_maxStep;
  }

  void setMaxStep(Scalar maxStep) {
	m_maxStep = maxStep;
  }

  template<class Solver, class SetType>
  inline static
  TolScalarType getEffectiveTolerance(Solver &solver, const SetType &s) {
	return capd::max(
		solver.getAbsoluteTolerance(),
		solver.getRelativeTolerance() * rightBound(solver.getCoeffNorm(0, s.degree()))
	);
  }

 protected:
  StepControlType m_stepControl;
  bool m_onOffStepControl;
  TolScalarType m_absoluteTolerance, m_relativeTolerance;
  Scalar m_maxStep;
};


}} // namespace capd::dynsys


#endif //  _CAPD_DYNSYS_MP_STEP_CONTROL_H_
/// @}
