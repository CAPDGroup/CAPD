/////////////////////////////////////////////////////////////////////////////
/// @file TimeMap.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_TIME_MAP_H_
#define _CAPD_POINCARE_TIME_MAP_H_

#include <string>
#include "capd/dynset/C1Set.h"
#include "capd/dynset/C2Set.h"
#include "capd/dynset/CnSet.h"
#include "capd/diffAlgebra/SolutionCurve.h"
#include "capd/poincare/SaveStepControl.h"

namespace capd{
namespace poincare{
/// @addtogroup poincare 
/// @{
/**
*  TimeMap class provides methods for transport of sets (or points) by a given flow over some time interval.
*  Template parameter is an abstract solver of dynamical system.
*  It is assumed that it provides methods for one-step transport (both for discrete and continuous case).
*/

template<typename SolverT>
class TimeMap
{
public:
  typedef SolverT Solver;
  typedef typename Solver::VectorFieldType VectorFieldType;
  typedef typename Solver::MatrixType MatrixType;
  typedef typename Solver::VectorType VectorType;
  typedef typename Solver::ScalarType ScalarType;

  typedef typename Solver::SolutionCurve CurveType;
  typedef capd::diffAlgebra::SolutionCurve<CurveType> SolutionCurve;
  typedef typename CurveType::HessianType HessianType;        ///< data structure for storing Hessians
  typedef typename CurveType::JetType JetType;                ///< data structure for storing jets (truncated Taylor series)
  typedef typename VectorType::size_type size_type;           ///< integral type used to index containers (vectors, matrices, etc)

  TimeMap(Solver& solver);

  /// For a vector v it computes its position after time time
  VectorType operator()(ScalarType time, VectorType& v);
  VectorType operator()(ScalarType time, VectorType& v, ScalarType& in_out_time);
  VectorType operator()(ScalarType time, VectorType& v, SolutionCurve& solution);

  /// For a vector v it computes its position after time 'time' and derivative with respect to initial conditions
  VectorType operator()(ScalarType time, VectorType& v, MatrixType& derivative);
  VectorType operator()(ScalarType time, VectorType& v, MatrixType& derivative, ScalarType& in_out_time);
  VectorType operator()(ScalarType time, VectorType& v, MatrixType& derivative, SolutionCurve& solution);
  VectorType operator()(ScalarType time, VectorType& v, const MatrixType& initMatrix, MatrixType& derivative);
  VectorType operator()(ScalarType time, VectorType& v, const MatrixType& initMatrix, MatrixType& derivative, ScalarType& in_out_time);
  VectorType operator()(ScalarType time, VectorType& v, const MatrixType& initMatrix, MatrixType& derivative, SolutionCurve& solution);

  /// Operators for integration of second order variational equations
  VectorType operator()(
        ScalarType time,
        VectorType& v, MatrixType& derivative, HessianType& hessian
      );
  VectorType operator()(
        ScalarType time,
        VectorType& v, MatrixType& derivative, HessianType& hessian,
        ScalarType& in_out_time
      );
  VectorType operator()(
        ScalarType time,
        VectorType& v, MatrixType& derivative, HessianType& hessian,
        SolutionCurve& solution
      );
  VectorType operator()(
        ScalarType time,
        VectorType& v, const MatrixType& initMatrix, const HessianType& initHessian,
        MatrixType& derivative, HessianType& hessian
      );
  VectorType operator()(
        ScalarType time,
        VectorType& v, const MatrixType& initMatrix, const HessianType& initHessian,
        MatrixType& derivative, HessianType& hessian,
        ScalarType& in_out_time
      );
  VectorType operator()(
        ScalarType time,
        VectorType& v, const MatrixType& initMatrix, const HessianType& initHessian,
        MatrixType& derivative, HessianType& hessian,
        SolutionCurve& solution
      );

  /// For a vector v it computes its position after time 'time' and higher order derivatives with respect to initial conditions
  VectorType operator()(ScalarType time, VectorType& v, JetType& jet);
  VectorType operator()(ScalarType time, VectorType& v, JetType& jet, ScalarType& in_out_time);
  VectorType operator()(ScalarType time, VectorType& v, JetType& jet, SolutionCurve& solution);
  VectorType operator()(ScalarType time, const JetType& initJet, JetType& jet);
  VectorType operator()(ScalarType time, const JetType& initJet, JetType& jet, ScalarType& in_out_time);
  VectorType operator()(ScalarType time, const JetType& initJet, JetType& jet, SolutionCurve& solution);

  template<class CnCoeffType>
  VectorType operator()(ScalarType time, VectorType& v, CnCoeffType& result);

  /// Integrates 'theSet' until time 'time'.
  /// @param [in,out] theSet - initial condition for integration. 'theSet.getCurrentTime()' should be an initial time for nonautonomous systems.
  /// On output 'theSet' represents solution to ODE at time 'time'
  /// @param [out] derivative monodromy matrix wrt to initial condition at time 'time'.
  template<typename SetType>
  VectorType operator()(ScalarType time, SetType &theSet, MatrixType& derivative);

  /// Integrates 'theSet' until time 'time'.
  /// @param [in,out] theSet - initial condition for integration. 'theSet.getCurrentTime()' should be an initial time for nonautonomous systems.
  /// On output 'theSet' represents solution to ODE at time 'time'
  template<typename SetType>
  VectorType operator()(ScalarType time, SetType &theSet);

  /// Integrates 'theSet' until time 'time'.
  /// @param [in,out] theSet - initial condition for integration. 'theSet.getCurrentTime()' should be an initial time for nonautonomous systems.
  /// On output 'theSet' represents solution to ODE at time 'time'
  /// @param[out] solution - after successful integration becomes a functional object representing entire solution to ODE that can be evaluated at any intermediate time.
  template<typename SetType>
  VectorType operator()(ScalarType time, SetType &theSet, SolutionCurve& solution);

  const Solver & getSolver() const;    ///< Returns read-only reference to solver used to integrate the system
  Solver & getSolver();                ///< Returns reference to solver used to integrate the system

  /// @deprecated
  const Solver& getDynamicalSystem() const;    ///< Returns read-only reference to solver used to integrate the system
  /// @deprecated
  Solver & getDynamicalSystem();                ///< Returns reference to solver used to integrate the system

  const VectorFieldType& getVectorField() const;          ///< Returns read-only reference to current vector field
  VectorFieldType& getVectorField();                      ///< Returns reference to current vector field

  size_type getOrder() const;                     ///< Returns order of numerical method
  void setOrder(size_type newOrder);              ///< sets new order of numerical method

  ScalarType getStep() const;               ///< Returns step of numerical method (does make sense only when step control is disabled)
  void setStep(const ScalarType& newStep);  ///< Sets step of numerical method (does make sense only when step control is disabled)

  void turnOnStepControl();                 ///< Disables automatic step control
  void turnOffStepControl();                ///< Enables automatic step control. Step control strategy is builded into the dynamical system.
  void onOffStepControl(bool sc);           ///< Disables or enables automatic step control

  void stopAfterStep(bool b);               ///< For dense output. If true, integration procedure returns after each successful step. The computation can be then resumed - see examples.
  bool completed() const;                   ///< For dense output. Returns true if the trajectory has been integrated to the requested time.
  const ScalarType& getCurrentTime() const; ///< For dense output. Returns current time during integration.

protected:
  template<class SetType>
  void initComputations(SetType& originalSet, void* current);

  template<class SetType>
  void moveSet(ScalarType time, SetType& originalSet, void* current);
/*
  void saveCurrentSet(capd::dynset::CnSet<MatrixType>& set){
    this->m_solver.setInitJet(set.currentSet());
  }
*/
  Solver& m_solver;           ///< an abstract solver of dynamical system

  void* m_currentSet;         ///< for dense output. Pointer to an object representing initial condition that is integrated.
  bool m_oneStep;             ///< for dense output. True if the procedure should return after each successful step.
  ScalarType m_currentTime;   ///< Holds current time during integration.
  bool m_completed;           ///< flag: true if trajectory has been integrated until requested time.
  SolutionCurve* m_solution;
}; // end of template TimeMap

/// @}
}} // namespace capd::poincare

#include "capd/poincare/TimeMap_inline.h"
#include "capd/poincare/TimeMap_template.h"

#endif // _CAPD_POINCARE_TIME_MAP_H_



