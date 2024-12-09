/// @addtogroup diffIncl
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusion.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Tomasz Kapela, 2007 */

#ifndef _CAPD_DIFFINCL_DIFFINCLUSION_H_
#define _CAPD_DIFFINCL_DIFFINCLUSION_H_

#include <string>
#include <map>
#include "capd/vectalg/Norm.h"
#include "capd/dynsys/OdeSolver.h"
#include "capd/dynsys/FirstOrderEnclosure.h"
#include "capd/dynsys/StepControl.h"
namespace capd{
/// A rigorous integration of the differential inclusions
namespace diffIncl{

/**
 * Base class for rigorous integration of differential inclusions.
 *
 * Template arguments:
 * - MapT     - MultiMap that stores RHS of the differential inclusion in the form : selection + 'error bounds'
 *              (we assume that it implements all methods that class capd::diffIncl::MultiMap has).
 * - DynSysT  - numerical method for ODE integration
 *
 * \see capd::diffIncl::DiffInclusionLN, \see capd::diffIncl::DiffInclusionCW
 */
//template<typename MapT, typename DynSysT = capd::dynsys::Solver< typename MapT::MapType > >
template<typename MapT, typename DynSysT = capd::dynsys::OdeSolver< typename MapT::MapType > >
class DiffInclusion
//: public capd::dynsys::StepControlInterface<typename DynSysT::StepControlType, typename MapT::ScalarType> {
{
// TODO : check if step control is doubled
// DW: yes, it is!
public:
  typedef MapT MultiMapType;
  typedef MultiMapType MapType;
  typedef typename MapT::MapType VectorFieldType;
  typedef typename MapT::FunctionType FunctionType;
  typedef typename MapType::MatrixType MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType, MatrixType> NormType;
  typedef DynSysT DynSysType;
  typedef typename DynSysType::SolutionCurve SolutionCurve;
  typedef typename DynSysT::StepControlType  StepControlType;

  DiffInclusion(   MultiMapType & diffInclusion,               // map of the form f(x)+ g(x,e)
      int order,                            // order for integration
      NormType const & norm                 // norm
  );

  virtual ~DiffInclusion();

  /// eclosure of solution of diff. inclusion during one time step starting at x
  virtual VectorType diffInclusionEnclosure(const ScalarType& t, const VectorType& x);

  /// eclosure of solution of selected ODE during one time step starting at x
  virtual VectorType dynamicalSystemEnclosure(const ScalarType& t, const VectorType& x) ;

  /// eclosure of solution of diff. inclusion during one time step starting at x
  virtual VectorType enclosure(const ScalarType& t,const VectorType & x);

  /// returns RHS of a diff. inclusion
  const MapType& getVectorField() const;
  /// returns RHS of a diff. inclusion
  MapType& getVectorField();

  /// returns order of numerical method
  int getOrder() const;
  /// sets order of numerical method
  void setOrder(int newOrder);

  /// returns current time step
  ScalarType getStep() const;
  /// sets currect time step
  void setStep(const ScalarType& newStep);

  ///< sets time step but does not change step control settings (compare setStep)
  void adjustTimeStep(const ScalarType& newStep);

  /// returns dynamical system (numerical ODE integrator)
  DynSysType & getDynamicalSystem();
  /// returns dynamical system (numerical ODE integrator)
  const DynSysType & getDynamicalSystem() const;

  int dimension() const {
    return m_dynamicalSystem.dimension();
  }
  const SolutionCurve& getCurve() {
    return  m_dynamicalSystem.getCurve();
  }
  virtual ScalarType getCoeffNorm(unsigned i, int degree) const{
    return  m_dynamicalSystem.getCoeffNorm(i,degree);
  }

  void clearCoefficients(){
    m_dynamicalSystem.clearCoefficients();
  }

  template <class SetType>
  ScalarType computeNextTimeStep(const SetType& x, const ScalarType& maxStep) {
    return m_dynamicalSystem.computeNextTimeStep(*this,x,maxStep);
  }

  template <class SetType>
  ScalarType getFirstTimeStep(const SetType& x, const ScalarType& maxStep) {
    return m_dynamicalSystem.getFirstTimeStep(*this, x, maxStep);
  }

  void turnOnStepControl() {
    m_dynamicalSystem.turnOnStepControl();
  }

  void turnOffStepControl() {
   m_dynamicalSystem.turnOffStepControl();
  }

  void onOffStepControl(bool _onOffStepControl) {
   m_dynamicalSystem.onOffStepControl(_onOffStepControl);
  }

  const StepControlType& getStepControl() const {
    return m_dynamicalSystem.getStepControl();
  }

  void setStepControl(const StepControlType& stepControl) {
   m_dynamicalSystem.setStepControl( stepControl);
  }

  bool isStepChangeAllowed() const {
    return m_dynamicalSystem.isStepChangeAllowed();
  }

  void setAbsoluteTolerance(double tol){
    m_dynamicalSystem.setAbsoluteTolerance(tol);
  }

  void setRelativeTolerance(double tol){
   m_dynamicalSystem.setRelativeTolerance(tol);
  }

  double getAbsoluteTolerance() const{
    return m_dynamicalSystem.getAbsoluteTolerance();
  }

  double getRelativeTolerance() const{
    return m_dynamicalSystem.getRelativeTolerance();
  }

  ScalarType getMaxStep() const{
    return m_dynamicalSystem.getMaxStep();
  }

  void setMaxStep(ScalarType maxStep){
    m_dynamicalSystem.setMaxStep(maxStep);
  }

  template <class Solver, class SetType>
  inline static
  double getEffectiveTolerance(Solver & solver, const SetType& s) {
    return DynSysType::getEffectiveTolerance(solver, s);
  }

protected:
  //  void operator=(const DiffInclusion & a_t){}
  // DiffInclusion(const DiffInclusion & diffIncl) : m_norm(m_norm.clone()), m_diffIncl(diffIncl.m_diffIncl){}

  NormType * m_norm;                  ///<  norm used in perturbation estimations
  MultiMapType & m_diffIncl;          ///<  RHS of differential inclusion
  DynSysType m_dynamicalSystem;       ///<  dynamical system of selected ODE (numerical integrator)
};



template<class T, class SetT
//  ,bool isSet =  capd::dynset::SetTraits<typename SetT::SetType>::isC0Set
//             or capd::dynset::SetTraits<typename SetT::SetType>::isC1Set
           >
struct DiffInclusionSetMove{
	  static void move(SetT& set, T& solver){
		  set.move(solver);
	  }
    static void move(SetT& set, SetT& result, T& solver){
		  set.move(solver, result);
	  }
};
// --------------- inline definitions -----------------


template <typename MapT, typename DynSysT>
inline const MapT& DiffInclusion<MapT, DynSysT>::getVectorField() const {
  return m_diffIncl;
}

template <typename MapT, typename DynSysT>
inline MapT & DiffInclusion<MapT, DynSysT>::getVectorField() {
  return m_diffIncl;
}


template <typename MapT, typename DynSysT>
inline int DiffInclusion<MapT, DynSysT>::getOrder() const {
  return m_dynamicalSystem.getOrder();
}

template <typename MapT, typename DynSysT>
inline void DiffInclusion<MapT, DynSysT>::setOrder(int newOrder) {
  m_dynamicalSystem.setOrder(newOrder);
}


template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::ScalarType DiffInclusion<MapT, DynSysT>::getStep() const {
  return m_dynamicalSystem.getStep();
}

template <typename MapT, typename DynSysT>
inline void DiffInclusion<MapT, DynSysT>::setStep(const ScalarType& newStep) {
  m_dynamicalSystem.setStep(newStep);
}

template <typename MapT, typename DynSysT>
inline void DiffInclusion<MapT, DynSysT>::adjustTimeStep(const ScalarType& newStep){ ///< sets time step but does not change step control settings (compare setStep)
  m_dynamicalSystem.adjustTimeStep(newStep);
}

/**
 * Computes enclosure of image of given set for differential inclusion during whole time step
 */
template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::VectorType DiffInclusion<MapT, DynSysT>::diffInclusionEnclosure(const ScalarType& t, const VectorType& x) {

  return capd::dynsys::FirstOrderEnclosure::enclosure(m_diffIncl, t, x, getStep());
}

template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::VectorType DiffInclusion<MapT, DynSysT>::dynamicalSystemEnclosure(const ScalarType& t, const VectorType & x)  {
  return m_dynamicalSystem.enclosure(t, x);
}

template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::VectorType DiffInclusion<MapT, DynSysT>::enclosure(const ScalarType& t, const VectorType& x) {

  return diffInclusionEnclosure(t, x);
}

template <typename MapT, typename DynSysT>
inline typename DiffInclusion<MapT, DynSysT>::DynSysType & DiffInclusion<MapT, DynSysT>::getDynamicalSystem(){
  return m_dynamicalSystem;
}

template <typename MapT, typename DynSysT>
inline const typename DiffInclusion<MapT, DynSysT>::DynSysType&  DiffInclusion<MapT, DynSysT>::getDynamicalSystem() const{
  return m_dynamicalSystem;
}


//###########################################################//


}} // namespace capd::diffIncl

#endif // _CAPD_DIFFINCL_DIFFINCLUSION_H_

/// @}
