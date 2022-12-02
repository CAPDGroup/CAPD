/// @addtogroup diffIncl
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file DiffInclusionCW.h
///
/// @author Tomasz Kapela
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Tomasz Kapela, 2007 */

#ifndef _CAPD_DIFFINCL_DIFFINCLUSIONCW_H_
#define _CAPD_DIFFINCL_DIFFINCLUSIONCW_H_

#include "capd/diffIncl/DiffInclusion.h"

namespace capd{
namespace diffIncl{

/**
 * Class for rigorous integration of differential inclusions.
 * It uses Component Wise method to bound perturbations.
 *
 * For more details on algorithm see: <br>
 *   T. Kapela, P. Zgliczyński, A Lohner-type algorithm for control systems and ordinary differential inclusions,
 *   Discrete and Continuous Dynamical Systems B, 11 (2009), 365-385.
 *
 * Template arguments:
 * - MapT     - MultiMap that stores RHS of the differential inclusion in the form : selection + 'error bounds'
 *              (we assume that it implements all methods that class capd::diffIncl::MultiMap has).
 * - DynSysT  - numerical method for ODE integration
 */
template<typename MapT, typename DynSysT = capd::dynsys::OdeSolver< typename MapT::MapType> >
class DiffInclusionCW : public capd::diffIncl::DiffInclusion<MapT, DynSysT>{

public:
  typedef DiffInclusion<MapT, DynSysT> BaseClass;
  typedef typename BaseClass::MultiMapType   MultiMapType;
  typedef typename BaseClass::MapType        MapType;
  typedef typename BaseClass::FunctionType   FunctionType;
  typedef typename BaseClass::MatrixType     MatrixType;
  typedef typename BaseClass::VectorType     VectorType;
  typedef typename BaseClass::ScalarType     ScalarType;
  typedef typename BaseClass::NormType       NormType;
  typedef typename MatrixType::size_type     size_type;

  DiffInclusionCW( MultiMapType & diffIncl,       // map of the form f(x)+ g(x,e)
      size_type order,                                  // order of Taylor method for integration
      NormType const & norm,                      // norm used in computation of perturbations
      ScalarType const & expErrorTolerance =  capd::TypeTraits<ScalarType>::epsilon()   // error tolerance in the exp(A) computations
  );

  /// Bounds for perturbation of solution of selected ODE after one time step
  VectorType perturbations(const ScalarType & time, const VectorType & x);

  using BaseClass::enclosure;
  using BaseClass::diffInclusionEnclosure;
  using BaseClass::dynamicalSystemEnclosure;
  using BaseClass::getStep;

  /// This operator computes image of the set (in given representation) using set.move function, see capd/dynsys/Move.h for details
  /// This template together with SetTraits prevent usage of various types of jets with incompatible solvers.
  /// The user will get an exception at runtime with clear message instead of unreadable compiler error.
  template <typename SetType>
  void operator()(SetType& set){
    DiffInclusionSetMove<DiffInclusionCW,SetType>::move(set,*this);
  }

  /// Computes image of the set (in set's representation) and stores it in the result set.
   /// @param[in]  set       C^0 or C^1  set representing initial conditions
   /// @param[out] result    on return contains image of the set
  template <typename SetType>
  void operator()(SetType& set, SetType& result){
//    this->saveCurrentSet(set);
	  DiffInclusionSetMove<DiffInclusionCW,SetType>::move(set, result, *this);
  }

protected:
  using BaseClass::m_norm;
  using BaseClass::m_diffIncl;

  ScalarType m_errorTolerance;   ///< error tolerance in the exp(A) computations

};



}} // namespace capd::diffIncl

#endif // _CAPD_DIFFINCL_DIFFINCLUSIONCW_H_

/// @}
