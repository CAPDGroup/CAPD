

/////////////////////////////////////////////////////////////////////////////
/// @file DynSys.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_DYNSYS_H_
#define _CAPD_DYNSYS_DYNSYS_H_

#include "capd/vectalg/Norm.h"
#include "capd/basicalg/TypeTraits.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

template<typename MatrixType>
class FlowballSet;

/**
 *  Class dynsys is an abstract class representing a discrete dynamical system.
 *  It is assumed that the generator Psi of the dynamical system decomposes into
 *                                Psi=Phi+Omega
 *  where Phi is a polynomial and Omega is the remainder.
 *  A realization of the class must implement three methods:
 *     Phi - the interval enclosure of the polynomial part
 *     JacPhi -  the interval enclosure of the jacobian of Phi
 *     Remainder - an enclosure of the remainder Omega.
 *  Although the class represents a discrete dynamical system, it may be also
 *  used with time t translations along trajectories of an ODE.
 *  This is done in the class OdeNum inheriting from DynSys and defined in the
 *  sequel.
 *
 */
template<typename MatrixT>
class DynSys{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef capd::vectalg::Norm<VectorType,MatrixType> NormType;

  /// Computes value of function (usually numerical scheme for an ODE) at time t and point iv
  virtual VectorType Phi(const ScalarType& t, const VectorType &iv) = 0;

  /// Computes derivative of function (usually numerical scheme for an ODE) at time t and point iv
  virtual MatrixType JacPhi(const ScalarType& t, const VectorType &iv) = 0;

  /// Computes and returns bound for local error of a function (for instance if Phi is finite Taylor series of a map then this method computes bound for Lagrange remainder).
  /// If DynSys is an ODE, then out_enc contains enclosure of trajectories over the time step.
  /// If the function cannot validate existence of solutions to ODE over the time step, out_enc might be in an inconsistent state.
  virtual VectorType Remainder(const ScalarType& t, const VectorType &iv, VectorType &out_enc) = 0;

  virtual ScalarType Lipschitz(const ScalarType& t, const VectorType &iv, NormType &n);

  /// Used for ODEs. It verifies the existence of solutions to IVP at time t and set of initial conditions x over the time step.
  /// If the function succeeds, a rigorous bound for the trajectories is returned.
  /// Otherwise, an exception is thrown.
  virtual VectorType enclosure(const ScalarType& t, const VectorType& x) = 0;

  /// For given set xx, time t and a point x from the set xx
  /// It simultaneously computes and returns enclosures for:
  /// - numerical scheme Phi(t,x),
  /// - derivative of numerical scheme D_x Phi(t,xx)
  /// - error of numerical scheme Remainder(t,xx)
  /// - an enclosure for the trajectories over the time step (ODEs only)
  /// If the function cannot compute any of the output results, an exception is thrown and the
  /// output parameters o_* might be in inconsistent state.
  virtual void encloseC0Map(
      const ScalarType& t,
      const VectorType& x,
      const VectorType& xx,
      VectorType& o_phi,
      VectorType& o_rem,
      VectorType& o_enc,
      MatrixType& o_jacPhi
  ) = 0;

  /// Returns time step of the dynamical system. By default it returns one - time step for discrete DS (maps).
  /// Shall be overridden in classes that implement numerical schemes for ODEs.
  virtual ScalarType getStep() const {
    return TypeTraits<ScalarType>::one();
  }

  virtual ~DynSys(){}
protected:
};

template<typename MatrixT>
class HOSolver{
public:
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;

  virtual void computePsiCoefficients(ScalarType t, const VectorType& x, const VectorType& xx, size_type order) = 0;

  virtual ~HOSolver(){}
};
/// @}
}} // namespace capd::dynsys

#endif // _CAPD_DYNSYS_DYNSYS_H_


