

/////////////////////////////////////////////////////////////////////////////
/// @file FirstOrderEnclosure.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

/* Author: Daniel Wilczak 2001-2017 */

#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/diffAlgebra/C2TimeJet.h"

#ifndef _CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_H_
#define _CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_H_

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{

class FirstOrderEnclosure{
public:
/**
  * Computes enclosure of solution of ODE during one time step  i.e \f$ \varphi([0,step],x) \f$
  *
  * @param vectorField vector field
  * @param currentTime (for nonautonomous ODEs)
  * @param x0          initial condition
  * @param timeStep    final time
  * @return  enclosure of solution of ODE
  */
template<typename MapType>
static typename MapType::VectorType enclosure(
                        MapType & vectorField,
                        typename MapType::ScalarType const & currentTime,
                        typename MapType::MatrixType::RowVectorType const & x,
                        typename MapType::ScalarType const & step
                        );

template<typename DS>
static typename DS::VectorType enclosure(
                        DS & ds,
                        typename DS::ScalarType const & currentTime,
                        typename DS::MatrixType::RowVectorType const & x
                        )
{
  return enclosure(ds.getVectorField(),currentTime,x,ds.getStep());
}

//###########################################################//

/**
   * Finds enclosure for Jacobian matrix (variational part) for whole time step
   *
   * @param vectorField  vector field
   * @param currentTime (for nonautonomous ODEs)
   * @param timeStep     time step
   * @param enclosure    enclosure of solution of ODE during whole time step
   * @param the_norm     logarithmic norm used to bound solution
   * @return matrix containing enclosure of Jacobian
   */
// source- "C^1-Lohner algorithm" by P. Zgliczynski
template<typename MapType, typename NormType>
static typename MapType::MatrixType jacEnclosure(
                        const MapType& vectorField,
                        const typename MapType::ScalarType& currentTime,
                        const typename MapType::ScalarType& step,
                        const typename MapType::VectorType& enc,
                        const NormType &the_norm,
                        typename MapType::ScalarType* o_logNormOfDerivative =0
                        );

//###########################################################//

/**
   * Finds enclosure for second order variational equations for whole time step
   *
   * @param vectorField         vector field
   * @param timeStep            time step
   * @param enclosure           enclosure of solution of ODE during whole time step
   * @param[out] jacEnclosure   computed enclosure of solution to first order variational equations during whole time step
   * @param[out] hessEnclosure  computed enclosure of solution to second order variational equations during whole time step
   * @return                    computed logarithmic norm of derivative
   */
// source- "C^r-Lohner algorithm" by D. Wilczak and P. Zgliczynski
template<typename MapType>
static typename MapType::ScalarType c2Enclosure(
          const MapType& vectorField,
          const typename MapType::ScalarType& step,
          const typename MapType::VectorType& enc,
          typename MapType::MatrixType& jacEnclosure,
          typename MapType::HessianType& hessEnclosure
      );

template<class DS>
inline static void computeEnclosureAndRemainder(DS& ds, const typename DS::ScalarType& t, const typename DS::VectorType& x, typename DS::VectorType& out_enc, typename DS::VectorType& out_rem){
  typedef typename DS::ScalarType Scalar;
  const static Scalar I(TypeTraits<Scalar>::zero().leftBound(),TypeTraits<Scalar>::one().rightBound());

  // we compute an enclosure for \varphi([0,timestep],iv)
  ds.setCurrentTime(t);
  out_enc = enclosure(ds.getVectorField(),t,x,ds.getStep());
  ds.computeRemainderCoefficients(t + I*ds.getStep(),out_enc);
  out_rem =  ds.getRemainderCoefficients()[ds.getOrder()+1]*power(ds.getStep(),ds.getOrder()+1);
}

template<class DS>
inline static void computeEnclosureAndRemainder(DS& ds, const typename DS::ScalarType& t, const typename DS::VectorType& x, capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_enc, capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_rem){
  typedef typename DS::ScalarType Scalar;
  const static Scalar I(TypeTraits<Scalar>::zero().leftBound(),TypeTraits<Scalar>::one().rightBound());

  // we compute an enclosure for \varphi([0,timestep],iv)
  out_enc.vector() = enclosure(ds.getVectorField(),t,x,ds.getStep());
  out_enc.matrix() = jacEnclosure(ds.getVectorField(),t,ds.getStep(),out_enc.vector(),capd::vectalg::EuclLNorm<typename DS::VectorType,typename DS::MatrixType>());

  Scalar fac = power(ds.getStep(),ds.getOrder()+1);
  ds.computeRemainderCoefficients(t+I*ds.getStep(),out_enc.vector(),out_enc.matrix());
  out_rem.vector()  = fac * ds.getRemainderCoefficients()[ds.getOrder()+1];
  out_rem.matrix()  = fac * ds.getMatrixRemainderCoefficients()[ds.getOrder()+1];
}

template<class DS>
inline static void computeEnclosureAndRemainder(DS& ds, const typename DS::ScalarType& t, const typename DS::VectorType& x, capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_enc, capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_rem){
  typedef typename DS::ScalarType Scalar;
  const static Scalar I(TypeTraits<Scalar>::zero().leftBound(),TypeTraits<Scalar>::one().rightBound());

  // we compute an enclosure for \varphi([0,timestep],iv)
  out_enc.vector() = enclosure(ds.getVectorField(),t,x,ds.getStep());
  c2Enclosure(ds.getVectorField(),ds.getStep(),out_enc.vector(),out_enc.matrix(),out_enc.hessian());

  Scalar fac = power(ds.getStep(),ds.getOrder()+1);
  ds.computeRemainderCoefficients(t+I*ds.getStep(),out_enc.vector(),out_enc.matrix(),out_enc.hessian());
  out_rem.vector()  = fac * ds.getRemainderCoefficients()[ds.getOrder()+1];
  out_rem.matrix()  = fac * ds.getMatrixRemainderCoefficients()[ds.getOrder()+1];
  out_rem.hessian() = fac * ds.getHessianRemainderCoefficients()[ds.getOrder()+1];
}

}; // class FirstOrderEnclosure
/// @}
}} //namespace capd::dynsys

#endif // _CAPD_DYNSYS_FIRST_ORDER_ENCLOSURE_H_


