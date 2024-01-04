
/////////////////////////////////////////////////////////////////////////////
/// @file AbstractSection.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_POINCARE_ABSTRACT_SECTION_H_
#define _CAPD_POINCARE_ABSTRACT_SECTION_H_

#include <string>
#include "SectionDerivativesEnclosure.h"
#include "capd/dynset/AbstractSet.h"
#include "capd/diffAlgebra/Hessian.h"
#include "capd/diffAlgebra/Jet.h"

namespace capd{
/// This namespace contains classes to compute Poincare Maps and Time Maps
namespace poincare{
/// @addtogroup poincare 
/// @{
/**
 *  PoicareSection is class that is default PoincareSection for class PoincareMap.
 *
 *  Internally it wraps a functional object \f$ \alpha: R^n\to R \f$. The Poincare section is defined as
 *  \f$ \Pi := \alpha^{-1}(0) \f$.
 *
 *  This class provides methods for checking the sign of \f$ \alpha \f$ evaluated at a point or a set.
 *  Efficient implementation of \f$ \alpha(set) \f$ depends on representation of the set.
 *
 *  This class provides also methods for recomputation of partial derivatives of the flow to partial derivatives of Poincare map.
 *
 *  Let
 *    - \f$ \varphi(t, x): R \times R^n -> R^n \f$ be a dynamical system generated
 *      by given vector field.
 *    - \f$ s: R^n \to R \f$ be a section function,
 *    - \f$ S = \{x \in R^n : s(x) = 0\} \f$
 *    - \f$ P: S \to S \f$ be a Poincare Map
 *    - for given point \f$ x \in S \f$ let T(x) be first return time (in given crossing direction)
 *      i.e.  \f$ P(x) = \varphi(T(x), x) \in S \f$
 *
 *  In the following we denote by
 *    - dP the derivative of Poincare Map P  : \f$ dP(x) = \frac{\partial P(x)}{\partial x} \f$
 *    - dT the derivative of T(x) : \f$ dT(x)  = \frac{\partial T(x)}{\partial x} \f$
 *    - dF the derivative of the flow : \f$ dF(x) = \frac{\partial \varphi(T(x), x)}{\partial x} \f$
 *  Then
 *    \f$ dP = dF + \frac{\partial \varphi}{\partial t} dT \f$
 */

template<class MatrixT>
class AbstractSection
{
public:
  virtual ~AbstractSection() = default;
  typedef MatrixT MatrixType;
  typedef typename MatrixType::RowVectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename VectorType::size_type size_type;             ///< integral type used to index containers (vectors, matrices, etc)
  typedef capd::dynset::AbstractSet<VectorType> Set;   ///< type of abstract base class for all sets
  typedef capd::diffAlgebra::Hessian<ScalarType,VectorType::csDim,VectorType::csDim> HessianType;
  typedef capd::diffAlgebra::Jet<MatrixT,0> JetType;
  typedef SectionDerivativesEnclosure<MatrixType> SectionDerivativesEnclosureType;

  virtual ScalarType operator() (const VectorType& v) const = 0;     ///< evaluates function at a given vector
  virtual VectorType gradient(const VectorType& u) const = 0;        ///< returns gradient of the function computed at vector u
  virtual ScalarType gradientByVector(const VectorType& x, const VectorType& u) const {
    return this->gradient(x)*u;
  }
  /** This is very important function.
      If it returns true, class PoincareMap delegates computation of value of section(set) to the section.
      Otherwise it is assumed that the set has more information to compute value of section(set) in most optimal way.
      This is quite natural as the set knows its own representation.
      This function returns true for instance if the PoincareSection is given by x_i=c, where x_i is i-th coordinate and c is constant.
  */
  virtual bool isSpecialSection() const{
    return false;
  }

  /** This function computes value of section function on a given set.*/
  virtual ScalarType evalAt(const Set& set) const = 0;

  /** computes gradient of return time
      @param[in] derivativeOfFlow - solution to first variational equation computed at return time
      @param[in] gradientOnPx - gradient of function that defines Poincare section evaluated at Px
      @param[in] denominator - scalar product of vector field evaluated at Px and gradientOnPx
      @param[out] result - computed gradient of return time
  */
  virtual void computeDT(
          const MatrixType& derivativeOfFlow,
          const VectorType& gradientOnPx,
          const ScalarType& denominator,
          VectorType& result
     ) const;

  /** Simultaneous computation of gradient of return time and derivative of Poincare Map dP.
      @param[in] Px - value of Poincare map
      @param[in] derivativeOfFlow - solution to first variational equation computed at return time
      @param[out] dT - computed gradient of return time
      @param[in] returnTime - enclosure of return time
      @note returnTime must be specified for nonautonomous flows only. Otherwise, default value 0 is valid.
  */
  virtual MatrixType computeDP(
          const VectorType& Px,
          const MatrixType& derivativeOfFlow,
          const VectorType& fieldOnPx,
          VectorType& dT
     ) const;

  /** Simultaneous computation of first and second Taylor coefficients of return time and Poincare map.
      @param[in] Px - value of Poincare map
      @param[in] derivativeOfFlow - solution to first variational equation computed at return time
      @param[in] hessianOfFlow - solution to first variational equation computed at return time
      @param[in] fieldOnPx - vector field evaluated at (t(Px),Px)
      @param[out] DP - computed derivative of Poincare map
      @param[out] D2P - computed second order Taylor coefficients of Poincare map
      @param[out] dT - computed gradient of return time
      @param[out] d2T - computed second order Taylor coefficients of return time
      @param[in] returnTime - return time to the section
      @note all input and output parameters are Taylor coefficients, not derivatives!
  */
  virtual void computeDP(
          const VectorType& Px,
          const MatrixType& derivativeOfFlow,
          const HessianType& hessianOfFlow,
          const VectorType& fieldOnPx,
          const VectorType& d2Phidt2,
          const MatrixType& derOfVectorFieldOnPx,
          MatrixType& DP,
          HessianType& D2P,
          VectorType& dT,
          MatrixType& d2T
      ) const;
/*
  /// This abstract function computes higher order jet of return time. The implementation can be highly optimized for affine sections and coordinate sections. Therefore the method is abstract.
  /// @param[in] Px - jet of flow, i.e \f$ \varphi(t,x) \f$ evaluated at t=return time
  /// @param[in] vfOnPx - jet of vector field evaluated at Px, i.e. \f$ (vectorField\circ \varpsi)(t,x) \f$ evaluated at t=return time
  /// @param[out] dT - computed jet of implicit function \f$ t(x) \f$
  virtual void computeDT(
          const JetType& Px,
          const JetType& vfOnPx,
          JetType& dT,
          size_type degree
      ) const = 0;

  /// This function computes higher order jet of Poincare map
  /// @param[in] Px - jet of flow, i.e \f$ \varphi(t,x) \f$ evaluated at t=return time
  /// @param[in] vfOnPx - jet of vector field evaluated at Px, i.e. \f$ (vectorField\circ \varpsi)(t,x) \f$ evaluated at t=return time
  /// @param[in] dT - jet of return time \f$ t(x) \f$
  /// @param[out] DP - computed jet of Poincare map, i.e \f$ P(x) = \varphi(t(x),x) \f$
  virtual void computeDP(
          const JetType& Px,
          const JetType& vfOnPx,
          const JetType& dT,
          JetType& DP,
          const size_type degree
      ) const;
*/
}; // end of template AbstractSection
/// @}
}} // namespace capd::poincare

#endif  /* _CAPD_POINCARE_ABSTRACT_SECTION_H_ */


