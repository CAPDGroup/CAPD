

/////////////////////////////////////////////////////////////////////////////
/// @file HighOrderEnclosure.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2017 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSYS_HIGH_ORDER_ENCLOSURE_CLASS_H_
#define _CAPD_DYNSYS_HIGH_ORDER_ENCLOSURE_CLASS_H_
#include "capd/intervals/Interval.hpp"
#include "capd/diffAlgebra/C1TimeJet.h"
#include "capd/diffAlgebra/C2TimeJet.h"
#include "capd/dynsys/SolverException.h"
#include "capd/dynsys/StepControl.h"
#include "capd/poincare/SaveStepControl.h"

namespace capd{
namespace dynsys{
/// @addtogroup dynsys
/// @{


/**
 * This file defines class for computation of [C0-C2] rough enclosure to
 * an ODE by high order Taylor method.
 * It is designed to be a template parameter for all rigorous solver classes.
 *
 */
class HighOrderEnclosure{
private:
///  auxiliary functions

/**
 * This function computes enclosure for the Taylor series up to some order of \varphi(t,x).
 * Then the coefficients for the remainder part are recomputed if they at least one of them is equal to zero.
 * This means that this is the first step of integration.
 */

//#define ENC_POLY

#ifdef ENC_POLY

  template<class DS>
  static void encloseRangeOfPolynomial(Object* coeffs, const Scalar& h, Object& result, int n)
  {
    capd::vectalg::evalPolynomial(coeffs,h,result,n);
  }

#else

  template<class Object, class Scalar>
  static void encloseRangeOfPolynomial(Object* coeffs, const Scalar& h, Object& result, int n)
  {
    typename Object::iterator b = result.begin(), e = result.end();
    for(int i=0; b!=e;++b,++i)
    {
      // compute derivative of polynomial
      // if monotone on [0,h] - take interval hull of P(0) and P(h)
      // otherwise - compute P([0,h])
      Scalar der = *(coeffs[n].begin()+i) * n;
      for(int j=n-1;j>0;--j)
        der = h*der + *(coeffs[j].begin()+i) * j;
      Scalar hh = isSingular(der) ? h : h.right();

      *b = *(coeffs[n].begin()+i);
      for(int j=n-1;j>=0;--j)
        *b = *b * hh + *(coeffs[j].begin()+i);
      *b = intervalHull(*b,*(coeffs[0].begin()+i));
    }
  }

#endif

  template<class It>
  static bool needToRecomputeRemainder(It b, It e, It i){
    for(; b!=e; ++b, ++i)
      if(b->leftBound() < i->leftBound() or b->rightBound() > i->rightBound())
        return true;
    return false;
  }

  template<class DS>
  static bool encloseRangeOfPolynomial(DS& ds,typename DS::VectorType& phi,const typename DS::ScalarType& h)
  {
    encloseRangeOfPolynomial(ds.getCoefficients(),h,phi,ds.getOrder());
    return needToRecomputeRemainder(phi.begin(),phi.end(),ds.beginRemainderCoefficients(0));
  }

  template<class DS>
  static bool encloseRangeOfPolynomial(
    DS& ds,
    capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& phi,
    const typename DS::ScalarType& h
    )
  {
    encloseRangeOfPolynomial(ds.getCoefficients(),h,phi.vector(),ds.getOrder());
    encloseRangeOfPolynomial(ds.getMatrixCoefficients(),h,phi.matrix(),ds.getOrder());
    return
        needToRecomputeRemainder(phi.vector().begin(),phi.vector().end(),ds.beginRemainderCoefficients(0)) or
        needToRecomputeRemainder(phi.matrix().begin(),phi.matrix().end(),ds.beginMatrixRemainderCoefficients(0));
  }

/**
 * This function computes enclosure for the Taylor series up to some order of \varphi(t,x).
 * Then the coefficients for the remainder part are recomputed if they at least one of them is equal to zero.
 * This usually means that this is the first step of integration.
 */
  template<class DS>
  static bool encloseRangeOfPolynomial(
      DS& ds,
      capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& phi,
      const typename DS::ScalarType& h
      )
  {
    encloseRangeOfPolynomial(ds.getCoefficients(),h,phi.vector(),ds.getOrder());
    encloseRangeOfPolynomial(ds.getMatrixCoefficients(),h,phi.matrix(),ds.getOrder());
    encloseRangeOfPolynomial(ds.getHessianCoefficients(),h,phi.hessian(),ds.getOrder());
    return
        needToRecomputeRemainder(phi.vector().begin(),phi.vector().end(),ds.beginRemainderCoefficients(0)) or
        needToRecomputeRemainder(phi.matrix().begin(),phi.matrix().end(),ds.beginMatrixRemainderCoefficients(0)) or
        needToRecomputeRemainder(phi.hessian().begin(),phi.hessian().end(),ds.beginHessianRemainderCoefficients(0));
  }

  template<class DS>
  static void highOrderEnclosureError(DS & ds, const typename DS::ScalarType& t, const char message[])
  {
    typename DS::VectorType problematicInitialCondition(ds.dimension());
    for(typename DS::size_type i=0;i<ds.dimension();++i)
      problematicInitialCondition[i] = ds.coefficient(i,0);
    throw SolverException<typename DS::VectorType>(message,t, problematicInitialCondition,ds.getStep());
  }

// ------------- computeRemainderCoefficients --------------------------------
  template<class DS>
  static inline void computeRemainderCoefficients(const typename DS::ScalarType& t, DS& ds, typename DS::VectorType& x){
    ds.computeRemainderCoefficients(t,x);
  }

  template<class DS>
  static inline void computeRemainderCoefficients(const typename DS::ScalarType& t, DS& ds, capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& x){
    ds.computeRemainderCoefficients(t,x.vector(),x.matrix());
  }

  template<class DS>
  static inline void computeRemainderCoefficients(const typename DS::ScalarType& t, DS& ds, capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& x){
    ds.computeRemainderCoefficients(t,x.vector(),x.matrix(),x.hessian());
  }

  // ----------------- predictNextEnclosure ------------------------------------
  template<class Scalar, class It1, class It2, class It3>
  static void predictNextEnclosure(
      const Scalar& stepToOrder,
      It1 phi, It1 phiEnd,
      It2 remCoeff,
      It3 out_enc, It3 out_remEnc)
  {
    const static typename Scalar::BoundType epsilon = 1.e-300;
    const static Scalar mulFactor(-2.,2.);
    for( ;phi!=phiEnd;++phi,++remCoeff,++out_enc,++out_remEnc)
    {
      *out_remEnc = mulFactor*(stepToOrder*(*remCoeff)+epsilon);
      *out_enc = *phi + *out_remEnc;
    }
  }

  template<class DS>
  static inline void predictNextEnclosure(
      DS& ds,
      const typename DS::ScalarType& stepToOrder,
      const typename DS::VectorType& phi,
      typename DS::VectorType& out_enc,
      typename DS::VectorType& out_remEnc
      )
  {
    predictNextEnclosure(stepToOrder,phi.begin(),phi.end(),ds.beginRemainderCoefficients(ds.getOrder()+1),out_enc.begin(),out_remEnc.begin());
  }

  template<class DS>
  static inline void predictNextEnclosure(
        DS& ds,
        const typename DS::ScalarType& stepToOrder,
        const capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& phi,
        capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_enc,
        capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_remEnc
        )
  {
    predictNextEnclosure(stepToOrder, phi.vector().begin(), phi.vector().end(), ds.beginRemainderCoefficients(ds.getOrder()+1), out_enc.vector().begin(), out_remEnc.vector().begin() );
    predictNextEnclosure(stepToOrder, phi.matrix().begin(), phi.matrix().end(), ds.beginMatrixRemainderCoefficients(ds.getOrder()+1),out_enc.matrix().begin(),out_remEnc.matrix().begin());
  }

  template<class DS>
  static inline void predictNextEnclosure(
        DS& ds,
        const typename DS::ScalarType& stepToOrder,
        const capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& phi,
        capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_enc,
        capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_remEnc
        )
  {
    predictNextEnclosure(stepToOrder,phi.vector().begin(),phi.vector().end(),ds.beginRemainderCoefficients(ds.getOrder()+1),out_enc.vector().begin(),out_remEnc.vector().begin());
    predictNextEnclosure(stepToOrder,phi.matrix().begin(),phi.matrix().end(),ds.beginMatrixRemainderCoefficients(ds.getOrder()+1),out_enc.matrix().begin(),out_remEnc.matrix().begin());
    predictNextEnclosure(stepToOrder,phi.hessian().begin(),phi.hessian().end(),ds.beginHessianRemainderCoefficients(ds.getOrder()+1),out_enc.hessian().begin(),out_remEnc.hessian().begin());
  }

  // ----------------- checkInclusion ------------------------------------
  template<class Scalar, class It1, class It2, class It3>
  static bool checkInclusion(
      const Scalar& stepToOrder,
      It1 remEnc, It1 remEncEnd,
      It2 remCoeff, It3 out_rem,
      Scalar& out_factor
    )
  {
    bool isSubset = true;
    for( ;remEnc!=remEncEnd;++remEnc,++remCoeff,++out_rem)
    {
      *out_rem = *remCoeff * stepToOrder;
      Scalar v = capd::max(capd::abs(out_rem->leftBound()),capd::abs(out_rem->rightBound()));
      if(!isSingular(v))
        out_factor = capd::min(out_factor,remEnc->right()/v);
      if(!subsetInterior(*out_rem,*remEnc))
        isSubset = false;
    }
    return isSubset;
  }

  template<class DS>
  static inline bool checkInclusion(
      DS& ds,
      const typename DS::VectorType& remEnc,
      const typename DS::ScalarType& stepToOrder,
      typename DS::VectorType& out_rem,
      typename DS::ScalarType& out_factor
      )
  {
    return checkInclusion(stepToOrder,remEnc.begin(),remEnc.end(),ds.beginRemainderCoefficients(ds.getOrder()+1),out_rem.begin(),out_factor);
  }

  template<class DS>
  static inline bool checkInclusion(
      DS& ds,
      const capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& remEnc,
      const typename DS::ScalarType& stepToOrder,
      capd::diffAlgebra::C1TimeJet<typename DS::MatrixType>& out_rem,
      typename DS::ScalarType& out_factor
      )
  {
    bool r0 = checkInclusion(stepToOrder,remEnc.vector().begin(),remEnc.vector().end(),ds.beginRemainderCoefficients(ds.getOrder()+1),out_rem.vector().begin(),out_factor);
    bool r1 = checkInclusion(stepToOrder,remEnc.matrix().begin(),remEnc.matrix().end(),ds.beginMatrixRemainderCoefficients(ds.getOrder()+1),out_rem.matrix().begin(),out_factor);
    return r0 and r1;
  }

  template<class DS>
  static inline bool checkInclusion(
      DS& ds,
      const capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& remEnc,
      const typename DS::ScalarType& stepToOrder,
      capd::diffAlgebra::C2TimeJet<typename DS::MatrixType>& out_rem,
      typename DS::ScalarType& out_factor
      )
  {
    bool r0 = checkInclusion(stepToOrder,remEnc.vector().begin(),remEnc.vector().end(),ds.beginRemainderCoefficients(ds.getOrder()+1),out_rem.vector().begin(),out_factor);
    bool r1 = checkInclusion(stepToOrder,remEnc.matrix().begin(),remEnc.matrix().end(),ds.beginMatrixRemainderCoefficients(ds.getOrder()+1),out_rem.matrix().begin(),out_factor);
    bool r2 = checkInclusion(stepToOrder,remEnc.hessian().begin(),remEnc.hessian().end(),ds.beginHessianRemainderCoefficients(ds.getOrder()+1),out_rem.hessian().begin(),out_factor);
    return r0 and r1 and r2;
  }

  // ----------------- refineEnclosure ------------------------------------

  template<class It1, class It2, class S>
  static void refineEnclosure(It1 phi, It1 e, It1 rem, It2 out_remEnc, It2 out_enc, S factor )
  {
    for(;phi!=e;++phi,++rem,++out_remEnc,++out_enc){
      if(rem->subsetInterior(*out_remEnc))
        continue;
      *out_remEnc *= factor;
      *out_enc = *phi + *out_remEnc;
    }
  }

  template<class VectorType>
  static inline void refineEnclosure( const VectorType& phi, const VectorType& rem, VectorType& out_remEnc, VectorType& out_enc )
  {
    typename TypeTraits<typename VectorType::ScalarType>::Real f(1.5);
    refineEnclosure(phi.begin(),phi.end(),rem.begin(),out_remEnc.begin(),out_enc.begin(),f);
  }

  template<class MatrixType>
  static inline void refineEnclosure(
      const capd::diffAlgebra::C1TimeJet<MatrixType>& phi,
      const capd::diffAlgebra::C1TimeJet<MatrixType>& rem,
      capd::diffAlgebra::C1TimeJet<MatrixType>& out_remEnc,
      capd::diffAlgebra::C1TimeJet<MatrixType>& out_enc
    )
  {
    typename TypeTraits<typename MatrixType::ScalarType>::Real f(1.5);
    refineEnclosure(phi.vector().begin(),phi.vector().end(),rem.vector().begin(),out_remEnc.vector().begin(),out_enc.vector().begin(),f);
    refineEnclosure(phi.matrix().begin(),phi.matrix().end(),rem.matrix().begin(),out_remEnc.matrix().begin(),out_enc.matrix().begin(),f);
  }

  template<class MatrixType>
  static inline void refineEnclosure(
      const capd::diffAlgebra::C2TimeJet<MatrixType>& phi,
      const capd::diffAlgebra::C2TimeJet<MatrixType>& rem,
      capd::diffAlgebra::C2TimeJet<MatrixType>& out_remEnc,
      capd::diffAlgebra::C2TimeJet<MatrixType>& out_enc
    )
  {
    typename TypeTraits<typename MatrixType::ScalarType>::Real f(1.5);
    refineEnclosure(phi.vector().begin(),phi.vector().end(),rem.vector().begin(),out_remEnc.vector().begin(),out_enc.vector().begin(),f);
    refineEnclosure(phi.matrix().begin(),phi.matrix().end(),rem.matrix().begin(),out_remEnc.matrix().begin(),out_enc.matrix().begin(),f);
    refineEnclosure(phi.hessian().begin(),phi.hessian().end(),rem.hessian().begin(),out_remEnc.hessian().begin(),out_enc.hessian().begin(),f);
  }

public:
  template<typename DS>
  static typename DS::VectorType enclosure(
                        DS & ds,
                        typename DS::ScalarType const & currentTime,
                        typename DS::MatrixType::RowVectorType const & x
                        )
  {
    capd::poincare::SaveStepControl<DS> saveStepControl(ds);
    ds.turnOffStepControl();
    typename DS::MatrixType::RowVectorType enc = x, rem = x;
    computeEnclosureAndRemainder(ds,currentTime,x,enc,rem);
    return enc;
  }
/**
  * This is the only public method required for enclosure policy.
  *
  * Computes enclosure of solution of ODE during one time step  i.e \f$ \varphi([0,h],x) \f$
  * for some \f$ h\in[0,hTrial] \f$, where \f$ hTrial \f$ is the predicted step
  * or for \f$ h=hTrial \f$ if changing of the time step is not allowed.
  *
  * If CoeffType is a C1Coff structure then the function computes enclosure
  * and the Lagrange remainder for the solutions to variational equations as well.
  *
  * @param ds             dynamical system
  * @param out_remainder  when success, contains the Lagrange remainder computed at enclosure
  * @param out_enclosure  enclosure of the set of solution over a time step
  */
  template<class DS, class CoeffType>
  static void computeEnclosureAndRemainder(DS& ds, const typename DS::ScalarType& t, const typename DS::VectorType&, CoeffType& out_enclosure, CoeffType& out_remainder)
  {
    typedef typename DS::ScalarType ScalarType;
    typedef typename ScalarType::BoundType Real;
    const typename DS::size_type dimension = ds.dimension();

    const static ScalarType I(TypeTraits<Real>::zero(),TypeTraits<Real>::one());
    ScalarType h = I*ds.getStep();
    ScalarType stepToOrder = power(h,ds.getOrder()+1);

    CoeffType phi(dimension), remEnc(dimension);
    if(encloseRangeOfPolynomial(ds,phi,h))
      computeRemainderCoefficients(t+h,ds,phi);

    // Now we are predicting the enclosure based on the remainder from the previous approved step
    // (enc and remEnc are returned by this function)
    predictNextEnclosure(ds,stepToOrder,phi,out_enclosure,remEnc);
    // and compute remainder coefficients on predicted enclosure.
    computeRemainderCoefficients(t+h,ds,out_enclosure);

    ScalarType factor(1.);
    // here we check if computed remainder coefficients are subsets of predicted enclosure of remainder coefficients
    // (factor can be decreased in this function)
    // out_remainder is returned by this function
    if (checkInclusion(ds,remEnc,stepToOrder,out_remainder,factor) )
      return;

    // we do not have inclusion, but maybe we can adjust the time step so that it will hold

    if(ds.isStepChangeAllowed())
    {
      if(isSingular(factor)
          or isInf(factor.leftBound()) or isInf(factor.rightBound())
          or isNaN(factor.leftBound()) or isNaN(factor.rightBound())
      ) highOrderEnclosureError(ds,t,"High Order Enclosure Error: cannot adjust time step. Cannot integrate.");

      factor = clearMantissaBits(exp(log(factor)/(ds.getOrder()+1)).leftBound());
      typename capd::TypeTraits<ScalarType>::Real newStep = leftBound(ds.getStep()*factor);

      if(capd::abs(newStep) < ds.getStepControl().getMinStepAllowed())
        highOrderEnclosureError(ds,t,"High Order Enclosure Error: minimal time step reached. Cannot integrate.");

      // success - we have enclosure with new (smaller) time step
      ds.adjustTimeStep(ScalarType(newStep));
      return;
    }

    // time step change is not allowed. We have to refine enclosure and recompute remainder
    int maxAttempts = 30;
    while(maxAttempts)
    {
      // In this function remEnc is enlarged and enc is set to phi + remEnc
      refineEnclosure(phi,out_remainder,remEnc,out_enclosure);
      computeRemainderCoefficients(t+h,ds,out_enclosure);

      // check the inclusion
      ScalarType factor(1.);
      if (checkInclusion(ds,remEnc,stepToOrder,out_remainder,factor) )
        return;
      maxAttempts--;
    } // endwhile

    highOrderEnclosureError(ds,t,"High Order Enclosure Error: cannot find enclosure guaranteeing bounds, loop limit exceeded");
  }
};

/// @}
}} // namespace capd::dynsys

#endif /* _CAPD_DYNSYS_HIGH_ORDER_ENCLOSURE_H_ */

