/////////////////////////////////////////////////////////////////////////////
///
/// @file C0HOTripletonSet.h
///
/// @author Daniel Wilczak
/// Created on: Aug 18, 2013
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2013 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0HOTRIPLETONSET_H_
#define _CAPD_DYNSET_C0HOTRIPLETONSET_H_

#include "capd/dynset/C0TripletonSet.h"

namespace capd{
namespace dynset{
// @addtogroup capd
/// @{

/**
 * Class \b C0HOTripletonSet \b represents a subset of \f$ R^n \f$ in the following form
 *
 * x + C*r0 + intersection(B*r , Q*q )
 *
 * where
 *  x is a point vector
 *  C,B,Q are point matrices, where B and Q are invertible and Q is close to orthogonal
 *  r0,r are interval vectors
 *
 *  Moreover it stores rigorous inverse matrices of B and Q
 *
 * The evaluation of the set by an ODE is realized by intersection of two methods:
 * the Taylor method and the Hermite-Obreshkov method.
 *
 * IMPORTANT: present implementation is valid only for orders of the Taylor method less or equal 32.
 * This is due to capacity of integer type used to store binomial coefficients.
 */
template <class MatrixT, class Policies>
class C0HOTripletonSet : public C0TripletonSet<MatrixT,Policies>{

public :
  typedef C0TripletonSet<MatrixT,Policies> BaseSet;
  typedef typename BaseSet::MatrixType MatrixType;
  typedef typename BaseSet::VectorType VectorType;
  typedef typename MatrixType::ScalarType ScalarType;
  typedef typename MatrixType::size_type size_type;
  typedef typename BaseSet::DynSysType DynSysType;

  C0HOTripletonSet(BaseSet & set);
  explicit C0HOTripletonSet(const VectorType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOTripletonSet(const VectorType& x, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOTripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOTripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t = TypeTraits<ScalarType>::zero());
  C0HOTripletonSet(const VectorType& x, const MatrixType& C,
      const VectorType& r0, const MatrixType& B,
      const VectorType& r,
      ScalarType t = TypeTraits<ScalarType>::zero()
    );

  /// This method computes value of functor f at interval vector represented by this set.
  /// It computes the value as an intersection of evalAt for two representations of this set: predictor and corrector.
  template<class Functional>
  ScalarType evalAt(const Functional& f) const {
    ScalarType r;
    if(!intersection(BaseSet::evalAt(f),predictor.evalAt(f),r)){
      throw std::logic_error("C0HOTripletonSet::evalAt - empty intersection. Report this error to CAPD developers!");
    }
    return r;
  }

  virtual ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& u) const{
    ScalarType r;
    if(!intersection(BaseSet::evalAffineFunctional(gradient,u),predictor.evalAffineFunctional(gradient,u),r)){
      throw std::logic_error("C0HODoubletonSet::evalAffineFunctional - empty intersection. Report this error to CAPD developers!");
    }
    return r;
  }

  template<class Solver>
  void move(Solver& solver){ this->move(solver,*this);  }

  template<class Solver>
  void move(Solver& solver, C0HOTripletonSet& result);

  std::string name() const { return "C0HOTripletonSet"; }

  BaseSet predictor;
};

// -------------------------------------------


template <class MatrixT, class Policies>
template<class Solver>
void C0HOTripletonSet<MatrixT,Policies>::move(Solver& solver, C0HOTripletonSet& result)
{
  const size_type dimension = this->m_x.dimension();
  // predictor step - Taylor method
  VectorType y(dimension), deltaY(dimension), rem(dimension), enc(dimension);
  MatrixType jacPhi(dimension, dimension);
  MatrixType B(dimension, dimension), Q(dimension,dimension);
  MatrixType deltaC(dimension, dimension);

  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly

  VectorType deltaX = this->m_currentSet - this->predictor.get_x();
  solver.encloseC0Map(this->getCurrentTime(),this->predictor.get_x(), this->m_currentSet, y, rem, enc, jacPhi);

  this->predictor.move(y,deltaX,jacPhi,rem,result.predictor);

  // corrector step
  VectorType psiPlus(dimension), psiMinus(dimension);
  MatrixType JPlus(dimension,dimension);
  MatrixType JMinus(dimension,dimension);

  size_type order = solver.getOrder();
  size_type q = order/2;
  size_type p = order - q;
  capd::dynset::computePsi(solver.getCurve(),p,q,solver.getStep(),psiPlus,JPlus);

  // we have to recompute coefficients at the image
  split((VectorType)(result.predictor),y,deltaY);
  ScalarType nextTime =  this->getCurrentTime() + solver.getStep();
  solver.computePsiCoefficients(nextTime,y,(VectorType)(result.predictor),q);
  capd::dynset::computePsi(solver.getPsiCurve(),q,p,-solver.getStep(),psiMinus,JMinus);

  // solve implicit equation
  int sign = q%2 ? -1 : 1;
  for(size_type i=0;i<dimension;++i)
    rem[i] = psiPlus[i] - psiMinus[i] + sign*rem[i]/newton(p+q,q);

  MatrixType midJMinusInverse = midMatrix(capd::matrixAlgorithms::gaussInverseMatrix(midMatrix(JMinus)));
  jacPhi = midJMinusInverse*JPlus;

  MatrixType R = - midJMinusInverse*JMinus;
  // add Identity
  for(size_type i=1;i<=dimension;++i) R(i,i) += TypeTraits<ScalarType>::one();

  VectorType eps = midJMinusInverse*rem + R*deltaY;
  VectorType z = y+eps;

  // first enclosure of the set after time step from corrector phase
  result.m_currentSet = z + jacPhi*deltaX;

  // second enclosure of the set after time step from the corrector phase and previous representation
  result.m_C = jacPhi*this->m_C;
  Q = result.m_Q = jacPhi*this->m_Q;
  B = result.m_B = jacPhi*this->m_B;
  intersection(result.m_currentSet,z + result.m_C*this->m_r0 + intersection(Q*this->m_q,B*this->m_r),result.m_currentSet);

  result.predictor.setCurrentTime(nextTime);
  result.predictor.reorganizeIfNeeded(result.predictor);

  if(subset((VectorType)result.predictor,result.m_currentSet)){
    result.BaseSet::operator=(result.predictor);
    result.setLastEnclosure(enc);
    return;
  }

  VectorType X = intersection((VectorType)result.predictor,result.m_currentSet);

  // computing representation for the next step
  split(result.m_C, deltaC);
  result.m_x = midVector(X);
  eps += y-result.m_x;
  VectorType deltaCr0 = deltaC*this->m_r0;
  rem = deltaCr0 + eps;
  this->Policies::computeBinvB(result.m_Q,result.m_invQ,this->m_q+rem);

  // we must enclose the following quantity
  // newq = (invQ*Q)*q + invQ*rem = (invQ*Q)*q + (invQ*deltaCr0) + (invQ*eps)
  result.m_q = (result.m_invQ * Q) * this->m_q
      + intersection(
          result.m_invQ*rem,
          (result.m_invQ*deltaC)*this->m_r0 + result.m_invQ*eps
        );

      // propagate by inverse matrix
  try{
    // Q is unnecessary now
    result.m_B = midMatrix(B);
    result.m_invB = capd::matrixAlgorithms::krawczykInverse(result.m_B);

    // we must enclose the following quantity
    // newR = (invB*B)*r + invB*rem = (invB*B)*r + (invB*deltaCr0) + (invB*eps)
    result.m_r = (result.m_invB * B) * this->m_r
        + intersection(
            result.m_invB*rem,
            (result.m_invB*deltaC)*this->m_r0 + result.m_invB*eps
          );

  }catch(...)
  {
    result.m_B = result.m_Q;
    result.m_r = result.m_q;
    result.m_invB = result.m_invQ;
  }

  if(&result != this)
    result.m_r0 = this->m_r0;

  this->reorganizeIfNeeded((BaseSet&)result);
  result.setCurrentTime(nextTime);
  result.setLastEnclosure(enc);

  if(subset(result.m_currentSet,(VectorType)result.predictor)){
    result.predictor=((BaseSet)result);
  } else{
    result.m_currentSet = X;
    result.predictor.setCurrentSet(X);
    if(!subset(result.predictor.get_x(),X))
    {
      result.predictor=((BaseSet)result);
    }
  }
}


/// @}
}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0HOTRIPLETONSET_H_

