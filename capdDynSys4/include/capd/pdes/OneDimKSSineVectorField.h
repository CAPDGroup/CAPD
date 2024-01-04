/// @addtogroup pdes
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file OneDimKSSineVectorField.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2008-2016 by the CAPD Group.
//
// Distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_PDES_OneDimKSSineVectorField_H_
#define _CAPD_PDES_OneDimKSSineVectorField_H_

#include <vector>
#include <iostream>
#include <future>
#include "capd/basicalg/minmax.h"
#include "capd/basicalg/power.h"
#include "capd/intervals/lib.h"
#include "capd/vectalg/lib.h"
#include "capd/pdes/DissipativeVectorField.h"
#include "capd/pdes/GeometricBound.h"

namespace capd {
namespace pdes {
/**
 * The class implements vector field of the one-dimensional real Kuramoto-Shivashinsky PDE under the following assumptions
 * 1 .The solutions are represented in the Fourier basis
 * 2. We impose periodic and odd solutions, thus only sine components are present.
 *    The functions are represented as \f$ u(t,x) = -2 \sum_{k=1}^\infty a_k(t) \sin (ikx) \f$
 * 3. Domain is restricted to some subset of analytic functions. We impose geometric-like decay of Fourier coefficients.
 *    |a_k| < C/(q^k k^s), for some real constants q>1 and s
 *
 * The implementation provides
 * - evaluation of vector field on a representable set of analytic functions
 * - computation of partial derivatives of the vector field with respect to finite number of variables
 * - automatic differentiation for d/dt^i, i-natural number
 * - automatic differentiation for d/da_k dt^i - -natural number k-bounded
*/
class OneDimKSSineVectorField: public DissipativeVectorField<capd::pdes::GeometricBound<capd::interval> >{
public:
  typedef capd::interval ScalarType;
  typedef capd::pdes::GeometricBound<capd::interval> VectorType;
  typedef capd::IMatrix MatrixType;

  typedef capd::IVector::size_type size_type;
  typedef std::vector<VectorType> VectorArray;
  typedef std::vector<VectorArray> MatrixArray;
  typedef std::vector<ScalarType> ScalarArray;
  typedef std::pair<VectorArray*, VectorArray*> C1Data;

  /// constructs vector field of KS-equation.
  /// @param[in] nu - viscosity parameter in the KS-equation
  OneDimKSSineVectorField(ScalarType nu, size_type dim, size_type firstDissipativeVariable);

  VectorType operator()(ScalarType h, const VectorType& v) {
    VectorArray a(2);
    a[0] = v;
    a[1] = v;
    this->computeODECoefficients(a,1);
    return a[1];
  }

  VectorType operator()(ScalarType h, const VectorType& v, MatrixType& A) {
    VectorArray a(2);
    MatrixArray J(A.numberOfColumns(),VectorArray(2,VectorType(this->dimension())));
    this->blockDerivative(v,a,J,A);
    return a[1];
  }

  /// this method computes finite-dimensional square block of the derivative of the vector field
  MatrixType derivative(ScalarType h, const VectorType& v){
    MatrixType A(this->dimension(),this->dimension());
    this->operator()(h,v,A);
    return A;
  }

  /// computes Taylor coefficients for C^0 part
  void computeODECoefficients(VectorArray& a, size_type order);
  /// given coefficients for C^0 part, it computes Taylor coefficient for variational equation (one column)
  /// with initial condition 'c'
  void computeODECoefficients(const VectorArray& a, VectorArray& c, size_type order);

  void makeSelfConsistentBound(VectorArray& a) { makeSelfConsistentBound(a,this->nonlinearPart, a[0].getGeometricDecay()); }
  void makeSelfConsistentBound(VectorArray& a, MatrixArray& J1, MatrixArray& J2, size_type numberOfColumns);

  size_type dimension() const { return m_dimension; }
  size_type firstDissipativeIndex() const { return m_firstDissipativeVariable; }
  void updateTail(VectorType& x, const VectorArray& enc, ScalarType h) const {
    this->updateTail(x,nonlinearPart,enc,h);
  }

  void updateTail(VectorArray& DyxId, VectorArray& Dyx, const MatrixArray& Enc, const MatrixArray& DyxEnc, ScalarType h) const {
    for(size_type i=0;i<Dyx.size();++i){
      this->updateTail(Dyx[i],dyxNonlinearPart[i],DyxEnc[i],h);
      this->updateTail(DyxId[i],jacNonlinearPart[i],Enc[i],h);
    }
  }

  /// This function should compute a matrix M such that
  /// M_ii is logarithmic norm of the diagonal block
  /// M_ij is a norm of ij block
  /// The infinite dimensional space is split onto m+1 blocks
  MatrixType blockNorms(const VectorType& a, size_type m) const{
    MatrixType result(m+1,m+1);
    dxxNorm(a,result);
    dxyNorm(a,result);
    dyxNorm(a,result);
    result(m+1,m+1) = dyyLogarithmicNorm(a,m);
    return result;
  }

  void setParameter(ScalarType nu) {
    this->nu = nu;
    for(size_type i=0;i<lambda.size();++i){
      size_type k2 = i*i;
      lambda[i] = (k2*(1.-nu*k2));
    }
  }

  std::tuple<ScalarType,ScalarType,ScalarType> computeD1D2DI(const VectorType& a){
    VectorArray tmp;
    tmp.push_back(a);
    return {computeD1(tmp,0),computeD2(tmp,0),computeConstantForInfinitePart(tmp,0)};
  }

  ScalarType getLambda(size_type k){
    for(size_type i=lambda.size();i<=k;++i){
      size_type k2 = i*i;
      lambda.push_back(k2*(1.-nu*k2));
    }
    return lambda[k];
  }

private:
  /// This function computes logarithmic norm on the diagonal and norm out of diagonal of Dxx(m \times m) block of the derivative of the vector field at a=(x,y).
  void dxxNorm(const VectorType& a, MatrixType& result) const;

  /// This function computes operator norm of each row in the block Dxy(m\times \infinity) of the derivative of the vector field at a=(x,y).
  void dxyNorm(const VectorType& a, MatrixType& result) const;

  /// This function computes operator norm of each column in the block Dyx(\infinity \times m) of the derivative of the vector field at a=(x,y).
  void dyxNorm(const VectorType& a, MatrixType& result) const;

  /// This function computes logarithmic norm of Dyy(\infinity \times\infinity) block of the derivative of the vector field at a=(x,y).
  ScalarType dyyLogarithmicNorm(const VectorType& a, size_type m) const;

  void blockDerivative(const VectorType& v, VectorArray& a, MatrixArray& J, MatrixType& A) {
      a[0] = v;
      a[1] = v;
      for(size_type j=0;j<A.numberOfColumns();++j){
        for(size_type c=1;c<=this->dimension();++c)
          J[j][0].setCoefficient(c, j+1==c ? 1. : 0.);
        J[j][0].setGeometricDecay(v.getGeometricDecay());
        J[j][0].setConstant(0.);
      }
      DissipativeVectorField<VectorType>::computeODECoefficients(a,J,1,A.numberOfColumns());
      for(size_type j=0;j<A.numberOfColumns();++j){
        for(size_type c=1;c<=this->dimension();++c)
          A(j+1,c) = J[j][1].getCoefficient(c);
      }
    }

  // auxiliary, type dependent functions to avoid redundant code
  // used in template function makeSelfConsistentBound
  static ScalarType getCoeff(const VectorArray& a, int k) { return a[0].getCoefficient(k); }
  static ScalarType getCoeff(const C1Data& c, int k) { return (*c.second)[0].getCoefficient(k); }
  static void setCoeff(VectorArray& a, int k, ScalarType s) { a[0].setCoefficient(k,s); }
  static void setCoeff(C1Data& c, int k, ScalarType s) { return (*c.second)[0].setCoefficient(k,s); }
  static ScalarType getConstant(const VectorArray& a) { return a[0].getConstant(); }
  static ScalarType getConstant(const C1Data& c) { return (*c.second)[0].getConstant(); }
  static void setConstant(VectorArray& a, ScalarType s) { a[0].setConstant(s); }
  static void setConstant(C1Data& c, ScalarType s) { (*c.second)[0].setConstant(s); }

  template<class V>
  void makeSelfConsistentBound(V& v, ScalarArray& nonlinearPart, ScalarType q);

  ScalarType computeExplicitCoefficient(const VectorArray& a, size_type i, size_type k);
  ScalarType computeExplicitCoefficient(const VectorArray& a, const VectorArray& J, size_type i, size_type k);
  ScalarType computeExplicitCoefficient(const C1Data c, size_type i, size_type k) { return computeExplicitCoefficient(*c.first,*c.second,i,k); }

  ScalarType computeConstantForInfinitePart(const VectorArray& a, size_type i);
  ScalarType computeConstantForInfinitePart(const VectorArray& a, const VectorArray& J, size_type i);
  ScalarType computeConstantForInfinitePart(const C1Data c, size_type i){ return computeConstantForInfinitePart(*c.first,*c.second,i); }

  ScalarType computeD1(const VectorArray& a, size_type i);
  ScalarType computeD1(const VectorArray& a, const VectorArray& J, size_type i);
  ScalarType computeD1(const C1Data c, size_type i) { return computeD1(*c.first,*c.second,i); }

  ScalarType computeD2(const VectorArray& a, size_type i);
  ScalarType computeD2(const VectorArray& a, const VectorArray& J, size_type i);
  ScalarType computeD2(const C1Data c, size_type i) { return computeD2(*c.first,*c.second,i); }

  void computeNextTail(VectorArray& a, size_type i, ScalarType DI);
  void updateTail(VectorType& x, const ScalarArray& nPart, const VectorArray& enc, ScalarType h) const;

  // members
  ScalarType nu;
  size_type m_dimension, m_firstDissipativeVariable;
  ScalarArray lambda;
  ScalarArray nonlinearPart;
  std::vector<ScalarArray> jacNonlinearPart;
  std::vector<ScalarArray> dyxNonlinearPart;
}; // end of class OneDimKSSineVectorField

// ***************************************************************************

OneDimKSSineVectorField::OneDimKSSineVectorField(ScalarType nu, size_type dim, size_type firstDissipativeVariable)
  : nu(nu), m_dimension(dim), m_firstDissipativeVariable(firstDissipativeVariable),
    nonlinearPart(dim+2), jacNonlinearPart(dim+2, ScalarArray(dim+2)), dyxNonlinearPart(dim+2, ScalarArray(dim+2))
{
  for(size_type k=1;k<=m_dimension+1;++k){
    this->getLambda(k);
  }
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeExplicitCoefficient(const VectorArray& a, size_type i, size_type k) {
  if(k>2*m_dimension)
    throw std::runtime_error("OneDimKSSineVectorField::computeExplicitCoefficient(a,i,k) incorrect index");

  ScalarType s = ScalarType(0.);
  for(size_type j=0;j<=i;++j){
    for(size_type n=1;n<=(k-1)/2;++n)
      s += a[j].getCoefficient(n)*( a[i-j].getCoefficient(n+k) - a[i-j].getCoefficient(k-n) );
    for(size_type n=(k+1)/2;n<=this->m_dimension;++n)
      s += a[j].getCoefficient(n)*a[i-j].getCoefficient(n+k);
  }
  if(k%2==0){
    k /=2;
    size_type j=0;
    for(;j<i/2;++j)
      s -= a[j].getCoefficient(k)*a[i-j].getCoefficient(k);
    if(i%2)
      s -= a[j].getCoefficient(k)*a[j+1].getCoefficient(k);
    else
      s -= 0.5*sqr(a[j].getCoefficient(k));
  }
  return 2.*s;
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeExplicitCoefficient(const VectorArray& a, const VectorArray& J, size_type i, size_type k) {
  ScalarType s = ScalarType(0.);
  if(k<=this->m_dimension)
  {
    for(size_type j=0;j<=i;++j){
      for(size_type n=1;n<k;++n)
        s += a[j].getCoefficient(n)*( J[i-j].getCoefficient(n+k) - J[i-j].getCoefficient(k-n) )
          + J[j].getCoefficient(n)*a[i-j].getCoefficient(n+k);
      for(size_type n=k;n<=this->m_dimension;++n){
        s += a[j].getCoefficient(n)*J[i-j].getCoefficient(n+k);
        s += J[j].getCoefficient(n)*a[i-j].getCoefficient(n+k);
      }
    }
  } else {
      for(size_type j=0;j<=i;++j){
        for(size_type n=1;n<=this->m_dimension;++n)
          s += a[j].getCoefficient(n)*( J[i-j].getCoefficient(n+k) - J[i-j].getCoefficient(k-n) )
            + J[j].getCoefficient(n)*( a[i-j].getCoefficient(n+k) - a[i-j].getCoefficient(k-n) );
        for(size_type n=this->m_dimension+1;n<k-this->m_dimension;++n)
          s -= a[j].getCoefficient(n)*J[i-j].getCoefficient(k-n);
      }
  }
  return 2.*s;
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeD1(const VectorArray& a, size_type i) {
  ScalarType s = 0.;
  for(size_type j=0;j<=i;++j){
    ScalarType q = 1;
    for(size_type n=1;n<=this->m_dimension;++n){
      q *= a[i-j].getGeometricDecay();
      s += (q+1./q)*a[i-j].getConstant() * abs(a[j].getCoefficient(n));
    }
  }
  s *= 2;
  s /= (2*this->m_dimension+1);

  for(size_type j=0;j<=i;++j)
    s += a[j].getConstant() * a[i-j].getConstant();

  return s.rightBound();
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeD1(const VectorArray& a, const VectorArray& J, size_type i) {
  ScalarType s = 0.;
  for(size_type j=0;j<=i;++j){
    ScalarType q = 1;
    for(size_type n=1;n<=this->m_dimension;++n){
      q *= J[i-j].getGeometricDecay();
      s += (q+1./q)*(a[i-j].getConstant() * abs(J[j].getCoefficient(n))+J[i-j].getConstant() * abs(a[j].getCoefficient(n)));
    }
  }
  s /= (2*this->m_dimension+1);

  for(size_type j=0;j<=i;++j)
    s += a[j].getConstant() * J[i-j].getConstant();

  s *= 2.;
  return s.rightBound();
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeD2(const VectorArray& a, size_type i) {
  const ScalarType q = a[i].getGeometricDecay();
  ScalarType q1 = power(q,this->m_dimension);
  double d = 0;
  for(size_type k=this->m_dimension+1;k<=2*this->m_dimension;++k){
    q1 *= q;
    double d1 = rightBound( abs(q1*computeExplicitCoefficient(a,i,k)/k) );
    if(d1>d) d = d1;
  }
  return d;
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeD2(const VectorArray& a, const VectorArray& J, size_type i) {
  const ScalarType q = a[i].getGeometricDecay();
  ScalarType q1 = power(q,this->m_dimension);
  double d = 0;
  for(size_type k=this->m_dimension+1;k<=2*this->m_dimension;++k){
    q1 *= q;
    double d1 = rightBound( abs(q1*computeExplicitCoefficient(a,J,i,k)/k) );
    if(d1>d) d = d1;
  }
  return d;
}


// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeConstantForInfinitePart(const VectorArray& a, size_type i) {
  ScalarType s = 0.;
  for(size_type j=0;j<=i;++j)
    s += a[j].getConstant() * a[i-j].getConstant();

  const ScalarType q = a[i].getGeometricDecay();
  return (s/(power(q,2*this->m_dimension)*(q*q-1.))).rightBound();
}

// ***************************************************************************

OneDimKSSineVectorField::ScalarType
OneDimKSSineVectorField::computeConstantForInfinitePart(const VectorArray& a, const VectorArray& J, size_type i) {
  ScalarType s = 0.;
  for(size_type j=0;j<=i;++j)
    s += a[j].getConstant() * J[i-j].getConstant();

  const ScalarType q = a[i].getGeometricDecay();
  return (2.*s/(power(q,2*this->m_dimension)*(q*q-1.))).rightBound();
}

// ***************************************************************************

void OneDimKSSineVectorField::makeSelfConsistentBound(VectorArray& a, MatrixArray& J1, MatrixArray& J2, size_type numberOfColumns) {
  makeSelfConsistentBound(a);
  for(int i=0;i<numberOfColumns;++i){
    C1Data c1 = make_pair(&a,&(J1[i]));
    makeSelfConsistentBound(c1,this->jacNonlinearPart[i],J1[i][0].getGeometricDecay());
    C1Data c2 = make_pair(&a,&(J2[i]));
    makeSelfConsistentBound(c2,this->dyxNonlinearPart[i],J2[i][0].getGeometricDecay());
  }
}

// ***************************************************************************

template<class V>
void OneDimKSSineVectorField::makeSelfConsistentBound(V& a, ScalarArray& nPart, ScalarType q) {
  size_type M1 = this->m_dimension+1;
  size_type M2 = M1*M1;
  size_type M3 = M1*M2;

  bool found;
  for(int maxIt=10000;maxIt>0;maxIt--){
    ScalarType DI = ScalarType(-2,2)*computeConstantForInfinitePart(a,0);
    found = true;
    for(size_type k=this->m_firstDissipativeVariable; k<=this->m_dimension; ++k){
      ScalarType coeff = getCoeff(a,k);
      double left= coeff.leftBound();
      double right=coeff.rightBound();
      nPart[k] = k*(computeExplicitCoefficient(a,0,k) + DI/power(q,k));
      if(! (right*lambda[k] + nPart[k] < 0.) ){
        if(right>0)
          right = (-1.01*nPart[k].rightBound()/lambda[k]).rightBound();
        else
          right = (-0.99*nPart[k].rightBound()/lambda[k]).rightBound();
        found = false;
      }
      if(! (left*lambda[k] + nPart[k] > 0.) ){
        if(left<0)
          left =(-1.01*nPart[k].leftBound()/lambda[k]).leftBound();
        else
          left =(-0.99*nPart[k].leftBound()/lambda[k]).leftBound();
        found = false;
      }
      if(left>right){
        throw std::runtime_error("OneDimKSSineVectorField::makeSelfConsistentBound - cannot make self-consistent bound. Inequality cannot be solved.");
      }
      setCoeff(a,k,ScalarType(left,right));
    }
    DI = ScalarType(-2,2)*computeConstantForInfinitePart(a,0);
    ScalarType D = ScalarType(-1,1)*max(computeD1(a,0),computeD2(a,0));
    nPart[M1] = DI/M1 + D;
    ScalarType C = abs(nonlinearPart[M1]/(M2*(nu-ScalarType(1.)/M2))).rightBound();
    if( getConstant(a) < C ){
      // otherwise enlarge and check again.
      setConstant(a,1.01*C.rightBound());
      found = false;
    }
    if(found) {
      for(size_type i=1;i<this->m_firstDissipativeVariable;++i)
        nPart[i] = i*(computeExplicitCoefficient(a,0,i) + DI/power(q,i));
      return;
    }
  }

  throw std::runtime_error("OneDimKSSineVectorField::makeSelfConsistentBound - loop limit exceeded.");
}

// ***************************************************************************

void OneDimKSSineVectorField::computeNextTail(VectorArray& a, size_type i, ScalarType DI) {
  // some heuristic to choose delta
  int M1 = this->m_dimension+1;
  ScalarType q = a[i].getGeometricDecay();
  ScalarType D = max(computeD1(a,i),computeD2(a,i));
  int M2 = M1*M1;
  int M3 = M2*M1;
  ScalarType C = (nu-ScalarType(1.)/M2)*a[i].getConstant() + DI/M3 + D/M2;

  double t = 3;
  ScalarType delta = (t-1+q)/t;
  ScalarType l = log(delta);
  ScalarType L = (this->m_dimension>4/l) ? power(M1,4)/power(delta,M1)
			     : 256./power(l*ScalarType::euler(),4);
  a[i+1].setGeometricDecay((q/delta).leftBound());

  // Fix constant for next level of derivative
  a[i+1].setConstant(L*C/(i+1));
}

// ***************************************************************************

void OneDimKSSineVectorField::computeODECoefficients(VectorArray& a, size_type p) {
  for(size_type i=0;i<p;++i){
    ScalarType DI = computeConstantForInfinitePart(a,i);
    ScalarType q = a[i].getGeometricDecay();
    ScalarType D = ScalarType(-2.,2.)*DI;

    for(size_type k=1; k<=this->m_dimension; ++k){
      ScalarType c = lambda[k]*a[i].getCoefficient(k) + k*(computeExplicitCoefficient(a,i,k) + D/q);
      a[i+1].setCoefficient(k,c/(i+1));
      q *= a[i].getGeometricDecay();
    }
    this->computeNextTail(a,i,DI);
  }
}

// ***************************************************************************

void OneDimKSSineVectorField::computeODECoefficients(const VectorArray& a, VectorArray& J, size_type p) {
  for(size_type i=0;i<p;++i){
    ScalarType DI = computeConstantForInfinitePart(a,J,i);
    ScalarType q = J[i].getGeometricDecay();
    ScalarType D = ScalarType(-2.,2.)*DI;

	  for(size_type k=1; k<=this->m_dimension; ++k){
      ScalarType c = lambda[k]*J[i].getCoefficient(k) + k*(computeExplicitCoefficient(a,J,i,k) + D/q);
      J[i+1].setCoefficient(k,c/(i+1));
      q *= J[i].getGeometricDecay();
    }
    this->computeNextTail(J,i,DI);
  }
}

// ***************************************************************************

void OneDimKSSineVectorField::updateTail(VectorType& x, const ScalarArray& nPart, const VectorArray& enc, ScalarType h) const {
  ScalarType E = enc[0].getConstant();

  for(size_type k=1;k<this->m_firstDissipativeVariable;++k){
    bool uppBoundOK = lambda[k]*getCoeff(enc,k).rightBound()+nPart[k]<0.;
    bool lowBoundOK = lambda[k]*getCoeff(enc,k).leftBound()+nPart[k]>0.;
    
    if(uppBoundOK and lowBoundOK) {
      ScalarType e = exp(lambda[k]*h);
      ScalarType N = -nPart[k]/lambda[k];
      ScalarType u = (enc[0].getCoefficient(k)-N)*e + N;
      if(! intersection(x.getCoefficient(k),u,u) )
        throw std::runtime_error("OneDimKSSineVectorField::updateTail Both - Intersection error\n");
      x.setCoefficient(k,u);
      continue;
    }
    if(uppBoundOK) {
      ScalarType e = exp(lambda[k]*h);
      ScalarType N = -nPart[k]/lambda[k];
      ScalarType u = (enc[0].getCoefficient(k).rightBound()-N)*e + N;
      if( ! (u.rightBound()>=x.getCoefficient(k).leftBound()) )
        throw std::runtime_error("OneDimKSSineVectorField::updateTail UP - Intersection error\n");
      u.setRightBound(capd::min(u.rightBound(),x.getCoefficient(k).rightBound()));
      x.setCoefficient(k,u);
      continue;
    }
    if(lowBoundOK) {
      ScalarType e = exp(lambda[k]*h);
      ScalarType N = -nPart[k]/lambda[k];
      ScalarType u = (enc[0].getCoefficient(k).leftBound()-N)*e + N;
      if( ! (u.leftBound()<=x.getCoefficient(k).rightBound()) )
        throw std::runtime_error("OneDimKSSineVectorField::updateTail UP - Intersection error\n");
      u.setLeftBound(capd::max(u.leftBound(),x.getCoefficient(k).leftBound()));
      x.setCoefficient(k,u);
      continue;
    }
  }

  for(size_type k=this->m_firstDissipativeVariable;k<=this->m_dimension;++k){
    ScalarType e = exp(lambda[k]*h);
    ScalarType N = -nPart[k]/lambda[k];
    ScalarType u = (enc[0].getCoefficient(k)-N)*e + N;
    if(! intersection(x.getCoefficient(k),u,u) )
      throw std::runtime_error("OneDimKSSineVectorField::updateTail - Intersection error\n");
    x.setCoefficient(k,u);
  }

  // now we update constant for far tail from linear approximation
  size_type M1 = this->m_dimension+1;
  size_type M2 = M1*M1;
  size_type M3 = M2*M1;
  ScalarType u = abs(nPart[M1]/(M2*(nu-ScalarType(1.)/M2))).rightBound();
  ScalarType K = (E-u)*exp(h*lambda[M1]) + u;
  x.setConstant(K.rightBound());
}

/// This function computes operator norm of Dxy(m \times\infinity) block of the derivative of the vector field at a=(x,y).
void OneDimKSSineVectorField::dxyNorm(const VectorType& a, MatrixType& result) const{
  // Assume l_\infty norm in the domain, so that we have to compute max of l_1 norms of corresponding rows.
  // The coefficients in the matrix in this block are given by
  // A_{k,c} = 2*k*(a_{k+c}+a_{c-k}) and here we have c>k
  int M = a.dimension();
  ScalarType q = a.getGeometricDecay();
  ScalarType C = 2.*a.getConstant()*power(q,-M)/(q-1);

  // numberOfColumns = m+1, where mxm is finite explicit block
  for(size_type k=1;k<result.numberOfColumns();++k){
    // split the sum into explicit coefficients and geometric series.
    // Geometric representation holds when c-k>M = dimension() => c> M+k
    // Start from a bound on geoemtric series
    ScalarType s = C;
    // case 1: both coefficients are explicit
    size_type c=result.numberOfColumns();
    for(;c<=M-k;++c)
      s += abs(a.getCoefficient(k+c)+a.getCoefficient(c-k));
    // case 2: a_{k+c} is represented uniformly and a_{c-k} is explicit
    for(;c<=M+k;++c)
      s += abs(a.getCoefficient(c-k));
    result(k,result.numberOfColumns()) = ((2*k)*s).rightBound();
  }
}

/// This function computes operator norm of Dyx(\infinity \times m) block of the derivative of the vector field at a=(x,y).
void OneDimKSSineVectorField::dyxNorm(const VectorType& a, MatrixType& result) const{
  // Assume l_\infty norm in the domain, so that we have to compute max of l_1 norms of corresponding rows.
  // The coefficients in the matrix in this block are given by
  // A_{k,c} = 2*k*(a_{k+c}-a_{k-c}) and here we have 1<=c<=m<k
  // If k>m+M both a_{k+c} and a_{k-c} are represented uniformly by Cq^{-(k+c)} and Cq^{-(k-c)}.
  // the function k-> k*q^{-k-c} takes its maximum at k = 1/Log[q].
  // Therefore the norm of all rows k>=k0:=max(1/Log(q),m+M) are less than the norm of k0 row.
  int M = a.dimension();
  int m = result.numberOfRows();
  ScalarType q = a.getGeometricDecay();
  size_type k0 = ceil(capd::max((1./log(q)).rightBound(),(double)(m-1+M)));

  for(size_type c=1;c<m;++c){
    // proceed c-column
    // Start from a bound in k0-th row
    result(m,c) = ((2*k0)*abs(a.getCoefficient(k0+c)-a.getCoefficient(k0-c))).rightBound();
    for(size_type k=m;k<k0;++k){
      result(m,c) = capd::max(result(m,c).rightBound(),((2*k)*abs(a.getCoefficient(k+c)-a.getCoefficient(k-c))).rightBound());
    }
  }
}

/// This function computes logarithmic norm of Dyy(\infinity \times\infinity) block of the derivative of the vector field at a=(x,y).
OneDimKSSineVectorField::ScalarType OneDimKSSineVectorField::dyyLogarithmicNorm(const VectorType& a, const size_type m) const{
  // Assume l_\infty norm in the domain, so that we have to compute max of l_1 norms of infinite number of rows, for k>m and c>m.
  // The coefficients in the matrix in this block are given by
  // A_{k,c} = 2*k*(a_{k+c}+a_{c-k})  if c>k
  // A_{k,c} = 2*k*(a_{k+c}-a_{k-c})  if c<k
  // A_{k,k} = L_k + 2*k*a_{2k}

  // MAX VALUE ON THE DIAGONAL:
  // Set N = L1Norm(a).
  // It turns out (some short computation) that sum of out-diagonal elements on the k-th row
  // is bounded by 4*k*L1Norm(a). Thus we have to find an upper bound of the quantity
  // h(k) = k^2(1-nu*k^2) + 4k*N
  // This function takes extremum if h'(k) = -4nu*k^3 + 2k + 4N =0
  // From the Cauchy estimate we have that all roots of this derivative are bounded by
  // k0 = 1 + max{0.5,N}/nu
  // Thus, we have to check max over rows m+1,...,\max{k0,m+1}, only.

  // First compute norm L1Norm(a). For k>M+1 we have a uniform geometric bound
  int M = a.dimension();
  ScalarType q = a.getGeometricDecay();
  /// bound for geometric series
  ScalarType N = a.getConstant()*power(q,-M)/(q-1);
  // add explicit coefficients
  for(int i=0;i<M;++i)
    N += abs(a[i]);
  int k0 = 1 + ceil((capd::max(N.rightBound(),0.5)/nu).rightBound());
  k0 = capd::max(k0,(int)m+1);

  // The lower bound for the result is the diagonal part at k0-row
  int k2 = k0*k0;
  ScalarType D = 4.*N;
  ScalarType result = (k2*(1.-nu*k2)+k0*D).rightBound();
  for(int k=m+1;k<k0;k++)
  {
    int k2 = k*k;
    result = capd::max((k2*(1.-nu*k2)+k*D).rightBound(),result.rightBound());
  }

  return result;
}

/// This function computes logarithmic norm of Dxx(m \times m) block of the derivative of the vector field at a=(x,y).
void OneDimKSSineVectorField::dxxNorm(const VectorType& a, MatrixType& result) const{
  int m = result.numberOfColumns()-1;
  for(int k=1;k<=m;k++)
  {
    int k2 = 2*k;
    result(k,k) = lambda[k] + k2*a.getCoefficient(k2); // diagonal part without abs
    for(int c=1;c<k;++c)
        result(k,c) = k2*abs(a.getCoefficient(k+c) - a.getCoefficient(k-c));
    for(int c=k+1;c<=m;++c)
        result(k,c) = k2*abs(a.getCoefficient(k+c) + a.getCoefficient(c-k));
  }
}

}} // namespace capd::pdes

#endif // _CAPD_PDES_OneDimKSSineVectorField_H_


/// @}
