/////////////////////////////////////////////////////////////////////////////
/// @file CnRect2Set.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_CNRECT2SET_HPP_
#define _CAPD_DYNSET_CNRECT2SET_HPP_

#include "capd/dynset/CnRect2Set.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/vectalg/iobject.hpp"
#include "capd/diffAlgebra/Jet.hpp"

namespace capd{
namespace dynset{

// ----------------------------- CONSTRUCTORS --------------------------------
template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>::CnRect2Set(const VectorType& x, size_type degree, ScalarType t)
  : SetType(x,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_B(x.dimension(),x.dimension()),
    m_invB(x.dimension(),x.dimension()),
    m_C(x.dimension(),x.dimension()),
    m_Bjac(x.dimension(),x.dimension()),
    m_invBjac(x.dimension(),x.dimension()),
    m_Cjac(x.dimension(),x.dimension())
{
  m_x() = x;
  split(m_x,m_r0);
  m_C.setToIdentity();
  m_B.setToIdentity();
  m_invB.setToIdentity();
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
  m_x.setMatrix(m_Cjac);
  this->m_currentSet.setMatrix(m_Cjac);
  this->m_lastEnclosure.setMatrix(m_Cjac);
}

template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>::CnRect2Set(const VectorType& x, const VectorType& r0, size_type degree, ScalarType t)
  : SetType(x+r0,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_B(x.dimension(),x.dimension()),
    m_invB(x.dimension(),x.dimension()),
    m_C(x.dimension(),x.dimension()),
    m_Bjac(x.dimension(),x.dimension()),
    m_invBjac(x.dimension(),x.dimension()),
    m_Cjac(x.dimension(),x.dimension())
{
  m_x() = x;
  m_r0() = r0;
  m_C.setToIdentity();
  m_B.setToIdentity();
  m_invB.setToIdentity();
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
  m_x.setMatrix(m_Cjac);
  this->m_currentSet.setMatrix(m_Cjac);
  this->m_lastEnclosure.setMatrix(m_Cjac);
}

template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>::CnRect2Set(const VectorType& x, const MatrixType& C, const VectorType& r0, size_type degree, ScalarType t)
  : SetType(x+C*r0,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_B(x.dimension(),x.dimension()),
    m_invB(x.dimension(),x.dimension()),
    m_C(x.dimension(),x.dimension()),
    m_Bjac(x.dimension(),x.dimension()),
    m_invBjac(x.dimension(),x.dimension()),
    m_Cjac(x.dimension(),x.dimension())
{
  m_x() = x;
  m_r0() = r0;
  m_C = C;
  m_B.setToIdentity();
  m_invB.setToIdentity();
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
  m_x.setMatrix(m_Cjac);
  this->m_currentSet.setMatrix(m_Cjac);
  this->m_lastEnclosure.setMatrix(m_Cjac);
}

template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>::CnRect2Set(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, size_type degree, ScalarType t)
  : SetType(x+C*r0+r,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_B(x.dimension(),x.dimension()),
    m_invB(x.dimension(),x.dimension()),
    m_C(x.dimension(),x.dimension()),
    m_Bjac(x.dimension(),x.dimension()),
    m_invBjac(x.dimension(),x.dimension()),
    m_Cjac(x.dimension(),x.dimension())
{
  m_x() = x;
  m_r0() = r0;
  m_r() = r;
  m_C = C;
  m_B.setToIdentity();
  m_invB.setToIdentity();
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
  m_x.setMatrix(m_Cjac);
  this->m_currentSet.setMatrix(m_Cjac);
  this->m_lastEnclosure.setMatrix(m_Cjac);
}

template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>::CnRect2Set(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r, size_type degree, ScalarType t)
  : SetType(x+C*r0+B*r,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_B(x.dimension(),x.dimension()),
    m_invB(x.dimension(),x.dimension()),
    m_C(x.dimension(),x.dimension()),
    m_Bjac(x.dimension(),x.dimension()),
    m_invBjac(x.dimension(),x.dimension()),
    m_Cjac(x.dimension(),x.dimension())
{
  m_x() = x;
  m_r0() = r0;
  m_r() = r;
  m_C = C;
  m_B = B;
  m_invB = capd::matrixAlgorithms::gaussInverseMatrix(B);
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
  m_x.setMatrix(m_Cjac);
  this->m_currentSet.setMatrix(m_Cjac);
  this->m_lastEnclosure.setMatrix(m_Cjac);
}

// ----------------------------- CONSTRUCTORS --------------------------------
template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>::CnRect2Set(const JetType& x, ScalarType t)
  : SetType(x,t),
    m_x(x),
    m_r(x.dimension(),x.degree()),
    m_r0(x.dimension(),x.degree()),
    phi(x.dimension(),x.degree()),
    rem(x.dimension(),x.degree()),
    enc(x.dimension(),x.degree()),
    m_B(x.dimension(),x.dimension()),
    m_invB(x.dimension(),x.dimension()),
    m_C(x.dimension(),x.dimension()),
    m_Bjac(x.dimension(),x.dimension()),
    m_invBjac(x.dimension(),x.dimension()),
    m_Cjac(x.dimension(),x.dimension())
{
  split(m_x,m_r0);
  m_C.setToIdentity();
  m_B.setToIdentity();
  m_invB.setToIdentity();
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
}

// ----------------------------------------------------------------------------

template<typename MatrixType, typename Policies, __size_type DEGREE>
CnRect2Set<MatrixType,Policies,DEGREE>& CnRect2Set<MatrixType,Policies,DEGREE>::operator=(const VectorType& v)
{
  if(v.dimension()!=this->dimension())
    throw std::runtime_error("CnRect2Set::operator=(Vector) - incorrect dimensions");
  VectorType temp(this->dimension());
  VectorType y = v;
  split(y,temp);
  m_x.clear();
  m_r.clear();
  m_r0.clear();
  this->m_currentSet.clear();

  m_C.setToIdentity();
  m_B.setToIdentity();
  m_invB.setToIdentity();
  m_Cjac.setToIdentity();
  m_Bjac.setToIdentity();
  m_invBjac.setToIdentity();
   
  m_x() = y;
  m_r0() = temp;
  this->m_currentSet() = v;
  this->m_lastEnclosure() = v;
  if(this->degree()>0)
  {
    for(size_type i=0;i<this->dimension();++i)
      this->m_lastEnclosure(i,i) = this->m_currentSet(i,i) = TypeTraits<ScalarType>::one();
  }
  return *this;
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_CNRECT2SET_HPP_
