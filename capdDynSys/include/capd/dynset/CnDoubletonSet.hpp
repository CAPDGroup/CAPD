/////////////////////////////////////////////////////////////////////////////
/// @file CnDoubletonSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library, 
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details. 

#ifndef _CAPD_DYNSET_CNDOUBLETON_HPP_
#define _CAPD_DYNSET_CNDOUBLETON_HPP_

#include "capd/dynset/CnDoubletonSet.h"
#include "capd/vectalg/iobject.hpp"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"

namespace capd{
namespace dynset{

// ----------------------------- CONSTRUCTORS --------------------------------

template<typename MatrixT, typename Policies, __size_type DEGREE>
CnDoubletonSet<MatrixT,Policies,DEGREE>::CnDoubletonSet(const VectorType& x, size_type degree, ScalarType t)
  : SetType(x,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_C(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_B(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_invB(1,x.dimension(),degree,MatrixT::Identity(x.dimension()))
{
  typename MatrixT::RefColumnVectorType refX = this->m_x();
  typename MatrixT::RefColumnVectorType refR0 = this->m_r0();
  refX = x;
  split(refX,refR0);

  // set Identity
  this->m_x.setMatrix(this->m_B(0));
  this->m_currentSet.setMatrix(this->m_B(0));
  this->m_lastEnclosure.setMatrix(this->m_B(0));
}

// -------------------------------------------------------------------------------------

template<typename MatrixT, typename Policies, __size_type DEGREE>
CnDoubletonSet<MatrixT,Policies,DEGREE>::CnDoubletonSet(const VectorType& x, const VectorType& r0, size_type degree, ScalarType t)
  : SetType(x+r0,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_C(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_B(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_invB(1,x.dimension(),degree,MatrixT::Identity(x.dimension()))
{
  typename MatrixT::RefColumnVectorType refX = this->m_x();
  typename MatrixT::RefColumnVectorType refR0 = this->m_r0();
  refX = x + r0;
  split(refX,refR0);

  // set Identity
  this->m_x.setMatrix(this->m_B(0));
  this->m_currentSet.setMatrix(this->m_B(0));
  this->m_lastEnclosure.setMatrix(this->m_B(0));
}

// -------------------------------------------------------------------------------------

template<typename MatrixT, typename Policies, __size_type DEGREE>
CnDoubletonSet<MatrixT,Policies,DEGREE>::CnDoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, size_type degree, ScalarType t)
  : SetType(x+C*r0,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_C(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_B(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_invB(1,x.dimension(),degree,MatrixT::Identity(x.dimension()))
{
  this->m_x() = x;
  this->m_r0() = r0;
  this->m_C(0) = C;

  typename MatrixT::RefColumnVectorType refX = this->m_x();
  typename MatrixT::RefColumnVectorType refR0 = this->m_r0();
  typename MatrixT::RefColumnVectorType refR = this->m_r();

  // center the set if necessary
  if(!subset(VectorType(x.dimension()),r0)){
    split(refR0,refR,refR0);
    this->m_x() = this->m_x() + C*refR;
  }
  split(refX,refR);
  // set Identity
  this->m_x.setMatrix(this->m_B(0));
  this->m_currentSet.setMatrix(this->m_B(0));
  this->m_lastEnclosure.setMatrix(this->m_B(0));
}

// -------------------------------------------------------------------------------------

template<typename MatrixT, typename Policies, __size_type DEGREE>
CnDoubletonSet<MatrixT,Policies,DEGREE>::CnDoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, size_type degree, ScalarType t)
  : SetType(x+C*r0+r,degree,t),
    m_x(x.dimension(),degree),
    m_r(x.dimension(),degree),
    m_r0(x.dimension(),degree),
    phi(x.dimension(),degree),
    rem(x.dimension(),degree),
    enc(x.dimension(),degree),
    m_C(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_B(1,x.dimension(),degree,MatrixT::Identity(x.dimension())),
    m_invB(1,x.dimension(),degree,MatrixT::Identity(x.dimension()))
{
  this->m_x() = x;
  this->m_r0() = r0;
  this->m_C(0) = C;

  typename MatrixT::RefColumnVectorType refX = this->m_x();
  typename MatrixT::RefColumnVectorType refR0 = this->m_r0();
  typename MatrixT::RefColumnVectorType refR = this->m_r();

  // center the set if necessary
  if(!subset(VectorType(x.dimension()),r0)){
    split(refR0,refR,refR0);
    this->m_x() = this->m_x() + C*refR;
  }
  this->m_x() = this->m_x() + r;
  split(refX,refR);
  // set Identity
  this->m_x.setMatrix(this->m_B(0));
  this->m_currentSet.setMatrix(this->m_B(0));
  this->m_lastEnclosure.setMatrix(this->m_B(0));
}


template<typename MatrixT, typename Policies, __size_type DEGREE>
CnDoubletonSet<MatrixT,Policies,DEGREE>::CnDoubletonSet(const JetType& x, ScalarType t)
  : SetType(x,t),
    m_x(x),
    m_r(x.dimension(),x.degree()),
    m_r0(x.dimension(),x.degree()),
    phi(x.dimension(),x.degree()),
    rem(x.dimension(),x.degree()),
    enc(x.dimension(),x.degree()),
    m_C(1,x.dimension(),x.degree(),MatrixT::Identity(x.dimension())),
    m_B(1,x.dimension(),x.degree(),MatrixT::Identity(x.dimension())),
    m_invB(1,x.dimension(),x.degree(),MatrixT::Identity(x.dimension()))
{
  split(m_x,m_r0);
  this->m_currentSet = x;
  this->m_lastEnclosure = x;
}


}} // namespace capd::dynset

#endif // _CAPD_DYNSET_CNDOUBLETON_HPP_
