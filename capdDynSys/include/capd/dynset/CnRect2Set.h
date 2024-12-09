/////////////////////////////////////////////////////////////////////////////
/// @file CnRect2Set.h
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_CNRECT2SET_H_
#define _CAPD_DYNSET_CNRECT2SET_H_

#include <stdexcept>
#include "capd/dynset/CnSet.h"
#include "capd/vectalg/algebraicOperations.hpp"

namespace capd{
namespace dynset{
/// @addtogroup dynset
/// @{

/**
 * Set that stores all derivatives to given order in doubleton form with reorganization moved by QR decomposition method.
 */
template<typename MatrixT, typename Policies , __size_type DEGREE=0>
class CnRect2Set : public Policies, public CnSet<MatrixT,DEGREE>{
public:
  typedef MatrixT MatrixType;
  typedef CnSet<MatrixT,DEGREE> SetType;
  typedef typename SetType::VectorType VectorType;
  typedef typename SetType::ScalarType ScalarType;
  typedef typename SetType::RefVectorType RefVectorType;
  typedef typename SetType::JetType JetType;
  typedef typename SetType::Multipointer Multipointer;
  typedef typename SetType::Multiindex Multiindex;
  typedef typename MatrixType::size_type size_type;

  // constructors
  CnRect2Set(const VectorType& x, size_type degree, ScalarType t = TypeTraits<ScalarType>::zero());
  CnRect2Set(const VectorType& x, const VectorType& r0, size_type degree, ScalarType t = TypeTraits<ScalarType>::zero());
  CnRect2Set(const VectorType& x, const MatrixType& C, const VectorType& r0, size_type degree, ScalarType t = TypeTraits<ScalarType>::zero());
  CnRect2Set(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, size_type degree, ScalarType t = TypeTraits<ScalarType>::zero());
  CnRect2Set(const VectorType& x, const MatrixType& C, const VectorType& r0, const MatrixType& B, const VectorType& r, size_type degree, ScalarType t = TypeTraits<ScalarType>::zero());

  CnRect2Set(const JetType& x, ScalarType t = TypeTraits<ScalarType>::zero());
  CnRect2Set(JetType x, const MatrixType& C, const VectorType& r0, ScalarType t = TypeTraits<ScalarType>::zero());

  CnRect2Set& operator=(const VectorType&);
  template<class DynSysT>
  void move(DynSysT& cndynsys) { this->move(cndynsys,*this); }

  template<class DynSysT>
  void move(DynSysT& cndynsys, CnRect2Set& result) const;

  using SetType::operator VectorType;
  using SetType::operator MatrixType;
  std::string show() const {
    throw std::runtime_error("CnRect2Set::show - not implemented");
  }

  ScalarType evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const {
    ScalarType r = TypeTraits<ScalarType>::zero();
    for(size_type i=0;i<x0.dimension();++i){
      r += gradient[i]*(this->m_x(i)-x0[i]);
      r += (gradient*this->m_B.column(i))*this->m_r(i);
      r += (gradient*this->m_C.column(i))*this->m_r0(i);
    }
    intersection(gradient*((VectorType)(*this)-x0),r,r);
    return r;
  }
private:
  JetType m_x,m_r,m_r0;
  JetType phi, rem, enc;
  MatrixType m_B, m_invB, m_C;
  MatrixType m_Bjac, m_invBjac, m_Cjac;
}; // end of class CnRect2Set


// ------------------------------ move ---------------------------------------

template<typename MatrixType, typename Policies, __size_type DEGREE>
template<class DynSysT>
void CnRect2Set<MatrixType,Policies,DEGREE>::move(DynSysT& cndynsys, CnRect2Set& result) const
{
  const size_type dimension = this->dimension();
  VectorType xx(*this);
  VectorType x(m_x), s(dimension);
  MatrixType S(dimension,dimension);

  VectorType y = cndynsys.encloseCnMap(this->getCurrentTime(),x,xx,result.phi,result.rem,result.enc);

  typename JetType::RefVectorType thisX = this->m_x();
  typename JetType::RefVectorType thisR = this->m_r();
  typename JetType::RefVectorType thisR0 = this->m_r0();

  typename JetType::RefVectorType resultX = result.m_x();
  typename JetType::RefVectorType resultR = result.m_r();
  typename JetType::RefVectorType resultR0 = result.m_r0();

  typename JetType::RefVectorType thisCurrentSet = this->m_currentSet();
  typename JetType::RefVectorType resultCurrentSet = result.m_currentSet();
  typename JetType::RefVectorType nonlinearPart = result.rem();

  if(this!=&result)
    result.m_r0 = this->m_r0;

  // C^0 part
  y += VectorType(result.rem());

  MatrixType A = MatrixType(result.phi);
  result.m_C = A*m_C;
  MatrixType B = A*m_B;
  result.m_B = B;
  resultCurrentSet = y + A*(xx-thisX);
  // here we compute enclosure of the image after one iteration of the map/flow
  if(!intersection(resultCurrentSet,y + (result.m_C)*thisR0 + B*thisR,resultCurrentSet))
    throw std::runtime_error("CnRect2Set - fatal error! Intersection of two different representations of a set is empty!");

  split(result.m_C,S);
  split(y,s);
  resultX = y;
  s += S * thisR0;
  this->Policies::computeBinvB(result.m_B,result.m_invB,thisR);
  resultR = (result.m_invB*B)*thisR + result.m_invB*s;
  if(!intersection(resultCurrentSet, resultX + (result.m_C)*resultR0 + (result.m_B)*resultR, resultCurrentSet))
    throw std::runtime_error("CnRect2Set - fatal error! Intersection of two different representations of a set is empty!");
  this->Policies::reorganizeIfNeeded(result.m_B,result.m_invB,resultR,result.m_C,resultR0);

//   Phi += Rem
  capd::vectalg::addAssignObjectObject(result.phi,result.rem);

  A = MatrixType(result.phi); // in fact A += Rem
  MatrixType C = result.m_Cjac = A * m_Cjac;
  split(result.m_Cjac,S);
  B = result.m_Bjac = A*m_Bjac;
  this->Policies::computeBinvB(result.m_Bjac,result.m_invBjac,(MatrixType)this->m_r);
  MatrixType J=result.m_invBjac*B;

// C^1 -- C^r part
  // Rem is now unnecessary. We will store in Rem a nonlinear part of composition
  if(this->degree()>0)
  {
    substitutionPowerSeries(result.phi,this->m_currentSet,result.rem,true);
  }

  size_type i=1;
  typename ScalarType::BoundType maxSizeR = 0., maxSizeR0 = 0.;
  for(
      thisX.next(), thisR.next(), thisR0.next(), resultX.next(), resultR.next(), resultR0.next(),
      nonlinearPart.next(), thisCurrentSet.next(), resultCurrentSet.next();

      i < binomial(dimension+this->degree(),this->degree());

      thisX.next(), thisR.next(), thisR0.next(), resultX.next(), resultR.next(), resultR0.next(),
      nonlinearPart.next(), thisCurrentSet.next(), resultCurrentSet.next(),++i
      )
  {
      VectorType temp = nonlinearPart + A*thisX;
      resultCurrentSet = temp + C*thisR0 + B*thisR;
      y = temp + S*thisR0;
      split(y,s);
      resultR = (result.m_invBjac*s)  + (J*thisR);
      resultX = y;
      maxSizeR = capd::max(maxSizeR,capd::vectalg::maxDiam(resultR).rightBound());
      maxSizeR0 = capd::max(maxSizeR0,capd::vectalg::maxDiam(resultR0).rightBound());
  }

  if(maxSizeR > 20.*maxSizeR0)
  {
    typename JetType::RefVectorType resultR = result.m_r();
    typename JetType::RefVectorType resultR0 = result.m_r0();
    B = (result.m_invBjac * result.m_Cjac);

    for(i=1, resultR.next(), resultR0.next();
        i < binomial(dimension+this->degree(),this->degree());
        resultR.next(), resultR0.next(),++i
        )
    {
      resultR0 = resultR + B * resultR0;
      resultR.clear();
    }

    result.m_Cjac = result.m_Bjac;
    result.m_Bjac.setToIdentity();
    result.m_invBjac.setToIdentity();
  }

  result.setCurrentTime(this->getCurrentTime() + cndynsys.getStep());

  // this might be seen as unnecessary step because we have
  // enc - in Cnrect2Set
  // m_lastEnclosure inherited from CnSet
  // enc, however, might be in an inconsistent state if enclosure could not be found.
  // m_lastEnclosure always contains valid enclosure copied here after successful integration step.
  result.setLastJetEnclosure(result.enc);
}

/// @}
}} // the end of the namespace capd::dynset

#endif // _CAPD_DYNSET_CNRECT2SET_H_
