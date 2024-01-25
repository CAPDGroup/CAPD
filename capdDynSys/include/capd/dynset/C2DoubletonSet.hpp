/////////////////////////////////////////////////////////////////////////////
/// @file C2DoubletonSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C2DOUBLETONSET_HPP_
#define _CAPD_DYNSET_C2DOUBLETONSET_HPP_

#include <sstream>
#include <stdexcept>
#include "capd/vectalg/iobject.hpp"
#include "capd/geomset/CenteredDoubletonSet.hpp"
#include "capd/geomset/MatrixDoubletonSet.hpp"
#include "capd/dynset/C2DoubletonSet.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/diffAlgebra/Hessian.hpp"
#include "capd/vectalg/algebraicOperations.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"
#include "capd/dynset/C1DoubletonSet.hpp"

namespace capd{
namespace dynset{

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, ScalarType t)
  : SetType(
      x,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t)
  : SetType(
      x+r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, r0),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t)
  : SetType(
      x+C*r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, C, r0),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t)
  : SetType(
      x+C*r0+r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, C, r0, r),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(
      const VectorType& x,
      const MatrixType& C, const VectorType& r0,
      const MatrixType& B, const VectorType& r,
      ScalarType t
   ) : SetType(
      x+C*r0+B*r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      HessianType(x.dimension()),
      HessianType(x.dimension()),
      t),
    C0BaseSet(x, C, r0, B, r),
    C1BaseSet(x.dimension()),
    Data(x.dimension()),
    m_Chess(x.dimension(),x.dimension()),
    m_Bhess(x.dimension(),x.dimension()),
    m_invBhess(x.dimension(),x.dimension()),
    m_alpha(x.dimension()),
    m_HR(x.dimension()),
    m_HR0(x.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t)
    : SetType(
       VectorType(c0part),
       VectorType(c0part.dimension()),
       MatrixType(c1part),
       MatrixType(c0part.dimension(),c0part.dimension()),
       HessianType(c0part.dimension()),
       HessianType(c0part.dimension()),
       t),
   C0BaseSet(c0part),
   C1BaseSet(c1part),
   Data(c0part.dimension()),
   m_Chess(c0part.dimension(),c0part.dimension()),
   m_Bhess(c0part.dimension(),c0part.dimension()),
   m_invBhess(c0part.dimension(),c0part.dimension()),
   m_alpha(c0part.dimension()),
   m_HR(c0part.dimension()),
   m_HR0(c0part.dimension())
{
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

template<typename MatrixType, typename Policies>
C2DoubletonSet<MatrixType,Policies>::C2DoubletonSet(
      const C0BaseSet & c0part,
      const C1BaseSet& c1part,
      const HessianType& H,
      const  ScalarType t
    )  : SetType(
       VectorType(c0part),
       VectorType(c0part.dimension()),
       MatrixType(c1part),
       MatrixType(c0part.dimension(),c0part.dimension()),
       HessianType(H),
       HessianType(c0part.dimension()),
       t),
   C0BaseSet(c0part),
   C1BaseSet(c1part),
   Data(c0part.dimension()),
   m_Chess(c0part.dimension(),c0part.dimension()),
   m_Bhess(c0part.dimension(),c0part.dimension()),
   m_invBhess(c0part.dimension(),c0part.dimension()),
   m_alpha(H),
   m_HR(c0part.dimension()),
   m_HR0(c0part.dimension())
{
  split(m_alpha,m_HR0);
  m_Chess.setToIdentity();
  m_Bhess.setToIdentity();
  m_invBhess.setToIdentity();
}

// ------------------------------------------------------------------------

template<typename MatrixType, class Policies>
void C2DoubletonSet<MatrixType,Policies>::move(DynSysType& c2dynsys, C2DoubletonSet& result)
{
  const size_type dim = this->m_currentSet.dimension();
  MatrixType B(dim,dim), deltaC(dim,dim);

  HessianType EH(dim), LH(dim), RH(dim);

  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  if(subset(this->m_x,this->m_currentSet)){
    result.x = this->m_x;
    subtractObjects( this->m_currentSet, this->m_x, result.deltaX );
  } else
    split(this->m_currentSet,result.x,result.deltaX);

  MatrixType V = this->m_currentMatrix;
  c2dynsys.encloseC2Map(this->getCurrentTime(),
                    result.x, this->m_currentSet,                 // input parameters
                    result.y, result.rem, result.enc,             // C^0 output
                    result.jacPhi, result.jacRem, result.jacEnc,  // C^1 output
                    LH, RH, EH                                    // output for C^2 part
                    );

  C0DoubletonSet<MatrixType,Policies>::move(*this,result,result.m_currentSet,result);

  result.jacPhi += result.jacRem;
  MatrixType J = result.jacPhi;
  result.m_currentMatrix = result.jacPhi*this->m_currentMatrix;
  C1DoubletonSet<MatrixType,Policies>::move(*this,result,result.m_currentMatrix,result);

  // ---------- C^2 part: enclosure -----------------

  // compute LH += RH;
  capd::vectalg::addAssignObjectObject(LH,RH);

  // first enclosure - store it in deltaAplha
  // newHessian = JacPhi*oldHessian + alpha
  HessianType alpha = LH*V;
  HessianType deltaAlpha = J*this->m_currentHessian + alpha;

  // second possible enclosure
  // newHessian = JacPhi*m_alpha + (JacPhi*Bhess)*HR + (JacPhi*Chess)*HR0 + alpha

  result.m_Chess = J*this->m_Chess;
  B = result.m_Bhess = J*this->m_Bhess;

  result.m_currentHessian = alpha;
  result.m_currentHessian += J*this->m_alpha;
  result.m_currentHessian += result.m_Chess*this->m_HR0;
  result.m_currentHessian += B*this->m_HR;

  if(! capd::vectalg::intersection(deltaAlpha,result.m_currentHessian,result.m_currentHessian))
    throw std::logic_error("C2Doubleton::move error: empty intersection of two enclosures of Hessian.");

  // computation of new representation
  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB
  this->Policies::computeBinvB(result.m_Bhess,result.m_invBhess,this->m_HR);

  result.m_alpha = J*this->m_alpha;
  result.m_alpha += alpha;
  capd::vectalg::split(result.m_alpha,deltaAlpha);
  capd::vectalg::split(result.m_Chess,deltaC);

  deltaAlpha += deltaC*m_HR0;
  result.m_HR = (result.m_invBhess*B)*this->m_HR;
  result.m_HR += result.m_invBhess * deltaAlpha;

  if(&result != this){
    result.m_HR0 = m_HR0;
  }

  // save enclosures and new time
  result.setCurrentTime(this->getCurrentTime()+c2dynsys.getStep());
  result.setLastEnclosure(result.enc);
  result.setLastMatrixEnclosure(result.jacEnc);
  result.setLastHessianEnclosure(EH);

  this->Policies::reorganizeIfNeeded(result.m_B,result.m_invB,result.m_r,result.m_C,result.m_r0);
  this->Policies::reorganizeC1IfNeeded(result.m_Bjac,result.m_invBjac,result.m_R,result.m_Cjac,result.m_R0);
  this->Policies::reorganizeC2IfNeeded(result.m_Bhess,result.m_invBhess,result.m_HR,result.m_Chess,result.m_HR0);
}

// -------------------------------------------------------------

template<typename MatrixType, class Policies>
std::string C2DoubletonSet<MatrixType,Policies>::show(void) const
{
  std::ostringstream descriptor;
  descriptor << name()
             << C0BaseSet::toString()
             << C1BaseSet::toString();
  return descriptor.str();
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C2RECT2SET_HPP_
