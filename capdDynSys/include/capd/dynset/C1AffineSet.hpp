/////////////////////////////////////////////////////////////////////////////
/// @file C1AffineSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C1AFFINESET_HPP_
#define _CAPD_DYNSET_C1AFFINESET_HPP_

#include <stdexcept>
#include "capd/vectalg/iobject.hpp"
#include "capd/dynset/C1AffineSet.h"
#include "capd/geomset/CenteredAffineSet.hpp"
#include "capd/geomset/MatrixAffineSet.hpp"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"

namespace capd {
namespace dynset {

template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies>::C1AffineSet(const VectorType& x, ScalarType t)
  : SetType(
      x,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x),
    C1BaseSet(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies>::C1AffineSet(const VectorType& x, const VectorType& r, ScalarType t)
  : SetType(
      x+r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x, r),
    C1BaseSet(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies>::C1AffineSet(const VectorType& x, const MatrixType& B, const VectorType& r, ScalarType t)
  : SetType(
      x+B*r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x, B, r),
    C1BaseSet(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies>::C1AffineSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t)
    : SetType(
       VectorType(c0part),
       VectorType(c0part.dimension()),
       MatrixType(c1part),
       MatrixType(c0part.dimension(),c0part.dimension()),
       t),
   C0BaseSet(c0part),
   C1BaseSet(c1part)
{}


template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies> C1AffineSet<MatrixType,Policies>::create(const VectorType& x, ScalarType t)
{
  return C1AffineSet(x,t);
}

template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies> C1AffineSet<MatrixType,Policies>::create(const VectorType& x, const VectorType& r, ScalarType t)
{
  return C1AffineSet(x,r,t);
}

template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies> C1AffineSet<MatrixType,Policies>::create(const VectorType& x, const MatrixType& B,
    const VectorType& r,ScalarType t
)
{
  return C1AffineSet(x,B,r,t);
}


template<typename MatrixType, typename Policies>
C1AffineSet<MatrixType,Policies> C1AffineSet<MatrixType,Policies>::create(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t)
{
  return C1AffineSet(c0part,c1part,t);
}

template <typename MatrixType, typename Policies>
void C1AffineSet<MatrixType, Policies>::move(DynSysType & dynsys) {
  move(dynsys, *this);
}

template <typename MatrixType, typename Policies>
void C1AffineSet<MatrixType, Policies>::move(DynSysType & dynsys, C1AffineSet& result) const
{
  // important: here we assume that m_r contains zero
  // this is assured by each constructor and each step of this algorithm

  const size_type dim = m_x.dimension();
  VectorType y(dim), rem(dim), enc(dim);
  MatrixType jacPhi(dim,dim), jacEnc(dim,dim), jacRem(dim,dim);
  MatrixType B(dim,dim), Q(dim,dim);

  VectorType xx = VectorType(*this);

  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  dynsys.encloseC1Map(this->getCurrentTime(),
                    this->m_x, xx,          // input parameters
                    y, rem, enc,            // C^0 output
                    jacPhi, jacRem, jacEnc  // C^1 output
                    );


  result.m_x = y + rem;
  B = result.m_B = jacPhi * m_B;

  // ---------- C^0 - part -----------------

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentSet = result.m_x + result.m_B*this->m_r;

  // here we compute representation for the new set
  // xx is unnecessary now
  split(result.m_x, xx);

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB and updates
  this->Policies::computeBinvB(result.m_B,result.m_invB,this->m_r);

  // eventually we compute new representation of r
  result.m_r = (result.m_invB * B) * m_r + result.m_invB * xx;

  // ---------- C^1 - part -----------------
  MatrixType J = jacPhi + jacRem;
  result.m_D = J*this->m_D;
  B = result.m_Bjac = J*this->m_Bjac;

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentMatrix = J*this->m_currentMatrix;
  intersection(result.m_currentMatrix,result.m_D + result.m_Bjac*m_R,result.m_currentMatrix);

  // here we compute representation for the new set
  // jacRem is unnecessary now
  split(result.m_D, jacRem);

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB and updates
  this->Policies::computeBinvB(result.m_Bjac,result.m_invBjac,this->m_R);

  // eventually we compute new representation of r
  result.m_R = (result.m_invBjac * B) * m_R + result.m_invBjac * jacRem;

  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
  result.setLastEnclosure(enc);
  result.setLastMatrixEnclosure(jacEnc);
}

template <typename MatrixType, typename Policies>
std::string C1AffineSet<MatrixType, Policies>::show(void) const {
  std::ostringstream descriptor;
  descriptor << name()
             << C0BaseSet::toString()
             << C1BaseSet::toString();
  return descriptor.str();
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C1AFFINESET_HPP_

