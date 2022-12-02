/////////////////////////////////////////////////////////////////////////////
/// @file C0AffineSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0AFFINETSET_HPP_
#define _CAPD_DYNSET_C0AFFINETSET_HPP_

#include <stdexcept>
#include "capd/vectalg/iobject.hpp"
#include "capd/dynset/C0AffineSet.h"
#include "capd/geomset/CenteredAffineSet.hpp"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"

namespace capd {
namespace dynset {

template<typename MatrixType, typename Policies>
C0AffineSet<MatrixType,Policies>::C0AffineSet(const VectorType& x, ScalarType t)
  : SetType(x,VectorType(x.dimension()),t),
    BaseSet(x)
{}

template<typename MatrixType, typename Policies>
C0AffineSet<MatrixType,Policies>::C0AffineSet(const VectorType& x, const VectorType& r, ScalarType t)
  : SetType(x+r,VectorType(x.dimension()),t),
    BaseSet(x,r)
{}

template<typename MatrixType, typename Policies>
C0AffineSet<MatrixType,Policies>::C0AffineSet(const VectorType& x, const MatrixType& B, const VectorType& r, ScalarType t)
  : SetType(x+B*r,VectorType(x.dimension()),t),
    BaseSet(x,B,r)
{}

template <typename MatrixType, typename Policies>
void C0AffineSet<MatrixType, Policies>::move(DynSysType & dynsys) {
  move(dynsys, *this);
}

template <typename MatrixType, typename Policies>
void C0AffineSet<MatrixType, Policies>::move(DynSysType & dynsys, C0AffineSet& result) const
{
  // important: here we assume that both m_r and m_r0 contains zero
  // this is assured by each constructor and each step of this algorithm

  VectorType y(m_x.dimension()), rem(m_x.dimension()), enc(m_x.dimension());
  MatrixType A(m_x.dimension(), m_x.dimension());

  VectorType xx = VectorType(*this);
  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  dynsys.encloseC0Map(this->getCurrentTime(),m_x, xx, y, rem, enc, A);

  result.m_x = y + rem;
  A = result.m_B = A * m_B;

  // here we compute enclosure of the image after one iteration of the map/flow
  result.m_currentSet = result.m_x + result.m_B*this->m_r;

  // here we compute representation for the new set
  // xx is unnecessary now
  split(result.m_x, xx);

  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB and updates
  this->Policies::computeBinvB(result.m_B,result.m_invB,this->m_r);

  // eventually we compute new representation of r
  result.m_r = (result.m_invB * A) * m_r + result.m_invB * xx;
  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
  result.setLastEnclosure(enc);
}

template <typename MatrixType, typename Policies>
std::string C0AffineSet<MatrixType, Policies>::show(void) const {
  std::ostringstream descriptor;
  descriptor << "C0AffineSet:\n" << BaseSet::toString();
  return descriptor.str();
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0AFFINETONSET_HPP_

