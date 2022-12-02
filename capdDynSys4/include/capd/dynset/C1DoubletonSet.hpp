/////////////////////////////////////////////////////////////////////////////
/// @file C1DoubletonSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C1DOUBLETONSET_HPP_
#define _CAPD_DYNSET_C1DOUBLETONSET_HPP_

#include <stdexcept>
#include "capd/vectalg/iobject.hpp"
#include "capd/dynset/C1DoubletonSet.h"
#include "capd/geomset/CenteredDoubletonSet.hpp"
#include "capd/geomset/MatrixDoubletonSet.hpp"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/dynset/C0DoubletonSet.hpp"

namespace capd {
namespace dynset {

template<typename MatrixType, typename Policies>
C1DoubletonSet<MatrixType,Policies>::C1DoubletonSet(const VectorType& x, ScalarType t)
  : SetType(
      x,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x),
    C1BaseSet(x.dimension()),
    Data(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1DoubletonSet<MatrixType,Policies>::C1DoubletonSet(const VectorType& x, const VectorType& r0, ScalarType t)
  : SetType(
      x+r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x, r0),
    C1BaseSet(x.dimension()),
    Data(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1DoubletonSet<MatrixType,Policies>::C1DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t)
  : SetType(
      x+C*r0,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x, C, r0),
    C1BaseSet(x.dimension()),
    Data(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1DoubletonSet<MatrixType,Policies>::C1DoubletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t)
  : SetType(
      x+C*r0+r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x, C, r0, r),
    C1BaseSet(x.dimension()),
    Data(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1DoubletonSet<MatrixType,Policies>::C1DoubletonSet(
      const VectorType& x,
      const MatrixType& C, const VectorType& r0,
      const MatrixType& B, const VectorType& r,
      ScalarType t
   ) : SetType(
      x+C*r0+B*r,
      VectorType(x.dimension()),
      MatrixType::Identity(x.dimension()),
      MatrixType(x.dimension(),x.dimension()),
      t),
    C0BaseSet(x, C, r0, B, r),
    C1BaseSet(x.dimension()),
    Data(x.dimension())
{}

template<typename MatrixType, typename Policies>
C1DoubletonSet<MatrixType,Policies>::C1DoubletonSet(const C0BaseSet & c0part, const C1BaseSet& c1part, ScalarType t)
    : SetType(
       VectorType(c0part),
       VectorType(c0part.dimension()),
       MatrixType(c1part),
       MatrixType(c0part.dimension(),c0part.dimension()),
       t),
   C0BaseSet(c0part),
   C1BaseSet(c1part),
   Data(c0part.get_x().dimension())
{}


template <typename MatrixType, typename Policies>
void C1DoubletonSet<MatrixType, Policies>::move(DynSysType & dynsys) {
  move(dynsys, *this);
}

// ------------------------------------------------------------------------------

template <typename MatrixType, typename Policies>
void C1DoubletonSet<MatrixType, Policies>::move(const C1BaseSet& set, C1BaseSet& result, MatrixType& bound, Data& data)
{
  // important: here we assume that both m_r and m_r0 contains zero
  // this is assured by each constructor and each step of this algorithm

  const size_type dim = set.m_D.numberOfRows();
  MatrixType deltaM(dim,dim);

  // another enclosure - intersect it with previous enclosure
  result.m_D = data.jacPhi*set.m_D;
  data.B = result.m_Bjac = data.jacPhi*set.m_Bjac;
  result.m_Cjac = data.jacPhi*set.m_Cjac;

  if( !capd::vectalg::intersection(bound, result.m_D+result.m_Cjac*set.m_R0 + result.m_Bjac*set.m_R, bound) )
    throw std::runtime_error("C1DoubletonSet - fatal error! Intersection of two different representations of matrix is empty!");

  // here we compute representation for the new set
  // jacRem is unnecessary now
  split(result.m_D, deltaM);
  split(result.m_Cjac, data.deltaC);
  deltaM += data.deltaC * set.m_R0;

  // we assume that Policies provides algorithms for computation
  Policies policies;
  policies.computeBinvB(result.m_Bjac,result.m_invBjac,set.m_R+deltaM);

  // eventually we compute new representation of r
  result.m_R = ( result.m_invBjac * data.B) * set.m_R + result.m_invBjac * deltaM ;

  if(&result != &set)
    result.m_R0 = set.m_R0;
}

// ------------------------------------------------------------------------------

template <typename MatrixType, typename Policies>
void C1DoubletonSet<MatrixType, Policies>::move(DynSysType & dynsys, C1DoubletonSet& result) const
{
  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  if(subset(this->m_x,this->m_currentSet)){
    result.x = this->m_x;
    subtractObjects( this->m_currentSet, this->m_x, result.deltaX );
  } else
    split(this->m_currentSet,result.x,result.deltaX);

  dynsys.encloseC1Map(this->getCurrentTime(),
                    result.x, this->m_currentSet,                // input parameters
                    result.y, result.rem, result.enc,            // C^0 output
                    result.jacPhi, result.jacRem, result.jacEnc  // C^1 output
                    );

  C0DoubletonSet<MatrixType,Policies>::move(*this,result,result.m_currentSet,result);

  // first enclosure of the image is given by simply interval evaluation
  result.jacPhi += result.jacRem;
  result.m_currentMatrix = result.jacPhi*this->m_currentMatrix;
  C1DoubletonSet::move(*this,result,result.m_currentMatrix,result);

  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
  result.setLastEnclosure(result.enc);
  result.setLastMatrixEnclosure(result.jacEnc);

  this->Policies::reorganizeIfNeeded(result.m_B,result.m_invB,result.m_r,result.m_C,result.m_r0);
  this->Policies::reorganizeC1IfNeeded(result.m_Bjac,result.m_invBjac,result.m_R,result.m_Cjac,result.m_R0);
}

// ------------------------------------------------------------------------------

template <typename MatrixType, typename Policies>
std::string C1DoubletonSet<MatrixType, Policies>::show(void) const {
  std::ostringstream descriptor;
  descriptor << name() << "\n"
             << C0BaseSet::toString() << "\n"
             << C1BaseSet::toString();
  return descriptor.str();
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C1DOUBLETONSET_HPP_

