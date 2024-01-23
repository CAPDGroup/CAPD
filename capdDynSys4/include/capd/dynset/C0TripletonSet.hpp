/////////////////////////////////////////////////////////////////////////////
/// @file C0TripletonSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2012 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0TRIPLETONSET_HPP_
#define _CAPD_DYNSET_C0TRIPLETONSET_HPP_

#include <stdexcept>
#include <sstream>
#include "capd/vectalg/iobject.hpp"
#include "capd/geomset/CenteredDoubletonSet.hpp"
#include "capd/geomset/CenteredTripletonSet.hpp"
#include "capd/dynset/C0TripletonSet.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/vectalg/algebraicOperations.hpp"
#include "capd/vectalg/Matrix_Interval.hpp"

namespace capd {
namespace dynset {

template<typename MatrixType, typename Policies>
C0TripletonSet<MatrixType,Policies>::C0TripletonSet(const VectorType& x, ScalarType t)
  : SetType(x,VectorType(x.dimension()),t),
    BaseSet(x), Data(x.dimension())
{}

// -----------------------------------------------------------------------------------

template<typename MatrixType, typename Policies>
C0TripletonSet<MatrixType,Policies>::C0TripletonSet(const VectorType& x, const VectorType& r, ScalarType t)
  : SetType(x+r,VectorType(x.dimension()),t),
    BaseSet(x,r), Data(x.dimension())
{}

// -----------------------------------------------------------------------------------

template<typename MatrixType, typename Policies>
C0TripletonSet<MatrixType,Policies>::C0TripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, ScalarType t)
  : SetType(x+C*r0,VectorType(x.dimension()),t),
    BaseSet(x,C,r0), Data(x.dimension())
{}

// -----------------------------------------------------------------------------------

template<typename MatrixType, typename Policies>
C0TripletonSet<MatrixType,Policies>::C0TripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r, ScalarType t)
  : SetType(x+C*r0+r,VectorType(x.dimension()),t),
    BaseSet(x,C,r0,r), Data(x.dimension())
{}

// -----------------------------------------------------------------------------------

template<typename MatrixType, typename Policies>
C0TripletonSet<MatrixType,Policies>::C0TripletonSet(const VectorType& x, const MatrixType& C,
    const VectorType& r0, const MatrixType& B,
    const VectorType& r, ScalarType t
)
  : SetType(x+C*r0+B*r,VectorType(x.dimension()),t),
    BaseSet(x,C,r0,B,r), Data(x.dimension())
{}

// -----------------------------------------------------------------------------------

template<typename MatrixType, typename Policies>
C0TripletonSet<MatrixType,Policies>::C0TripletonSet(
                const VectorType& x,
                const MatrixType& C, const VectorType& r0,
                const MatrixType& B, const VectorType& r,
                const MatrixType& Q, const VectorType& q,
                ScalarType t)
  : SetType(x+C*r0+intersection(B*r,Q*q),VectorType(x.dimension()),t),
    BaseSet(x,C,r0,B,r), Data(x.dimension())
{}

// -----------------------------------------------------------------------------------

template <typename MatrixType, typename Policies>
void C0TripletonSet<MatrixType, Policies>::move(DynSysType & dynsys) {
  move(dynsys, *this);
}

// -----------------------------------------------------------------------------------

template <typename MatrixType, typename Policies>
void C0TripletonSet<MatrixType, Policies>::move(DynSysType & dynsys, C0TripletonSet& result) const
{
  // the following function can throw an exception leaving output parameters in an inconsistent state
  // do not overwrite parameters of the set until we are sure that they are computed correctly
  if(subset(this->m_x,this->m_currentSet)){
    result.x = this->m_x;
    result.deltaX = this->m_currentSet-this->m_x;
  }else
    split(this->m_currentSet, result.x, result.deltaX);

  dynsys.encloseC0Map(this->getCurrentTime(),result.x,this->m_currentSet, result.y, result.rem, result.enc, result.jacPhi);

  C0TripletonSet::move(*this,result,result.m_currentSet,result);

  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
  result.setLastEnclosure(result.enc);
  if(this->isReorganizationNeeded(result)){
    this->reorganize(result);
    result.m_Q = result.m_B;
    result.m_q = result.m_r;
  }
}

// -------------------------------------------------------------

template <typename MatrixType, typename Policies>
void C0TripletonSet<MatrixType, Policies>::move(const BaseSet& set, BaseSet& result, VectorType& bound, Data& data)
{
  // first enclosure of the image is given by simply interval evaluation
  subtractObjects( set.m_x, data.x, data.deltaY );
  addObjects( data.y, data.rem, result.m_x );
  addObjects( result.m_x, data.jacPhi*data.deltaX, bound );

  // another enclosure - intersect it with previous enclosure
  addAssignObjectObject( result.m_x, data.jacPhi*data.deltaY );
  result.m_C = data.jacPhi * set.m_C;
  matrixByMatrix( data.jacPhi, set.m_B, data.B );
  matrixByMatrix( data.jacPhi, set.m_Q, data.Q );
  result.m_B = data.B;

  // first we enclose intersection(B*r,Q*qr)
  intersection( data.B*set.m_r, data.Q*set.m_q, data.qr );

  // another enclosure is phi + C*r0 + qr
  if( !intersection(bound,result.m_x + result.m_C*set.m_r0 + data.qr, bound) ){
      std::ostringstream out;
      out << "C0TripletonSet: empty intersection of two representations. Report this error to CAPD developers!\n";
      out << "currentTime=" << dynamic_cast<C0TripletonSet&>(result).getCurrentTime() << std::endl;
      out << "bound=" << bound << std::endl;
      out << "\nbound=" << bound;
      out << "\nresult.m_x=" << result.m_x;
      out << "\nresult.m_r=" << result.m_r;
      out << "\nresult.m_r0=" << result.m_r0;
      out << "\nresult.m_C=" << result.m_C;
      out << "\nresult.m_B=" << result.m_B;
      out << "\ndata.qr=" << data.qr << std::endl ;
      out << "\nresult.m_x + result.m_C*set.m_r0 + data.qr=" << result.m_x + result.m_C*set.m_r0 + data.qr << std::endl ;
      throw std::logic_error(out.str());
  }
  // here we compute representation for the new set
  split(result.m_C, data.deltaC);
  split(result.m_x, data.y);
  addAssignMatrixByVector(data.y,data.deltaC,set.m_r0);
  // we assume that Policies provides algorithms for computation
  // of B, its inverse invB and updates
  Policies policies;
  policies.computeBinvB(result.m_B,result.m_invB,set.m_r+data.y);

  // eventually we compute new representation of r
  matrixByMatrix(result.m_invB , data.B, data.deltaC);
  result.m_r = data.deltaC * set.m_r;
  addAssignMatrixByVector(result.m_r,result.m_invB,data.y);

  // propagate by inverse matrix
  try{
    mid( data.Q, result.m_Q );
    result.m_invQ = capd::matrixAlgorithms::krawczykInverse(result.m_Q);
    result.m_q = (result.m_invQ*data.Q)*set.m_q + result.m_invQ*data.y;
  }catch(...)
  {
    result.m_Q = result.m_B;
    result.m_q = result.m_r;
    result.m_invQ = result.m_invB;
  }

  // copy r0 if necessary
  if(&result != &set)
    result.m_r0 = set.m_r0;
}

template <typename MatrixType, typename Policies>
std::string C0TripletonSet<MatrixType, Policies>::show(void) const {
  std::ostringstream descriptor;
  descriptor << name() << ":\n";
  descriptor << " x=" << this->m_x << "\n B=" << this->m_B << "\n r=" << this->m_r << " ";
  descriptor << " m_Q = " << this->m_Q << "\n" << "m_q=" << this->m_q << " "
             << "\n C="  << this->m_C
             << "\n r0=" << this->m_r0 << " ";
  return descriptor.str();
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0TRIPLETONSET_HPP_
