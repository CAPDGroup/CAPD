/////////////////////////////////////////////////////////////////////////////
/// @file CenteredTripletonSet.hpp
///
/// @author Daniel Wilczak
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2014 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_GEOMSET_CENTEREDTRIPLETONSET_HPP_
#define _CAPD_GEOMSET_CENTEREDTRIPLETONSET_HPP_

#include <stdexcept>
#include <sstream>
#include "capd/vectalg/iobject.hpp"
#include "capd/geomset/CenteredDoubletonSet.hpp"
#include "capd/geomset/CenteredTripletonSet.h"

namespace capd {
namespace geomset {
/// @addtogroup geomset
/// @{

template<typename MatrixType>
void CenteredTripletonSet<MatrixType>::copyQ(){
  this->m_q = this->m_r;
  this->m_Q = this->m_B;
  this->m_invQ = this->m_invB;
}

// -------------------------------------------------------------

template<typename MatrixType>
CenteredTripletonSet<MatrixType>::CenteredTripletonSet(const VectorType& x)
  : BaseSet(x),
    m_q(x.dimension()),
    m_Q(x.dimension(),x.dimension()),
    m_invQ(x.dimension(),x.dimension())
{
  copyQ();
}

// -------------------------------------------------------------

template<typename MatrixType>
CenteredTripletonSet<MatrixType>::CenteredTripletonSet(const VectorType& x, const VectorType& r)
  : BaseSet(x,r),
    m_q(x.dimension()),
    m_Q(x.dimension(),x.dimension()),
    m_invQ(x.dimension(),x.dimension())
{
  copyQ();
}

// -------------------------------------------------------------

template<typename MatrixType>
CenteredTripletonSet<MatrixType>::CenteredTripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0)
  : BaseSet(x,C,r0),
    m_q(x.dimension()),
    m_Q(x.dimension(),x.dimension()),
    m_invQ(x.dimension(),x.dimension())
{
  copyQ();
}

// -------------------------------------------------------------

template<typename MatrixType>
CenteredTripletonSet<MatrixType>::CenteredTripletonSet(const VectorType& x, const MatrixType& C, const VectorType& r0, const VectorType& r)
  : BaseSet(x,C,r0,r),
    m_q(x.dimension()),
    m_Q(x.dimension(),x.dimension()),
    m_invQ(x.dimension(),x.dimension())
{
  copyQ();
}

// -------------------------------------------------------------

template<typename MatrixType>
CenteredTripletonSet<MatrixType>::CenteredTripletonSet(const VectorType& x, const MatrixType& C,
    const VectorType& r0, const MatrixType& B,
    const VectorType& r
)
  :BaseSet(x,C,r0,B,r),
    m_q(x.dimension()),
    m_Q(x.dimension(),x.dimension()),
    m_invQ(x.dimension(),x.dimension())
{
  copyQ();
}

// -------------------------------------------------------------

template<typename MatrixType>
CenteredTripletonSet<MatrixType>::CenteredTripletonSet(
                const VectorType& x,
                const MatrixType& C, const VectorType& r0,
                const MatrixType& B, const VectorType& r,
                const MatrixType& Q, const VectorType& q
                )
  : BaseSet(x,C,r0,B,r),
    m_q(q),
    m_Q(Q),
    m_invQ(x.dimension(),x.dimension())
{
  m_invQ = capd::matrixAlgorithms::krawczykInverse(Q);
}

// -------------------------------------------------------------

template <typename MatrixType>
typename CenteredTripletonSet<MatrixType>::VectorType
CenteredTripletonSet<MatrixType>::affineTransformation(
      const MatrixType& A_M, const VectorType& A_C
  ) const {
  return A_M*(this->m_x-A_C) + (A_M*this->m_C)*this->m_r0 + intersection((A_M*this->m_B)*this->m_r,(A_M*this->m_Q)*this->m_q);
}

// -------------------------------------------------------------

template <typename MatrixType>
typename CenteredTripletonSet<MatrixType>::ScalarType
CenteredTripletonSet<MatrixType>::evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const {
  ScalarType r = TypeTraits<ScalarType>::zero();
  ScalarType rB = TypeTraits<ScalarType>::zero();
  ScalarType rQ = TypeTraits<ScalarType>::zero();

  for(size_type i=0;i<this->m_x.dimension();++i){
    rB += (gradient*this->m_B.column(i))*this->m_r[i];
    rQ += (gradient*this->m_Q.column(i))*this->m_q[i];
    r += (gradient*this->m_C.column(i))*this->m_r0[i];
    r += gradient[i]*(this->m_x[i] - x0[i]);
  }
  ScalarType y;
  if(!intersection(rB,rQ,y)){
    std::ostringstream message;
    message << "CenteredTripletonSet::evalAffineFunctional - empty intersection of rB and rQ. Report this error to CAPD developers!\n";
    message << "rB=" << rB << std::endl;
    message << "rQ=" << rQ << std::endl;
    throw std::logic_error(message.str());
  }
  return r+y;
}
/// @}
}} // namespace capd::geomset

#endif // _CAPD_GEOMSET_CENTEREDTRIPLETONSET_HPP_

