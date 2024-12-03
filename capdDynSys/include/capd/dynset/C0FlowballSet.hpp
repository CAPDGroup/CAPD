/////////////////////////////////////////////////////////////////////////////
/// @file C0FlowballSet.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_DYNSET_C0FLOWBALLSET_HPP_
#define _CAPD_DYNSET_C0FLOWBALLSET_HPP_

#include <stdexcept>
#include "capd/dynset/C0FlowballSet.h"
#include "capd/vectalg/Norm.hpp"
#include "capd/vectalg/Vector_Interval.hpp"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"

namespace capd{
namespace dynset{

template<typename MatrixType>
C0FlowballSet<MatrixType>::C0FlowballSet(const VectorType& x, const ScalarType& r, const NormType& aNorm, ScalarType t)
  : C0Set<MatrixType>(intervalBall(x,r),VectorType(x.dimension()),t),
    m_x(x),m_r(r),m_n(aNorm.clone())
{}

template<typename MatrixType>
void C0FlowballSet<MatrixType>::move(capd::dynsys::DynSys<MatrixType> &dynsys)
{
  this->move(dynsys,*this);
}

template<typename MatrixType>
void C0FlowballSet<MatrixType>::move(capd::dynsys::DynSys<MatrixType> &dynsys, C0FlowballSet &result) const
{
  VectorType s(m_x.dimension());
  VectorType vi = VectorType(*this);
  result.m_r = m_r*dynsys.Lipschitz(this->getCurrentTime(),vi,*m_n);

  VectorType y = dynsys.Phi(this->getCurrentTime(),m_x);
  split(y,s);
  s += dynsys.Remainder(this->getCurrentTime(),m_x,vi);
  result.m_x = y;
  result.m_r += (*m_n)(s);
  result.setCurrentTime(this->getCurrentTime()+dynsys.getStep());
  result.setLastEnclosure(vi);
}

template<typename MatrixType>
std::string C0FlowballSet<MatrixType>::show(void) const
{
  std::ostringstream descriptor;
  descriptor << "C0FlowballSet(" << m_n->show() << "): x=";
  descriptor << m_x << " r=";
  descriptor << m_r << " ";
  return descriptor.str();
}

template<typename MatrixType>
typename C0FlowballSet<MatrixType>::VectorType C0FlowballSet<MatrixType>::affineTransformation(
    const MatrixType& A_M, const VectorType& A_C
  ) const
{
  return A_M*(intervalBall(m_x,m_r)-A_C);
}

}} // namespace capd::dynset

#endif // _CAPD_DYNSET_C0FLOWBALLSET_HPP_
