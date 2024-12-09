
/////////////////////////////////////////////////////////////////////////////
/// @file AffineSet.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.ii.uj.edu.pl/ for details.

#ifndef _CAPD_GEOMSET_AFFINESET_HPP_
#define _CAPD_GEOMSET_AFFINESET_HPP_

#include "capd/vectalg/Vector_Interval.hpp"
#include "capd/vectalg/algebraicOperations.hpp"
#include "capd/vectalg/iobject.hpp"
#include "capd/geomset/AffineSet.h"
#include "capd/matrixAlgorithms/floatMatrixAlgorithms.hpp"
#include "capd/vectalg/Matrix_Interval.hpp"

#include <sstream>

namespace capd {
namespace geomset {
/// @addtogroup geomset
/// @{


/// sets x:=0, r:=0, B:=Id
template<typename MatrixType>
AffineSet<MatrixType>::AffineSet(size_type dim) :
  m_x(dim), m_r(dim), m_B(MatrixType::Identity(dim)), m_invB(MatrixType::Identity(dim)) {
}

/// sets x:=mid(v),  B:=Id, r:=[-radius(v), radius(v)]
template<typename MatrixType>
AffineSet<MatrixType>::AffineSet(const VectorType& v) :
  m_x(v), m_r(v.dimension(), false), m_B(MatrixType::Identity(v.dimension())), m_invB(MatrixType::Identity(v.dimension())) {
  split(m_x, m_r);
}

/// sets x:=x,  B:=Id, r:=0
template<typename MatrixType>
AffineSet<MatrixType>::AffineSet(const VectorType& x, bool) :
  m_x(x), m_r(x.dimension()), m_B(MatrixType::Identity(x.dimension())), m_invB(MatrixType::Identity(x.dimension())) {
}

/// sets x:=x r:=r B:=Id
template<typename MatrixType>
AffineSet<MatrixType>::AffineSet(const VectorType& x, const VectorType& r) :
  m_x(x), m_r(r), m_B(MatrixType::Identity(x.dimension())), m_invB(MatrixType::Identity(x.dimension())) {
}
/// sets x:=x r:=r B:=B
template<typename MatrixType>
AffineSet<MatrixType>::AffineSet(const VectorType& x,
                                 const MatrixType& B, const VectorType& r)
  :  m_x(x), m_r(r), m_B(B)
{
  m_invB = capd::matrixAlgorithms::krawczykInverse(B);
}

template<typename MatrixType>
std::string AffineSet<MatrixType>::toString(void) const {
  std::ostringstream descriptor;
  descriptor << " x=" << m_x << "\n B=" << m_B << "\n r=" << m_r << " ";
  return descriptor.str();
}

/**
 *  makes affine transformation of the set with matrix M and center c
 *
 *   M * (set - c)
 */
template<typename MatrixType>
typename AffineSet<MatrixType>::VectorType AffineSet<MatrixType>::affineTransformation(const MatrixType& A_M,
                                                                                       const VectorType& A_C) const {
  return A_M * (m_x - A_C) + (A_M * m_B) * m_r;
}

template<typename MatrixType>
typename AffineSet<MatrixType>::ScalarType
AffineSet<MatrixType>::evalAffineFunctional(const VectorType& gradient, const VectorType& x0) const{
  ScalarType r = TypeTraits<ScalarType>::zero();
  for(size_type i=0;i<this->m_x.dimension();++i){
    r += gradient[i]*(this->m_x[i]-x0[i]);
    r += (gradient*this->m_B.column(i))*this->m_r[i];
  }
  return r;
}

/// @}
}
} // namespace capd::dynset

#endif // _CAPD_GEOMSET_AFFINESET_HPP_
