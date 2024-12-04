/// @addtogroup vectalg
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file Matrix_Interval.hpp
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_VECTALG_IMATRIX_HPP_
#define _CAPD_VECTALG_IMATRIX_HPP_

#include "capd/vectalg/iobject.hpp"
#include "capd/vectalg/Vector_Interval.hpp"

namespace capd{
namespace vectalg{

template<typename IMatrixType>
IMatrixType midMatrix(const IMatrixType& v)
{
  IMatrixType result(v.numberOfRows(),v.numberOfColumns(),true);
  mid(v,result); // defined in iobject.hpp
  return result;
}


}} // namespace capd::vectalg


namespace capd{
namespace matrixAlgorithms{

template<class MatrixType>
void krawczykCorrection(const MatrixType& A, MatrixType& invA)
{
  typedef typename MatrixType::size_type size_type;
  MatrixType C = midMatrix(invA);
  // compute T = C*A-Id
  MatrixType T = C*A;
  for(size_type i=1;i<=T.numberOfRows();++i)
    T(i,i) -= TypeTraits<typename MatrixType::ScalarType>::one();
  invA = intersection(C - T*invA,invA);
}

template<class MatrixType>
MatrixType krawczykInverse(const MatrixType& A)
{
  MatrixType invA = gaussInverseMatrix(A);
  krawczykCorrection(A,invA);
  return invA;
}

}
}
#endif // _CAPD_VECTALG_IMATRIX_HPP_

/// @}
