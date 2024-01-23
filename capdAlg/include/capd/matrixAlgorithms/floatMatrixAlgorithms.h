/// @addtogroup matrixAlgorithms
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file floatMatrixAlgorithms.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_MATRIXALGORITHMS_FLOATMATRIXALGORITHMS_H_
#define _CAPD_MATRIXALGORITHMS_FLOATMATRIXALGORITHMS_H_

#include <vector>
#include <stdexcept>
#include <iostream>
#include "capd/basicalg/TypeTraits.h"

namespace capd{
/// Matrix algorithms: Gauss elimination, orthonormalization, QR decomposition etc.
namespace matrixAlgorithms{

// ------------- Gauss elimination - auxiliary function ------------------------------- //

template<typename MatrixType, typename ResultType>
void gauss(MatrixType a, ResultType b, ResultType& result);

/**
 *  Gauss elimination
 *
 * this function solves equaiton A*x=b for x
 * where A is nonsingular matrix
 */
template<typename MatrixType, typename VectorType>
VectorType gauss(const MatrixType& A, const VectorType& b);


// -------------------------------- orthonormalize -------------------------

// Gramm-Schmit column orthonormalization
// with permutation of columns
template<typename MatrixType>
void orthonormalize(MatrixType& Q);


// Gramm-Schmit column orthonormalization
// with permutation of columns depending also on a vector v
template<typename MatrixType>
void orthonormalize(MatrixType& Q, const typename MatrixType::RowVectorType& v);

template<typename MatrixType>
void orthonormalize(MatrixType& Q, const MatrixType& R);

// -------------------------------- QR_decompose -------------------------

template<typename MatrixType>
void QR_decompose(const MatrixType& A, MatrixType& Q, MatrixType& R);

// -------------------------------- diagonalize -------------------------

template<typename MatrixType>
int symMatrixDiagonalize(
  const MatrixType& A,
  MatrixType& D,
  typename MatrixType::ScalarType diagonalizingRelTolerance = capd::TypeTraits<typename MatrixType::ScalarType>::epsilon()
);


/// this function computes bound for spectral radius of a symmetric matrix
/// first it computes matrix which has the same eigenvalues and which is close to diagonal,
/// next bound is computed from Gerschgorin theorem
template<typename MatrixType>
typename MatrixType::ScalarType spectralRadiusOfSymMatrix(
  const MatrixType & A,
  typename MatrixType::ScalarType diagonalizingRelTolerance = capd::TypeTraits<typename MatrixType::ScalarType>::epsilon()
);


/// this function computes bound for maximal eigenvalue of a symmetric matrix
/// first it computes matrix which has the same eigenvalues and which is close to diagonal,
/// next bound is computed from Gerschgorin theorem
template<typename MatrixType>
typename MatrixType::ScalarType maxEigenValueOfSymMatrix(
  const MatrixType &A,
  typename MatrixType::ScalarType diagonalizingRelTolerance = capd::TypeTraits<typename MatrixType::ScalarType>::epsilon()
);

template<typename MatrixType>
MatrixType matrixExp(
  const MatrixType& M,
  typename MatrixType::ScalarType tolerance = capd::TypeTraits<typename MatrixType::ScalarType>::epsilon()
);

/**
  Crout Decomposition of a matrix
  As a result matrix D is a lower triangle
  and G is an upper triangle with 1 on diagonal
*/
template<typename MatrixType>
void croutDecomposition(const MatrixType &A, MatrixType &D, MatrixType &G);

template<typename MatrixType>
MatrixType invLowerTriangleMatrix(const MatrixType &A);

template<typename MatrixType>
MatrixType invUpperTriangleMatrix(const MatrixType &A);


template<typename MatrixType>
MatrixType inverseMatrix(const MatrixType &A);

template<typename MatrixType>
MatrixType gaussInverseMatrix(const MatrixType& A);

}} // namespace capd::matrixAlgorithms

#endif // _CAPD_MATRIXALGORITHMS_FLOATMATRIXALGORITHMS_H_

/// @}
