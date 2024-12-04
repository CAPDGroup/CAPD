/// @addtogroup matrixAlgorithms
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file intMatrixAlgorithms.hpp
///
/// @author Marian Mrozek
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

#ifndef _CAPD_MATRIXALGORITHMS_INTMATRIXALGORITHMS_HPP_
#define _CAPD_MATRIXALGORITHMS_INTMATRIXALGORITHMS_HPP_

#include <capd/auxil/Logger.h>

#include "Invert.h"
#include "MatrixOp.h"
#include "capd/auxil/debuger.h"
#include "capd/vectalg/MatrixSlice.h"
#include "capd/vectalg/Matrix.h"
#include <cmath>
#include "CAPDIntMatrixAlgorithms.h"

namespace capd{
  namespace matrixAlgorithms{

            // *** partial reduction *** //
    typedef CAPDIntMatrixAlgorithms DefaultIntMatrixAlgorithms;
//typedef PARIIntMatrixAlgorithms DefaultIntMatrixAlgorithms;

/* ------------------------  ------------------------ */
template<class matrix>
void kernelImage(matrix& B,matrix& kernel,matrix& image){
  int m=B.numberOfRows();
  int n=B.numberOfColumns();
  matrix R=matrix::Identity(n);
  matrix Rinv = emptyMatrix<matrix>();
  int k;
  columnEchelon(B,R,Rinv,k);
  kernel=MatrixSlice<matrix>(R,1,n,k+1,n);
  image=MatrixSlice<matrix>(B,1,m,1,k);
}

/* ------------------------  ------------------------ */

template<class matrix, class sqMatrix1, class sqMatrix2>
void smithForm(matrix& B,
               sqMatrix1& Q, sqMatrix1& Qinv,
               sqMatrix2& R, sqMatrix2& Rinv,
               int &s,int &t)
{
  DefaultIntMatrixAlgorithms intMatrixAlgorithms;
  typedef typename DefaultIntMatrixAlgorithms::template SmithForm<matrix>::type SmithFormType;
  SmithFormType smithForm = intMatrixAlgorithms.smithForm(B, Q.numberOfRows() > 0, Qinv.numberOfRows() > 0,
                                                          R.numberOfRows() > 0, Rinv.numberOfRows() > 0);

  smithForm();

  Q = smithForm.getQ();
  Qinv = smithForm.getQinv();
  R = smithForm.getR();
  Rinv = smithForm.getRinv();

  s = smithForm.getS();
  t = smithForm.getT();
}

/* ------------------------  ------------------------ */

template<class matrix, class vector, class colVector>
bool solveLinearEquation(const matrix& A,const colVector& b,vector& x)
{
  const int m=A.numberOfRows();
  const int n=A.numberOfColumns();

  CAPD_TRACE("solveLinearEquation for matrix size: " << m << "x" << n);
  const bool trace = (n <= 100 && m <= 100);

  if (trace) {
    CAPD_TRACE("solveLinearEquation args: " << cppReprezentation(A, "A", "TYPE") << " "
	       << cppReprezentation(b, "b", "TYPE"));
  }

  DefaultIntMatrixAlgorithms intMatrixAlgorithms;
  typename DefaultIntMatrixAlgorithms::template SolveLinearEquation<const matrix>::type solver
    = intMatrixAlgorithms.solveLinearEquation(A);

  const bool result = solver(b, x);

  if (result) {
    if (trace) {
      CAPD_TRACE("solveLinearEquation result: " << cppReprezentation(x, "x", "TYPE"));
    }
  } else {
    CAPD_TRACE("solveLinearEquation result: " << cppReprezentation(vector((unsigned)0), "x", "TYPE"));
  }

  return result;
}


/* ------------------------  ------------------------ */

template<class matrix>
bool invert(const matrix& A,matrix& Ainv)
{
  Invert<matrix> invertOp(A);

  return invertOp(Ainv);
}


/* ------------------------  ------------------------ */
template<class matrix,class IntVector>
void quotientBaseMatrix(
  const matrix& A_W,                       // input: basis of G \subset Z^p as columns of A_W,
  const matrix& A_V,                       // input: basis of H\subset G \subset Z^p as columns of A_V,
  matrix& A_U,                             // output: pseudobasis of G/H as columns of A_U,
  IntVector& A_orders                      // output: IntVector of finite orders
                                           // of elements of pseudobasis A_U
){
  DefaultIntMatrixAlgorithms intMatrixAlgorithms;
  typedef typename DefaultIntMatrixAlgorithms::template QuotientBaseMatrix<const matrix>::type QuotientBaseMatrix;
  QuotientBaseMatrix quotientBaseMatrix = intMatrixAlgorithms.quotientBaseMatrix(A_W, A_V);
  quotientBaseMatrix();
  A_U = quotientBaseMatrix.pseudoBasis();
  quotientBaseMatrix.getOrders(A_orders);
}
/* ------------------------  ------------------------ */
template<class matrix>
void copy(const matrix &A, matrix &result, int row, int col){ // copy matrix A to matrix result starting at position
																														  // result[row][col]
	for(int i=0;i<A.numberOfRows();++i)
		for(int j=0;j<A.numberOfColumns();++j)
			result[row+i][col+j] = A[i][j];
}

/* ------------------------  ------------------------ */
template<class matrix>
void spaceIntersection(const matrix &A,const matrix &B,matrix &C){ // input: basis of space A as columns of A,
																																	 // input: basis of space B as columns of B,
																																	 // output: basis of A intersection B as columns of C
	int A_n = A.numberOfColumns();
	int A_m = A.numberOfRows();
	int B_n = B.numberOfColumns();
	//int B_m = B.numberOfRows();
	int D_n = A_n + B_n, D_m = A_m;
	matrix D(D_m, D_n);
	copy(A,D,0,0);
	copy(B,D,0,A_n);

	matrix kernel,image,Dcopy = D;
	kernelImage(Dcopy,kernel,image);
	D = MatrixSlice<matrix>(D,1,A_m,1,A_n);
	kernel = MatrixSlice<matrix>(kernel,1,A_n,1,kernel.numberOfColumns());
	C = D * kernel;
}

/* ------------------------  ------------------------ */

  } // end of namespace matrixAlgorithms

} // end of namespace capd;
#endif // _CAPD_MATRIXALGORITHMS_INTMATRIXALGORITHMS_HPP_

/// @}
