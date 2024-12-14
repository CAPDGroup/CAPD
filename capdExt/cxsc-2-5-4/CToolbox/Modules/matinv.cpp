/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: matinv.cpp,v 1.14 2014/01/30 17:49:27 cxsc Exp $ */

//============================================================================
//
//                              Program/Module
//                                   from
//                 C++ TOOLBOX FOR VERIFIED COMPUTING I
//                         Basic Numerical Problems
//
//      Copyright (c) 1995   Rolf Hammer, Matthias Hocks, Dietmar Ratz
//
// For details on theory, algorithms, and programs, see the book
//
//  R. Hammer, M. Hocks, U. Kulisch, D. Ratz:  C++ Toolbox for
//  Verified Computing I - Basic Numerical Problems. Springer-Verlag,
//  Heidelberg, New York, 1995.
//
//============================================================================
//----------------------------------------------------------------------------
// File: matinv (implementation)
// Purpose: Computation of an approximate inverse of a real square matrix.
// Method: LU decomposition applying Crout's algorithm.
// Global functions:
//   MatInv()      : matrix inversion
//   MatInvErrMsg(): to get an error message text
//----------------------------------------------------------------------------
// In the comments below the notations #*(...) and ##(...) are used to indi-
// cate the evaluation of the expression specified in round brackets using
// the exact scalar product. I.e. each component of the result is computed
// with only one final rounding. The symbol #* holds for rounding to nearest
// whereas ## holds for rounding to the smallest enclosing interval. An exact
// scalar product may be implemented using dotprecision accumulators.
//----------------------------------------------------------------------------
#include <string.h>       // String handling
#include <matinv.hpp>

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,     // No error occurred
  NotSquare    = 1,     // Matrix to be inverted is not square
  Singular     = 2;     // Matrix to be inverted is probably singular

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* MatInvErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
      case NotSquare:
        sprintf(Hlp,"Matrix to be inverted is not square");           break;
      case Singular:
        sprintf(Hlp,"Inversion failed, matrix is probably singular"); break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // MatInvErrMsg

//----------------------------------------------------------------------------
// Purpose: The function 'MatInv()' may be used for the computation of an
//    approximate inverse of a real square matrix.
// Parameters:
//    In : 'A'  : matrix to be inverted
//    Out: 'R'  : approximate inverse
//         'Err': error code
// Description:
//    Inversion of a regular matrix A stored in 'A' using LU decomposition.
//    For LU decomposition, formally a permutation matrix P is determined so
//    that the product P*A may be decomposed into a lower-triangular matrix L
//    and an upper-triangular matrix U with P*A = L*U. The diagonal elements
//    of L are set to 1. Using Crout's algorithm, the elements of both matri-
//    ces L and U are stored by temporary overwriting a copy of the input
//    matrix 'A'. The permutation matrix P is not explicitly generated. The
//    indices of row interchanges are stored in the index vector 'p' instead.
//    The i-th element of P*b may be accessed indirectly using the p[i]-th
//    entry of 'b'. Accurate expressions are used to avoid cancellation
//    errors. The k-th column of the inverse R of P*A is computed by
//    forward/backward substitution with the k-th unit vector e_k as the
//    right-hand side of the system: U*y = P*e_k ==> y, L*x = y ==> x. For
//    error codes, see above.
//----------------------------------------------------------------------------
void MatInv ( rmatrix A, rmatrix& R, int& Err )
{
  const real
    Tiny = 1E-200;          // A divisor less than 'Tiny' is handled like zero

  dotprecision Accu;                       // Long accumulator
  real    Max, Temp;                       // Help variables
  int     p1 = Lb(A,1), q1 = Lb(A,2);      // Lower bounds of 'A'.
  int     pn = Ub(A,1), qm = Ub(A,2);      // Upper bounds of 'A'.
  int     n = pn-p1+1;                     // Length of the rows of 'A'
  int     m = qm-q1+1;                     // Length of the columns of 'A'
  int     i, j, k, l, kk;                  // For loops

  Err = NoError;                           // Initialization

  if (n != m)                              // Error: 'A' is not square
    { Err = NotSquare; return; }

  // Normalize index range of 'A' to standard range (1..n,1..n)
  //-----------------------------------------------------------
  SetLb(A,ROW,1);  SetLb(A,COL,1);
  Resize(R,p1,pn,q1,qm);                   // Resize 'R' to same shape as 'A'

  if (n == 2) {                            // Special case: (2,2)-matrix
    Accu = 0.0;                            // Compute the determinant of 'A'
    accumulate(Accu, A[1][1],A[2][2]);
    accumulate(Accu,-A[2][1],A[1][2]);
    Temp = rnd(Accu);

    if (abs(Temp) < Tiny)                  // Error: 'A' is probably singular
      { Err = Singular; return; }

    R[p1][q1] =  A[2][2] / Temp;  R[p1][qm] = -A[1][2] / Temp;
    R[pn][q1] = -A[2][1] / Temp;  R[pn][qm] =  A[1][1] / Temp;
    return;
  }

  // Start usual case: Dimension of 'A' > 2
  //---------------------------------------
  rvector v(n), x(n);          // Help vectors
  int*    p = new int [n+1];   // Dynamic array of type integer
                               // Note: p[0] not used !

  for (i = 1; i <= n; i++)  p[i] = i; // Initializing permutation vector

  // Start LU factorization
  //-----------------------
  i = 1;
  while ( (Err == NoError) && (i <= n) ) {
    // Compute the numerators of those elements of the i-th column
    // of L which are not updated so far and store them in 'v'.
    //------------------------------------------------------------
    for (k = i; k <= n; k++) {   // v_k = #*(A_ki - sum(j=1(1)i-1,A_kj*A_ji))
      Accu = A[k][i];            //------------------------------------------
      for (j = 1; j < i; j++) accumulate(Accu,-A[k][j],A[j][i]);
      v[k] = rnd(Accu);
    }

    // Look for the column pivot
    //--------------------------
    j = i;  Max = abs(v[i]);
    for (k = i+1; k <= n; k++)
      if ( (Temp = abs(v[k])) > Max ) { j = k; Max = Temp; }

    // Swap rows of 'A' and 'v', store the permutation in 'p'
    //-------------------------------------------------------
    if (j != i) {
      x = A[i];     A[i] = A[j];  A[j] = x;
      k = p[i];     p[i] = p[j];  p[j] = k;
      Temp = v[i];  v[i] = v[j];  v[j] = Temp;
    }

    if (Max < Tiny)                   // Pivot element < 'Tiny', inversion
      { Err = Singular; return; }     // failed, matrix 'A' assumed to be
                                      // singular
    Temp = v[i];
    A[i][i] = Temp;    // Update the diagonal element of U

    // Update U's row and L's column elements
    //---------------------------------------
    for (k = i+1; k <= n; k++) { // A_ik = #*(A_ik - sum(j=1(1)i-1,A_ij*A_jk))
      Accu = A[i][k];            //-------------------------------------------
      for (j = 1; j < i; j++) accumulate(Accu,-A[i][j],A[j][k]);
      A[i][k] = rnd(Accu);
      A[k][i] = v[k] / Temp;
    }
    i++;
  } // while

  // Now 'A' is overwritten with the subdiagonal elements of L in its
  // lower left triangle and with the elements of U in its diagonal and
  // its upper right triangle. The elements of the inverse matrix are
  // computed column by column using forward/backward substitution.
  //-------------------------------------------------------------------
  for (k = 1; k <= n; k++) {
    // Forward substitution: L*x = P*e_k, where e_k is the k-th unit
    // vector. Note: If P*e_k has m leading zeros, this results in
    // x_i = 0 for 1,..,l-1 and x_l = 1. Thus, forward substitution
    // starts at index l+1.
    //--------------------------------------------------------------
    l = 1;
    while (p[l] != k) { x[l] = 0.0; l++; }
    x[l] = 1.0;
    for (i = l+1; i <= n; i++) {    // x_i = - #*(sum(j=1(1)i-1,A_ij*x_j))
      Accu = 0.0;                   //------------------------------------
      for (j = l; j < i; j++) accumulate(Accu,A[i][j],x[j]);
      x[i] = -rnd(Accu);
    }

    // Backward substitution: U * x = x, where the right-hand side is
    // the result of the forward substitution. It will be overwritten
    // by the solution of the system, the k-th column of the inverse
    // matrix.
    //---------------------------------------------------------------
    kk = q1+k-1;                 // Index offset for the k-th column of R
    for (i = n; i >= 1; i--) {   // x_i = #*(x_i-sum(j=i+1(1)n,A_ij*x_j)/A_ii)
      Accu = x[i];               //-------------------------------------------
      for (j = i+1; j <= n; j++) accumulate(Accu,-A[i][j],x[j]);
      x[i] = rnd(Accu) / A[i][i];
      R[p1+i-1][kk] = x[i];      // Remember index offset !
    }
  } // for (k = 1,...

  delete [] p;   // Free dynamically allocated memory
} // MatInv




