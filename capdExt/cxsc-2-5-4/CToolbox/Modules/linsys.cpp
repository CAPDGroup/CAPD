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

/* CVS $Id: linsys.cpp,v 1.14 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: linsys (implementation)
// Purpose: Computation of a verified solution of a square linear system of
//    equations A*x = b with full real matrix A and real right-hand side b.
// Method: Transformation of A*x = b to fixed-point form and applying an
//    interval residual iteration.
// Global functions:
//    LinSolve()      : to get a verified enclosure of the solution (two
//                      versions)
//    LinSolveErrMsg(): to get an error message text
//----------------------------------------------------------------------------
// In the comments below the notations #*(...) and ##(...) are used to indi-
// cate the evaluation of the expression specified in round brackets using
// the exact scalar product. I.e. each component of the result is computed
// with only one final rounding. The symbol #* holds for rounding to nearest
// whereas ## holds for rounding to the smallest enclosing interval. An exact
// scalar product may be implemented using dotprecision accumulators.
//----------------------------------------------------------------------------
#include <string.h>         // String handling
#include <imatrix.hpp>      // Interval matrix/vector arithmetic
#include <mvi_util.hpp>     // Interval matrix/vector utility functions
#include <matinv.hpp>       // Matrix inversion
#include <linsys.hpp>

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,   // No error occurred
  NotSquare    = 1,   // System to be solved is not square
  DimensionErr = 2,   // Dimensions of A and b are not compatible
  InvFailed    = 3,   // System is probably singular
  VerivFailed  = 4;   // Verification failed, system is probably
                      // ill-conditioned

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* LinSolveErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
      case NotSquare:
        strcpy(Hlp,"System to be solved is not square"); break;
      case DimensionErr:
        strcpy(Hlp,"Dimensions of A and b are not compatible"); break;
      case InvFailed:
        strcpy(Hlp,"System is probably singular"); break;
      case VerivFailed:
        strcpy(Hlp,"Verification failed, system is probably ill-conditioned");
        break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // LinSolveErrMsg

//----------------------------------------------------------------------------
// Computes an upper bound for the maximum norm of a real matrix 'M'.
//----------------------------------------------------------------------------
static real MaxNorm ( const rmatrix& M )
{
  int          i, j;
  real         Max, Tmp;
  dotprecision Accu;

  Max = 0.0;
  for (i = Lb(M,1); i <= Ub(M,1); i++) {
    // Compute Tmp = #>( sum(j=1(1)n,|M_ij|) )
    //----------------------------------------
    Accu = 0.0;
    for (j = Lb(M,2); j <= Ub(M,2); j++)
      Accu += abs(M[i][j]);
    Tmp = rnd(Accu,RND_UP);
    if (Tmp > Max) Max = Tmp;
  }
  return Max;
} // MaxNorm

//----------------------------------------------------------------------------
// The vectors x and y are successive approximations for the solution of a
// linear system of equations computed by iterative refinement. If a compo-
// nent of y is diminished by more than 'Factor', it is a good candidate for
// a zero entry. Thus, it is set to zero.
//----------------------------------------------------------------------------
static void CheckForZeros ( rvector& x, rvector& y )
{
  const real Factor = 1E+5;
  int        i;

  for (i = Lb(y); i <= Ub(y); i++)
    if ( abs(y[i])*Factor < abs(x[i]) ) y[i] = 0.0;
}

//----------------------------------------------------------------------------
// The vectors x and y are successive iterates. The function returns true if
// the relative error of all components x_i and y_i is <= 10^(-12), i.e. y_i
// has about 12 correct decimals. If x_i or y_i vanishes, the relative error
// is implicitly set to zero.
//----------------------------------------------------------------------------
static int Accurate ( rvector& x, rvector& y )
{
  const real Delta = 1E-12;   // Relative error bound
  int        i, n, ok;
  real       abs_yi;

  i = Lb(y); n = Ub(y);
  do {
    if (sign(x[i])*sign(y[i]) <= 0)       // Relative error set to 0
      ok = 1;
    else {
     abs_yi = abs(y[i]);                // Relative error > Delta ?
     ok = (abs(y[i] - x[i]) <= Delta * abs_yi );
    }
    i++;
  } while (ok && (i <= n));
  return ok;
} // Accurate

//----------------------------------------------------------------------------
// This function 'VerificationStep()' performs the iteration
// [y] = Blow([x],Eps), [x] = [z] + [C]*[y] for k = 1,2,... until the new
// iterate [x] lies in the interior of [y] or until the maximum number of
// iterations is exceeded. The flag 'IsVerified' is set if an inclusion in
// the interior could be established.
//----------------------------------------------------------------------------
static void VerificationStep ( ivector& xx, ivector& zz,
                               imatrix& C,  int& IsVerified )
{
  const int  MaxIter = 3;          // Maximal number of iteration steps
  const real Epsilon = 1000.0;     // Factor for the epsilon inflation

  int     p;
  ivector yy(Lb(xx),Ub(xx));

  xx = zz; p = 0;                      // Initialize:  [x] = [z]
  do {
    yy = Blow(xx,Epsilon);             // Epsilon inflation
    xx = zz + C*yy;                    // New iterate: [x] = [z]+[C]*[y]
    IsVerified = in(xx,yy);            // Inclusion in the interior?
    p++;
  } while (!IsVerified && (p < MaxIter));
}

//----------------------------------------------------------------------------
// Purpose: The function 'LinSolveMain()' computes a verified solution of a
//    square linear system of equations A*x=b.
// Parameters:
//    In : 'A'          : matrix of the system, passed as reference
//                        parameter to save time for copying it
//         'b'          : right-hand side of the system, passed as
//                        reference parameter to save time for copying it
//         'ComputeCond': flag signalling whether a condition number
//                        estimate is to be computed
//    Out: 'xx'         : enclosure of the unique solution, resized to
//                        standard index range with lower index bound 1
//         'Cond'       : condition number estimate
//         'Err'        : error code
// Description: An approximate inverse 'R' of 'A' is computed by calling
//   function 'MatInv()'. Then an approximate solution 'x' is computed
//   applying a conventional real residual iteration. For the final verifica-
//   tion, an interval residual iteration is performed. An enclosure of the
//   unique solution is returned in the interval vector 'xx'. The function
//   also returns a condition number estimate 'Cond' if the flag 'ComputeCond'
//   is set. 'Cond' is initialized by -1. A negative value for 'Cond' signals
//   that an estimate could not be computed.
//----------------------------------------------------------------------------
static void LinSolveMain ( rmatrix  A,
                           rvector  b,
                           ivector&  xx,
                           real&     Cond,
                           int       ComputeCond,
                           int&      Err )
{
  const int
    MaxResCorr = 10;    // Maximal number of real residual corrections

  int     n = ColLen(A);                // Length of the columns of 'A'
  int     m = RowLen(A);                // Length of the rows of 'A'
  int     IsVerified, i, j, k;
  int     bLow, ARowLow, AColLow;       // Lower index bounds od 'A' and 'b'
  rmatrix R;                            // To store the inverse of 'A'

  Cond = -1.0;                          // Initial value for the condition

  if (n != m)                           // Error: 'A' is not square
    { Err = NotSquare; return; }

  if (n != VecLen(b))                   // Error: Dimensions of 'A' and 'b'
    { Err = DimensionErr; return; }     //        are not compatible

  MatInv(A,R,Err);
  if (Err != NoError)                   // Error: Inversion failed
    { Err = InvFailed; return; }

  // Start algorithm
  //----------------
  rvector       x(n), y(n), d(n);       // Allocate dynamic arrays
  ivector       yy(n), zz(n);
  imatrix       C(n,n);
  dotprecision  Accu;                   // Allocate accumulators
  idotprecision IAccu;

  if (ComputeCond)                      // Compute condition number
    Cond = mulu(MaxNorm(A),MaxNorm(R)); //-------------------------

  // Normalize index range of 'A', 'R' and 'b' to standard range 1..n
  // and resize 'xx' to length 'n'. Save lower bounds of 'b' and 'A'
  // to restore them before leaving 'LinSolveMain()'.
  //-----------------------------------------------------------------
  Resize(xx,n);
  bLow = Lb(b); SetLb(b,1);
  ARowLow = Lb(A,ROW); SetLb(A,ROW,1);
  AColLow = Lb(A,COL); SetLb(A,COL,1);
  SetLb(R,ROW,1); SetLb(R,COL,1);

  x = R*b; k = 0;   // Real residual iteration
  do {              //------------------------
    y = x;
    for (i = 1; i <= n; i++) {   // Compute: d = #*(b-A*y)
      Accu = b[i];               //-----------------------
      accumulate(Accu,-A[i],y);
      d[i] = rnd(Accu);
    }
    for (i = 1; i <= n; i++) {   // Compute: x = #*(y+R*d)
      Accu = y[i];               //-----------------------
      accumulate(Accu,R[i],d);
      x[i] = rnd(Accu);
    }
    CheckForZeros(y,x);
    k++;
  } while (!Accurate(y,x) && (k < MaxResCorr));

  // Prepare verification step, i.e. compute enclosures [C]
  // and [z] of C = (I - R*A) and z = R*(b - A*x).
  //-------------------------------------------------------
  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++) {            // Compute [C] = ##(I-R*A)
      Accu = (i == j) ? 1.0 : 0.0;        //------------------------
      accumulate(Accu,-R[i],A[Col(j)]);
      rnd(Accu,C[i][j]);
    }
  for (i = 1; i <= n; i++) {              // Compute: d = #*(b-A*x)
    Accu = b[i];                          //-----------------------
    accumulate(Accu,-A[i],x);
    d[i] = rnd(Accu);
  }
  for (i = 1; i <= n; i++) {              // Compute: [y] = ##(b-A*x-d)
    Accu = b[i];                          //---------------------------
    accumulate(Accu,-A[i],x);
    Accu -= d[i];
    rnd(Accu,yy[i]);
  }

  // Now b-A*x lies in d+[y]. Thus R*(b-A*x) lies in [z] = ##(R*d+R*[y])
  //--------------------------------------------------------------------
  for (i = 1; i <= n; i++) {
    IAccu = 0.0;
    accumulate(IAccu,R[i],d);
    accumulate(IAccu,R[i],yy);
    zz[i] = rnd(IAccu);
  }

  // If R*(b-A*x) = 0 then x is an exact solution else try to find a
  // verified enclosure [x] for the absolute error of x.
  //----------------------------------------------------------------
  if (Zero(zz))
    xx = x;
  else {
    VerificationStep(xx,zz,C,IsVerified);   // Attempt to compute [x]
    if ( !IsVerified )
      Err = VerivFailed;
    else
      xx = x + xx;                          // The exact solution lies x+[x]
  }

  // Restore lower index bounds of 'A' and 'b'
  //------------------------------------------
  SetLb(b,bLow);
  SetLb(A,ROW,ARowLow); SetLb(A,COL,AColLow);
} //LinSolveMain

//----------------------------------------------------------------------------
// Purpose: The function 'LinSolve()' computes a verified solution of a
//    square linear system of equations A*x=b without returning a condition
//    number estimate.
// Parameters:
//    In : 'A'  : matrix of the system, passed as reference
//                parameter to save time for copying it
//         'b'  : right-hand side of the system, passed as
//                reference parameter to save time for copying it
//    Out: 'xx' : enclosure of the unique solution
//         'Err': error code
// Description: Calls 'LinSolveMain()' for solving the linear system with the
//    flag 'ComputeCond' not set.
//----------------------------------------------------------------------------
void LinSolve ( const rmatrix& A, const rvector& b, ivector& xx, int& Err )
{
  real DummyCond;         // Dummy parameter for call of 'LinSolveMain()'

  LinSolveMain(A,b,xx,DummyCond,0,Err);
}

//----------------------------------------------------------------------------
// Purpose: The function 'LinSolve()' computes a verified solution of a
//    square linear system of equations A*x=b and returns a condition
//    number estimate.
// Parameters:
//    In : 'A'   : matrix of the system, passed as reference
//                 parameter to save time for copying it
//         'b'   : right-hand side of the system, passed as
//                 reference parameter to save time for copying it
//    Out: 'xx'  : enclosure of the unique solution
//         'Cond': condition number estimate
//         'Err' : error code
// Description: Calls 'LinSolveMain()' for solving the linear system with the
//    flag 'ComputeCond' set.
//----------------------------------------------------------------------------
void LinSolve ( const rmatrix& A, const rvector& b, ivector& xx, real& Cond, int& Err )
{
  LinSolveMain(A,b,xx,Cond,1,Err);
}




