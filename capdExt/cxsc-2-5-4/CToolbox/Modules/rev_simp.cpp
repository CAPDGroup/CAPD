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

/* CVS $Id: rev_simp.cpp,v 1.16 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: rev_simp (implementation)
// Purpose: Determine an optimal value 'z', an optimal basic index set 'v',
//    and an optimal solution vector 'x' for a linear programming problem
//    P = (A,b,c) given in the standard form:
//                   ( z = c^t * x = max! )
//         (LP)      (     A * x = b      )
//                   (       x >= 0       ).
// Method: An initial basic index set is determined using the 2-phases method.
//    Then the revised simplex method described by Syslo, Deo, and Kowalic in
//    'Discrete Optimization Algorithms with Pascal' [p.14], Prentice Hall,
//    New Jersey (1983) is used.
// Global function:
//    RevSimplex()      : determines the values z, v, x for the linear
//                        programming problem P = (A,b,c)
//    RevSimplexErrMsg(): to get an error message text
//----------------------------------------------------------------------------
#include <rev_simp.hpp>
#include <cstring>       // strcpy

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError           = 0,   // No error occurred
  WrongDimension    = 1,   // Wrong dimension of problem
  NoOptimalSolution = 2,   // No initial optimal solution found
  NoInitialVertex   = 3;   // No initial vertex found

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* RevSimplexErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
      case WrongDimension:
        strcpy(Hlp,"Wrong dimension of problem (e.g. m >= n)"); break;
      case NoOptimalSolution:
        strcpy(Hlp,"No initial optimal solution found"); break;
      case NoInitialVertex:
        strcpy(Hlp,"No initial vertex found"); break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // RevSimplexErrMsg

//----------------------------------------------------------------------------
// Purpose: Determination of the solutions of a linear optimization problem.
// Parameters:
//   In : 'A'  : matrix of constraints
//        'b'  : right-hand side of constraints
//        'c'  : objective function by its coefficients
//   Out: 'x'  : optimal solution vector
//        'v'  : optimal basic index set
//        'z'  : optimal value
//        'Err': error flag
// Description: For a detailed description of the revised simplex method see
//   Syslo, Deo, and Kowalic.
//   The result vectors 'x' and 'v' are automatically resized to standard
//   index range with lower index bounds 1.
//----------------------------------------------------------------------------
void RevSimplex ( rmatrix A, rvector b, rvector c,
                  rvector& x, intvector& v, real& z, int& Err )
{
  // Check for correct dimensions of the input data
  //-----------------------------------------------
  int  m = ColLen(A), n = RowLen(A);

  if ( m >= n || VecLen(b) != m || VecLen(c) != n )
    { Err = WrongDimension;  return; }

  SetLb(A,1,1); SetLb(A,2,1);           // Resize input data to standard index
  SetLb(b,1); SetLb(c,1);               // range with lower index bound 1
                                        //------------------------------------
  Resize(x,n);                                 // Resize the result vectors to
  Resize(v,m);                                 // appropriate length
                                               //-----------------------------
  // Start algorithm
  //----------------
  const real
    eps = 1E-30;    // if (|a| < eps), it is handled as (a == 0)

  real      min, s;
  rmatrix   T(m+2,n+2), U(m+2,m+2);
  rvector   w(n+2), y(n+2);
  int       ex, SolutionFound;
  int       i, j, k, l, p, q, phase;

  Err = NoError; phase = 1; l=0;                             // Initialization
  p = m + 2; q = m + 2; k = m + 1;                           //---------------
  SolutionFound = 0;

  T(m,n) = A;                                               // Compute tableau
                                                            //----------------
  for (j = 1; j <= n; j++) {
    T[k][j] = -c[j];
    s = 0.0;
    for (i = 1; i <= m; i++) s -= T[i][j];
    T[p][j] = s;
  }
  U = Id(U);

  s = 0.0; w = 0.0;                    // Determine initial basic index set v,
  for (i = 1; i <= m; i++) {           // initial optimal value s, and initial
    v[i] = n - m + i;        /*!!!*/   // optimal solution x
    w[i] = b[i];                       //-------------------------------------
    s -= b[i];
  }
  w[k] = 0.0; w[p] = s;

  do {
    if (w[p] >= -eps*10.0 && phase == 1)         // Only phase 1:
      { phase = 2; q = m+1; }                    // Solution for phase 1 found
                                                 //---------------------------
    min = 0.0;                             // Phase 1 and phase 2:
    for (j = 1; j <= n; j++) {             // Determine index k that minimizes
      s = Row(U,q) * Col(T,j);             // U[q][*] * T[*][j] (j == 1,..,n)
      if (s < min ) { min = s; k = j; }    //---------------------------------
    }

    if (-eps < min) {                      // Vector of reduced costs vanishes
      if (phase == 1)                      //---------------------------------
        { Err = NoInitialVertex; return; }
      SolutionFound = 1;                                         // phase == 2
      z = w[q];                                      // Optimal solution found
    }                                                //-----------------------
    else {                                        // Determine candidate l for
      for (i = 1; i <= q; i++)                    // exchange of indices
        y[i] = Row(U,i) * Col(T,k);               //--------------------------
      ex = 1;
      for (i = 1; i <= m; i++)
        if (y[i] >= eps) {
          s = w[i]/y[i];
          if (ex || s < min) { min = s; l = i; }
          ex = 0;
        }

      if (ex) { Err = NoOptimalSolution; return; }       // No candidate found
                                                         //-------------------
      v[l] = k;                               // Determine new basic index set
      s = 1.0 / y[l];                         // and compute new tableau
      for (j = 1; j <= m; j++) U[l][j] *= s;  //------------------------------
      i = (l == 1) ? 2 : 1;

      do {
        s = y[i];
        w[i] = w[i] - min * s;
        for (j = 1; j <= m; j++) U[i][j] -= U[l][j] * s;
        i = (i == l-1) ? i+2 : i+1;
      } while (i <= q);

      w[l] = min;
    } // Determine candidate l ...

  } while (!SolutionFound);

  x = 0.0;                                  // Return optimal value z, optimal
  for (i = 1; i <= m; i++)                  // basic index set v, and optimal
    x[v[i]] = w[i];                         // solution vector x
                                            //--------------------------------
} // RevSimplex




