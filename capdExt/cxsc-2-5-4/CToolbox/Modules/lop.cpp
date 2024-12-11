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

/* CVS $Id: lop.cpp,v 1.16 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: lop (implementation)
// Purpose: Determine enclosures for the optimal value 'z_opt', for the set
//    of optimal basic index sets 'V_opt', and for the set of solution
//    vectors 'X_opt' for a linear programming problem P = (A,b,c) given in
//    the standard form:
//               ( z = c^t * x = max! )
//         (LP)  (     A * x = b      )
//               (      x >= 0        )
//    with an initial optimal basic index set.
// Method: Starting from an initial basic index set, all neighboring index
//    sets are examined to determine all optimal solutions. To validate the
//    optimality of a solution, all linear systems of equations are solved
//    using a verifying linear system solver.
// Global functions:
//    LinOpt()      : determines the enclosures of z_opt, V_opt, and X_opt
//                    for LP P = (A,b,c)
//    LinOptErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#include <string.h>       // String handling
#include <set_ari.hpp>    // Integer set arithmetic
#include <linsys.hpp>     // Verifying linear system solver
#include <lop_ari.hpp>    // Linear linked list for index sets
#include <lop.hpp>

using namespace cxsc;
using namespace std;

const int
  kmax = 100;             // Maximum number of iterations

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError                 = 0,  // No error occurred
  WrongDimension          = 1,  // Wrong dimension of problem
  SolutionSetIsEmpty      = 2,  // Set of feasible solutions is empty
  FunctionUnbounded       = 3,  // Objective function is unbounded
  SubmatrixSingular       = 4,  // Submatrix A_B is singular
  StartIndexSetNotOptimal = 5,  // Initial basic index set not optimal
  IterationError          = 6;  // Maximum number of iteration exceeded

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* LinOptErrMsg ( int Err )
{
  static char Msg[100] = "";

  if (Err != NoError) {
    char Hlp[80];

    switch (Err) {
      case WrongDimension:
        strcpy(Hlp,"Wrong dimension of problem (e.g. m >= n)");       break;
      case SolutionSetIsEmpty:
        strcpy(Hlp,"Set of feasible solutions is probably empty");    break;
      case FunctionUnbounded:
        strcpy(Hlp,"Objective function is probably unbounded");       break;
      case SubmatrixSingular:
        strcpy(Hlp,"Submatrix A_B in [A_B] is probably singular");    break;
      case StartIndexSetNotOptimal:
        strcpy(Hlp,"Initial basic index set not optimal (i.e. ");
        strcat(Hlp,"B_start not in V_opt)");                          break;
      case IterationError:
        sprintf(Hlp,"Maximum number of iterations (=%1d) exceeded",
                    kmax);                                            break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // LinOptErrMsg

static ivector Solution ( imatrix  TT,            // Determine solution vector
                          IndexSet Base,          // corresponding to Base
                          IndexSet NonBase )      //--------------------------
{
  int      Index, i, sB = Size(Base);
  ivector  res(sB+Size(NonBase));
  interval Rplus = interval(0.0,MaxReal);

  res = 0.0;
  for (i = 1; i <= sB; i++) {
    Index = Base[i];
    res[Index] = TT[i+1][1] & Rplus;
  }
  return res;
}

//----------------------------------------------------------------------------
// Compute the interval tableau
//                            ([z]   | [d]^t    )
//                      [T] = (------+----------)
//                            ([x]_B | [H]      ).
//----------------------------------------------------------------------------
static void ComputeTableau ( rmatrix   A,
                             rvector   b,
                             rvector   c,
                             IndexSet  Base,
                             IndexSet  NonBase,
                             imatrix&  TT,
                             int&      Err )
{
  const int
    m = Ub(A,1), n = Ub(A,2), sB = Size(Base), sN = Size(NonBase);

  rmatrix   A_B(m,sB), A_N(m,sN);
  rvector   c_B(sB), c_N(sN);
  ivector   xx_B(sB), yy(sB), dd(sN);
  imatrix   HH(sB,sN);
  interval  zz;
  int       i, j, Index, local_Error;

  // Determine submatrices and subvectors according to index sets Base and
  // NonBase
  A_B = extract(A, Base); A_N = extract(A, NonBase);  // Determine submatrices
  c_B = extract(c, Base); c_N = extract(c, NonBase);  // and subvectors accor-
                                                      // ding to B and N
  // Solve linear systems of equations                //----------------------
  LinSolve(A_B, b, xx_B, local_Error);

  if (local_Error == NoError)
    LinSolve(transp(A_B), c_B, yy, local_Error);

  if (local_Error == NoError)
    for (i = 1; i <= n-m; i++) {
      Index = NonBase[i];
      if (local_Error == NoError) {
        // Temporary variables colHH and colA needed since automatic
        // resizing of Col(HH,i) and Col(A,Index) is not supported by C++
        //---------------------------------------------------------------
        ivector colHH = Col(HH,i);
        rvector colA = Col(A,Index);
        LinSolve(A_B, colA, colHH, local_Error);
        Col(HH,i) = colHH;
      }
    }

  if (local_Error == NoError) {
    dd = transp(A_N) * yy - c_N;             // Compute components of interval
    zz = (c_B * xx_B) & (b * yy);            // tableau
    TT[1][1] = zz;                           //-------------------------------
    for (j = 2; j <= n-m+1; j++) TT[1][j] = dd[j-1];
    for (i = 2; i <= m+1; i++) TT[i][1] = xx_B[i-1];
    for (i = 2; i <= m+1; i++)
      for (j = 2; j <= n-m+1; j++)
        TT[i][j] = HH[i-1][j-1];
  }
  else
    Err = SubmatrixSingular;
} // ComputeTableau

//----------------------------------------------------------------------------
// Determine whether the interval tableau [T] represents a basis stable
// solution, i.e. (Inf([x]_B) > 0) and (Inf([d]) > 0).
//----------------------------------------------------------------------------
static int BasisStable ( imatrix& TT )
{
  int IsStable = 1, i = 1, j = 1;

  while (i < Ub(TT,1) && IsStable) {
    i++;                                    // [x]_B = [T_i][1] (i = 2 .. m+1)
    if (Inf(TT[i][1]) <= 0.0) IsStable = 0;
  }
  while (j < Ub(TT,2) && IsStable) {
    j++;                                    // [d] = [T_1][j] (j = 2 .. n-m+1)
    if (Inf(TT[1][j]) <= 0.0) IsStable = 0;
  }
  return IsStable;
}

//----------------------------------------------------------------------------
// Determine whether the interval tableau [T] represents a possibly optimal
// solution, i.e (Sup([x]_B) >= 0) and (Sup([d]) >= 0).
//----------------------------------------------------------------------------
static int PossiblyOptimalSolution ( imatrix& TT )
{
  int IsOptimal = 1, i = 1, j = 1;

  while (i < Ub(TT,1) && IsOptimal) {
    i++;                                    // [x]_B = [T_i][1] (i = 2 .. m+1)
    if (Sup(TT[i][1]) < 0.0) IsOptimal = 0;
  }
  while (j < Ub(TT,2) && IsOptimal) {
    j++;                                    // [d] = [T_1][j] (j = 2 .. n-m+1)
    if (Sup(TT[1][j]) < 0.0) IsOptimal = 0;
  }
 return IsOptimal;
}

//----------------------------------------------------------------------------
// Check whether the set of feasible solutions is empty, i.e. not for all(beta
// in B) with (0 in [x]_beta) there exists a (nu in N) with ([H]_beta,nu < 0)
//----------------------------------------------------------------------------
static int EmptySolutionSet ( imatrix& TT )
{
  int IsEmpty = 0, i = 1, j;

  while (!IsEmpty && i < Ub(TT,1)) {
    i++;
    if ( in(0.0,TT[i][1]) ) {
      IsEmpty = 1;
      j = 1;
      while (IsEmpty && j < Ub(TT,2)) {
        j++;
        if (Sup(TT[i][j]) < 0.0) IsEmpty = 0;
      }
    }
  }
  return IsEmpty;
}

//----------------------------------------------------------------------------
// Check whether the objective function is unbounded, i.e. not for all (nu
// in N) with (0 in [d]_nu) there exists a (beta in B) with ([H]_beta,nu > 0)
//----------------------------------------------------------------------------
static int Unbounded ( imatrix& TT )
{
  int IsUnbounded = 0, i, j = 1;

  while (!IsUnbounded && j < Ub(TT,2)) {
    j++;
    if ( in(0.0,TT[1][j]) ) {
      IsUnbounded = 1;
      i = 1;
      while (IsUnbounded && i < Ub(TT,1)) {
        i++;
        if (Inf(TT[i][j]) > 0.0) IsUnbounded = 0;
      }
    }
  }
  return IsUnbounded;
}

//----------------------------------------------------------------------------
// Determine list of neighboring basic index sets L for index set Base that
// are good candidates for being optimal basic index sets.
//----------------------------------------------------------------------------
static BaseList NeighboringList ( imatrix&   TT,
                                  IndexSet&  Base,
                                  IndexSet&  NonBase )
{
  BaseList  L = NULL;
  IndexSet  NewBase;
  int       i, j, m = Size(Base),n = m+Size(NonBase);
  int       beta, nu;
  real      colmin, rowmax;
  interval  xx, dd, HH;

  colmin = MaxReal;                              // Search for candidates (nu
  // for all (nu in N)                           // in N) and (beta in B) by
  for (j = 1; j <= n-m; j++) {                   // determining of the minimum
    dd = TT[1][j+1];                             // of the quotients
    if ( in(0.0,dd) ) {                          // [x]_beta / [H]_beta,nu
      // for all (beta in B)                     //---------------------------
      for (i = 1; i <= m; i++) {                     //  Determine minimum of
        xx = TT[i+1][1];                             // [x]_beta / [H]_beta,nu
        HH = TT[i+1][j+1];                           // for column nu of [T]
                                                     //-----------------------
        if (Inf(HH) > 0.0 && divd(Sup(xx),Inf(HH)) < colmin)
          colmin = divd(Sup(xx),Inf(HH));
      } // (beta in B)

      // for all (beta in B)
      for (i = 1; i <= m; i++) {                  // Determine candidate (beta
        xx = TT[i+1][1];                          // in B) for exchange of
        HH = TT[i+1][j+1];                        // indices
                                                  //--------------------------
        if (Sup(HH) > 0.0 && divd(Inf(xx),Sup(HH)) <= colmin) {
          beta = Base[i];                       // Determine new index set and
          nu = NonBase[j];                      // add it to L
          NewBase = Base - beta + nu;           //----------------------------
          insert(L, NewBase);
        }
      } // (beta in B)
    }
  } // (nu in N)

  rowmax = -MaxReal;                             // Search for candidates beta
  // for all (beta in B)                         // in B) and (nu in N) by
  for (i = 1; i <= m; i++) {                     // determining the maximum
    xx = TT[i+1][1];                             // of the quotients
    if ( in(0.0,xx) ) {                          // [d]_nu / [H]_beta,nu
      // for all (nu in N)                       //---------------------------
      for (j = 1;j <= n-m; j++) {                      // Determine maximum of
        dd = TT[1][j+1];                               // [d]_nu/[H]_beta,nu
        HH = TT[i+1][j+1];                             // for row beta of [T]
                                                       //---------------------
        if (Sup(HH) < 0.0 && divu(Inf(dd),Sup(HH)) > rowmax)
          rowmax = divu(Sup(dd),Sup(HH));
      } // (nu in N)

      // for all (nu in N)
      for (j = 1;j <= n-m; j++) {                    // Determine candidate
        dd = TT[1][j+1];                             // (nu in N) for exchange
        HH = TT[i+1][j+1];                           // of indices
                                                     //-----------------------
        if (Inf(HH) < 0.0 && divu(Inf(dd),Inf(HH)) >= rowmax) {
          beta = Base[i];                       // Determine new index set and
          nu = NonBase[j];                      // add it to L
          NewBase = Base - beta + nu;           //----------------------------
          insert(L, NewBase);
        }
      } // (nu in N)
    }
  } // (beta in B)

  return L;
} // NeighboringList

//----------------------------------------------------------------------------
// Purpose: Determination of enclosures for the solutions of a linear opti-
//    mization problem:      ( z = c^t * x = max! )
//                     (LP)  (     A * x = b      )
//                           (      x >= 0        )
//    with an initial optimal basic index set.
// Parameters:
//    In : 'A'             : matrix of constraints
//         'b'             : right-hand side of constraints
//         'c'             : represents the objective function by its
//                           coefficients
//         'B_Start_Vector': initial basic index set
//    Out: 'z'             : enclosure of the optimal value 'z_opt'
//         'V'             : superset of all optimal basic index sets 'V_opt'
//         'X'             : superset of all optimal solution vectors 'X_opt'
//         'No'            : number of solutions computed
//         'Err'           : error code
// Description: Starting from an initial basic index set 'B_start', a true
//    superset of all neighboring index sets 'L' that are candidates for being
//    optimal index sets are determined and stored in a list of candidates
//    'CList'. As long as there are candidates left, they are examined using
//    a verifying linear system solver. The linear system solver returns
//    enclosures for the solutions of the square linear systems A_B [x] = b,
//    where A_B denotes a submatrix of the constraint matrix A according to a
//    basic index set B. The result arrays 'V' and 'X' are automatically
//    resized to standard index range with lower index bounds 1.
//----------------------------------------------------------------------------
void LinOpt ( rmatrix     A,
              rvector     b,
              rvector     c,
              intvector   B_Start_Vector,
              interval&   z,
              intmatrix&  V,
              imatrix&    X,
              int&        No,
              int&        Err )
{
  const int
    MinSize = 10;   // Minimum size (number of rows) of the result matrices

  // Check for correct dimensions of the input data
  //-----------------------------------------------
  int  m = ColLen(A), n = RowLen(A);

  if ( m >= n || VecLen(b) != m
              || VecLen(c) != n
              || VecLen(B_Start_Vector) != m )
    { Err = WrongDimension;  return; }

  SetLb(A,1,1); SetLb(A,2,1);           // Resize input data to standard index
  SetLb(b,1); SetLb(c,1);               // range with lower index bound 1
  SetLb(B_Start_Vector,1);              //------------------------------------

  Resize(V,MinSize,m);                        // Resize the result matrices to
  Resize(X,MinSize,n);                        // 'MinSize' rows
                                              //------------------------------

  IndexSet  B_Start(n), Base(n,'\0'), NonBase(n,'\0');
  BaseList  L, CList, E;
  imatrix   TT(m+1,n-m+1);
  int       k, i, TempSize = MinSize;


  for (i = 1; i <= m; i++) B_Start = B_Start + B_Start_Vector[i];

  V = 0; X = 0.0; k = 0; No = 0; Err = NoError;
  E = NULL; CList = NULL; insert(CList, B_Start);

  do { // while (CList != EmptyList ... )
    Base = select(CList);                                // Select Base out of
    NonBase = Complement(Base);                          // CList and update
    del(CList, Base); insert(E, Base);                   // examined list E
                                                         //-------------------
    ComputeTableau(A, b, c, Base, NonBase, TT, Err);

    if (Err == NoError) {
      if (BasisStable(TT)) {
        No = 1;                               // Store unique optimal solution
        z = TT[1][1];                         //------------------------------

        Row(X,No) = Solution(TT, Base, NonBase);
        SetToVector(Base,Row(V,No));
      }
      else { // not basis stable
        if (EmptySolutionSet(TT))
          Err = SolutionSetIsEmpty;
        else if (Unbounded(TT))
          Err = FunctionUnbounded;
        else if (PossiblyOptimalSolution(TT) && Err == NoError) {
          No++;                                      // Store optimal solution
          if (No == 1)                               //-----------------------
            z = TT[1][1];
          else
            z = TT[1][1] | z;

          if (No > TempSize)                 // Double size of result matrices
            { DoubleSize(V); DoubleSize(X); TempSize *= 2; }

          Row(X,No) = Solution(TT, Base, NonBase);
          SetToVector(Base,Row(V,No));

          L = NeighboringList(TT,Base,NonBase);       // Determine list of
          remove(L, E); remove(L, CList);             // neighboring basic
                                                      // index sets
                                                      //------------------
          append(CList, L);              // Compute new list of candidates
        }                                //-------------------------------
      } // not basis stable
    }
    
    k++;

  } while (CList != EmptyList && Err == NoError && k < kmax);

  if (Err == NoError) {                                // Determine error code
    if (CList != EmptyList && k == kmax)               //---------------------
      Err = IterationError;
    else if (No == 0)
      Err = StartIndexSetNotOptimal;
  }
  
  // Resize the result matrices to appropriate number of rows
  // if no error has occurred. Otherwise, resize them to minimum
  // dimensions with undefined value.
  //------------------------------------------------------------
  if (Err == NoError) {
    Resize(V,No,m);
    Resize(X,No,n);
  }
  else {
    Resize(V,1,1);     // Values undefined
    Resize(X,1,1);     //-----------------
  }

  FreeAll(CList); FreeAll(E);
} // LinOpt




