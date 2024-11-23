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

/* CVS $Id: gop.cpp,v 1.16 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: gop (implementation)
// Purpose: Computing enclosures for all global minimizers and for the global
//    minimum value of a twice continuously differentiable multi-dimensional,
//    scalar valued function, assuming that the global minimum is a stationa-
//    ry point. If it is a boundary point of the search area with gradient of
//    the function being different from zero, the method fails in its form
//    presented here.
// Method: Bisection method combined with midpoint, monotonicity, concavity
//    test and extended interval Newton step.
// Global functions:
//    AllGOp()      : computes enclosures for all zeros
//    AllGOpErrMsg(): delivers an error message text
// Updates:
//    04.03.1996 : 'p = 0' eliminated in 'return' statements of 'NewtonStep'
//----------------------------------------------------------------------------
#include <string.h>          // String handling
#include <xi_ari.hpp>        // Extended interval arithmetic
#include <lst_ari.hpp>       // List arithmetic
#include <mv_util.hpp>       // Real matrix/vector utility functions
#include <mvi_util.hpp>      // Interval matrix/vector utility functions
#include <matinv.hpp>        // Inversion of real matrices
#include <hess_ari.hpp>      // Hessian differentiation arithmetic
#include <gop.hpp>

using namespace cxsc;
using namespace std;

static int MaxOptiNo;        // Maximum number of optimizers to be computed

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,  // No error occurred
  IllMaxOptiNo = 1,  // Illegal value for the maximum number of optimizers
  NotAllOptis  = 2,  // Not all optimizers computed
  NoStatOpti   = 3;  // No stationary point is global optimizer

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* AllGOpErrMsg ( int Err )
{
  static char Msg[160] = "";

  switch (Err) {
    case NoError:
      break;  // Blank message
    case IllMaxOptiNo:
      sprintf(Msg,"Error: Parameter for maximum number of optimizers "
                  "must lie in 1,...,%1d!",MaxCountGOp); break;
    case NotAllOptis:
      sprintf(Msg,"Warning: Not all optimizers found due to the user "
                           "limit of %1d optimizer(s).\n"
                  "         The enclosure of the global minimum value "
                           "could not be optimal!",MaxOptiNo); break;
    case NoStatOpti:
      strcpy(Msg,"Warning: No global optimizer (stationary point!) found!");
      break;
    default:
      strcpy(Msg,"Error: Code not defined!");
  }
  return(Msg);
}

static int MaxDiamComp ( ivector& iv )       // Determine the component 'mc'
{                                            // with maximum interval diameter
  int mc, k;                                 // of 'iv[mc]'.
  rvector d(Lb(iv), Ub(iv));                 //-------------------------------

  d = diam(iv);  mc = Lb(iv);
  for (k = Lb(iv)+1; k <= Ub(iv); k++)
    if (d[k] > d[mc]) mc = k;
  return mc;
}

static int MonotonicityTest ( ivector& gradY )              // returns true if
{                                                           // 'f' is monotone
  int i, n, Delete;                                         //----------------

  Delete = 0; i = 1;  n = Ub(gradY);
  while (i <= n && !Delete) {
    if (0.0 < Inf(gradY[i]) || 0.0 > Sup(gradY[i]) )   // 'f' is monotone
      Delete = 1;                                      //----------------
    i++;
  }
  return Delete;
}

static int ConcavityTest ( imatrix& HessY )               // returns true if
{                                                         // 'f' is not convex
  int i, n, Delete;                                       //------------------

  Delete = 0;  i = 1;  n = Ub(HessY,1);
  while (i <= n && !Delete) {
    if (Sup(HessY[i][i]) < 0.0)     // 'f' is not convex
      Delete = 1;                   //-------------------
    i++;
  }
  return Delete;
}

//----------------------------------------------------------------------------
// Purpose: Execution of one single interval Newton Gauss-Seidel step for the
//    interval vector 'Y' and the gradient of 'f'.
// Parameters:
//    In : 'f'    : must be declared for the type 'HessType' to enable the
//                  internal use of the differentiation arithmetic 'hess_ari'
//         'Y'    : specifies the starting interval
//         'HessY': Hessian matrix of 'f(Y)', already computed outside of
//                  'NewtonStep'
//    Out: 'V'    : stores the enclosures 'V[i]' for the splittings generated
//                  by the Newton step
//         'p'    : number of non-empty interval vectors 'V[i]'
// Description:
//    'NewtonStep' executes the extended interval Newton Gauss-Seidel step
//    for 'Y' with result interval vector(s) 'V[i]' which can be empty.
//    'p' gives the number of non-empty interval vectors stored as vectors
//    'V[1]',...,'V[p]'.
//----------------------------------------------------------------------------
static void NewtonStep ( HTscalar_FctPtr f,
                         ivector         Y,
                         imatrix&        HessY,
                         imatrix&        V,
                         int&            p )
{
  const int     n = Ub(Y);

  rvector       c(n);
  ivector       GradC(n), b(n), Y_minus_c(n);
  rmatrix       R(n,n);
  imatrix       A(n,n);
  int           i, i0, InvErr, j;
  interval      h, fC;
  ivector       z(2);
  idotprecision dd;

  c = mid(Y);                                  // Midpoint gradient evaluation
  fgEvalH(f,_ivector(c),fC,GradC);             //-----------------------------

  MatInv(mid(HessY),R,InvErr);                   // Invert the midpoint matrix
  if (InvErr) R = Id(R);                         //---------------------------

  A = R * HessY; b = R * GradC;          // Compute data for Gauss-Seidel step
  Y_minus_c = Y - c;                     //-----------------------------------

  p = 0;  i0 = 0;                   // Initializations, A[i0,i0] contains zero
                                    //----------------------------------------
  for (i = 1; i <= n; i++ ) {                // Interval-Gauss-Seidel-step for
                                             // non-zero A[i,i] elements
    if ( in(0.0,A[i][i]) ) {                 //-------------------------------
      i0 = i;       // Largest i with 0 in A[i,i]
      continue;     // Next i
    }

    dd = b[i];
    for (j = 1; j < i; j++) accumulate(dd,A[i][j],Y_minus_c[j]);
    for (j = i+1; j <= n; j++) accumulate(dd,A[i][j],Y_minus_c[j]);
    rnd(dd,h);
    h = c[i] - h / A[i][i];

    if ( Disjoint(Y[i],h) ) { return; }                // No solution possible
                                                       //---------------------
    Y[i] = Y[i] & h;  Y_minus_c[i] = Y[i] - c[i];
  } // for (i ...

  for (i = 1; i <= i0; i++ ) {               // Interval-Gauss-Seidel-step for
                                             // zero A[i,i] elements
                                             //-------------------------------
    if ( !in(0.0,A[i][i]) ) continue;                                // Next i
                                                                     //-------
    dd = b[i];
    for (j = 1;j < i; j++) accumulate(dd,A[i][j],Y_minus_c[j]);
    for (j = i+1;j <= n; j++) accumulate(dd,A[i][j],Y_minus_c[j]);
    rnd(dd,h);
    z = Y[i] & (c[i] - h % A[i][i]);             // Extended interval division
                                                 //---------------------------
    if (z[1] == EmptyIntval())                // z[1] == z[2] == EmptyIntval()
      { return; }                             // i.e. no solution possible
                                              //------------------------------
                                                            // Compute new 'Y'
    Y[i] = z[1];  Y_minus_c[i] = Y[i] - c[i];               //----------------
    if (z[2] != EmptyIntval())                     // Store further splittings
      { p++; V[p] = Y; V[p][i] = z[2]; }           //-------------------------
  } // for (i ...

  p++; V[p] = Y;
} // NewtonStep

//----------------------------------------------------------------------------
// Purpose: Execution of the global optimization method including a bisection
//    method, midpoint test, monotonicity test, concavity test, and extended
//    interval Newton steps.
// Parameters:
//    In : 'f'         : must be declared for type 'HTvector' to enable
//                       the internal use of the differentiation arithmetic
//                       'hess_ari'
//         'Start:     : specifies the starting interval
//         'Epsilon'   : specifies the desired relative accuracy
//    Out: 'ResultList': stores the candidates for enclosure of a global
//                       minimizer
//         'ListLength': length of the result list, that is the number of
//                       candidates in the result list
//         'Minimum'   : stores the enclosure of the global minimum value
// Description:
//    The function manages the list 'L' of pending subintervals that may
//    contain global minimizers. Subintervals are removed from the list and
//    placed in the accepted list 'ResultList' when they satisfy relative
//    error acceptance criteria. Subintervals are also removed from the list
//    by the midpoint, monotonicity, concavity tests, or by the interval
//    Newton steps. Subintervals are added to the pending list when an element
//    from the list is bisected or when the extended interval Newton step
//    yields two candidate intervals. 'ResultList' returns the list of
//    enclosures of the global minimizers. The maximum number of enclosures
//    to be computed is given by the value of the global variable 'MaxOptiNo'.
//    'Minimum' returns the enclosure of the global minimum value.
//----------------------------------------------------------------------------
static void GlobalOptimize ( HTscalar_FctPtr f,
                             ivector&        Start,
                             real            Epsilon,
                             PairPtr&        ResultList,
                             int&            ListLength,
                             interval&       Minimum )
{
  const int  n = Ub(Start);

  Pair       PairY;                    // For the pair (Y, inf(f(Y))
  ivector    Y(n);
  imatrix    U(2,n);                   // Subboxes of Y
  imatrix    V(n+1,n);                 // Subboxes of U
  interval   fY, fU, fV, fC, fCV;      // Function evaluations of f
  ivector    GradU(n), GradV(n);       // Gradient evaluations of f
  imatrix    HessU(n,n);               // Hessian evaluation of f
  real       fmax;                     // Upper bound for minimum
  rvector    c(n), cV(n);              // Midpoints of Y and V
  PairPtr    WorkList;                 // List of pairs
  int        i, j, p, k;               // Control variables
  int        Bisect;                   // Flag for iteration


ivector Vjhelp(n);  

  
  c = mid(Start);                           // Compute upper bound for minimum
  fEvalH(f,_ivector(c),fC);                 //--------------------------------
  fmax = Sup(fC);

  ResultList = EmptyList;
  ListLength = 0;

  if ( !UlpAcc(Start,1) ) {                                 // Start iteration
    Y = Start; WorkList = EmptyList;                        //----------------
    do {
      k = MaxDiamComp(Y);                           // d(Y[k] == max_i d(Y[i])
      U[1] = Y;  U[2] = Y;                          // Bisect 'Y' with respect
      SetSup(U[1][k],c[k]); SetInf(U[2][k],c[k]);   // to component 'k'.
      for (i = 1; i <= 2; i++) {                    //------------------------
        fgEvalH(f,U[i],fU,GradU);

        if (MonotonicityTest(GradU)) continue;   // Next i
                                                 // Try centered form to get
        fU = (fC + GradU*(U[i] - c)) & fU;       // better enclosure of 'f(U)'
                                                 //---------------------------
        if (fmax < Inf(fU)) continue;            // Next i

        fghEvalH(f,U[i],fU,GradU,HessU);           // Compute interval Hessian
                                                   //-------------------------
        if (ConcavityTest(HessU)) continue;        // Next i

        NewtonStep(f,U[i],HessU,V,p);         // Extended interval Newton step
                                              //------------------------------
        for (j = 1; j <= p; j++) {
          fgEvalH(f,V[j],fV,GradV);               // Compute interval gradient
                                                  //--------------------------
          if (MonotonicityTest(GradV)) continue;          // Next j

Vjhelp = V[j];
cV = mid(Vjhelp);     
//        cV = mid(V[j]);                                 // Try centered form
          fEvalH(f,_ivector(cV),fCV);                     // to get better en-
          fV = (fCV + GradV*(V[j] - cV)) & fV;            // closure of 'f(U)'
                                                          //------------------
          if (fmax >= Inf(fV)) {                                    // Store V
            PairY = _Pair(V[j],Inf(fV));                            //--------
            WorkList = WorkList + PairY;
          }
        } // for (j = ...
      } // for (i = ...

      Bisect = 0;                                           // Get next 'Y' of
      while (WorkList != EmptyList && !Bisect ) {           // the work list
        PairY = Head(WorkList);  DelHead(WorkList);         //----------------

        Y = GetInt(PairY);  c = mid(Y);
        fEvalH(f,_ivector(c),fC);                              // Compute f(c)
                                                               //-------------
        if (Sup(fC) < fmax) fmax = Sup(fC);
        MultiDelete(WorkList,fmax);

        Minimum = interval(GetFyi(PairY),fmax);

        if (RelDiam(Minimum) <= Epsilon || MaxRelDiam(Y) <= Epsilon) {
          ListLength++;                    // Check termination criteria. Stop
                                           // due to user specified maximum
          if (ListLength > MaxOptiNo) {    // for the number of optimizers.
            FreeAll(WorkList);             //---------------------------------
            Minimum = interval(GetFyi(Head(ResultList)),fmax);
            return;
          }

          ResultList = ResultList + PairY;                    // Store 'PairY'
        }                                                     //--------------
        else
          Bisect = 1;
      } // while

    } while(Bisect);
  } // if ( !UlpAcc ...
  else {
    fEvalH(f,Start,fY);                             // Store starting interval
    ResultList = EmptyList + _Pair(Start,Inf(fY));  // and interval evaluation
  }                                                 //------------------------
                                                // Compute good enclosure of
  if (ResultList != EmptyList)                  // of the global minimum value
    Minimum = interval(GetFyi(Head(ResultList)),fmax);//---------------------
} // GlobalOptimize

//----------------------------------------------------------------------------
// 'MaxNorm' delivers an upper bound of the maximum norm of a symmetric in-
// terval matrix 'H', i.e. the row sum norm or infinity norm.
//----------------------------------------------------------------------------
static real MaxNorm ( imatrix& H )
{
  real  Nm, MaxNm;
  int   i, j;

  MaxNm = 0.0;
  for (i = Lb(H,1); i <= Ub(H,1); i++) {
    Nm = 0.0;
    for (j = Lb(H,1); j <= Ub(H,1); j++) Nm = addu(Nm,Sup(abs(H[i][j])));
    if (Nm > MaxNm) MaxNm = Nm;
  }
  return MaxNm;
}

//----------------------------------------------------------------------------
// 'PosDef' delivers true if it is guaranteed that all real symmetric matri-
// ces enclosed in 'H' are positive definite.
//----------------------------------------------------------------------------
static int PosDef ( imatrix& H )
{
  const int  kmax = 10;      // Maximum number of iterations

  int      pd;
  real     kappa, eps;
  imatrix  S(Lb(H,1),Ub(H,1),Lb(H,2),Ub(H,2));
  ivector  Z(Lb(H,1),Ub(H,1)),U(Lb(H,1),Ub(H,1));
  int      k;

  kappa = MaxNorm(H);  S = Id(H) - H/kappa;
  for (k = Lb(Z); k <= Ub(Z); k++) Z[k] = interval(-1.0,1.0);
  k = 0;  eps = 0.25;
  do {
    U = Blow(Z,eps);  Z = S * U;  pd = in(Z,U);
    k++;  eps = 8.0*eps;
  } while( !(pd || k == kmax) );
  return pd;
}

//----------------------------------------------------------------------------
// Purpose: Execution of a verification step including the use of an epsilon
//    inflation.
// Parameters:
//    In    : 'f'      : function of 'HessType'
//    Out   : 'yUnique': returns 'true' if the verification is successful
//    In/Out: 'y'      : interval enclosure to be verified
// Description: This function checks the uniqueness of the local minimizer
//    enclosed in the interval variable 'y' by a verification step including
//    the use of an epsilon inflation of the iterates.
//----------------------------------------------------------------------------
static void VerificationStep ( HTscalar_FctPtr f, const imatrix_subv& y, int& yUnique )
{
  const int
    kmax = 10,        // Maximum number of iterations
    n = Ub(y);        // Length of vectors in use

  interval  fY;
  ivector   yIn(n), yOld(n), GradY(n);
  imatrix   HessY(n,n);
  imatrix   yp(n+1,n);
  int       k, p;
  real      eps;

  yUnique = (Inf(y) == Sup(y));                // y is a point interval vector
  if (yUnique) return;                         //-----------------------------

  yIn = y;  yOld = y;  k = 0;  eps = 0.25;                  // Initializations
                                                            //----------------
  while (!yUnique && k < kmax) {       // Do 'kmax' loops to achieve inclusion
                                       //-------------------------------------
    yOld = Blow(y,eps);                            // Epsilon inflation of 'y'
                                                   //-------------------------
    // Perform interval Newton step
    //------------------------------
    k++;
    fghEvalH(f,yOld,fY,GradY,HessY);
    NewtonStep(f,yOld,HessY,yp,p);

    if (p != 1) break;                             // No verification possible
                                                   //-------------------------
    if (yp[1] == yOld)
      eps = eps * 8.0;                                       // Increase 'eps'
    else {                                                   //---------------
	for (int ii=Lb(y); ii<=Ub(y); ii++)
	    y[ii] = yp[1][ii];
      yUnique = in(y,yOld);                 // Inner inclusion ===> uniqueness
    }                                       // of a stationary point
  }   // while                              //--------------------------------

  if (yUnique) yUnique = PosDef(HessY);        // Positive definite ===> local
                                               // minimizer is verified
  if (!yUnique)
      for (int ii=Lb(y); ii<=Ub(y); ii++)
          y[ii] = yIn[ii];      
} // VerificationStep

//----------------------------------------------------------------------------
// Purpose: Computation of enclosures for all global minimizers and for the
//    global minimum value of a twice continuously differentiable multi-di-
//    mensional, scalar-valued function.
// Parameters:
//    In : 'f'               : objective function, must be declared for the
//                             'HessType' to enable the internal use of
//                             the differentiation arithmetic 'hess_ari'
//         'Start',          : specifies the starting interval vector
//         'Epsilon'         : specifies the desired relative accuracy
//                             (interval diameter) of the result boxes
//         'MaxNumberOfOptis': maximum number of optimizers to be computed
//                             (default value: 'MaxCountGOp' = compute all)
//    Out: 'OptiVector'      : stores (row by row) and returns the boxes
//                             (enclosures) for the global optimizers of
//                             'f'. 'OptiVector' is a vector of interval
//                             vectors (ie. an interval matrix).
//         'InfoVector'      : stores the corresponding information on the
//                             uniqueness of the local optimizers in these
//                             enclosures
//         'NumberOfOptis'   : number of enclosures computed
//         'Minimum'     '   : enclosure for the minimum value
//         'Err'             : error code
// Description:
//    The enclosures for the global minimizers of 'f' are computed by calling
//    function 'GlobalOptimize'. Then a verification step is applied.
//    The enclosures (boxes) for the global minimizers of 'f' are stored in
//    the interval matrix 'OptiVector' row by row, the corresponding
//    information on the uniqueness of the local minimizers in these
//    enclosures is stored in the integer vector 'InfoVector'. The number of
//    enclosures computed is returned in the integer variable 'NumberOfOptis'.
//    The enclosure for the global minimum value is returned in the variable
//    'Minimum'. If an error occurs, the value of 'Err' is different from 0.
//    The result arrays 'OptiVector' and 'InfoVector' are automatically
//    resized to standard index range with lower index bounds 1. The maximum
//    number of optimizers which can be computed is given by
//    'MaxNumberOfOptis'. It is bounded by 'MaxCountGOp'. If your system can
//    handle more than 'MaxCountGOp' optimizers, you can increase its value.
//----------------------------------------------------------------------------
void AllGOp ( HTscalar_FctPtr f,
              ivector         Start,
              real            Epsilon,
              imatrix&        OptiVector,
              intvector&      InfoVector,
              int&            NumberOfOptis,
              interval&       Minimum,
              int&            Err,
              int             MaxNumberOfOptis )
{
  // Check range for the maximum number of optimizers
  //-------------------------------------------------
  if (MaxNumberOfOptis < 1 || MaxCountGOp < MaxNumberOfOptis)
    { Err = IllMaxOptiNo;  NumberOfOptis = 0;  return; }

  int     i;
  PairPtr ResultList;
  real    MinEpsilon;

  SetLb(Start,1);                     // Resize start vector to standard index
                                      // range with lower index bound 1
                                      //--------------------------------------
  Err = NoError;
  MaxOptiNo = MaxNumberOfOptis;           // Initialize global static variable
  MinEpsilon = succ(1.0) - 1.0;           // Relative machine accuracy (1 ulp)
                                          //----------------------------------
  if (Epsilon < MinEpsilon) Epsilon = MinEpsilon;          // Set 'Epsilon' to
                                                           // 1 ulp accuracy
  GlobalOptimize(f,Start,Epsilon,ResultList,NumberOfOptis,Minimum);

  if (NumberOfOptis > MaxNumberOfOptis)
    { Err = NotAllOptis;  NumberOfOptis = MaxNumberOfOptis; }

  if (NumberOfOptis == 0)
    Err = NoStatOpti;
  else                                                    // Store results
    ListToMatrix(ResultList,OptiVector,InfoVector);       // in the matrix
                                                          //--------------
  // Verification step for the enclosure intervals
  //----------------------------------------------
  for (i = 1; i <= NumberOfOptis; i++)
    VerificationStep(f,OptiVector[i],InfoVector[i]);

  FreeAll(ResultList);
}




