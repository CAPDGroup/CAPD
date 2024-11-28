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

/* CVS $Id: gop1.cpp,v 1.16 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: gop1 (implementation)
// Purpose: Computing enclosures for all global minimizers and for the global
//    minimum value of a twice continuously differentiable one-dimensional,
//    scalar valued function.
// Method: Bisection method combined with midpoint, monotonicity, concavity
//    test and extended interval Newton step.
// Global functions:
//    AllGOp1()      : computes enclosures for all global optimizers
//    AllGOp1ErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#include <string.h>         // String handling
#include <i_util.hpp>       // Interval utility functions
#include <xi_ari.hpp>       // Extended interval arithmetic
#include <lst1_ari.hpp>     // List arithmetic
#include <gop1.hpp>

using namespace cxsc;
using namespace std;

static int MaxOptiNo;       // Maximum number of optimizers to be computed

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,  // No error occurred
  IllMaxOptiNo = 1,  // Illegal value for the maximum number of optimizers
  NotAllOptis  = 2;  // Not all optimizers computed

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* AllGOp1ErrMsg ( int Err )
{
  static char Msg[160] = "";

  switch (Err) {
    case NoError:
      break;  // Blank message
    case IllMaxOptiNo:
      sprintf(Msg,"Error: Parameter for maximum number of optimizers "
                  "must lie in 1,...,%1d!",MaxCountGOp1); break;
    case NotAllOptis:
      sprintf(Msg,"Warning: Not all optimizers found due to the user "
                           "limit of %1d optimizer(s).\n"
                  "         The enclosure of the global minimum value "
                           "could not be optimal!",MaxOptiNo); break;
    default:
      strcpy(Msg,"Error: Code not defined!");
  }
  return(Msg);
}

//----------------------------------------------------------------------------
// Purpose: Execution of one extended interval Newton step for the derivative
//    of function 'f'.
// Parameters:
//    In : 'f'   : must be declared for the 'DerivType' to enable the
//                 internal use of the differentiation arithmetic 'ddf_ari'
//         'Y'   : specifies the starting interval
//         'ddfY': f''(Y), already computed outside of 'NewtonStep'
//    Out: 'V'   : Pair of intervals V[1] and V[2]
//         'p'   : Number of valid intervals in V (0, 1, or 2)
// Description:
//    One extended interval Newton step for 'Y' is executed resulting in 'p'
//    interval(s) 'V[1]' and 'V[2]' which can be empty.
//----------------------------------------------------------------------------
static void NewtonStep ( ddf_FctPtr f,
                         interval   Y,
                         interval   ddfY,
                         ivector&   V,
                         int&       p )
{
  real     c;
  interval fC, dfC;

  c = mid(Y); dfEval(f,_interval(c),fC,dfC);     // Midpoint evaluation of 'f'
  V = Y & (c - dfC % ddfY)     ;                 // Execution of Newton step
  if (V[1] != EmptyIntval()) p = 1; else p = 0;  // Fix number of non-empty
  if (V[2] != EmptyIntval()) p++;                // intervals
}                                                //---------------------------

//----------------------------------------------------------------------------
// Purpose: Execution of the global optimization method including a bisection
//    method, midpoint test, monotonicity test, concavity test, and extended
//    interval Newton steps.
// Parameters:
//    In : 'f'         : must be declared for the 'DerivType' to enable the
//                       internal use of the differentiation arithmetic
//                       'ddf_ari'
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
static void GlobalOptimize ( ddf_FctPtr f,
                             interval   Start,
                             real       Epsilon,
                             Pair1Ptr&  ResultList,
                             int&       ListLength,
                             interval&  Minimum )
{
  Pair1     PairY;                        // For the pair (Y, inf(f(Y))
  interval  Y, fY, BdP;        // Initial box, f(Y), f'(Y), f''(Y)
  ivector   U(2), V(2);                   // Subboxes of Y and V
  interval  fU, fV, fC, fCV, fBdP;        // Function evaluations of f
  interval  dfU, dfV, ddfU;               // Derivative evaluations of f
  real      fmax;                         // Upper bound for minimum value
  real      c, cV;                        // Midpoint of interval Y and V
  Pair1Ptr  WorkList;                     // List of pairs
  int       i, j, p;                      // Control variables
  int       Bisect;                       // Flag for iteration and algorithm

  c = mid(Start); fEval(f,_interval(c),fC); // Compute upper bound for minimum
  fmax = Sup(fC);                           //--------------------------------

  if ( !UlpAcc(Start,1) ) {
    Y = Start;                                                 // Start method
    WorkList = EmptyList;                                      //-------------
    ResultList = EmptyList;
    ListLength = 0;

    BdP = Inf(Y);                          // Treat boundary points separately
    for (i = 1; i <= 2; i++) {             //---------------------------------
      fEval(f,BdP,fBdP);
      if (Sup(fBdP) < fmax) fmax = Sup(fBdP);
      if (fmax >= Inf(fBdP))
        WorkList = WorkList + _Pair1(BdP,Inf(fBdP));
      BdP = Sup(Y);
    }

    do {                                                    // Start iteration
      U[1] = interval(Inf(Y),c);                            // Bisect 'Y'
      U[2] = interval(c,Sup(Y));                            //----------------
      for (i = 1; i <= 2; i++) {
        dfEval(f,U[i],fU,dfU);                          // Compute dfu = f'(U)
                                                        //--------------------
        if ( !in(0.0,dfU) ) continue;     // Monot. test: if not 0 in dfU stop
                                          //----------------------------------
        fU = (fC + dfU * (U[i] - c)) & fU;             // Centered form 'f(U)'
                                                       //---------------------
        if (fmax < Inf(fU)) continue;

        ddfEval(f,U[i],fU,dfU,ddfU);                    // Compute dfu = f'(U)
                                                        //--------------------
        if (Sup(ddfU) < 0.0) continue;     // Concavity test: if 0 > ddfU stop
                                           //---------------------------------
        NewtonStep(f,U[i],ddfU,V,p);          // Extended interval Newton step
                                              //------------------------------
        for (j = 1; j <= p; j++) {
          dfEval(f,V[j],fV,dfV);
          if ( !in(0.0,dfV) ) continue;   // Monot. test: if not 0 in dfV stop
                                          //----------------------------------
          cV = mid(V[j]);                                 // Try centered form
          fEval(f,_interval(cV),fCV);                     // to get better en-
          fV = (fCV + dfV * (V[j] - cV)) & fV;            // closure of 'f(U)'
                                                          //------------------
          if (fmax >= Inf(fV))                                      // Store V
            WorkList = WorkList + _Pair1(V[j],Inf(fV));             //--------
        } // for j ...
      } // for i ...

      Bisect = 0;                                           // Get next 'Y' of
      while ( (WorkList != EmptyList) && !Bisect ) {        // the work list
        PairY = Head(WorkList); DelHead(WorkList);          //----------------

        Y = GetInt(PairY); c = mid(Y); fEval(f,_interval(c),fC);
        if (Sup(fC) < fmax) fmax = Sup(fC);
        MultiDelete(WorkList,fmax);

        Minimum = interval(GetFyi(PairY),fmax);
        if ( (RelDiam(Minimum) <= Epsilon) || (RelDiam(Y) <= Epsilon) ) {

          ListLength++;                // Checking termination criterion. Stop
                                       // due to user specified maximum for
          if (ListLength > MaxOptiNo) {// the number of optimizers.
            FreeAll(WorkList);         //-------------------------------------
            Minimum = interval(GetFyi(Head(ResultList)),fmax);
            return;
          }

          ResultList = ResultList + PairY;
        }
        else
          Bisect = 1;
      } // while

    } while (Bisect);
  } // if ( !UlpAcc(Start,1) )
  else {                                            // Store starting interval
    fEval(f,Start,fY);                              //------------------------
    ResultList = EmptyList + _Pair1(Start,Inf(fY));
    ListLength = 1;
  }
                                                     // Compute good enclosure
  Minimum = interval(GetFyi(Head(ResultList)),fmax); // of the global minimum
} // GlobalOptimize                                  //-----------------------

//----------------------------------------------------------------------------
// Purpose: Execution of a verification step including the use of an epsilon
//    inflation.
// Parameters:
//    In    : 'f'      : function of 'DerivType'
//    Out   : 'yUnique': returns 'true' if the verification is successful
//    In/Out: 'y'      : interval enclosure to be verified
// Description: This function checks the uniqueness of the local minimizer
//    enclosed in the interval variable 'y' by a verification step including
//    the use of an epsilon inflation of the iterates.
//----------------------------------------------------------------------------
static void VerificationStep ( ddf_FctPtr f, interval& y, int& yUnique )
{
  const
    int kmax = 10;  // Maximum number of iterations

  interval  yIn, yOld, fm, dfm, fy, dfy, ddfy;
  real      m, eps;
  int       k;

  yUnique = (Inf(y) == Sup(y));                       // y is a point interval
  if (yUnique) return;                                //----------------------

  yIn = y; k = 0; eps = 0.25;                               // Initializations
                                                            //----------------
  while ( !yUnique && (k < kmax) ) {     // Do kmax loops to achieve inclusion
                                         //-----------------------------------
    yOld = Blow(y,eps);                            // Epsilon inflation of 'y'
    ddfEval(f,yOld,fy,dfy,ddfy);                   //-------------------------

    if (Inf(ddfy) <= 0.0) break;      // No verification of a minimum possible
                                      //--------------------------------------
                                               // Perform interval Newton step
    k++; m = mid(yOld);                        //-----------------------------
    dfEval(f,_interval(m),fm,dfm);
    y = (m - dfm / ddfy);

    if (Disjoint(y,yOld)) break;                   // No verification possible
                                                   //-------------------------
    yUnique = in(y,yOld);            // Inner inclusion ===> unique zero of f'
    y = y & yOld;                    // Intersection with old value
    if (y == yOld)                   //---------------------------------------
      eps = eps * 8.0;                          // Increase the value of 'eps'
                                                //----------------------------
  } // while

  if (!yUnique) y = yIn;
} // VerificationStep

//----------------------------------------------------------------------------
// Purpose: Computation of enclosures for all global minimizers and for the
//    global minimum value of a twice continuously differentiable one-dimen-
//    sional, scalar-valued function.
// Parameters:
//    In : 'f'               : objective function, must be declared for the
//                             'DerivType' to enable the internal use of
//                             the differentiation arithmetic 'ddf_ari'
//         'Start',          : specifies the starting interval
//         'Epsilon'         : specifies the desired relative accuracy
//                             (interval diameter) of the result intervals
//         'MaxNumberOfOptis': maximum number of optimizers to be computed
//                             (default value: 'MaxCountGOp1' = compute all)
//    Out: 'OptiVector'      : stores and returns the enclosures for the
//                             global optimizers of 'f'
//         'InfoVector'      : stores the corresponding information on the
//                             uniqueness of the local optimizers in these
//                             enclosures
//         'NumberOfOptis'   : returns the number of enclosures computed
//         'Minimum'     '   : returns the enclosure for the minimum value
//         'Err'             : returns an error code
// Description:
//    The enclosures for the global minimizers of 'f' are computed by calling
//    function 'GlobalOptimize'. Then a verification step is applied.
//    The enclosures for the global minimizers of 'f' are stored in the
//    interval vector 'OptiVector', the corresponding information on the
//    uniqueness of the local minimizers in these enclosures is stored in the
//    integer vector 'InfoVector'. The number of enclosures computed is
//    returned in the integer variable 'NumberOfOptis'. The enclosure for the
//    global minimum value is returned in the variable 'Minimum'. If an error
//    occurs, the value of 'Err' is different from 0. The result vectors
//    'OptiVector' and 'InfoVector' are automatically resized to standard
//    index range with lower index bounds 1. The maximum number of optimizers
//    which can be computed is given by 'MaxNumberOfOptis'. It is bounded by
//    'MaxCountGOp1'. If your system can handle more than 'MaxCountGOp1' optimizers,
//    you can increase its value.
//----------------------------------------------------------------------------
void AllGOp1 ( ddf_FctPtr  f,
               interval    Start,
               real        Epsilon,
               ivector&    OptiVector,
               intvector&  InfoVector,
               int&        NumberOfOptis,
               interval&   Minimum,
               int&        Err,
               int         MaxNumberOfOptis )
{
  int       i;
  real      MinEpsilon;
  Pair1Ptr  ResultList;

  if (1 <= MaxNumberOfOptis && MaxNumberOfOptis <= MaxCountGOp1) {
    // Start global optimization method
    //---------------------------------
    Err = NoError;
    MaxOptiNo = MaxNumberOfOptis;         // Initialize global static variable

    MinEpsilon = succ(1.0) - 1.0;         // Relative machine accuracy (1 ulp)
    if (Epsilon < MinEpsilon)
      Epsilon = MinEpsilon;               // Set 'Epsilon' to 1 ulp accuracy

    GlobalOptimize(f,Start,Epsilon,ResultList,NumberOfOptis,Minimum);

    if (NumberOfOptis > MaxNumberOfOptis)
      { Err = NotAllOptis;  NumberOfOptis = MaxNumberOfOptis; }

    ListToVector(ResultList,OptiVector,InfoVector);

    for (i = 1; i <= NumberOfOptis; i++)                  // Verification step
      VerificationStep(f,OptiVector[i],InfoVector[i]);    // for the enclosure
                                                          // intervals
    FreeAll(ResultList);                                  //------------------
  }
  else {
    Err = IllMaxOptiNo;
    NumberOfOptis = 0;
  }
}




