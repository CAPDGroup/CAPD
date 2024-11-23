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

/* CVS $Id: nlfzero.cpp,v 1.15 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: nlfzero (implementation)
// Purpose: Computing enclosures for all zeros of a continuously
//    differentiable one-dimensional, real-valued function.
// Method: Extended interval Newton's method.
// Global functions:
//    AllZeros()      : computes enclosures for all zeros
//    AllZerosErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#include <string.h>         // String handling
#include <interval.hpp>     // Interval arithmetic
#include <i_util.hpp>       // Interval utility functions
#include <xi_ari.hpp>       // Extended interval arithmetic
#include <nlfzero.hpp>

using namespace cxsc;
using namespace std;

static int MaxZeroNo;       // Maximum number of zeros to be computed

//----------------------------------------------------------------------------
// Error codes used in this module.
//----------------------------------------------------------------------------
const int
  NoError      = 0,  // No error occurred
  IllMaxZeroNo = 1,  // Illegal value for the maximum number of zeros
  NotAllZeros  = 2;  // Not all zeros computed

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* AllZerosErrMsg ( int Err )
{
  static char Msg[80] = "";

  switch (Err) {
    case NoError:
      break;  // Blank message
    case IllMaxZeroNo:
      sprintf(Msg,"Error: Parameter for maximum number of zeros "
                  "must lie in 1,...,%1d!",MaxCount); break;
    case NotAllZeros:
      sprintf(Msg,"Warning: Not all zeros found due to the user "
                  "limit of %1d zero(s)!",MaxZeroNo); break;
    default:
      strcpy(Msg,"Error: Code not defined!");
  }
  return(Msg);
}

//----------------------------------------------------------------------------
// Purpose: Recursive function for the execution of the extended interval
//    Newton's method for the function 'f'.
// Parameters:
//    In    : 'f'         : must be declared for the 'DerivType' to enable
//                          the internal use of the differentiation
//                          arithmetic 'ddf_ari'
//            'y'         : specifies the starting interval
//            'Epsilon'   : specifies the desired relative accuracy
//                          (interval diameter) of the result intervals
//            'yUnique'   : signals whether it is already verified that the
//                          actual interval 'y' contains a unique zero
//    Out   : 'ZeroVector': stores the enclosures for the zeros of 'f'
//            'InfoVector': stores the corresponding information on the
//                          uniqueness of the zero in each enclosure
//    In/Out: 'ZeroNo'    : represents the total number of enclosures
//                          computed at the end of the recursion
// Description:
//    The function 'XINewton' is recursively called whenever the extended
//    interval division results in a bisection of the actual interval 'y' in
//    two intervals 'y1' and 'y2', and the tolerance condition is not ful-
//    filled yet. Otherwise, the enclosures for the zeros of 'f' are stored
//    in the interval vector 'ZeroVector', the corresponding information on
//    the uniqueness of the zero in each enclosure is stored in the integer
//    vector 'InfoVector'. The maximum number of enclosures to be computed
//    is given by the value of the global variable 'MaxZeroNo'.
//----------------------------------------------------------------------------
static void XINewton ( ddf_FctPtr f,
                       interval   y,
                       real       Epsilon,
                       int        yUnique,
                       ivector&   ZeroVector,
                       intvector& InfoVector,
                       int&       ZeroNo )
{
  if (ZeroNo > MaxZeroNo) return;       // Stop due to user specified limit
                                        //---------------------------------
  xinterval  z;
  interval   fc, fy, dfy;
  ivector    V(2);
  real       c;
  int        i;

  dfEval(f,y,fy,dfy);                   // Compute f(y) and f'(y)
  if ( !in(0.0,fy) ) return;            // Start if 0 in f(y), else do nothing
                                        //------------------------------------
  c  = mid(y);
  fEval(f,_interval(c),fc);    // Compute f(c) and f'(c).
  z  = c - fc % dfy;           // Extended interval Newton step.
  V = y & z;                   // Intersect interval y and extended interval z
                               // resulting in two intervals V[1] and V[2].
                               //---------------------------------------------
  if (V[1] == y) {                                       // Stagnation, so 'y'
    V[1] = interval(Inf(y),c);                           // must be bisected.
    V[2] = interval(c,Sup(y));                           //-------------------
  }

  if ( (V[1] != EmptyIntval()) && (V[2] == EmptyIntval()) )
    yUnique = yUnique || in(V[1],y);        // Inner inclusion ===> uniqueness
  else                                      //--------------------------------
    yUnique = 0;

  for (i = 1; i <= 2; i++) {
    if (V[i] == EmptyIntval()) continue;   // Next i

    if (RelDiam(V[i]) <= Epsilon) {        // No more Newton steps

        fEval(f,V[i],fy);                       // Compute f(V[i])
        if ( in(0.0,fy) ) {                     // Store enclosure and info
          ZeroNo++;                             //-------------------------

          if (ZeroNo > MaxZeroNo) return;       // Stop due to user
                                                // specified limit
          if (ZeroNo > Ub(ZeroVector))          //-----------------
            { DoubleSize(ZeroVector); DoubleSize(InfoVector); }

          ZeroVector[ZeroNo] = V[i];            // Store enclosure of the zero
          InfoVector[ZeroNo] = yUnique;         // Store uniqueness info
        }                                       //----------------------------
      }
    else  // Recursive call of 'XINewton' for interval 'V[i]'
      XINewton(f,V[i],Epsilon,yUnique,ZeroVector,InfoVector,ZeroNo);
  } // for (i ...
} // XINewton

//----------------------------------------------------------------------------
// Purpose: Execution of a verification step including the use of an epsilon
//    inflation.
// Parameters:
//    In    : 'f'      : function of 'DerivType'
//    Out   : 'yUnique': returns 'true' if the verification is successful
//    In/Out: 'y'      : interval enclosure to be verified
// Description: This function checks the uniqueness of the zero enclosed in
//    the variable 'y' by an additional verification step including the use
//    of an epsilon inflation of the iterates.
//----------------------------------------------------------------------------
static void VerificationStep ( ddf_FctPtr f, interval& y, int& unique )
{
  const int kmax = 10;

  interval yIn, yOld, fc, fy, dfy;
  real     c, eps;
  int      k;

  k = 0;  yIn = y;  eps = 0.25; unique = 0;                 // Initializations
                                                            //----------------
  while ( !unique && (k < kmax) ) {                // Do kmax loops to achieve
                                                   // inclusion
    yOld  = Blow(y,eps);                           // Epsilon inflation of 'y'
    dfEval(f,yOld,fy,dfy);                         //-------------------------

    if ( in(0.0,dfy) ) break;                      // No verification possible
                                                   //-------------------------
    k++;  c = mid(yOld);                       // Perform interval Newton step
    fEval(f,_interval(c),fc);                  //-----------------------------
    y = c - fc / dfy;

    if (Disjoint(y,yOld)) break;                   // No verification possible
                                                   //-------------------------
    unique = in(y,yOld);                    // Inner inclusion ===> uniqueness
    y      = y & yOld;                      // Intersection with old value
    if (y == yOld)                          //--------------------------------
      eps = eps * 8.0;                      // Increase the value of 'eps'
  } // while

  if (!unique) y = yIn;
} // VerificationStep

//----------------------------------------------------------------------------
// Purpose: Computation of enclosures for all zeros of a continuously
//    differentiable one-dimensional, real-valued function.
// Parameters:
//    In : 'f'                : objective function, must be declared for the
//                              'DerivType' to enable the internal use of
//                              the differentiation arithmetic 'ddf_ari'
//         'Start',           : specifies the starting interval
//         'Epsilon'          : specifies the desired relative accuracy
//                              (interval diameter) of the result intervals
//    Out: 'ZeroVector'       : stores and returns the enclosures for the
//                              zeros of 'f'
//         'InfoVector'       : stores the corresponding information on the
//                              uniqueness of the zeros in these enclosures
//         'NumberOfZeros'    : returns the number of enclosures computed
//         'Err'              : returns an error code
//         'MaxNumberOfZeros' : maximum number of zeros to be computed
//                              (default value: 'MaxCount' = compute all)
// Description:
//    The enclosures for the zeros of 'f' are computed by calling function
//    'XINewton'. Then an additional verification step is applied to those
//    enclosures which have not been verified. If an error occurs the value
//    of 'Err' is different from 0. The result vectors 'ZeroVector' and
//    'InfoVector' are automatically resized to standard index range with
//    lower index bounds 1. The maximum number of zeros which can be computed
//    is given by 'MaxNumberOfZeros'. It is bounded by 'MaxCount'. If your
//    system can handle more than 'MaxCount' zeros, you can increase its
//    value.
//----------------------------------------------------------------------------
void AllZeros ( ddf_FctPtr f,
                interval   Start,
                real       Epsilon,
                ivector&   ZeroVector,
                intvector& InfoVector,
                int&       NumberOfZeros,
                int&       Err,
                int        MaxNumberOfZeros )
{
  int       i, MinSize;
  real      MinEpsilon;

  if (1 <= MaxNumberOfZeros && MaxNumberOfZeros <= MaxCount) {
    MaxZeroNo = MaxNumberOfZeros;   // Initialize global static variable

    // Resize the result vectors to 'MinSize' components with standard index
    // range. Use 'MaxNumberOfZeros' if specified and less than 'MaxCount'.
    //----------------------------------------------------------------------
    MinSize = (MaxNumberOfZeros == MaxCount) ? 10 : MaxNumberOfZeros;

    Resize(ZeroVector,MinSize);
    Resize(InfoVector,MinSize);

    // Start extended interval Newton method
    //--------------------------------------
    Err = NoError;  NumberOfZeros = 0;

    MinEpsilon = succ(1.0) - 1.0;       // Relative machine accuracy (1 ulp)
    if (Epsilon < MinEpsilon)
      Epsilon = MinEpsilon;             // Set 'Epsilon' to 1 ulp accuracy

    XINewton(f,Start,Epsilon,0,ZeroVector,InfoVector,NumberOfZeros);

    if (NumberOfZeros > MaxNumberOfZeros)
      { Err = NotAllZeros;  NumberOfZeros = MaxNumberOfZeros; }
  }
  else {
    Err = IllMaxZeroNo;
    NumberOfZeros = 0;
  }

  // Resize the result vectors to appropriate length if no error
  // has occurred. Otherwise, resize them to minimum dimensions
  // with undefined value.
  //------------------------------------------------------------
  if (NumberOfZeros == 0) {
    Resize(ZeroVector,1);
    Resize(InfoVector,1);
  }
  else {
    Resize(ZeroVector,NumberOfZeros);     // Values undefined
    Resize(InfoVector,NumberOfZeros);     //-----------------
  }

  for (i = 1; i <= NumberOfZeros; i++)                   // Verification step
    if (InfoVector[i] == 0)                              // for the enclosures
      VerificationStep(f,ZeroVector[i],InfoVector[i]);   //-------------------
} // AllZeros




