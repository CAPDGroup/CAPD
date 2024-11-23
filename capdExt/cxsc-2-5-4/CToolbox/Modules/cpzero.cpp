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


// Purpose: Determination and enclosure of a root of a complex polynomial,
//    and of the coefficients of the deflated polynomial.
// Method: The root of a complex polynomial matches with the eigenvalue of
//    its companion matrix and the coefficients of the deflated polynomial
//    match to the components of the corresponding eigenvector. The eigen-
//    vector and eigenvalue are determined using the simplified Newton
//    iteration with iterative refinement.
// Global functions:
//    CPolyZero()      : computes an enclosure for a root and for the
//                       deflated complex polynomial
//    CPolyZeroErrMsg(): delivers an error message text
//----------------------------------------------------------------------------
#include <string.h>        // String Handling
#include <r_util.hpp>      // Real utility functions
#include <ci_util.hpp>     // Complex interval utility functions
#include <cpzero.hpp>

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// Error codes used in this module. In the comments below p[0],..., p[n] are
// the coefficients of a polynomial p.
//----------------------------------------------------------------------------
const int
  NoError   = 0,   // No Error occurred
  ZeroPoly  = 1,   // Zero polynomial, i.e. n = 0 and p[0] = 0
  ConstPoly = 2,   // Constant polynomial, i.e. n = 0  and p[0] != 0
  InvFailed = 3,   // Inversion of the Jacobian failed
  VerFailed = 4,   // Verified inversion of the Jacobian failed
  IncFailed = 5;   // Inclusion failed

//----------------------------------------------------------------------------
// Error messages depending on the error code.
//----------------------------------------------------------------------------
char* CPolyZeroErrMsg ( int Err )
{
  static char Msg[80] = "";

  if (Err != NoError) {
    char Hlp[60];

    switch (Err) {
      case ZeroPoly:
        sprintf(Hlp,"Zero polynomial occurred");                  break;
      case ConstPoly:
        sprintf(Hlp,"Constant polynomial != 0 occurred");         break;
      case InvFailed:
        sprintf(Hlp,"Inversion of the Jacobian failed");          break;
      case IncFailed:
        sprintf(Hlp,"Inclusion failed");                          break;
      case VerFailed:
        sprintf(Hlp,"Verified inversion of the Jacobian failed"); break;
      default:
        strcpy(Hlp,"Code not defined");
    }
    sprintf(Msg,"Error: %s!",Hlp);
  }
  return(Msg);
} // CPolyZeroErrMsg

//----------------------------------------------------------------------------
// Function used to get the maximum norm of 'p' where 'p' is interpreted as a
// complex vector.
//----------------------------------------------------------------------------
static real MaxNorm ( CPolynomial p )
{
  int  i;
  real max, tmp;

  max = abs(p[0]);
  for (i = 1; i <= Deg(p); i++) {
    tmp = abs(p[i]);
    if (tmp > max)  max = tmp;
  }
  return max;
}

//----------------------------------------------------------------------------
// Purpose: Determination of a floating-point approximation for a root of
//    polynomial p and for the deflated polynomial.
// Parameters:
//    In    : 'p'  : represents a complex polynomial by its coefficients
//    In/Out: 'z'  : starting approximation (in) for a root and is returned
//                   as an improved floating-point approximation
//    Out   : 'q'  : returns an approximation (for the coefficients of) the
//                   deflated polynomial
//            'Err': error flag
// Description: For a given starting approximation 'z', a floating-point
//    Newton iteration is performed to get improved approximations for a root
//    'z' of polynomial 'p', as well as for the coefficients 'q' of the
//    deflated polynomial.
//----------------------------------------------------------------------------
static void Approximation ( CPolynomial  p, complex& z,
                            CPolynomial& q, int&     Err )
{
  const int      kmax  = 50;              // Maximum number of iteration steps
  const real     eps   = 1E-9;            // Relative error of approximation

  int            stop, k, i, n = Deg(p);
  complex        t;
  cdotprecision  accu;
  CPolynomial    delta(n-1), d(n-1), w(n-1), qHlp(n-2);

  Err = NoError;                            // Initialization

  q[n-1] = p[n];                            // Determination of an approximate
  for (i = n-2; i >= 0; i--)                // eigenvector q
    q[i] = q[i+1]*z + p[i+1];               //--------------------------------


  k = 0;      // Floating-point iteration of z and q
  do {        //------------------------------------
    accu = -p[0];                             // Computation of the defect
    accumulate(accu,-z,q[0]);                 //     d = (A - z*I)*q
    d[0] = rnd(accu);                         //--------------------------
    for (i = 1; i <= n-1; i++) {
      accu = q[i-1];
      accumulate(accu,-z,q[i]);
      accu -= p[i];
      d[i] = rnd(accu);
    }


    w[n-1] = q[n-1];                          // Computation of R, the inverse
    for (i = n-2; i >= 0; i--) {              // of the Jacobian J
      accu = q[i];                            //------------------------------
      accumulate(accu,z,w[i+1]);
      w[i] = rnd(accu);
    }

    if (w[0] == 0.0)
      Err = InvFailed;
    else {
      delta[n-1] = d[n-1];                  // Computation of delta^(k+1) =
      for (i = n-2; i >= 0; i--)            //  g(delta^(k))
        delta[i] = d[i] + z*delta[i+1];     //-----------------------------
      t = delta[0]/w[0];
      for (i = 0; i <= n-2; i++)
        delta[i] = -delta[i+1] + t*w[i+1];
      delta[n-1] = t;

      for (i = 0; i <= n-2; i++) {          // Determine the new iterates of
        q[i] = q[i] + delta[i];             // q and z and store the first
        qHlp[i] = q[i];                     // n-1 components of q in qHlp
      }                                     //------------------------------
      z = z + delta[n-1];

      k++;
    } //else

    // Stop iteration if the relative round-off error of the coefficients
    // and of the zero is approximately <= a given epsilon. Stop also if
    // the maximum number of iteration steps is exceeded, or if an error
    // occurred. Some C++ compilers (e.g. HP C++ V2.0) do not support
    // temporaries of class variables in an || expression. Thus, 'stop'
    // is used to avoid this error.
    //------------------------------------------------------------------
    // stop = in(ddelta,ddelta_old) || (k == kmax);
    // stop = (MaxNorm(delta) / Max(MaxNorm(qHlp),abs(z)) <= eps) ||
    //        (Err != NoError) || (k == kmax);
    stop = (MaxNorm(delta) / Max(MaxNorm(qHlp),abs(z)) <= eps);
    stop = stop || (Err != NoError) || (k == kmax);

  } while ( !stop );

} //Approximation

//----------------------------------------------------------------------------
// Purpose: Determination of enclosures for a root of polynomial p and for
//    the deflated polynomial.
// Parameters:
//    In : 'p'  : represents a complex polynomial by its coefficients
//         'q'  : approximation of (the coefficients of) the deflated
//                polynomial
//         'z'  : floating-point approximation of a root
//    Out: 'qq' : returns an enclosure for (the coefficients of) the
//                deflated polynomial
//         'zz' : returns an enclosure of the root
//         'Err': error flag
// Description: For starting floating-point approximations 'z' of a root
//    and 'q' of the deflated polynomial an interval residual iteration is
//    performed to get enclosures 'zz' for the root of polynomial 'p', as
//    well as for the coefficients of the deflated polynomial.
//----------------------------------------------------------------------------
static void IntervalIteration ( CPolynomial p,  CPolynomial   q,
                                complex     z,  CIPolynomial& qq,
                                cinterval&  zz, int&    Err       )
{
  const int       kmax = 10;              // Maximum number of iteration steps

  int             stop, i, k, n = Deg(p);
  real            eps;
  cinterval       tt;
  cidotprecision  accu;
  CIPolynomial    dd(n-1), vv(n-1), ww(n-1);
  CIPolynomial    ddelta(n-1), ddelta_0(n-1), ddelta_old(n-1);

  Err = NoError;                            // Initialization

  accu = -p[0];                             // Compute an inclusion of the
  accumulate(accu,-z,q[0]);                 // defect [d] = (A - z*I)*q
  rnd(accu,dd[0]);                          //----------------------------
  for (i = 1; i <= n-1; i++) {
    accu = q[i-1];
    accumulate(accu,-z,q[i]);
    accumulate(accu,-1.0,p[i]);
    dd[i] = rnd(accu);
  }

  ww[n-1] = q[n-1];                        // Compute an inclusion of [R], the
  for (i = n-2; i >=0; i--) {              // inverse of the Jacobian J
    accu = q[i];                           //---------------------------------
    accumulate(accu,z,ww[i+1]);
    ww[i] = rnd(accu);
  }

  if (_cinterval(0.0) < ww[0])
    Err = VerFailed;
  else {
    ddelta_0[n-1] = dd[n-1];               // Compute a starting
    for (i = n-2; i >= 0; i--) {           // inclusion [y]^(0)
      accu = dd[i];                        //-------------------
      accumulate(accu,z,ddelta_0[i+1]);
      ddelta_0[i] = rnd(accu);
    }
    tt = ddelta_0[0]/ww[0];
    for (i = 0; i <= n-2; i++) {
      accu = -ddelta_0[i+1];
      accumulate(accu,tt,ww[i+1]);
      ddelta_0[i] = rnd(accu);
    }
    ddelta_0[n-1] = tt;

    k = 0;  ddelta = ddelta_0;             // Interval iteration
    do {                                   //-------------------
      ddelta_old = ddelta;

      if (k <= 3)                          // Slightly enlarge the
        eps = 0.125;                       // inclusion interval
      else if (k <= 6)                     //---------------------
        eps = 0.5;
      else
        eps = 5.0;
      ddelta_old = Blow(ddelta,eps);

      vv[n-1] = 0.0;                       // Determine a new
      for (i = n-2; i >= 0; i--) {         // inclusion interval
        accu = 0.0;                        //-------------------
        accumulate(accu,ddelta_old[n-1],ddelta_old[i]);
        accumulate(accu,z,vv[i+1]);
        vv[i] = rnd(accu);
      }
      vv[0] = vv[0] / ww[0];
      for (i = 0; i <= n-2; i++) {
        accu = ddelta_0[i];
        accu += vv[i+1];
        accumulate(accu,-vv[0],ww[i+1]);
        ddelta[i] = rnd(accu);
      }
      ddelta[n-1] = ddelta_0[n-1] - vv[0];

      k++;

      // Some C++ compilers (e.g. HP C++ V2.0) do not support temporaries
      // of class variables in an || expression. Thus, 'stop' is used to
      // avoid this error.
      //-----------------------------------------------------------------
      // stop = in(ddelta,ddelta_old) || (k == kmax);
      stop = in(ddelta,ddelta_old);
      stop = stop || (k == kmax);

    } while ( !stop );

  } //else

  if ( in(ddelta,ddelta_old) ) {           // Verification of the result
    for (i = 0; i <= n-2; i++) {           //---------------------------
      accu = q[i];
      accu += ddelta[i];
      qq[i] = rnd(accu);
    }
    qq[n-1] = p[n];
    accu = z;
    accu += ddelta[n-1];
    zz = rnd(accu);
  }
  else
    Err = IncFailed;

} //IntervalIteration

//----------------------------------------------------------------------------
// Purpose: Determination of enclosures for a root of polynomial p and for
//    the deflated polynomial
// Parameters:
//    In : 'p'  : represents a complex polynomial by its coefficients
//         'z'  : floating-point approximation of a root
//    Out: 'qq' : returns an enclosure for (the coefficients of) the
//                deflated polynomial
//         'zz' : returns an enclosure of the root
//         'Err': error flag
// Description: For starting approximations 'z' of a root, enclosures 'zz'
//    for a root of polynomial 'p', as well as 'qq' for the coefficients of
//    the deflated polynomial are computed. Since the root of a complex
//    polynomial matches with an eigenvalue of its companion matrix and the
//    coefficients of the deflated polynomial match to the components of the
//    corresponding eigenvector, the eigenvector and eigenvalue are determined
//    by iterative refinement. First, good floating-point approximations are
//    computed, and then verified enclosures for a root and for a deflated
//    polynomial are returned.
//----------------------------------------------------------------------------
void CPolyZero ( CPolynomial   p,  complex    z,
                 CIPolynomial& qq, cinterval& zz, int& Err )
{
  int  n = Deg(p);

  Err = NoError;                                 // Initialization

  if ( (n == 0) || (n == 1) )
    if (n == 0)                                  // Polynomial of degree n = 0
                                                 //---------------------------
      if (p[0] == 0.0) Err = ZeroPoly;           // Zero polynomial
      else             Err = ConstPoly;          // Constant polynomial

    else                                         // Polynomial of degree n = 1
      if (p[1] == 0.0)                           //---------------------------
        if (p[0] == 0.0) Err = ZeroPoly;         // Zero polynomial
        else             Err = ConstPoly;        // Constant polynomial
      else {
        zz = -_cinterval(p[0]) / p[1];
        qq[0] = p[1];
      }
  else {                                         // Polynomial of degree n > 1
    CPolynomial  q(n-1);                         //---------------------------

    Approximation(p,z,q,Err);
    if (Err == NoError) IntervalIteration(p,q,z,qq,zz,Err);
  }
} //CPolyZero




