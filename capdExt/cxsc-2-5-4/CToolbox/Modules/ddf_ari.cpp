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

/* CVS $Id: ddf_ari.cpp,v 1.17 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: ddf_ari (implementation)
// Purpose: Definition of an interval differentiation arithmetic which allows
//    function evaluation with automatic differentiation up to second order.
// Method: Overloading of operators and elementary functions for operations
//    of data type 'DerivType'.
// Class DerivType:
//    DerivType()             : constructors
//    operators +, -, *, /    : operators of diff. arithmetic
//    operator =              : assignment operator
//    DerivConst()
//    DerivVar()              : to define derivative constants/variables
//    fValue()
//    dfValue()
//    ddfValue()              : to get function and derivative values
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEval()                 : to compute function value only
//    dfEval()                : to compute function and first derivative
//                              value
//    ddfEval()               : to compute function, first, and second
//                              derivative value
//----------------------------------------------------------------------------
#include <imath.hpp>       // Interval mathematical functions
#include <ddf_ari.hpp>

using namespace cxsc;
using namespace std;

// The local variable 'DerivOrder' is used to select the highest order of
// derivative which is computed. Its default value is 2, and normally the
// first and the second derivatives are computed.
//----------------------------------------------------------------------------
#ifdef _WIN32
static __declspec(thread) int DerivOrder = 2;
#elif __APPLE__ && !CXSC_FORCE_TLS
static int DerivOrder = 2;
#else
static __thread int DerivOrder = 2;
#endif

//----------------------------------------------------------------------------
// Constructors and assignment operator
//----------------------------------------------------------------------------
DerivType::DerivType ( )
{
  f   = 0.0;
  df  = 0.0;
  ddf = 0.0;
}

DerivType::DerivType ( const interval& x_f, const interval& x_df, const interval& x_ddf )
{
  f   = x_f;
  df  = x_df;
  ddf = x_ddf;
}

DerivType::DerivType ( const DerivType& u )
{
  f   = u.f;
  df  = u.df;
  ddf = u.ddf;
}

DerivType& DerivType::operator= ( const DerivType& u )
{
  f   = u.f;
  df  = u.df;
  ddf = u.ddf;
  return *this;
}

//----------------------------------------------------------------------------
// Transfer functions for constants and variables
//----------------------------------------------------------------------------
DerivType DerivConst ( const real& c )                    // Generate constant
{                                                         //------------------
  DerivType u;

  u.f   = c;
  u.df  = 0.0;
  u.ddf = 0.0;
  return u;
}

DerivType DerivConst ( const interval& c )                // Generate constant
{                                                         //------------------
  DerivType u;

  u.f   = c;
  u.df  = 0.0;
  u.ddf = 0.0;
  return u;
}

DerivType DerivVar ( const real& v )                      // Generate variable
{                                                         //------------------
  DerivType u;

  u.f   = v;
  u.df  = 1.0;
  u.ddf = 0.0;
  return u;
}

DerivType DerivVar ( const interval& v )                  // Generate variable
{                                                         //------------------
  DerivType u;

  u.f   = v;
  u.df  = 1.0;
  u.ddf = 0.0;
  return u;
}

//----------------------------------------------------------------------------
// Unary operators + and - for DerivType operands
//----------------------------------------------------------------------------
inline DerivType operator+ ( DerivType& u )
  { return u; }

DerivType operator- ( const DerivType& u )
{
  DerivType res;

  res.f = -u.f;
  if (DerivOrder > 0) {
    res.df = -u.df;
    if (DerivOrder > 1) res.ddf = -u.ddf;
  }
  return res;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for two DerivType operands
//----------------------------------------------------------------------------
DerivType operator+ ( const DerivType& u, const DerivType& v )
{
  DerivType res;

  res.f = u.f + v.f;
  if (DerivOrder > 0) {
    res.df = u.df + v.df;
    if (DerivOrder > 1) res.ddf = u.ddf + v.ddf;
  }
  return res;
}

DerivType operator- ( const DerivType& u, const DerivType& v )
{
  DerivType res;

  res.f = u.f - v.f;
  if (DerivOrder > 0) {
    res.df = u.df - v.df;
    if (DerivOrder > 1) res.ddf = u.ddf - v.ddf;
  }
  return res;
}

DerivType operator* ( const DerivType& u, const DerivType& v )
{
  DerivType res;

  res.f = u.f * v.f;
  if (DerivOrder > 0) {
    res.df = u.df*v.f + u.f*v.df;
    if (DerivOrder > 1) res.ddf = u.ddf*v.f + 2.0*u.df*v.df + u.f*v.ddf;
  }
  return res;
}

DerivType operator/ ( const DerivType& u, const DerivType& v )
{
  DerivType res;
  interval  h1, h2;

  h1 = u.f / v.f;   // Can propagate 'division by zero' error
  res.f = h1;       //---------------------------------------
  if (DerivOrder > 0) {
    h2 = (u.df - h1*v.df) / v.f;
    res.df = h2;
    if (DerivOrder > 1) res.ddf = (u.ddf - h1*v.ddf - 2.0*h2*v.df)/v.f;
  }
  return res;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for one interval and one DerivType operand
//----------------------------------------------------------------------------
DerivType operator+ ( const interval& u, const DerivType& v )
{
  DerivType res;

  res.f = u + v.f;
  if (DerivOrder > 0) {
    res.df = v.df;
    if (DerivOrder > 1) res.ddf = v.ddf;
  }
  return res;
}

DerivType operator- ( const interval& u, const DerivType& v )
{
  DerivType res;

  res.f = u - v.f;
  if (DerivOrder > 0) {
    res.df = -v.df;
    if (DerivOrder > 1) res.ddf = - v.ddf;
  }
  return res;
}

DerivType operator* ( const interval& u, const DerivType& v )
{
  DerivType res;

  res.f = u*v.f;
  if (DerivOrder > 0) {
    res.df = u * v.df;
    if (DerivOrder > 1) res.ddf = u*v.ddf;
  }
  return res;
}

DerivType operator/ ( const interval& u, const DerivType& v )
{
  DerivType res;
  interval  h1, h2;

  h1  = u / v.f;    // Can propagate 'division by zero' error
  res.f = h1;       //---------------------------------------
  if (DerivOrder > 0) {
    h2   = -h1 * v.df / v.f;
    res.df = h2;
    if (DerivOrder > 1) res.ddf = (-h1*v.ddf - 2.0*h2*v.df)/v.f;
  }
  return res;
}

DerivType operator+ ( const DerivType& u, const interval& v )
{
  DerivType res;

  res.f = u.f + v;
  if (DerivOrder > 0) {
    res.df = u.df;
    if (DerivOrder > 1) res.ddf = u.ddf;
  }
  return res;
}

DerivType operator- ( const DerivType& u, const interval& v )
{
  DerivType res;

  res.f = u.f - v;
  if (DerivOrder > 0) {
    res.df = u.df;
    if (DerivOrder > 1) res.ddf = u.ddf;
  }
  return res;
}

DerivType operator* ( const DerivType& u, const interval& v )
{
  DerivType res;

  res.f = u.f * v;
  if (DerivOrder > 0) {
    res.df = u.df * v;
    if (DerivOrder > 1) res.ddf = u.ddf * v;
  }
  return res;
}

DerivType operator/ ( const DerivType& u, const interval& v )
{
  DerivType res;

  res.f = u.f / v;       // Can propagate 'division by zero' error
  if (DerivOrder > 0) {  //---------------------------------------
    res.df = u.df / v;
    if (DerivOrder > 1) res.ddf = u.ddf / v;
  }
  return res;
}

//----------------------------------------------------------------------------
// Operands for +, -, *, and / for one real and one DerivType operand
//----------------------------------------------------------------------------
DerivType operator+ ( const real& u, const DerivType& v )
  { return( interval(u) + v ); }

DerivType operator- ( const real& u, const DerivType& v )
  { return( interval(u) - v ); }

DerivType operator* ( const real& u, const DerivType& v )
  { return( interval(u) * v ); }

DerivType operator/ ( const real& u, const DerivType& v )
  { return( interval(u) / v ); }    // Can propagate 'division by zero' error
                                    //---------------------------------------
DerivType operator+ ( const DerivType& u, const real& v )
  { return( u + interval(v) ); }

DerivType operator- ( const DerivType& u, const real& v )
  { return( u - interval(v) ); }

DerivType operator* ( const DerivType& u, const real& v )
  { return( u * interval(v) ); }

DerivType operator/ ( const DerivType& u, const real& v )
  { return( u / interval(v) ); }    // Can propagate 'division by zero' error
                                    //---------------------------------------

//----------------------------------------------------------------------------
// Elementary functions for DerivType arguments
//----------------------------------------------------------------------------
DerivType sqr ( const DerivType& u )
{
  DerivType res;

  res.f = Power(u.f,2);
  if (DerivOrder > 0) {
    res.df = 2.0*u.f*u.df;
    if (DerivOrder > 1) res.ddf = 2.0 * (Power(u.df,2) + u.f*u.ddf);
  }
  return res;
}

DerivType power ( const DerivType& u, int k )
{
  DerivType res;
  interval  h1;

  if (k == 0)
     res = DerivConst(1.0);
  else if (k == 1)
    res = u;
  else {
    res.f = Power(u.f, k);
    if (DerivOrder > 0) {
      h1 = double(k) * Power(u.f, k-1);
      res.df = h1 * u.df;
      if (DerivOrder > 1)
        res.ddf = h1 * u.ddf + double(k*(k-1))*Power(u.f,k-2)*Power(u.df,2);
    }
  }
  return res;
}

DerivType sqrt ( const DerivType& u )
{
  DerivType res;
  interval  h1, h2;

  h1 = sqrt(u.f);             // Can propagate domain error
  res.f = h1;                 //---------------------------
  if (DerivOrder > 0) {
    h1 = 0.5/h1;
    h2 = u.df*h1;
    res.df = h2;
    if (DerivOrder > 1) res.ddf = u.ddf*h1 - 0.5*u.df/u.f*h2;
  }
  return res;
}

DerivType exp ( const DerivType& u )
{
  DerivType res;
  interval  h1, h2;

  h1 = exp(u.f);
  res.f = h1;
  if (DerivOrder > 0) {
    h2 = h1*u.df;
    res.df = h2;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*u.df;
  }
  return res;
}

DerivType ln ( const DerivType& u )
{
  DerivType res;
  interval  h;

  res.f = ln(u.f);            // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h = u.df/u.f;
    res.df = h;
    if (DerivOrder > 1) res.ddf = (u.ddf - h*u.df) / u.f;
  }
  return res;
}

DerivType sin ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1;

  h0 = sin(u.f);
  res.f = h0;
  if (DerivOrder > 0) {
    h1 = cos(u.f);
    res.df = h1*u.df;
    if (DerivOrder > 1) res.ddf = h1*u.ddf - h0*sqr(u.df);
  }
  return res;
}

DerivType cos ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1;

  h0 = cos(u.f);
  res.f = h0;
  if (DerivOrder > 0) {
    h1 = -sin(u.f);
    res.df = h1*u.df;
    if (DerivOrder > 1) res.ddf = h1*u.ddf - h0*sqr(u.df);
  }
  return res;
}

DerivType tan ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1, h2;

  h0 = tan(u.f);              // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (DerivOrder > 0) {                     // The subdistributive law implies
    h1 = sqr(h0)+1.0;                       //   h0 * (h0^2 + 1) <= h0^3 + h0.
    h2 = 2.0*h0*h1;                         // So, the first form is used.
    res.df = h1*u.df;                       //--------------------------------
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType cot ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1, h2;

  h0 = cot(u.f);              // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (DerivOrder > 0) {                     // The subdistributive law implies
    h1 = -(sqr(h0)+1.0);                    //   h0 * (h0^2 + 1) <= h0^3 + h0.
    h2   = -2.0*h0*h1;                      // So, the first form is used.
    res.df = h1*u.df;                       //--------------------------------
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType asin ( const DerivType& u )
{
  DerivType res;
  interval  h, h1, h2;

  res.f = asin(u.f);          // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h  = 1.0 - sqr(u.f);
    h1 = 1.0/sqrt(h);
    res.df = h1*u.df;
    h2 = u.f*h1/h;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType acos ( const DerivType& u )
{
  DerivType res;
  interval  h, h1, h2;

  res.f = acos(u.f);          // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h  = 1.0 - sqr(u.f);
    h1 = -1.0/sqrt(h);
    res.df = h1*u.df;
    h2 = u.f*h1/h;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType atan ( const DerivType& u )
{
  DerivType res;
  interval  h1, h2;

  res.f = atan(u.f);
  if (DerivOrder > 0) {
    h1 = 1.0 / (1.0 + sqr(u.f));
    res.df = h1*u.df;
    h2 = -2.0*u.f*sqr(h1);
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType acot ( const DerivType& u )
{
  DerivType res;
  interval  h1, h2;

  res.f = acot(u.f);          // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h1 = -1.0 / (1.0 + sqr(u.f));
    res.df = h1*u.df;
    h2 = 2.0*u.f*sqr(h1);
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType sinh ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1;

  h0 = sinh(u.f);
  res.f = h0;
  if (DerivOrder > 0) {
    h1 = cosh(u.f);
    res.df = h1*u.df;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h0*sqr(u.df);
  }
  return res;
}

DerivType cosh ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1;

  h0 = cosh(u.f);
  res.f = h0;
  if (DerivOrder > 0) {
    h1 = sinh(u.f);
    res.df = h1*u.df;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h0*sqr(u.df);
  }
  return res;
}

DerivType tanh ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1, h2;

  h0 = tanh(u.f);
  res.f = h0;
  if (DerivOrder > 0) {                     // The subdistributive law implies
    h1 = 1.0 - sqr(h0);                     //   h0 * (h0^2 - 1) <= h0^3 - h0.
    h2 = -2.0*h0*h1;                        // So, the first form is used.
    res.df = h1*u.df;                       //--------------------------------
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType coth ( const DerivType& u )
{
  DerivType res;
  interval  h0, h1, h2;

  h0 = coth(u.f);             // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (DerivOrder > 0) {                     // The subdistributive law implies
    h1 = 1.0 - sqr(h0);                     //   h0 * (h0^2 - 1) <= h0^3 - h0.
    h2 = -2.0*h0*h1;                        // So, the first form is used.
    res.df = h1*u.df;                       //--------------------------------
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType asinh ( const DerivType& u )
{
  DerivType res;
  interval  h, h1, h2;

  res.f = asinh(u.f);         // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h  = 1.0 + sqr(u.f);
    h1 = 1.0/sqrt(h);
    res.df = h1*u.df;
    h2 = -u.f*h1/h;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType acosh ( const DerivType& u )
{
  DerivType res;
  interval  h, h1, h2;

  res.f = acosh(u.f);         // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h  = sqr(u.f) - 1.0;
    h1 = 1.0/sqrt(h);
    res.df = h1*u.df;
    h2 = -u.f*h1/h;
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType atanh ( const DerivType& u )
{
  DerivType res;
  interval  h1, h2;

  res.f = atanh(u.f);         // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h1 = 1.0 / (1.0 - sqr(u.f));
    res.df = h1*u.df;
    h2 = 2.0*u.f*sqr(h1);
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

DerivType acoth ( const DerivType& u )
{
  DerivType res;
  interval  h1, h2;

  res.f = acoth(u.f);         // Can propagate domain error
  if (DerivOrder > 0) {       //---------------------------
    h1 = 1.0 / (1.0 - sqr(u.f));
    res.df = h1*u.df;
    h2 = 2.0*u.f*sqr(h1);
    if (DerivOrder > 1) res.ddf = h1*u.ddf + h2*sqr(u.df);
  }
  return res;
}

//----------------------------------------------------------------------------
// Predefined routines for evaluation of DerivType-functions
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of 'DerivType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
// Description: This function sets 'DerivOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
void fEval ( ddf_FctPtr f, interval x, interval& fx )
{
  DerivType xD, fxD;

  DerivOrder = 0;
  xD         = DerivVar(x);
  fxD        = f(xD);
  fx         = fxD.f;
  DerivOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the value of the first
//    derivative.
// Parameters:
//    In : 'f'  : function of 'DerivType'
//         'x'  : argument for evaluation of 'f'
//    Out: 'fx' : returns the function value 'f(x)'
//         'dfx': returns the first derivative value 'f'(x)'
// Description: This function sets 'DerivOrder' to 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the first derivative.
//----------------------------------------------------------------------------
void dfEval ( ddf_FctPtr f, interval x, interval& fx, interval& dfx )
{
  DerivType xD, fxD;

  DerivOrder = 1;
  xD         = DerivVar (x);
  fxD        = f(xD);
  fx         = fxD.f;
  dfx        = fxD.df;
  DerivOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value, the value of the first, and
//    the value of the second derivative.
// Parameters:
//    In : 'f'   : function of 'DerivType'
//         'x'   : argument for evaluation of 'f'
//    Out: 'fx'  : returns the function value 'f(x)'
//         'dfx' : returns the value of the first derivative 'f'(x)'
//         'ddfx': returns the value of the second derivative 'f''(x)'
// Description: This function keeps 'DerivOrder' = 2, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value, the value
//    of the first, and the value of the second derivative.
//----------------------------------------------------------------------------
void ddfEval ( ddf_FctPtr f, interval x, interval& fx, interval& dfx,
                                                       interval& ddfx )
{
  DerivType xD, fxD;

  xD   = DerivVar(x);
  fxD  = f(xD);
  fx   = fxD.f;
  dfx  = fxD.df;
  ddfx = fxD.ddf;
}




