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

/* CVS $Id: grad_ari.cpp,v 1.15 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: grad_ari (implementation)
// Purpose: Definition of a multi-dimensional interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to first order (i.e. gradient or Jacobian matrix).
// Method: Overloading of operators and elementary functions for operations
//    of data types 'GradType' and 'GTvector'.
// Class GradType:
//    GradType()              : constructors
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for arguments of types
//                              GradType, interval, and real
//    GradVar()               : to define GradType variables
//    fValue(), gradValue()   : to get function and gradient values
//    operators +, -, *, /    : operators of diff. arithmetic
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEvalG()                : to compute function value only
//    fgEvalG()               : to compute function and first derivative
//                              value (gradient)
// Class GTvector:
//    GTvector()              : constructors
//    ~GTvector()             : destructor
//    Dim()                   : to get the actual dimension
//    operator =              : assignment operator
//    operator []             : component access
//    fValue(), JacValue()    : to get function and derivative values
//    fEvalJ()                : to compute function value only
//    fJEvalJ()               : to compute function and first derivative
//                              value (Jacobian matrix)
//----------------------------------------------------------------------------
#include <cstdlib>          // For function 'exit()'
#include <iostream>         // I/O Handling
#include <imath.hpp>        // Interval mathematical functions
#include <i_util.hpp>       // Interval utility functions
#include <grad_ari.hpp>

using namespace cxsc;
using namespace std;

// The local variable 'GradOrder' is used to select the highest order of
// derivative which is computed. Its default value is 1, and normally
// the gradient or the Jacobian matrix are computed.
//----------------------------------------------------------------------------
static int GradOrder = 1;

//----------------------------------------------------------------------------
// Constructors and destructors
//----------------------------------------------------------------------------
GradType::GradType ( )                                          // Constructor
  { nmax = -1; }                                                //------------

GradType::GradType ( int n )                                  // Constructor
{                                                             //------------
  nmax = n;
  Resize(g,0,n);
}

GradType::GradType ( const GradType& A )                   // Copy constructor
{                                                          //-----------------
  nmax = A.nmax;
  g    = A.g;
}

GTvector::GTvector ( int n )                                  // Constructor
{                                                             //------------
  nmax = n;
  gt = new GradType[nmax];
  for (int i = 0; i < nmax; i++) Resize(gt[i],nmax);
}

GTvector::GTvector ( const GTvector& v )                   // Copy constructor
{                                                          //-----------------
  nmax = v.nmax;
  gt = new GradType[nmax];
  for (int i = 0; i < nmax; i++)
    { Resize(gt[i],nmax); gt[i] = v.gt[i];}
}

GTvector::~GTvector ( )                                          // Destructor
{                                                                //-----------
  nmax = 0;
  delete [] gt;
}

//----------------------------------------------------------------------------
// Functions for resizing and range check
//----------------------------------------------------------------------------
void Resize ( GradType& A, int n )
{
  A.nmax = n;
  Resize(A.g,0,n);
}

void TestSize ( const GradType& A, const GradType& b, const char *fname )
{
  if (A.nmax != b.nmax) {
    cout << "Parameters must be of same size in '" << fname
         << "'!" << endl;
    exit(-1);
  }
}

void TestSize ( const GTvector& v, const GTvector& w, const char *fname)
{
  if (v.nmax != w.nmax) {
    cout << "Parameters must be of same size in '" << fname
         << "'!" << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------------
// Operators for component access
//----------------------------------------------------------------------------
interval& GradType::operator[] ( int i )               // Private [] operator!
{                                                      //---------------------
  if ( (i < 0) || (i > nmax) ) {
    cout << "Index out of range in 'interval& GradType::operator[] ( int )'!"
         << endl;
    exit(-1);
  }
  return g[i];
}

GradType& GTvector::operator[] ( int ind )            // Public [] operator!
{                                                       //--------------------
  int n = ind;
  if ( (n < 1) || (n > nmax) ) {
    cout << "Index out of range in "
         << "'GradType& GTvector::operator[] ( index )'!" << endl;
    exit(-1);
  }
  return gt[n-1];
}

const GradType& GTvector::operator[] ( int ind ) const  // Public [] operator!
{                                                       //--------------------
  int n = ind;
  if ( (n < 1) || (n > nmax) ) {
    cout << "Index out of range in "
         << "'GradType& GTvector::operator[] ( index )'!" << endl;
    exit(-1);
  }
  return gt[n-1];
}

//----------------------------------------------------------------------------
// Assignment operators for 'GradType' and 'GTvector' variables with the
// right-hand side of type ...
//----------------------------------------------------------------------------
GradType& GradType::operator= ( const GradType& A )            // ... GradType
{                                                              //-------------
  TestSize(*this,A,"operator= ( GradType&, GradType& )");
  g = A.g;
  return *this;
}

GradType& GradType::operator= ( const interval& x )            // ... interval
{                                                              //-------------
  g = 0.0;
  g[0] = x;
  return *this;
}

GradType& GradType::operator= ( const real& r )                    // ... real
{                                                                  //---------
  *this = interval(r);
  return *this;
}

GTvector& GTvector::operator= ( const GTvector& v )            // ... GTvector
{                                                              //-------------
  // Check if left-hand side is already allocated!
  //----------------------------------------------
  TestSize(*this,v,"operator= ( GTvector&, GTvector& )");
  for (int i = 0; i < nmax; i++) gt[i] = v.gt[i];
  return *this;
}

//----------------------------------------------------------------------------
// Transfer functions for variables
//----------------------------------------------------------------------------
GTvector GradVar ( const ivector& x )                     // Generate variable
{                                                         //------------------
  int       d = -Lb(x)+1;
  int       ubd = Ub(x)+d;
  GTvector  hht(ubd);

  for (int i = 1; i <= ubd; i++) {
    hht[i][0] = x[i-d];     // GradVar[][0] = x;
    for (int k = 1; k <= ubd; k++)
      if (i == k) hht[i][k] = 1.0; else hht[i][k] = 0.0;
  }
  return hht;
}

GTvector GradVar ( const rvector& v )                     // Generate variable
{                                                         //------------------
  int     i, n = Lb(v), m = Ub(v);
  ivector u(n,m);

  for (i = n; i <= m; i++) u[i] = v[i];
  return GradVar(u);
}

//----------------------------------------------------------------------------
// Access functions for the function value and the gradient
//----------------------------------------------------------------------------
interval fValue ( const GradType& f )                    // Get function value
  { return f.g[0]; }                                     //-------------------

ivector gradValue ( const GradType& f )                  // Get gradient value
{                                                        //-------------------
  ivector hv(f.nmax);
  for (int i = 1; i <= f.nmax; i++) hv[i] = f.g[i];
  return hv;
}

//----------------------------------------------------------------------------
// Access functions for vector-valued functions (function value or Jacobian)
//----------------------------------------------------------------------------
ivector fValue ( const GTvector& f )                        // Get function value of
{                                                     // n-dimensional vector-
  ivector hv(f.nmax);                                 // valued function
                                                      //----------------------
  for (int i = 1; i <= f.nmax; i++)
    hv[i] = f.gt[i-1][0];
  return hv;
}

imatrix JacValue ( const GTvector& f )                      // Get Jacobian value of
{                                                     // n-dimensional scalar-
  imatrix hm(f.nmax,f.nmax);                          // valued function
                                                      //----------------------
  for (int i = 1; i <= f.nmax; i++)
    for (int j = 1; j <= f.nmax; j++)
      hm[i][j] = f.gt[i-1][j];
  return hm;
}

//----------------------------------------------------------------------------
// Unary operators + and - for GradType operands
//----------------------------------------------------------------------------
GradType operator+ ( GradType& u )
  { return u; }

GradType operator- ( const GradType& u )
{
  GradType umin(u.nmax);

  umin[0] = -u.g[0];
  if (GradOrder > 0)
    for (int i = 1; i <= u.nmax; i++)
      umin[i] = -u.g[i];
  return umin;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for two GradType operands
//----------------------------------------------------------------------------
GradType operator+ ( const GradType& u, const GradType& v)
{
  GradType add(u.nmax);

  TestSize(u,v,"operator+ ( GradType&, GradType& )");
  add.g[0] = u.g[0] + v.g[0];
  if (GradOrder > 0)
    for (int i = 1; i <= u.nmax; i++)
      add.g[i] = u.g[i] + v.g[i];
  return add;
}

GradType operator- ( const GradType& u, const GradType& v )
{
  GradType sub(u.nmax);

  TestSize(u,v,"operator- ( GradType&, GradType& )");
  sub.g[0] = u.g[0] - v.g[0];
  if (GradOrder > 0)
    for (int i=1; i <= u.nmax; i++)
      sub.g[i] = u.g[i] - v.g[i];
  return sub;
}

GradType operator* ( const GradType& u, const GradType& v )
{
  GradType mul(u.nmax);

  TestSize(u,v,"operator* ( GradType&, GradType& )");
  mul.g[0] = u.g[0] * v.g[0];
  if (GradOrder > 0)
    for (int i = 1; i <= u.nmax; i++)
      mul.g[i] = v.g[0]*u.g[i] + u.g[0]*v.g[i];
  return mul;
}

GradType operator/ ( const GradType& u, const GradType& v )
{
  GradType div(u.nmax);

  TestSize(u,v,"operator/ ( GradType&, GradType& )");
  div.g[0] = u.g[0]/v.g[0];    // Can propagate 'division by zero' error
  if (GradOrder > 0)           //---------------------------------------
    for (int i = 1; i <= u.nmax; i++)
      div.g[i] = (u.g[i] - div.g[0]*v.g[i]) / v.g[0];
  return div;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for one interval and one GradType operand
//----------------------------------------------------------------------------
GradType operator+ ( const GradType& u, const interval& b )
{
  GradType add(u.nmax);

  add = u; add.g[0] += b;
  return add;
}

GradType operator- ( const GradType& u, const interval& b )
{
  GradType sub(u.nmax);

  sub = u; sub.g[0] -= b;
  return sub;
}

GradType operator* ( const GradType& u, const interval& b )
{
  GradType mul(u.nmax);

  mul.g[0] = u.g[0]*b;
  if (GradOrder > 0)
    for (int i = 1; i <= u.nmax; i++)
      mul.g[i] = b*u.g[i];
  return mul;
}

GradType operator/ ( const GradType& u, const interval& b )
{
  GradType div(u.nmax);

  div.g[0] = u.g[0]/b;      // Can propagate 'division by zero' error
  if (GradOrder > 0)        //---------------------------------------
    for (int i = 1; i <= u.nmax; i++)
      div.g[i] = u.g[i]/b;
  return div;
}

GradType operator+ ( const interval& a, const GradType& v )
{
  GradType add(v.nmax);

  add = v; add.g[0] = a + v.g[0];
  return add;
}

GradType operator- ( const interval& a, const GradType& v )
{
  GradType sub(v.nmax);

  sub.g[0] = a - v.g[0];
  if (GradOrder > 0)
    for (int i = 1; i <= v.nmax; i++)
      sub.g[i] = -v.g[i];
  return sub;
}

GradType operator* ( const interval& a, const GradType& v )
{
  GradType mul(v.nmax);

  mul.g[0] = a*v.g[0];
  if (GradOrder > 0)
    for (int i = 1; i <= v.nmax; i++)
      mul.g[i] = a*v.g[i];
  return mul;
}

GradType operator/ ( const interval& a, const GradType& v )
{
  GradType div(v.nmax);
  interval p;

  div.g[0] = a/v.g[0];      // Can propagate 'division by zero' error
  if (GradOrder > 0) {      //---------------------------------------
    p = -div.g[0]/v.g[0];
    for (int i = 1; i <= v.nmax; i++)
      div.g[i] = p*v.g[i];
  }
  return div;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for one real and one GradType operand
//----------------------------------------------------------------------------
GradType operator+ ( const GradType& u, const real& b )
  { return u + interval(b); }

GradType operator- ( const GradType& u, const real& b )
  { return u - interval(b); }

GradType operator* ( const GradType& u, const real& b )
  { return u * interval(b); }

GradType operator/ ( const GradType& u, const real& b )
  { return u / interval(b); }    // Can propagate 'division by zero' error
                                  //---------------------------------------
GradType operator+ ( const real& a, const GradType& v )
  { return interval(a) + v; }

GradType operator- ( const real& a, const GradType& v )
  { return interval(a) - v; }

GradType operator* ( const real& a, const GradType& v )
  { return interval(a) * v; }

GradType operator/ ( const real& a, const GradType& v )
  { return interval(a) / v; }     // Can propagate 'division by zero' error
                                  //---------------------------------------

//----------------------------------------------------------------------------
// Elementary function for GradType arguments
//----------------------------------------------------------------------------
GradType sqr ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = sqr(u.g[0]);   // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h1 = 2.0*u.g[0];
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType power ( const GradType& u, const int k )
{
  GradType res(u.nmax);
  interval h1;

  if (k == 0)
    { res = 1.0; return res; }
  if (k == 1)
    return u;
  if (k == 2)
    return sqr(u);

  res.g[0] = Power(u.g[0],k);
  if (GradOrder > 0) {
    h1 = (double)k * Power(u.g[0],k-1);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1 * u.g[i];
  }
  return res;
}

GradType sqrt ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0, h1;

  h0 = sqrt(u.g[0]);
  res.g[0] = h0;
  if (GradOrder > 0) {
    h1 = 0.5/h0;
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType exp ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0;

  h0 = exp(u.g[0]);
  res.g[0] = h0;
  if (GradOrder > 0)
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h0*u.g[i];
  return res;
}

GradType ln ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = ln(u.g[0]);    // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h1 = 1.0/u.g[0];
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType sin ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res[0] = sin(u.g[0]);
  if (GradOrder > 0) {
    h1 = cos(u.g[0]);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType cos ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = cos(u.g[0]);
  if (GradOrder > 0) {
    h1 = -sin(u.g[0]);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType tan ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0, h1;

  h0 = tan(u.g[0]);            // Can propagate domain error
  res.g[0] = h0;                 //---------------------------
  if (GradOrder > 0) {
    h1 = sqr(h0)+1.0;                     // The subdistributive law implies
                                          //   h0 * (h0^2 + 1) <= h0^3 + h0.
    for (int i = 1; i <= u.nmax; i++)     // So, the first form is used.
      res.g[i] = h1*u.g[i];               //--------------------------------
  }
  return res;
}

GradType cot ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0, h1;

  h0 = cot(u.g[0]);              // Can propagate domain error
  res.g[0] = h0;                 //---------------------------
  if (GradOrder > 0) {
    h1 = -(sqr(h0)+1.0);                  // The subdistributive law implies
                                          //   h0 * (h0^2 + 1) <= h0^3 + h0.
    for (int i = 1; i <= u.nmax; i++)     // So, the first form is used.
      res.g[i] = h1*u.g[i];               //--------------------------------
  }
  return res;
}

GradType asin ( const GradType& u )
{
  GradType res(u.nmax);
  interval h, h1;

  res.g[0] = asin(u.g[0]);  // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h  = 1.0-sqr(u.g[0]);
    h1 = 1.0/sqrt(h);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType acos ( const GradType& u )
{
  GradType res(u.nmax);
  interval h, h1;

  res.g[0] = acos(u.g[0]);  // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h  = 1.0-sqr(u.g[0]);
    h1 = -1.0/sqrt(h);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType atan ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = atan(u.g[0]);
  if (GradOrder > 0) {
    h1 = 1.0/(1.0+sqr(u.g[0]));
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType acot ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = acot(u.g[0]);  // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h1 = -1.0/(1.0+sqr(u.g[0]));
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType sinh ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0, h1;

  h0 = sinh(u.g[0]);
  res.g[0] = h0;
  if (GradOrder > 0) {
    h1 = cosh(u.g[0]);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType cosh ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = cosh(u.g[0]);
  if (GradOrder > 0) {
    h1 = sinh(u.g[0]);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i]  = h1*u.g[i];
  }
  return res;
}

GradType tanh ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0, h1;

  h0 = tanh(u.g[0]);
  res.g[0] = h0;
  if (GradOrder > 0) {
    h1 = 1.0-sqr(h0);                     // The subdistributive law implies
                                          //   h0 * (h0^2 - 1) <= h0^3 - h0.
    for (int i = 1; i <= u.nmax; i++)     // So, the first form is used.
      res.g[i]  = h1*u.g[i];                  //--------------------------------
  }
  return res;
}

GradType coth ( const GradType& u )
{
  GradType res(u.nmax);
  interval h0, h1;

  h0 = coth(u.g[0]);             // Can propagate domain error
  res.g[0] = h0;                 //---------------------------
  if (GradOrder > 0) {
    h1 = 1.0-sqr(h0);                     // The subdistributive law implies
                                          //   h0 * (h0^2 - 1) <= h0^3 - h0.
    for (int i = 1; i <= u.nmax; i++)     // So, the first form is used.
      res.g[i] = h1*u.g[i];                   //--------------------------------
  }
  return res;
}

GradType asinh ( const GradType& u )
{
  GradType res(u.nmax);
  interval h, h1;

  res.g[0] = asinh(u.g[0]); // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h  = 1.0+sqr(u.g[0]);
    h1 = 1.0/sqrt(h);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType acosh ( const GradType& u )
{
  GradType res(u.nmax);
  interval h, h1;

  res.g[0] = acosh(u.g[0]); // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h  = sqr(u.g[0])-1.0;
    h1 = 1.0/sqrt(h);
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType atanh ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = atanh(u.g[0]); // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h1 = 1.0/(1.0-sqr(u.g[0]));
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

GradType acoth ( const GradType& u )
{
  GradType res(u.nmax);
  interval h1;

  res.g[0] = acoth(u.g[0]); // Can propagate domain error
  if (GradOrder > 0) {      //---------------------------
    h1 = 1.0/(1.0-sqr(u.g[0]));
    for (int i = 1; i <= u.nmax; i++)
      res.g[i] = h1*u.g[i];
  }
  return res;
}

//----------------------------------------------------------------------------
// Predefined routines for evaluation of GradType-functions
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of 'GradType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
// Description: This function sets 'GradOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
void fEvalG ( GTscalar_FctPtr f, ivector x, interval& fx )
{
  GradOrder = 0;
  fx        = fValue(f(GradVar(x)));
  GradOrder = 1;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the gradient value.
// Parameters:
//    In : 'f' : function of 'GradType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'gx': returns the gradient value 'grad f(x)'
// Description: This function keeps 'GradOrder' = 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the gradient.
//----------------------------------------------------------------------------
void fgEvalG ( GTscalar_FctPtr f, ivector x, interval& fx, ivector& gx )
{
  // GradType fxG(Ub(x));     // not implemented in HP C++ V2.0!
  int      n = Ub(x);
  GradType fxG(n);

  fxG       = f(GradVar(x));
  fx        = fValue(fxG);
  gx        = gradValue(fxG);
}

//----------------------------------------------------------------------------
// Predefined routines for evaluation of GTvector-functions (with Jacobians)
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of type 'GTvector'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the vector function value 'f(x)'
// Description: This function sets 'GradOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
void fEvalJ ( GTvector_FctPtr f, ivector x, ivector& fx )
{
  GradOrder = 0;
  fx        = fValue(f(GradVar(x)));
  GradOrder = 1;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the Jacobian matrix value.
// Parameters:
//    In : 'f' : function of type 'GTvector'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'Jx': returns the Jacobian value 'Jf(x)'
// Description: This function keeps 'GradOrder' = 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the Jacobian matrix.
//----------------------------------------------------------------------------
void fJEvalJ ( GTvector_FctPtr f, ivector x, ivector& fx, imatrix& Jx )
{
  // GTvector fxGTv(Ub(x));     // not implemented in HP C++ V2.0!
  int      n = Ub(x);
  GTvector fxGTv(n);

  fxGTv     = f(GradVar(x));
  fx        = fValue(fxGTv);
  Jx        = JacValue(fxGTv);
}






















