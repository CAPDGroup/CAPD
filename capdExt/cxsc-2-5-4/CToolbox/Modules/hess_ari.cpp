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

/* CVS $Id: hess_ari.cpp,v 1.19 2014/01/30 18:13:37 cxsc Exp $ */

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
// File: hess_ari (implementation)
// Purpose: Definition of a multi-dimensional interval differentiation
//    arithmetic which allows function evaluation with automatic differen-
//    tiation up to second order (i.e. Hessian matrix).
// Method: Overloading of operators and elementary functions for operations
//    of data types 'HessType' and 'HTvector'. To represent a Hessian matrix
//    by a lower triangular matrix, the type 'LowTriMatrix' is supplied.
// Types:
//    HTscalar_FctPtr         : pointer for a scalar valued function
//    HTvector_FctPtr         : pointer for a vector valued function
// Class LowTriMatrix:
//    LowTriMatrix()          : constructors
//    ~LowTriMatrix()         : destructor
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for triangular matrices
//                              and real values
//    operator []             : component access
// Class HessType:
//    HessType()              : constructors
//    Resize()                : for resizing to a fixed dimension
//    operator =              : assignment operators for arguments of types
//                              HessType, interval, and real
//    Dim()                   : returns the dimension of the Hessian matrix
//    HessVar()               : to define HessType variables
//    fValue(), gradValue(),
//    hessValue()             : to get function and derivative values
//    operators +, -, *, /    : operators of diff. arithmetic
//    sqr(), sqrt(), power(),
//    exp(), sin(), cos(), ...: elementary functions of diff. arithmetic
//    fEvalH()                : to compute function value only
//    fgEvalH()               : to compute function and first derivative
//                              value (gradient)
//    fghEvalH()              : to compute function value, gradient, and
//                              Hessian matrix value
// Class HTvector:
//    HTvector()              : constructors
//    ~HTvector()             : destructor
//    operator =              : assignment operator
//    operator []             : component access
//    Dim()                   : to get the actual dimension
//    fValue(), JacValue()    : to get function and derivative values
//    fEvalJ()                : to compute function value only
//    fJEvalJ()               : to compute function and first derivative
//                              value (Jacobian)
//----------------------------------------------------------------------------
#include <cstdlib>          // For function 'exit()'
#include <iostream>         // I/O Handling
#include <imath.hpp>        // Interval mathematical functions
#include <i_util.hpp>       // Interval utility functions
#include <hess_ari.hpp>

using namespace cxsc;
using namespace std;

// The local variable 'HessOrder' is used to select the highest order of
// derivative which is computed. Its default value is 2, and normally
// the gradient and the Hessian matrix are computed.
//----------------------------------------------------------------------------
#ifdef _WIN32
static __declspec(thread) int HessOrder = 2;
#elif __APPLE__ && !CXSC_FORCE_TLS
static int HessOrder = 2;
#else
static __thread int HessOrder = 2;
#endif

//----------------------------------------------------------------------------
// Constructors, destructor, and resize function for the class 'LowTriMatrix'.
//----------------------------------------------------------------------------
void LowTriMatrix::init ( int n )         // Initializes all components of the
{                                         // lower triangular matrix with 0.0
  if (n <= 0)                             //----------------------------------
    { dim = 0; rows = NULL; return; }
  dim = n;
  rows = new ivector[n];
  for (int i = 0; i < n; i++)                  // i-th row of the matrix takes
    { Resize(rows[i],i+1); rows[i] = 0.0; }    // i components, i = 1,2,...,n
}                                              //-----------------------------

LowTriMatrix::LowTriMatrix ( )                      // Default: no matrix rows
  { dim = 0; rows = NULL; }                         //------------------------


LowTriMatrix::LowTriMatrix ( int n )             // Constructor: dimension n
  { init(n); }                              //-------------------------

LowTriMatrix::~LowTriMatrix ( )                                  // Destructor
  { delete [] rows; }                                            //-----------

void Resize ( LowTriMatrix& a, int n )           // Resizes matrix to n rows
  { delete [] a.rows; a.init(n); }          //-------------------------

//----------------------------------------------------------------------------
// Assignment operators and copy constructor of the class 'LowTriMatrix'.
//----------------------------------------------------------------------------
LowTriMatrix& LowTriMatrix::operator= ( const LowTriMatrix& a )   // matrix = matrix
{                                                           //----------------
  if (this == &a) return *this;
  delete [] rows;
  dim = a.dim;
  if (dim == 0) { rows = NULL; return *this; }
  rows = new ivector[dim];
  for (int i = 0; i < dim; i++) rows[i] = a.rows[i];
  return *this;
}

LowTriMatrix& LowTriMatrix::operator= ( const real& r )    // matrix = real
{                                                          //--------------
  for (int i = 0; i < dim; i++) rows[i] = r;
  return *this;
}

LowTriMatrix::LowTriMatrix ( const LowTriMatrix& a )       // Copy constructor
{                                                          //-----------------
  rows = NULL;
  *this = a;   // use previously defined assignment operator
}

//----------------------------------------------------------------------------
// Component access for the type 'LowTriMatrix'
//----------------------------------------------------------------------------
ivector& LowTriMatrix::operator[] ( int i )
{ 
  if ( (i < 1) || (i > dim) ) {
    cout << "Index out of range in "
         << "'ivector& LowTriMatrix::operator[] ( int )'!" << endl;
    exit(-1);
  }
  return rows[i-1];      // Index range starts at 1
}

const ivector& LowTriMatrix::operator[] ( int i ) const
{ 
  if ( (i < 1) || (i > dim) ) {
    cout << "Index out of range in "
         << "'ivector& LowTriMatrix::operator[] ( int )'!" << endl;
    exit(-1);
  }
  return rows[i-1];      // Index range starts at 1
}
////----------------------------------------------------------------------------
// Constructors, destructor, and resize function for the classes 'HessType'
// and 'HTvector'.
//----------------------------------------------------------------------------
HessType::HessType ( )                               // Default: no components
  { nmax = 0; }                                      //-----------------------

HessType::HessType ( int n )                // Constructor: resizes internal
{                                             // arrays to dimension 'nmax'
  nmax = n;                            //------------------------------
  if (nmax <= 0) { nmax = 0; return; }
  if (HessOrder > 0) Resize(g,nmax);
  if (HessOrder > 1) Resize(h,nmax);
}

void Resize ( HessType& u, int n )
{
  int dim = n;
  u.nmax = dim;
  if (HessOrder > 0) Resize(u.g,dim);
  if (HessOrder > 1) Resize(u.h,dim);
}

HTvector::HTvector ( int n )                   // Constructor: allocates 'n'
{                                                // 'HessType' components of
  nmax = n;                               // dimension 'n'
  if (nmax <= 0)                                 //---------------------------
    { nmax = 0; ht = NULL; return; }
  ht = new HessType[nmax];
  for (int i = 0; i < nmax; i++) Resize(ht[i],nmax);
}

HTvector::~HTvector ( )                                          // Destructor
  { delete [] ht; }                                              //-----------

//----------------------------------------------------------------------------
// Assignment operators and copy constructors of the classes 'HessType' and
// 'HTvector'.
//----------------------------------------------------------------------------
HessType& HessType::operator= ( const HessType& u )     // HessType = HessType
{                                                       //--------------------
  if (this == &u) return *this;
  nmax = u.nmax;
  if (nmax > 0) { 
    f = u.f;  
    if (HessOrder > 0) g = u.g;  
    if (HessOrder > 1) h = u.h; 
  }
  return *this;
}

HessType::HessType ( const HessType& u )                   // Copy constructor
{                                                          //-----------------
  nmax = u.nmax;
  *this = u;   // use previously defined assignment operator
}

HessType& HessType::operator= ( const interval& x )    // HessType =  interval
{                                                      //---------------------
  f = x;  
  if (HessOrder > 0) g = 0.0;  
  if (HessOrder > 1) h = 0.0;
  return *this;
}

HessType& HessType::operator= ( const real& r )            // HessType =  real
{                                                          //-----------------
  *this = interval(r);
  return *this;
}

HTvector& HTvector::operator= ( const HTvector& v )     // HTvector = HTvector
{                                                       //--------------------
  if (this == &v) return *this;
  delete [] ht;
  nmax = v.nmax;
  if (nmax == 0) { ht = NULL; return *this; }
  ht = new HessType[nmax];
  for (int i = 0; i < nmax; i++) ht[i] = v.ht[i];
  return *this;
}

HTvector::HTvector ( const HTvector& v )                   // Copy constructor
{                                                          //-----------------
  ht = NULL;
  *this = v;   // use previously defined assignment operator
}

//----------------------------------------------------------------------------
// Component access for the type 'HTvector'.
//----------------------------------------------------------------------------
HessType& HTvector::operator[] ( int ind )
{
  int n = ind;
  if ( (n < 1) || (n > nmax) ) {
    cout << "Index out of range in "
         << "'HessType& HTvector::operator[] ( index )'!" << endl;
    exit(-1);
  }
  return ht[n-1];        // Index range starts at 1
}

const HessType& HTvector::operator[] ( int ind ) const
{
  int n = ind;
  if ( (n < 1) || (n > nmax) ) {
    cout << "Index out of range in "
         << "'HessType& HTvector::operator[] ( index )'!" << endl;
    exit(-1);
  }
  return ht[n-1];        // Index range starts at 1
}

//----------------------------------------------------------------------------
// Functions for range check of 'HessType' and 'HTvector' objects.
//----------------------------------------------------------------------------
void TestSize ( const HessType& a, const HessType& b, const char *fname )
{
  if (a.nmax != b.nmax) {
    cout << "Parameters must be of same size in '" << fname
         << "'!" << endl;
    exit(-1);
  }
}

void TestSize ( const HTvector& v, const HTvector& w, const char *fname)
{
  if (v.nmax != w.nmax) {
    cout << "Parameters must be of same size in '" << fname
         << "'!" << endl;
    exit(-1);
  }
}

//----------------------------------------------------------------------------
// Transfer functions for variables
//----------------------------------------------------------------------------
HTvector HessVar ( const ivector& x )                     // Generate variable
{                                                         //------------------
  int       d = -Lb(x)+1;
  int       ubd = Ub(x)+d;
  HTvector  hht(ubd);

  for (int i = 1; i <= ubd; i++) {
    hht[i].f = x[i-d];
    if (HessOrder > 0)
      for (int k = 1; k <= ubd; k++)
        hht[i].g[k] = (i == k) ? 1.0 : 0.0;
    if (HessOrder > 1) hht[i].h = 0.0;
  }
  return hht;
}

HTvector HessVar ( const rvector& v )                     // Generate variable
{                                                         //------------------
  int     i, n = Lb(v), m = Ub(v);
  ivector u(n,m);

  for (i = n; i <= m; i++) u[i] = v[i];
  return HessVar(u);
}

//----------------------------------------------------------------------------
// Access functions for the function value, the gradient, or the Hessian
//----------------------------------------------------------------------------
interval fValue ( const HessType& u )                    // Get function value
  { return u.f; }                                        //-------------------

ivector gradValue ( const HessType& u )                  // Get gradient value
  { return u.g; }                                        //-------------------

imatrix hessValue ( const HessType& u )                  // Get Hessian value
{                                                        //------------------
  imatrix hm(u.nmax,u.nmax);

  for (int i = 1; i <= u.nmax; i++) {
    hm[i][i] = u.h[i][i];
    for (int j=1; j <= i-1; j++) {
      hm[i][j] = u.h[i][j];
      hm[j][i] = u.h[i][j];
    }
  }
  return hm;
}

//----------------------------------------------------------------------------
// Access functions for vector-valued functions (function value or Jacobian)
//----------------------------------------------------------------------------
ivector fValue ( const HTvector& u )                  // Get function value of
{                                                     // n-dimensional vector-
  ivector hv(u.nmax);                                 // valued function
                                                      //----------------------
  for (int i = 1; i <= u.nmax; i++) hv[i] = fValue(u[i]);
  return hv;
}

imatrix JacValue ( const HTvector& u )                // Get Jacobian value of
{                                                     // n-dimensional scalar-
  imatrix hm(u.nmax,u.nmax);                          // valued function
                                                      //----------------------
  for (int i = 1; i <= u.nmax; i++) hm[i] = gradValue(u[i]);
  return hm;
}

//----------------------------------------------------------------------------
// Unary operators + and - for HessType operands
//----------------------------------------------------------------------------
HessType operator+ ( HessType& u )
  { return u; }

HessType operator- ( const HessType& u )
{
  HessType umin(u.nmax);

  umin.f = -u.f;
  if (HessOrder > 0)
    for (int i = 1; i <= u.nmax; i++) {
      umin.g[i] = -u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) (umin.h)[i][j] = -(u.h)[i][j];
    }
  return umin;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for two HessType operands
//----------------------------------------------------------------------------
HessType operator+ ( const HessType& u, const HessType& v )
{
  HessType add(u.nmax);

  TestSize(u,v,"operator+ ( HessType&, HessType& )");
  add.f = u.f + v.f;
  if (HessOrder > 0)
    for (int i = 1; i <= u.nmax; i++) {
      add.g[i] = u.g[i] + v.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) add.h[i][j] = u.h[i][j] + v.h[i][j];
    }
  return add;
}

HessType operator- ( const HessType& u, const HessType& v )
{
  HessType sub(u.nmax);

  TestSize(u,v,"operator- ( HessType&, HessType& )");
  sub.f = u.f - v.f;
  if (HessOrder > 0)
    for (int i=1; i <= u.nmax; i++) {
      sub.g[i] = u.g[i] - v.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) sub.h[i][j] = u.h[i][j] - v.h[i][j];
    }
  return sub;
}

HessType operator* ( const HessType& u, const HessType& v )
{
  HessType mul(u.nmax);

  TestSize(u,v,"operator* ( HessType&, HessType& )");
  mul.f = u.f * v.f;
  if (HessOrder > 0)
    for (int i = 1; i <= u.nmax; i++) {
      mul.g[i] = v.f*u.g[i] + u.f*v.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          mul.h[i][j] = v.f*u.h[i][j] + u.g[i]*v.g[j] + v.g[i]*u.g[j]
                                                      + u.f*v.h[i][j];
    }
  return mul;
}

HessType operator/ ( const HessType& u,const HessType& v )
{
  HessType div(u.nmax);

  TestSize(u,v,"operator/ ( HessType&, HessType& )");
  div.f = u.f / v.f;              // Can propagate 'division by zero' error
  if (HessOrder > 0)              //---------------------------------------
    for (int i = 1; i <= u.nmax; i++) {
      div.g[i] = (u.g[i] - div.f*v.g[i]) / v.f;
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          div.h[i][j] = (u.h[i][j] - div.g[i]*v.g[j] - v.g[i]*div.g[j]
                                                     - div.f*v.h[i][j]) / v.f;
    }
  return div;
}

//----------------------------------------------------------------------------
// Operators +, -, *, and / for one interval and one HessType operand
//----------------------------------------------------------------------------
HessType operator+ ( const HessType& u, const interval& b )
{
  HessType add(u.nmax);

  add = u; add.f += b;
  return add;
}

HessType operator- ( const HessType& u, const interval& b )
{
  HessType sub(u.nmax);

  sub = u; sub.f -= b;
  return sub;
}

HessType operator* ( const HessType& u, const interval& b )
{
  HessType mul(u.nmax);

  mul.f = u.f * b;
  if (HessOrder > 0)
    for (int i = 1; i <= u.nmax; i++) {
      mul.g[i] = b * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) mul.h[i][j] = b * u.h[i][j];
    }
  return mul;
}

HessType operator/ ( const HessType& u, const interval& b )
{
  HessType div(u.nmax);

  div.f = u.f / b;                // Can propagate 'division by zero' error
  if (HessOrder > 0)              //---------------------------------------
    for (int i = 1; i <= u.nmax; i++) {
      div.g[i] = u.g[i] / b;
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) div.h[i][j] = u.h[i][j] / b;
    }
  return div;
}

HessType operator+ ( const interval& a, const HessType& v )
{
  HessType add(v.nmax);

  add = v; add.f = a + v.f;
  return add;
}

HessType operator- ( const interval& a, const HessType& v )
{
  HessType sub(v.nmax);

  sub.f = a - v.f;
  if (HessOrder > 0)
    for (int i = 1; i <= v.nmax; i++) {
      sub.g[i] = -v.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) sub.h[i][j] = -v.h[i][j];
    }
  return sub;
}

HessType operator* ( const interval& a, const HessType& v )
{
  HessType mul(v.nmax);

  mul.f = a * v.f;
  if (HessOrder > 0)
    for (int i = 1; i <= v.nmax; i++) {
      mul.g[i] = a * v.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++) mul.h[i][j] = a * v.h[i][j];
    }
  return mul;
}

HessType operator/ ( const interval& a, const HessType& v )
{
  HessType div(v.nmax);
  interval p, q;

  div.f = a / v.f;                // Can propagate 'division by zero' error
  if (HessOrder > 0) {            //---------------------------------------
    p = -div.f / v.f;
    q = (-2.0*p) / v.f;
    for (int i = 1; i <= v.nmax; i++) {
      div.g[i] = p * v.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          div.h[i][j] = p*v.h[i][j] + q*v.g[i]*v.g[j];
    }
  }
  return div;
}

//----------------------------------------------------------------------------
// Operators +, -, * and / for one real and one HessType operand
//----------------------------------------------------------------------------
HessType operator+ ( const HessType& u, const real& b )
  { return u + interval(b); }

HessType operator- ( const HessType& u, const real& b )
  { return u - interval(b); }

HessType operator* ( const HessType& u, const real& b )
  { return u * interval(b); }

HessType operator/ ( const HessType& u, const real& b )
  { return u / interval(b); }    // Can propagate 'division by zero' error
                                  //---------------------------------------
HessType operator+ ( const real& a, const HessType& v )
  { return interval(a) + v; }

HessType operator- ( const real& a, const HessType& v )
  { return interval(a) - v; }

HessType operator* ( const real& a, const HessType& v )
  { return interval(a) * v; }

HessType operator/ ( const real& a, const HessType& v )
  { return interval(a) / v; }    // Can propagate 'division by zero' error
                                  //---------------------------------------

//----------------------------------------------------------------------------
// Elementary function for HessType arguments
//----------------------------------------------------------------------------
HessType sqr ( const HessType& u )
{
  HessType res(u.nmax);
  interval h1;

  res.f = sqr(u.f);               // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 2.0 * u.f;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + 2.0*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType power ( const HessType& u, const int k )
{
  HessType res(u.nmax);
  interval h1;

  if (k == 0)
    { res = 1.0; return res; }
  if (k == 1)
    return u;
  if (k == 2)
    return sqr(u);

  res.f = Power(u.f,k);
  if (HessOrder > 0) {
    h1 = double(k) * Power(u.f,k-1);
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] +
                        double(k*(k-1))*Power(u.f,k-2)*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType sqrt ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = sqrt(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = 0.5 / h0;
    h2 = -0.5*h1 / u.f;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType exp ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0;

  h0 = exp(u.f);
  res.f = h0;
  if (HessOrder > 0)
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h0 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h0 * (u.h[i][j] + u.g[i]*u.g[j]);
    }
  return res;
}

HessType ln ( const HessType& u )
{
  HessType res(u.nmax);
  interval h1,h2;

  res.f = ln(u.f);                // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 1.0/u.f;
    h2 = -sqr(h1);
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType sin ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = sin(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = cos(u.f);
    h2 = -h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType cos ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = cos(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = -sin(u.f);
    h2 = -h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType tan ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = tan(u.f);              // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (HessOrder > 0) {
    h1 = sqr(h0) + 1.0;                   // The subdistributive law implies
    h2 = 2.0*h0*h1;                       //   h0 * (h0^2 + 1) <= h0^3 + h0.
    for (int i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g[i] = h1 * u.g[i];             //--------------------------------
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType cot ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = cot(u.f);              // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (HessOrder > 0) {
    h1 = -(sqr(h0) + 1.0);                // The subdistributive law implies
    h2 = -2.0*h0*h1;                      //   h0 * (h0^2 + 1) <= h0^3 + h0.
    for (int i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g[i] = h1 * u.g[i];             //--------------------------------
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType asin ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  res.f = asin(u.f);              // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = 1.0-sqr(u.f);
    h1 = 1.0/sqrt(h0);
    h2 = u.f*h1/h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType acos ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  res.f = acos(u.f);              // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = 1.0-sqr(u.f);
    h1 = -1.0/sqrt(h0);
    h2 = u.f*h1/h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType atan ( const HessType& u )
{
  HessType res(u.nmax);
  interval h1, h2;

  res.f = atan(u.f);
  if (HessOrder > 0) {
    h1 = 1.0/(1.0+sqr(u.f));
    h2 = -2.0*u.f*sqr(h1);
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType acot ( const HessType& u )
{
  HessType res(u.nmax);
  interval h1, h2;

  res.f = acot(u.f);              // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = -1.0/(1.0+sqr(u.f));
    h2 = 2.0*u.f*sqr(h1);
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType sinh ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = sinh(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = cosh(u.f);
    h2 = h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType cosh ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = cosh(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = sinh(u.f);
    h2 = h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i]  = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType tanh ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = tanh(u.f);
  res.f = h0;
  if (HessOrder > 0) {
    h1 = 1.0-sqr(h0);                     // The subdistributive law implies
    h2 = -2.0*h0*h1;                      //   h0 * (h0^2 - 1) <= h0^3 - h0.
    for (int i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g[i]  = h1 * u.g[i];            //--------------------------------
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType coth ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  h0 = coth(u.f);             // Can propagate domain error
  res.f = h0;                 //---------------------------
  if (HessOrder > 0) {
    h1 = 1.0-sqr(h0);                     // The subdistributive law implies
    h2 = -2.0*h0*h1;                      //   h0 * (h0^2 - 1) <= h0^3 - h0.
    for (int i = 1; i <= u.nmax; i++) {   // So, the first form is used.
      res.g[i] = h1 * u.g[i];             //--------------------------------
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType asinh ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  res.f = asinh(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = 1.0+sqr(u.f);
    h1 = 1.0/sqrt(h0);
    h2 = -u.f*h1/h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType acosh ( const HessType& u )
{
  HessType res(u.nmax);
  interval h0, h1, h2;

  res.f = acosh(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h0 = sqr(u.f)-1.0;
    h1 = 1.0/sqrt(h0);
    h2 = -u.f*h1/h0;
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType atanh ( const HessType& u )
{
  HessType res(u.nmax);
  interval h1, h2;

  res.f = atanh(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 1.0/(1.0-sqr(u.f));
    h2 = 2.0*u.f*sqr(h1);
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

HessType acoth ( const HessType& u )
{
  HessType res(u.nmax);
  interval h1,h2;

  res.f = acoth(u.f);             // Can propagate domain error
  if (HessOrder > 0) {            //---------------------------
    h1 = 1.0/(1.0-sqr(u.f));
    h2 = 2.0*u.f*sqr(h1);
    for (int i = 1; i <= u.nmax; i++) {
      res.g[i] = h1 * u.g[i];
      if (HessOrder > 1)
        for (int j = 1; j <= i; j++)
          res.h[i][j] = h1*u.h[i][j] + h2*u.g[i]*u.g[j];
    }
  }
  return res;
}

//----------------------------------------------------------------------------
// Predefined routines for evaluation of HessType-functions
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of 'HessType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
// Description: This function sets 'HessOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
void fEvalH ( HTscalar_FctPtr f, ivector x, interval& fx )
{
  HessOrder = 0;
  fx        = fValue(f(HessVar(x)));
  HessOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the gradient value.
// Parameters:
//    In : 'f' : function of 'HessType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'gx': returns the gradient value 'grad f(x)'
// Description: This function sets 'HessOrder' to 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the gradient.
//----------------------------------------------------------------------------
void fgEvalH ( HTscalar_FctPtr f, ivector x, interval& fx, ivector& gx )
{
  // HessType fxH(Ub(x));     // not implemented in HP C++ V2.0!
  int      n = Ub(x);
  HessType fxH(n);

  HessOrder = 1;
  fxH       = f(HessVar(x));
  fx        = fValue(fxH);
  gx        = gradValue(fxH);
  HessOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value, the gradient value, and the
//    Hessian matrix value.
// Parameters:
//    In : 'f' : function of 'HessType'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'gx': returns the gradient value 'grad f(x)'
//         'hx': returns the Hessian matrix value 'hess f(x)'
// Description: This function keeps 'HessOrder' = 2, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value, the
//    value of the gradient, and the value of the Hessian matrix.
//----------------------------------------------------------------------------
void fghEvalH ( HTscalar_FctPtr f, ivector x, interval& fx, ivector& gx,
                                                            imatrix& hx )
{
  // HessType fxH(Ub(x));     // not implemented in HP C++ V2.0!
  int      n = Ub(x);
  HessType fxH(n);

  fxH = f(HessVar(x));
  fx  = fValue(fxH);
  gx  = gradValue(fxH);
  hx  = hessValue(fxH);
}

//----------------------------------------------------------------------------
// Predefined routines for evaluation of HTvector-functions (with Jacobians)
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing only the function value.
// Parameters:
//    In : 'f' : function of type 'HTvector'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the vector function value 'f(x)'
// Description: This function sets 'HessOrder' to 0, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value only.
//----------------------------------------------------------------------------
void fEvalJ ( HTvector_FctPtr f, ivector x, ivector& fx )
{
  HessOrder = 0;
  fx        = fValue(f(HessVar(x)));
  HessOrder = 2;
}

//----------------------------------------------------------------------------
// Purpose: Evaluation of function 'f' for argument 'x' in differentiation
//    arithmetic computing the function value and the Jacobian matrix value.
// Parameters:
//    In : 'f' : function of type 'HTvector'
//         'x' : argument for evaluation of 'f'
//    Out: 'fx': returns the function value 'f(x)'
//         'Jx': returns the Jacobian value 'Jf(x)'
// Description: This function sets 'HessOrder' to 1, evaluates 'f(x)' in
//    differentiation arithmetic, and returns the function value and the
//    value of the Jacobian matrix.
//----------------------------------------------------------------------------
void fJEvalJ ( HTvector_FctPtr f, ivector x, ivector& fx, imatrix& Jx )
{
  // HTvector fxGTv(Ub(x));     // not implemented in HP C++ V2.0!
  int      n = Ub(x);
  HTvector fxGTv(n);

  HessOrder = 1;
  fxGTv     = f(HessVar(x));
  fx        = fValue(fxGTv);
  Jx        = JacValue(fxGTv);
  HessOrder = 2;
}




