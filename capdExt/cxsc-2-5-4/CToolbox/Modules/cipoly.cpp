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

/* CVS $Id: cipoly.cpp,v 1.15 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: cipoly (implementation)
// Purpose: Definition of the class for the representation of a complex
//    interval polynomial by its coefficients.
// Class CIPolynomial:
//    Deg()         : to get the degree of the polynomial
//    CIPolynomial(): constructors
//    operator []   : component access
//    in()          : Contained-in-the-interior relation for two
//                    complex interval polynomials
//    Blow()        : Epsilon inflation
//    operator >>   : input operator for data of type CIPolynomial
//    operator <<   : output operator for data of type CIPolynomial
//----------------------------------------------------------------------------
#include <ci_util.hpp>     // Complex interval utility functions
#include <cipoly.hpp>

using namespace cxsc;
using namespace std;

int Deg ( const CIPolynomial& p ) { return Ub(p.coeff); }

CIPolynomial::CIPolynomial( int n )
{
  Resize(coeff,0,n);
  coeff = 0.0;
}

CIPolynomial::CIPolynomial( const CIPolynomial& p )
{
  Resize(coeff,0,Deg(p));
  coeff = p.coeff;
}

int in ( const CIPolynomial& p, const CIPolynomial& q )               // Contained-in-the-
{                                                         // interior relation
  int i, incl = 1;                                        //------------------

  for (i = 0; (incl == 1) && (i <= Deg(p)); i++)
    incl = in(p[i],q[i]);
  return incl;
}

CIPolynomial Blow ( CIPolynomial p, const real eps )      // Epsilon inflation
{                                                         //------------------
  for (int i = 0; i <= Deg(p); i++)
    p[i] = Blow(p[i],eps);
  return p;
}

istream& operator>> ( istream& in, CIPolynomial& p )
{
  cout << "  x^0 * ";  in >> p[0];
  for (int i = 1; i <= Deg(p); i++)
    { cout << "+ x^" << i << " * ";  in >> p[i]; }
  return in;
}

ostream& operator<< ( ostream& out, CIPolynomial p )
{
  int PolyIsZero, n = Deg(p);

  PolyIsZero = 1;
  for (int i = 0; i <= n; i++) {
    if (p[i] == 0.0) continue;
    if (PolyIsZero)
      out << "  ";
    else
      out << "+ ";
    out << p[i] << " * x^" << i << endl;
    PolyIsZero = 0;
  }
  if (PolyIsZero) out << "  0 (= zero polynomial)" << endl;
  return out;
}




