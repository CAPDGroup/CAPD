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

/* CVS $Id: cipoly.hpp,v 1.16 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: cipoly (header)
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
#ifndef __CIPOLY_HPP
#define __CIPOLY_HPP

#include <civector.hpp>     // Complex interval vector arithmetic

using namespace cxsc;
using namespace std;

class CIPolynomial {
  private:
    civector coeff;
  public:
    CIPolynomial ( int );
    CIPolynomial ( const CIPolynomial& );
    cinterval& operator[] ( int i ) { return coeff[i]; }
    const cinterval& operator[] ( int i ) const { return coeff[i]; }

    friend int Deg ( const CIPolynomial& );
    friend int in ( const CIPolynomial&, const CIPolynomial& );
    friend CIPolynomial Blow ( CIPolynomial, const real );
    friend istream& operator>> ( istream&, CIPolynomial& );
    friend ostream& operator<< ( ostream&, CIPolynomial );
};

int Deg ( const CIPolynomial& );
int in ( const CIPolynomial&, const CIPolynomial& );
CIPolynomial Blow ( CIPolynomial, const real );
istream& operator>> ( istream&, CIPolynomial& );
ostream& operator<< ( ostream&, CIPolynomial );
#endif




