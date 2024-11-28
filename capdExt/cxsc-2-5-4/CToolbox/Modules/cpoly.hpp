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

/* CVS $Id: cpoly.hpp,v 1.15 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: cpoly (header)
// Purpose: Definition of the class for the representation of a complex
//    polynomial by its coefficients.
// Class CPolynomial:
//    Deg()        : to get the degree of the polynomial
//    CPolynomial(): constructors
//    operator []  : component access
//    operator >>  : input operator for data of type CPolynomial
//    operator <<  : output operator for data of type CPolynomial
//----------------------------------------------------------------------------
#ifndef __CPOLY_HPP
#define __CPOLY_HPP

#include <cvector.hpp>     // Complex vector arithmetic

using namespace cxsc;
using namespace std;

class CPolynomial {
  private:
    cvector coeff;
  public:
    CPolynomial ( int );
    CPolynomial ( const CPolynomial& );
    complex& operator[] ( int i ) { return coeff[i]; }

    friend int Deg ( const CPolynomial& );
    friend istream& operator>> ( istream&, CPolynomial& );
    friend ostream& operator<< ( ostream&, CPolynomial );
};

int Deg ( const CPolynomial& );
istream& operator>> ( istream&, CPolynomial& );
ostream& operator<< ( ostream&, CPolynomial );
#endif




