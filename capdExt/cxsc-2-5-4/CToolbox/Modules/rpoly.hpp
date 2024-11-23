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

/* CVS $Id: rpoly.hpp,v 1.15 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: rpoly (header)
// Purpose: Definition of the class for the representation of a real
//    polynomial by its coefficients.
// Class RPolynomial:
//    Deg()        : to get the degree of the polynomial
//    RPolynomial(): constructors
//    operator []  : component access
//    operator >>  : input operator for data of type RPolynomial
//    operator <<  : output operator for data of type RPolynomial
//----------------------------------------------------------------------------
#ifndef __RPOLY_HPP
#define __RPOLY_HPP

#include <rvector.hpp>     // Real vector arithmetic

using namespace cxsc;
using namespace std;

class RPolynomial {
  private:
    rvector coeff;
  public:
    RPolynomial ( int );
    RPolynomial ( const RPolynomial& );
    real& operator[] ( int i ) { return coeff[i]; }

    friend int Deg ( const RPolynomial& );
    friend istream& operator>> ( istream&, RPolynomial& );
    friend ostream& operator<< ( ostream&, RPolynomial );
};

int Deg ( const RPolynomial& );
istream& operator>> ( istream&, RPolynomial& );
ostream& operator<< ( ostream&, RPolynomial );
#endif




