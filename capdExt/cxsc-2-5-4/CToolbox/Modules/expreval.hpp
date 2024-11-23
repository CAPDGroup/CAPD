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

/* CVS $Id: expreval.hpp,v 1.15 2014/01/30 17:49:26 cxsc Exp $ */

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
// File: expreval (header)
// Purpose: Computation of an enclosure of the value of a real arithmetic
//    expression composed of the operations +, -, *, /, and ^, where ^
//    denotes exponentiation by an integer.
// Type:
//    Stagg_FctPtr        : pointer for a function of type 'Staggered'
// Class Staggered (staggered data representation):
//    Staggered()         : constructors
//    operators =         : assignment of a real argument
//    operators +, -, *, /: both operands of type 'Staggered' or one of
//                          type 'Staggered' and one of type 'real'
//    Power()             : argument of type 'Staggered', exponent of type
//                          'integer' (Exponentiation by an integer)
//    Eval()              : main function for evaluation of an expression
//    EvalErrMsg()        : to get an error message text
// Class StaggArray (arrays of type 'Staggered'):
//    StaggArray()        : constructors
//    ~StaggArray()       : destructor
//    operator []         : component access
//----------------------------------------------------------------------------
#ifndef __EXPREVAL_HPP
#define __EXPREVAL_HPP

#include <interval.hpp>     // Interval arithmetic
#include <i_util.hpp>       // Interval utility functions
#include <rvector.hpp>      // Real vector arithmetic

using namespace cxsc;
using namespace std;

extern char* EvalErrMsg ( int );

class Staggered;
class StaggArray;

typedef Staggered  (*Stagg_FctPtr)(StaggArray&);

class Staggered {
  private:
    rvector   Val;
    interval  Err;

    friend void InitEntry       ( real );
    friend void UpdateError     ( interval );
    friend void UpdateStaggComp ( int );

  public:
    Staggered ( );
    Staggered ( const Staggered& );

    Staggered& operator= ( const real& );
    Staggered& operator= ( const Staggered& );

    friend Staggered  operator+ ( const Staggered& , const Staggered& );
    friend Staggered  operator- ( const Staggered& , const Staggered& );
    friend Staggered  operator* ( const Staggered& , const Staggered& );
    friend Staggered  operator/ ( const Staggered& , const Staggered& );

    friend Staggered  operator+ ( const Staggered& , const real& );
    friend Staggered  operator- ( const Staggered& , const real& );
    friend Staggered  operator* ( const Staggered& , const real& );
    friend Staggered  operator/ ( const Staggered& , const real& );

    friend Staggered  operator+ ( const real& , const Staggered& );
    friend Staggered  operator- ( const real& , const Staggered& );
    friend Staggered  operator* ( const real& , const Staggered& );
    friend Staggered  operator/ ( const real& , const Staggered& );

    friend Staggered Power ( const Staggered& , const int );

    friend void Eval ( Stagg_FctPtr, rvector, real, real&, interval&, int&, int& );

}; // class Staggered


void InitEntry       ( real );
void UpdateError     ( interval );
void UpdateStaggComp ( int );
Staggered  operator+ ( const Staggered& , const Staggered& );
Staggered  operator- ( const Staggered& , const Staggered& );
Staggered  operator* ( const Staggered& , const Staggered& );
Staggered  operator/ ( const Staggered& , const Staggered& );

Staggered  operator+ ( const Staggered& , const real& );
Staggered  operator- ( const Staggered& , const real& );
Staggered  operator* ( const Staggered& , const real& );
Staggered  operator/ ( const Staggered& , const real& );

Staggered  operator+ ( const real& , const Staggered& );
Staggered  operator- ( const real& , const Staggered& );
Staggered  operator* ( const real& , const Staggered& );
Staggered  operator/ ( const real& , const Staggered& );

Staggered Power ( const Staggered& , const int );

void Eval ( Stagg_FctPtr, rvector, real, real&, interval&, int&, int& );



class StaggArray {
  private:
    Staggered *SA;
    int       Dim;

  public:
    StaggArray  ( );
    StaggArray  ( int );
    ~StaggArray ( );
    StaggArray  ( const StaggArray& );

    Staggered& operator[] ( int );
}; // class StaggArray
#endif







