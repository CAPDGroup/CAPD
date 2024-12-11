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

/* CVS $Id: lst1_ari.hpp,v 1.16 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lst1_ari (header)
// Purpose: Definition of a list arithmetic used in connection with an
//    interval bisection method in global optimization for storing pairs of
//    an interval and a real value.
// Class Pair1:
//    Pair1()       : constructors
//    ~Pair1()      : destructor
//    operator =    : assignment operator
//    _Pair1()      : for explicit generating a pair
//    GetInt()      : get the interval component
//    GetFyi()      : get the function value
// Global type:
//    Pair1Ptr      : list of pairs
// Global constant:
//    EmptyList     : the empty list (NULL pointer)
// Global functions and operators:
//    operator +    : adding a new element to a list
//    Head()        : access functions for the head of a list
//    Length()      : access function to length of list
//    FreeAll()     : free complete list
//    MultiDelete() : deletes several elements in a list
//    DelHead()     : deletes the first element of a list
//    ListToVector(): transfer a list into a dynamic interval vector
//----------------------------------------------------------------------------
#ifndef __LST1_ARI_HPP
#define __LST1_ARI_HPP

#include <intvector.hpp>     // Integer vector type
#include <ivector.hpp>      // Interval vector arithmetic

using namespace cxsc;
using namespace std;

// #define EmptyList NULL      // The empty list
#define EmptyList 0         // The empty list
                            //---------------

struct                      // Structure used for a list of pairs
  Pair1Elmt;                //-----------------------------------

typedef                     // Pointer to a list of pairs
  Pair1Elmt* Pair1Ptr;      //---------------------------

class Pair1 {               // Pair of an interval 'intv' and
  private:                  // the real value 'inf(f(intv))'.
    interval intv;          //-------------------------------
    real     fyi;
  public:
    Pair1 ( ) { };
    Pair1 ( const Pair1& );

    Pair1& operator= ( const Pair1& );

    friend Pair1    _Pair1 ( const interval&, const real& );
    friend interval GetInt ( const Pair1& );
    friend real     GetFyi ( const Pair1& );
};

Pair1    _Pair1 ( const interval&, const real& );
interval GetInt ( const Pair1& );
real     GetFyi ( const Pair1& );

extern void     MultiDelete  ( Pair1Ptr&, const real& );
extern void     ListToVector ( Pair1Ptr, ivector&, intvector& );
extern Pair1Ptr operator+    ( Pair1Ptr, Pair1 );

extern void     FreeAll      ( Pair1Ptr& );
extern int      Length       ( Pair1Ptr );
extern Pair1    Head         ( Pair1Ptr );
extern void     DelHead      ( Pair1Ptr& );
#endif




