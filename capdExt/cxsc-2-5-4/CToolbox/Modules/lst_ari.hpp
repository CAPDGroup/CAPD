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

/* CVS $Id: lst_ari.hpp,v 1.17 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lst_ari (header)
// Purpose: Definition of a list arithmetic used in connection with an
//    interval bisection method in global optimization for storing pairs of
//    an interval vector and a real value.
// Class Pair:
//    Pair()        : constructors
//    ~Pair()       : destructor
//    operator =    : assignment operator
//    _Pair()       : for explicit generating a pair
//    GetInt()      : get the interval vector component
//    GetFyi()      : get the function value
// Global type:
//    PairPtr       : list of pairs
// Global constant:
//    EmptyList     : the empty list (NULL pointer)
// Global functions and operators:
//    operator +    : adding a new element to a list
//    Next(), Head(): access functions for lists
//    Length()      : access function to length of list
//    FreeAll()     : free complete list
//    MultiDelete() : deletes several elements in a list
//    DelHead()     : deletes the first element of a list
//    ListToMatrix(): transfer a list into a dynamic interval matrix
//----------------------------------------------------------------------------
#ifndef __LST_ARI_HPP
#define __LST_ARI_HPP

#include <intvector.hpp>     // Integer vector type
#include <imatrix.hpp>      // Interval matrix/vector arithmetic

using namespace cxsc;
using namespace std;

// #define EmptyList NULL      // The empty list
#define EmptyList 0         // The empty list
                            //---------------

struct                      // Structure used for a list of pairs
  PairElmt;                 //-----------------------------------

typedef                     // Pointer to a list of pairs
  PairElmt* PairPtr;        //---------------------------

class Pair {                // Pair of an interval vector 'intv'
  private:                  // and the real value 'inf(f(intv))'
    ivector intv;           //----------------------------------
    real    fyi;
  public:
    Pair ( ) { };
    Pair ( const Pair& );

    Pair& operator= ( const Pair& );

    friend Pair    _Pair  ( const ivector&, real );
    friend Pair    _Pair  ( const ivector_slice&, real );
    friend Pair    _Pair  ( const imatrix_subv&, real );
    friend ivector GetInt ( const Pair& );
    friend real    GetFyi ( const Pair& );
};

Pair    _Pair  ( const ivector&, real );
Pair    _Pair  ( const ivector_slice&, real );
Pair    _Pair  ( const imatrix_subv&, real );
ivector GetInt ( const Pair& );
real    GetFyi ( const Pair& );

extern void    MultiDelete  ( PairPtr&, const real& );
extern void    ListToMatrix ( PairPtr, imatrix&, intvector& );
extern PairPtr operator+    ( PairPtr, const Pair& );

extern void    FreeAll      ( PairPtr& );
extern int     Length       ( PairPtr );
extern Pair    Head         ( PairPtr );
extern void    DelHead      ( PairPtr& );
#endif




