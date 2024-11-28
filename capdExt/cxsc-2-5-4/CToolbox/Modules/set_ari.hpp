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

/* CVS $Id: set_ari.hpp,v 1.16 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: set_ari (header)
// Purpose: Implementation of basic operations for sets of integers (also
//    called indices).
// Class IndexSet:
//    IndexSet()   : constructors
//    ~IndexSet()  : destructor
//    operator =   : assignment operator for sets
//    operator +   : unification of two sets or of a set and an index
//    operator -   : difference of two sets or of a set and an index
//    operator []  : returns an element of the set
//    operator ==  : test for equality of two sets
// Global functions:
//    Size()       : returns the actual number of elements in the set
//    Complement() : returns the complementary set
//    SetToVector(): assignment of an index set to an integer vector
//----------------------------------------------------------------------------
#ifndef __SET_ARI_HPP
#define __SET_ARI_HPP

#include <intvector.hpp>     // Integer vector type


using namespace cxsc;
using namespace std;


class IndexSet {
  int   MaxSize;
  char* Index;

  public:
    ~IndexSet ( );
    IndexSet  ( int = 0, char = '\0' );
    IndexSet  ( const IndexSet& );
    IndexSet& operator= ( const IndexSet& );
    IndexSet  operator+ ( const IndexSet& ) const;
    IndexSet  operator- ( const IndexSet& ) const;
    IndexSet  operator+ ( int ) const;
    IndexSet  operator- ( int ) const;
    int       operator[] ( int ) const;
    int       operator== ( const IndexSet& );

  friend int      Size        ( const IndexSet& );
  friend IndexSet Complement  ( const IndexSet& );
  friend void     SetToVector ( const IndexSet&, const intmatrix_subv& );
};

int      Size        ( const IndexSet& );
IndexSet Complement  ( const IndexSet& );
void     SetToVector ( const IndexSet&, const intmatrix_subv& );
#endif




