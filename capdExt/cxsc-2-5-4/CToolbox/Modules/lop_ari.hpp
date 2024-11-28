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

/* CVS $Id: lop_ari.hpp,v 1.14 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lop_ari (header)
// Purpose: Definition of a linearly linked list as an abstract data type for
//    the representation of a list of integer sets.
// Global type:
//    BaseList : representation of a linearly linked list of index sets
// Global constant:
//    EmptyList: the empty list (NULL pointer)
// Global functions:
//    FreeAll(): free complete list
//    extract(): extracts submatrix/subvector depending on the
//               actual index set
//    in()     : returns TRUE, if index set is in list
//    select() : selects first element of a list
//    insert() : inserts index set at the head of a list
//    del()    : deletes index set from list
//    append() : appends 2nd list to end of 1st list
//    remove() : removes elements of 2nd list from 1st list
//----------------------------------------------------------------------------
#ifndef __LOP_ARI_HPP
#define __LOP_ARI_HPP

#include <set_ari.hpp>     // Index set handling
#include <rmatrix.hpp>     // Real matrix/vector arithmetic

using namespace cxsc;
using namespace std;

#define EmptyList NULL           // The empty list
                                 //---------------

struct                           // Structure used for a list of sets
  BaseListElement;               //----------------------------------

typedef                          // Pointer to a list of integer sets
  BaseListElement* BaseList;     //----------------------------------

extern void     FreeAll ( BaseList& );
extern rmatrix  extract ( rmatrix&, const IndexSet& );
extern rvector  extract ( rvector&, const IndexSet& );
extern int      in      ( const IndexSet&, BaseList );
extern IndexSet select  ( BaseList );
extern void     insert  ( BaseList&, const IndexSet& );
extern void     del     ( BaseList&, const IndexSet& );
extern void     append  ( BaseList&, BaseList& );
extern void     remove  ( BaseList&, BaseList );
#endif




