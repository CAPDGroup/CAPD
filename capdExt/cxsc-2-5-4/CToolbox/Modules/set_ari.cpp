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

/* CVS $Id: set_ari.cpp,v 1.14 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: set_ari (implementation)
// Purpose: Implementation of basic operations for sets of integers (also
//    called indices).
// Method : A set of integers is represented by an array of characters INDEX.
//    An index N is element of the set if INDEX[N] != '\0'. For example
//
//       INDEX =  ('\0','1','1','\0','1','\0','1')
//
//    represents the set {2,3,5,7}. The maximum number of set elements is
//    restricted to MAXSIZE elements.
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
#include <cstdlib>         // For error handling with 'exit()'
#include <cstring>         // For memory copying with 'memcpy()'
#include <iostream>        // I/O handling with 'cerr'
#include <mv_util.hpp>     // Integer and real matrix/vector utility functions
#include <set_ari.hpp>

using namespace cxsc;
using namespace std;

//----------------------------------------------------------------------------
// Fatal errors for set handling
//----------------------------------------------------------------------------
const int
  IllSetRange = 1,  // Illegal set range
  DiffSetSize = 2;  // Different maximum sizes of two sets

static void FatalError ( int Err )
{
  switch (Err) {
    case IllSetRange : cerr << "Illegal set range, index can not be stored!";
                       break;
    case DiffSetSize : cerr << "Sets have different maximum size!"; break;
  }
  cerr << endl;
  exit(-1);
}

//----------------------------------------------------------------------------
// Destructor and constructors
//----------------------------------------------------------------------------
IndexSet::~IndexSet ( )                                          // Destructor
  { delete[] Index; }                                            //-----------

//----------------------------------------------------------------------------
// By default all components of the index array are initialized by '\0'.
// Thus, the default initialization for a set is the empty set {}.
//----------------------------------------------------------------------------
IndexSet::IndexSet ( int NewSize, char FillChar )
{
  MaxSize = NewSize;
  if (MaxSize == 0)
    { Index = NULL; return; }
  Index = new char[MaxSize];
  for (int i = 0; i < MaxSize; i++) Index[i] = FillChar;
}

IndexSet::IndexSet ( const IndexSet& s )                   // Copy constructor
{                                                          //-----------------
  MaxSize = s.MaxSize;
  if (MaxSize == 0)
    { Index = NULL; return; }
  Index = new char[MaxSize];
  memcpy(Index, s.Index, MaxSize*sizeof(char));
}

//----------------------------------------------------------------------------
// Assignment operator
//----------------------------------------------------------------------------
IndexSet& IndexSet::operator= ( const IndexSet& s )
{
  MaxSize = s.MaxSize;
  Index = new char[MaxSize];
  memcpy(Index, s.Index, MaxSize*sizeof(char));
  return *this;
}

//----------------------------------------------------------------------------
// Unification and difference of two sets. Both sets must have the same
// maximum size.
//----------------------------------------------------------------------------
IndexSet IndexSet::operator+ ( const IndexSet& s ) const
{
  // Sets to be unified must have the same maximum size!
  //----------------------------------------------------
  if (MaxSize != s.MaxSize) FatalError(DiffSetSize);

  IndexSet res(MaxSize);
  for (int i = 0; i < MaxSize; i++) res.Index[i] = Index[i] | s.Index[i];
  return res;
}

IndexSet IndexSet::operator- ( const IndexSet& s ) const
{
  // Sets to be subtracted must have the same maximum size!
  //-------------------------------------------------------
  if (MaxSize != s.MaxSize) FatalError(DiffSetSize);

  IndexSet res(MaxSize);
  for (int i = 0; i < MaxSize; i++)
    res.Index[i] = (s.Index[i]) ? '\0' : Index[i];
  return res;
}

//----------------------------------------------------------------------------
// Union and difference of a set and a single element. The integer value must
// lie within the range of the set.
//----------------------------------------------------------------------------
IndexSet IndexSet::operator+ ( int i ) const
{
  // 'i' must lie in the range of the set!
  //--------------------------------------
  if (i < 1 || MaxSize < i) FatalError(IllSetRange);

  IndexSet res = *this;
  res.Index[i-1] = '\1';
  return res;
}

IndexSet IndexSet::operator- ( int i ) const
{
  // 'i' must lie in the range of the set!
  //--------------------------------------
  if (i < 1 || MaxSize < i) FatalError(IllSetRange);

  IndexSet res = *this;
  res.Index[i-1] = '\0';
  return res;
}

//----------------------------------------------------------------------------
// Returns the n-th element of the set. The elements of the set are
// assumed to be sorted in increasing order.
//----------------------------------------------------------------------------
int IndexSet::operator[] ( int n ) const
{
  int i = 0, cnt = 0;

  while (i < MaxSize && cnt != n) {
    if (Index[i] != '\0') cnt++;
    i++;
  }
  return i;
}

int IndexSet::operator== ( const IndexSet& s )            // Test for equality
{                                                         //------------------
  // Sets to be compared must have the same size!
  //---------------------------------------------
  if (MaxSize != s.MaxSize) FatalError(DiffSetSize);

  int res = 1, i = 0;

  while (res == 1 && i < MaxSize) {
    res = (Index[i] == s.Index[i]);
    i++;
  }
  return res;
}

int Size ( const IndexSet& s )        // Returns the number of elements in 's'
{                                     //--------------------------------------
  int i, cnt = 0;

  for (i = 0; i < s.MaxSize; i++)
    if (s.Index[i] != '\0') cnt++;
  return cnt;
}

//----------------------------------------------------------------------------
// Returns the complementary set of 's', i.e. {1,2,..,MaxSize} - s
//----------------------------------------------------------------------------
IndexSet Complement ( const IndexSet& s )
{
  IndexSet t(s.MaxSize,'\1');
  return t - s;
}

void SetToVector ( const IndexSet& B, const intmatrix_subv& x )    // Assign index set to
{                                                       // integer vector
                                                        //--------------------
  // Lower index bound 1 assumed!
  //-----------------------------
  for (int i = 1; i <= Ub(x); i++) x[i] = B[i];
}




