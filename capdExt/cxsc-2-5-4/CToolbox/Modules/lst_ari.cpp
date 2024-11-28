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

/* CVS $Id: lst_ari.cpp,v 1.15 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lst_ari (implementation)
// Purpose: Definition of a list arithmetic used in connection with an
//    interval bisection method in global optimization for storing pairs of
//    an interval vector and a real value.
// Method: Overloading of functions and operators of class type 'Pair'
//    (list element) and data type 'PairPtr' (list).
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
//    Head()        : access functions for the head of a list
//    Length()      : access function to length of list
//    FreeAll()     : free complete list
//    MultiDelete() : deletes several elements in a list
//    DelHead()     : deletes the first element of a list
//    ListToMatrix(): transfer a list into a dynamic interval matrix
//----------------------------------------------------------------------------
#include <lst_ari.hpp>

using namespace cxsc;
using namespace std;

struct PairElmt {         // Type for a list of pairs
  Pair    P;              //-------------------------
  PairPtr N;
};

static PairPtr            // Local variable storing list of free
  FreeList = EmptyList;   // elements (automatic garbage recycling)
                          //---------------------------------------

//----------------------------------------------------------------------------
// Constructor and assignment operator
//----------------------------------------------------------------------------
Pair::Pair ( const Pair& p )
  { intv = p.intv; fyi  = p.fyi; }

Pair& Pair::operator= ( const Pair& p )
{
  intv = p.intv;
  fyi  = p.fyi;
  return *this;
}

//----------------------------------------------------------------------------
// Transfer and access functions for pairs
//----------------------------------------------------------------------------
Pair _Pair ( const ivector& i, real r )    // Generate a pair
{                                    //----------------
  Pair p;
  p.intv = i;
  p.fyi  = r;
  return p;
}

Pair _Pair ( const ivector_slice& i, real r )    // Generate a pair
{                                    //----------------
  Pair p;
  p.intv = i;
  p.fyi  = r;
  return p;
}

Pair _Pair ( const imatrix_subv& i, real r )    // Generate a pair
{                                    //----------------
  Pair p;
  p.intv = i;
  p.fyi  = r;
  return p;
}

ivector GetInt ( const Pair& p )           // Get int-component
  { return p.intv; }                 //-------------------

real GetFyi ( const Pair& p )        // Get fyi-component
  { return p.fyi; }                  //-------------------

//----------------------------------------------------------------------------
// Functions for generating and freeing of list elements (pairs)
//----------------------------------------------------------------------------
static void NewPP ( PairPtr& pp )      // 'NewPP' generates a new list element
{                                      // or gets one from 'FreeList'.
  if (FreeList == EmptyList)           //-------------------------------------
    { pp = new PairElmt; pp->N = NULL; }
  else
    { pp = FreeList; FreeList = FreeList->N; pp->N = NULL; }
}

static void Free ( PairPtr& pp )             // 'Free' enters one element of a
{                                            // list in the 'FreeList'.
  if (pp != NULL)                            //-------------------------------
    { pp->N = FreeList; FreeList = pp; pp = NULL; }
}

void FreeAll ( PairPtr& List )                // 'FreeAll' enters all elements
{                                             // of 'List' in the 'FreeList'
  PairPtr H;                                  //------------------------------

  if (List != EmptyList) {
    H = List;
    while (H->N != NULL) H = H->N;
    H->N = FreeList; FreeList = List; List = NULL;
  }
}

//----------------------------------------------------------------------------
// Operator + enters the pair 'p' on the list 'List' in such a way that after
// entering, one of the four following condition holds:
//   1) o.fyi <= p.fyi < q.fyi,
//   2)          p.fyi < q.fyi  and  'p' is the first element of 'List',
//   3) o.fyi <= p.fyi          and  'p' is the last  element of 'List',
//   4)                              'p' is the only  element of 'List',
// where 'o' is the preceding and 'q' is the succeeding element of 'p' in
// the resulting list.
//----------------------------------------------------------------------------
PairPtr operator+ ( PairPtr List, const Pair& p )
{
  PairPtr H, HN;
  int     ready, alreadyIn;

  if (List == EmptyList) {                          // List is empty, so new
    NewPP(H); H->P = p;                             // list contains only 'p'.
    return H;                                       //------------------------
  }
  else if (GetFyi(List->P) > GetFyi(p)) {             // 'p' becomes new first
    NewPP(H); H->P = p; H->N = List;                  // element of the list.
    return H;                                         //----------------------
  }
  else {
    H = List; HN = H->N; ready = 0;
    alreadyIn = ( GetInt(H->P) == GetInt(p) );

    while ( !(ready || alreadyIn) ) {                // Search for the right
      if (HN == NULL)                                // position to enter 'p'.
        ready = 1;                                   //-----------------------
      else if (GetFyi(HN->P) > GetFyi(p))
        ready = 1;
      else {
        H = HN; HN = H->N;
        alreadyIn = ( GetInt(H->P) == GetInt(p) );
      }
    } // while

    if (!alreadyIn)                                       // Enter 'p' between
      { NewPP(H->N); H = H->N; H->P = p; H->N = HN; }     // H and HN.
                                                          //------------------
    return List;
  }
} // operator+

//----------------------------------------------------------------------------
// 'MultiDelete' deletes all elements 'p' in 'List' for which the condition
// 'p.fyi > fmax' holds.  This function assumes that the 'fyi' components of
// the list elements are sorted in increasing order (see operator +).
//----------------------------------------------------------------------------
void MultiDelete ( PairPtr& List, const real& fmax )
{
  PairPtr  DelPrev, Del;
  int      ready;

  if (List != EmptyList) {
    if (GetFyi(List->P) > fmax)                   // All list elements fulfill
      { Del  = List; List = EmptyList; }          // 'p.fyi > fmax'.
    else {                                        //--------------------------
      DelPrev = List; Del = DelPrev->N; ready = (Del == NULL);

      while (!ready) {
        if (Del == NULL)
          ready = 1;
        else if (GetFyi(Del->P) > fmax)
          { ready = 1; DelPrev->N = NULL; }
        else
          { DelPrev = Del; Del = Del->N; }
      }
    }
    FreeAll(Del);
  }
} // Multidelete

Pair Head ( PairPtr List )                       // Delivers first pair of the
  { return List->P; }                            // list, i.e. the pair P with
                                                 // the smallest value P.fyi.
                                                 //---------------------------

void DelHead ( PairPtr& List )                   // Deletes the first pair of
{                                                // the List
  PairPtr Del;                                   //---------------------------

  Del = List; List = List->N; Free(Del);
}


int Length ( PairPtr List )                    // 'Length' delivers the number
{                                              // of elements in list 'List'
  int     i = 0;                               //-----------------------------

  while (List != EmptyList)
    { i++; List = List->N; }
  return i;
}

//----------------------------------------------------------------------------
// Copies the interval vector components of every pair of the 'List' to the
// interval matrix 'Mat' and allocates the intvector 'bv'. For an empty list
// matrix 'Mat' and vector 'bv' are sized to unit length with undefined
// components.
//----------------------------------------------------------------------------
void ListToMatrix (PairPtr List, imatrix& Mat, intvector& bv)
{
  int      i, n, m;
  PairPtr  H;

  m = Length(List);                              // Resize the arrays:
  if (m == 0)                                    //   - to minimum length with
    { Resize(Mat,1,1); Resize(bv,1); }           //     undefined components
  else {                                         //
    n = Ub(GetInt(List->P));                     // or
    Resize(bv,m);                                //   - to correct dimensions
    Resize(Mat,m,n);                             //---------------------------
  }

  H = List;                                      // Copy vector components of
  for (i = 1; i <= m; i++)                       // the list to matrix 'Mat'
    { Mat[i] = GetInt(H->P); H = H->N; }         //--------------------------
}




