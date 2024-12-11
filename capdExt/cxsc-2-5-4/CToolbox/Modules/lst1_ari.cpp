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

/* CVS $Id: lst1_ari.cpp,v 1.15 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lst1_ari (implementation)
// Purpose: Definition of a list arithmetic used in connection with an
//    interval bisection method in global optimization for storing pairs of
//    an interval and a real value.
// Method: Overloading of functions and operators of class type 'Pair1' (list
//    element) and data type 'Pair1Ptr' (list). See module 'lst_ari' for the
//    definition of a class 'Pair' for pairs with a vector component 'fyi'.
//    The suffix 1 for the type 'Pair1' is needed to avoid a duplicate name
//    error when linking all toolbox modules to a library.
// Class Pair1:
//    Pair1()       : constructors
//    ~Pair1()      : destructor
//    operator =    : assignment operator
//    _Pair1()      : for explicit generating a pair
//    GetInt()      : get the interval component
//    GetFyi()      : get the function value
// Global type:
//    Pair1Ptr      : list of pairs
// Global functions and operators:
//    operator +    : adding a new element to a list
//    Head()        : access functions for the head of a list
//    Length()      : access function to length of list
//    FreeAll()     : free complete list
//    MultiDelete() : deletes several elements in a list
//    DelHead()     : deletes the first element of a list
//    ListToVector(): transfer a list into a dynamic interval vector
//----------------------------------------------------------------------------
#include <interval.hpp>     // Interval arithmetic
#include <intvector.hpp>    // Integer vector type
#include <lst1_ari.hpp>

using namespace cxsc;
using namespace std;

struct Pair1Elmt {        // Type for a list of pairs
  Pair1    P;             //-------------------------
  Pair1Ptr N;
};

static Pair1Ptr           // Local variable storing list of free
  FreeList = EmptyList;   // elements (automatic garbage recycling)
                          //---------------------------------------

//----------------------------------------------------------------------------
// Constructors and assignment operator
//----------------------------------------------------------------------------
Pair1::Pair1 ( const Pair1& a )
  { intv = a.intv; fyi  = a.fyi; }

Pair1& Pair1::operator= ( const Pair1& a )
{
  intv = a.intv;
  fyi  = a.fyi;
  return *this;
}

//----------------------------------------------------------------------------
// Transfer and access functions for pairs
//----------------------------------------------------------------------------
Pair1 _Pair1 ( const interval& i, const real& r )
{
  Pair1 p;
  p.intv = i;
  p.fyi  = r;
  return p;
}

interval GetInt ( const Pair1& p )         // Get int-component
  { return p.intv; }                 //-------------------

real GetFyi ( const Pair1& p )             // Get fyi-component
  { return p.fyi; }                  //-------------------

//----------------------------------------------------------------------------
// Functions for generating and freeing of list elements (pairs)
//----------------------------------------------------------------------------
static void NewPP ( Pair1Ptr& pp )     // 'NewPP' generates a new list element
{                                      // or gets one from 'FreeList'.
  if (FreeList == EmptyList)           //-------------------------------------
    { pp = new Pair1Elmt; pp->N = NULL; }
  else
    { pp = FreeList; FreeList = FreeList->N; pp->N = NULL; }
}

static void Free ( Pair1Ptr& pp )            // 'Free' enters one element of a
{                                            // list in the 'FreeList'.
  if (pp != NULL)                            //-------------------------------
    { pp->N = FreeList; FreeList = pp; pp = NULL; }
}

void FreeAll ( Pair1Ptr& List )               // 'FreeAll' enters all elements
{                                             // of 'List' in the 'FreeList'
  Pair1Ptr H;                                 //------------------------------

  if (List != EmptyList) {
    H = List;
    while (H->N != NULL) H = H->N;
    H->N = FreeList; FreeList = List; List = NULL;
  }
}

//----------------------------------------------------------------------------
// Operator + enters the pair 'P' on the list 'List' in such a way that after
// entering, one of the four following condition holds:
//   1) O.fyi <= P.fyi < Q.fyi,
//   2)          P.fyi < Q.fyi  and  'P' is the first element of 'List',
//   3) O.fyi <= P.fyi          and  'P' is the last  element of 'List',
//   4)                              'P' is the only  element of 'List',
// where 'O' is the preceding and 'Q' is the succeeding element of 'P' in
// the resulting list.
//----------------------------------------------------------------------------
Pair1Ptr operator+ ( Pair1Ptr List, Pair1 p )
{
  Pair1Ptr H, HN;
  int     ready, alreadyIn;

  if (List == EmptyList) {                            // List is empty, so new
    NewPP(H); H->P = p; H->N = NULL;                  // list contains only P.
    return H;                                         //----------------------
  }
  else if (GetFyi(List->P) > GetFyi(p)) {             // p becomes new first
    NewPP(H); H->P = p; H->N = List;                  // element of the list.
    return H;                                         //----------------------
  }
  else {
    H = List; HN = H->N; ready = 0;
    alreadyIn = (GetInt(H->P) == GetInt(p));

    while ( !(ready || alreadyIn) ) {                  // Search for the right
      if (HN == NULL)                                  // position to enter P.
        ready = 1;                                     //---------------------
      else if (GetFyi(HN->P) > GetFyi(p))
        ready = 1;
      else {
        H = HN; HN = H->N;
        alreadyIn = (GetInt(H->P) == GetInt(p));
      }
    } // while

    if (!alreadyIn)                                       // Enter p between H
      { NewPP(H->N); H = H->N; H->P = p; H->N = HN; }     // and HN.
                                                          //------------------
    return List;
  }
} // operator+

//----------------------------------------------------------------------------
// 'MultiDelete' deletes all elements 'P' in 'List' for which the condition
// 'P.fyi > fmax' holds.  This function assumes that the 'fyi' components of
// the list elements are sorted in increasing order (see operator +).
//----------------------------------------------------------------------------
void MultiDelete ( Pair1Ptr& List, const real& fmax )
{
  Pair1Ptr DelPrev, Del;
  int     ready;

  if (List != EmptyList) {
    if (GetFyi(List->P) > fmax)                   // All list elements fulfill
      { Del  = List; List = EmptyList; }          // 'P.fyi > fmax'.
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

Pair1 Head ( Pair1Ptr List )                     // Delivers first pair of the
  { return List->P; }                            // list, i.e. the pair P with
                                                 // the smallest value P.fyi.
                                                 //---------------------------

void DelHead ( Pair1Ptr& List )                   // Deletes the first pair of
{                                                 // the List
  Pair1Ptr Del;                                   //--------------------------

  Del = List; List = List->N; Free(Del);
}


int Length ( Pair1Ptr List )                   // 'Length' delivers the number
{                                              // of elements in list 'List'
  int     i = 0;                               //-----------------------------

  while (List != EmptyList)
    { i++; List = List->N; }
  return i;
}

//----------------------------------------------------------------------------
// Copies the interval components of every pair of the 'List' to the interval
// vector 'Vec' and allocates the intvector 'bv'. For an empty list the
// vectors 'Vec' and 'bv' are sized to unit length with undefined components.
//----------------------------------------------------------------------------
void ListToVector ( Pair1Ptr List, ivector& Vec, intvector& bv )
{
  int     Len, i;
  Pair1Ptr H;

  Len = Length(List);
                                                 // Resize the vectors:
  if (Len == 0)                                  //   - to minimum length with
    { Resize(Vec,1); Resize(bv,1); }             //     undefined components
  else                                           // or
    { Resize(Vec,Len); Resize(bv,Len); }         //   - to correct length
                                                 //---------------------------
  H = List;
  for (i = 1; i <= Len; i++)
    { Vec[i] = GetInt(H->P); H = H->N; }
}




