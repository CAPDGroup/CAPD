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

/* CVS $Id: lop_ari.cpp,v 1.15 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: lop_ari (implementation)
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
#include <lop_ari.hpp>

using namespace cxsc;
using namespace std;

struct BaseListElement {          // Type for a list of integer sets
  IndexSet Member;                //--------------------------------
  BaseList Next;
};

static BaseList                   // Local variable storing list of free
  FreeList = EmptyList;           // elements (automatic garbage recycling)
                                  //---------------------------------------

//----------------------------------------------------------------------------
// Functions for generating and freeing of list elements
//----------------------------------------------------------------------------
static void NewBL ( BaseList& P )     // 'NewBL' generates a new list element
{                                     // or gets one from 'FreeList'.
  if (FreeList == EmptyList)          //--------------------------------------
    { P = new BaseListElement;  P->Next = NULL; }
  else
    { P = FreeList;  FreeList = FreeList->Next;  P->Next= NULL; }
}

static void Free ( BaseList& P )            // 'Free' enters one element of a
{                                           // list in the 'FreeList'.
  if (P != NULL)                            //--------------------------------
    { P->Next = FreeList;  FreeList = P;  P = NULL; }
}

void FreeAll ( BaseList& List )               // 'FreeAll' enters all elements
{                                             // of 'List' in the 'FreeList'.
  BaseList H;                                 //------------------------------

  if (List != EmptyList) {
    H = List;
    while (H->Next != NULL) H = H->Next;
    H->Next = FreeList;  FreeList = List;  List = NULL;
  }
}

rmatrix extract ( rmatrix& A, const IndexSet& B )     // Extract submatrix of
{                                                     // of columns A[*,Index]
  int     k, n = Size(B);                             // with (Index in B)
  rmatrix res(Ub(A,1),n);                             //----------------------

  for (k = 1; k <= n; k++) Col(res,k) = Col(A,B[k]);
  return res;
}

rvector extract ( rvector& x, const IndexSet& B )      // Extract subvector
{                                                      // x[Index] with (Index
  int     k, n = Size(B);                              // in B)
  rvector res(n);                                      //---------------------

  for (k = 1; k <= n; k++) res[k] = x[B[k]];
  return res;
}

//----------------------------------------------------------------------------
// Global functions for lists of index sets
//----------------------------------------------------------------------------
int in ( const IndexSet& B, BaseList L )                    // Check if B in L
{                                                           //----------------
  while (L != EmptyList)
    if (L->Member == B)
      return 1;
    else
      L = L->Next;
  return 0;
}

IndexSet select ( BaseList L )                // Select first index set from L
  { return L->Member; }                       //------------------------------

void insert ( BaseList& L, const IndexSet& B )    // Insert index set B to the
{                                                 // head of the list of index
  BaseList P;                                     // sets L
                                                  //--------------------------
  if ( !in(B,L) ) {
    P = L;
    NewBL(L);  L->Member = B;  L->Next = P;
  }
}

void del ( BaseList& L, const IndexSet& B )    // Delete index set B from list
{                                              // of index sets L
  BaseList  P, D;                              //-----------------------------

  if (L != EmptyList) {
    if (L->Member == B)                        // B is 1st element
      { D = L; L = L->Next; Free(D); }
    else {                                     // B is not 1st element
      P = L;
      while (P->Next != NULL)
        if (P->Next->Member == B)
          { D = P->Next;  P->Next = P->Next->Next;  Free(D); }
        else
          P = P->Next;
    }
  }
}

void append ( BaseList& Res, BaseList& Add )  // Append list of index sets Add
{                                             // to list Res
  BaseList P;                                 //------------------------------

  if (Res == EmptyList)
    Res = Add;
  else {
    P = Res;
    while (P->Next != NULL) P = P->Next;
    P->Next = Add;
  }
  Add = EmptyList;
}

void remove ( BaseList& Res, BaseList Sub )   // Remove list of index sets Sub
{                                             // from Res
  IndexSet B;                                 //------------------------------

  while (Sub != EmptyList) {
    B = Sub->Member;
    del(Res,B);
    Sub = Sub->Next;
  }
}




