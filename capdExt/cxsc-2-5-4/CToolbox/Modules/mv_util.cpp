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

/* CVS $Id: mv_util.cpp,v 1.14 2014/01/30 17:49:27 cxsc Exp $ */

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
// File: mv_util (implementation)
// Purpose: Utilities of type 'intvector', 'intmatrix', 'rvector', and
//    'rmatrix'.
// Global functions:
//    VecLen()    : length of an integer or real vector
//    RowLen()    : length of the rows of an integer or real matrix
//    ColLen()    : length of the columns of an integer or real matrix
//    Id()        : identity real matrix
//    transp()    : transposed of a real matrix
//    DoubleSize(): for doubling the size of an integer vector or matrix
//    operator << : output of integer vector
//----------------------------------------------------------------------------
#include <mv_util.hpp>
/*
using namespace cxsc;
using namespace std;

int VecLen ( const intvector& v )          // Length of an integer vector
  { return Ub(v)-Lb(v)+1; }                //----------------------------

int RowLen ( const intmatrix& A )    // Length of the rows of a integer matrix
  { return Ub(A,2)-Lb(A,2)+1; }      //---------------------------------------

int ColLen ( const intmatrix& A )// Length of the columns of an integer matrix
  { return Ub(A,1)-Lb(A,1)+1; }  //-------------------------------------------

int VecLen ( const rvector& v )                     // Length of a real vector
  { return Ub(v)-Lb(v)+1; }                         //------------------------

int RowLen ( const rmatrix& A )         // Length of the rows of a real matrix
  { return Ub(A,2)-Lb(A,2)+1; }         //------------------------------------

int ColLen ( const rmatrix& A )      // Length of the columns of a real matrix
  { return Ub(A,1)-Lb(A,1)+1; }      //---------------------------------------

rmatrix Id ( const rmatrix& A )                        // Real identity matrix
{                                                      //---------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  rmatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = (i==j) ? 1.0 : 0.0;
  return B;
}

rmatrix transp ( rmatrix& A )                       // Transposed matrix
{                                                         //------------------
  int      n;
  rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));

  for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
  return res;
}

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------
void DoubleSize ( intvector& x )
{
  int n = Lb(x);
  Resize(x,n,2*Ub(x)-n+1);
}

void DoubleSize ( intmatrix& A )
{
  int n = Lb(A,1);
  Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
}

ostream& operator<< ( ostream& os, intvector& v )         // Output of integer
{                                                         // vectors
  int i, newline = (Ub(v)-Lb(v) > 15);                    //------------------

  for (i = Lb(v); i <= Ub(v); i++) {
    os << v[i] << ' ';
    if (newline) os << endl;
  }
  if (!newline) os << endl;
  return os;
}
*/




