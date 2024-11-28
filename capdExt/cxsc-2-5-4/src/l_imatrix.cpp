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

/* CVS $Id: l_imatrix.cpp,v 1.18 2014/01/30 17:23:46 cxsc Exp $ */

#define _CXSC_CPP

#include "l_imatrix.hpp"
#include "matrix.inl"
#include "l_imatrix.inl"
#include "livecrmat.inl"
#include "liveclrmat.inl"
#include "livecimat.inl"

// they should also be included here
#include "lrvecimat.inl"
#include "iveclrmat.inl"
#include "lrmatimat.inl"

namespace cxsc {
l_imatrix Id ( const l_imatrix& A )                // l_interval identity matrix
{                                                  //---------------------------
  int i,j;
  int lbi = Lb(A,1), ubi = Ub(A,1);
  int lbj = Lb(A,2), ubj = Ub(A,2);
  l_imatrix B(lbi,ubi,lbj,ubj);

  for (i = lbi; i <= ubi; i++)
    for (j = lbj; j <= ubj; j++)
      B[i][j] = l_interval( (i==j) ? 1.0 : 0.0 );
  return B;
}

l_imatrix transp ( const l_imatrix& A )              // Transposed matrix
{                                                    //------------------
  int      n;
  l_imatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));
      
  for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
  return res;
}

l_real MaxRelDiam ( const l_imatrix_subv& v )    // Maximum relative diameter
{                                                //--------------------------
  l_real r;
  int  i, l=Lb(v), u=Ub(v);

  r = 0.0;
  for (i = l; i <= u; i++)
    if (RelDiam(v[i]) > r) r = RelDiam(v[i]);
  return r;
}

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------

void DoubleSize ( l_imatrix& A )
{
  int n = Lb(A,1);
  Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
}

} // namespace cxsc

