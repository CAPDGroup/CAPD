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

/* CVS $Id: l_rmatrix.cpp,v 1.16 2014/01/30 17:23:46 cxsc Exp $ */

#define _CXSC_CPP

#include "l_rmatrix.hpp"
#include "matrix.inl"
#include "l_rmatrix.inl"
#include "lrvecrmat.inl"

namespace cxsc {
      l_rmatrix Id ( const l_rmatrix& A )           // l_real identity matrix
      {                                             //-----------------------
        int i,j;
        int lbi = Lb(A,1), ubi = Ub(A,1);
        int lbj = Lb(A,2), ubj = Ub(A,2);
        l_rmatrix B(lbi,ubi,lbj,ubj);
      
        for (i = lbi; i <= ubi; i++)
          for (j = lbj; j <= ubj; j++)
            B[i][j] = (i==j) ? 1.0 : 0.0;
        return B;
      }
      
      l_rmatrix transp ( const l_rmatrix& A )            // Transposed matrix
      {                                                  //------------------
        int      n;
        l_rmatrix  res(Lb(A,2),Ub(A,2),Lb(A,1),Ub(A,1));
      
        for (n = Lb(A,1); n <= Ub(A,1); n++) Col(res,n) = Row(A,n);
        return res;
      }

      void DoubleSize ( l_rmatrix& A )
      {
        int n = Lb(A,1);
        Resize(A,n,2*Ub(A,1)-n+1,Lb(A,2),Ub(A,2));
      }

} // namespace cxsc

