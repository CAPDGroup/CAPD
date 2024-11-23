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

/* CVS $Id: intvector.cpp,v 1.15 2014/01/30 17:23:45 cxsc Exp $ */

#define _CXSC_CPP

#include "intvector.hpp"
#include "vector.inl"
#include "intvector.inl"

namespace cxsc {

int abs(int a)
{ return a<0?-a:a; }

// The 'DoubleSize' functions double the number of rows of a matrix
// or double the length of a vector preserving existing components.
//------------------------------------------------------------------
void DoubleSize ( intvector& x )
{
  int n = Lb(x);
  Resize(x,n,2*Ub(x)-n+1);
}

std::ostream& operator<< ( std::ostream& os, intvector& v ) // Output of integer
{                                                          // vectors
  int i, newline = (Ub(v)-Lb(v) > 15);                     //------------------

  for (i = Lb(v); i <= Ub(v); i++) {
    os << v[i] << ' ';
    if (newline) os << std::endl;
  }
  if (!newline) os << std::endl;
  return os;
}

} // namespace cxsc

