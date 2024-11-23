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

/* CVS $Id: cxscmatr.cpp,v 1.22 2014/01/30 17:23:44 cxsc Exp $ */




#include "cxscmatr.hpp"

namespace cxsc {

cxscmatrix<T>::cxscmatrix<T>() throw()
{
	dat=NULL;
	lb1=lb2=1;
	ub1=ub2=xsize=ysize=0;
}

cxscmatrix<T>::cxscmatrix<T>(const int &m, const int &n) throw()
{
	lb1=lb2=1;
	ub1=ysize=m;
	ub2=xsize=n;
	dat=new T[m*n];
}

cxscmatrix<T>::cxscmatrix<T>(const int &m1, const int &n1, const int &m2, const int &n2)
{
	lb1=m1; ub1=m2;
	lb2=n1; ub2=n2;
	xsize=n2-n1+1;
	ysize=m2-m2+1;
	dat=new T[xsize*ysize];
}

} // namespace cxsc

