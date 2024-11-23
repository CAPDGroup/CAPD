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

/* CVS $Id: matrtest.cpp,v 1.24 2014/01/30 17:23:47 cxsc Exp $ */

#define CXSC_INDEX_CHECK 1
#include <iostream>

#include "imatrix.hpp"
#include "lrvector.hpp"

namespace cxsc {

using namespace std;

int main()
{
	imatrix A(6,6);
	l_rvector v(6);
			  
	for(int i=1;i<=6;i++)
	{
		v[i]=2;
		for(int j=1;j<=6;j++)
			A[i][j]=(i-1)*6+j-1;
	}
	cout<<A<<endl;
	cout<<A[2]<<endl;
	cout<<v<<endl;
	cout<<"A[2]*v="<<(A[2]*v)<<endl;
}

} // namespace cxsc

