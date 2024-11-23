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

/* CVS $Id: matrix.cpp,v 1.23 2014/01/30 17:23:47 cxsc Exp $ */

#include <iostream>

namespace cxsc {

	template <class M>
	std::ostream &_mout(std::ostream &s,const M &r) throw()
	{
		int i,j;
		for (i=0;i<r.ysize;i++)
		{
			for (j=0;j<r.xsize;j++)
			{
				s << r.dat[i*r.xsize+j]<<" ";
			}
			s<<std::endl;
		}
		return s;
	}

	template <class M>
	std::istream &_min(std::istream &s,M &r) throw()
	{
		int i,j;
		for (i=0;i<r.ysize;i++)
		{
			for (j=0;j<r.xsize;j++)
			{
				s >> r.dat[i*r.xsize+j];
			}
		}
		return s;
	}

	template <class MS>
	std::ostream &_msout(std::ostream &s,const MS &r) throw()
	{
		int i,j;
		for (i=0;i<r.sysize;i++)
		{
			for (j=0;j<r.sxsize;j++)
			{
				s << r.dat[(i+r.offset1)*r.mxsize+j+r.offset2]<<" ";
			}
			s<<std::endl;
		}
		return s;
	}

	template <class MS>
	std::istream &_msin(std::istream &s,MS &r) throw()
	{
		int i,j;
		for (i=0;i<r.sysize;i++)
		{
			for (j=0;j<r.sxsize;j++)
			{
				s >> r.dat[(i+r.offset1)*r.mxsize+j+r.offset2];
			}
		}
		return s;
	}

} // namespace cxsc

