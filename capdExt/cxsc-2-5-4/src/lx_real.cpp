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

/* CVS $Id: lx_real.cpp,v 1.8 2014/01/30 17:23:47 cxsc Exp $ */

/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

#include "lx_real.hpp"
#include "lx_interval.hpp"

namespace cxsc {

inline std::ostream& operator << (std::ostream& s, const lx_real& b) throw()
// A value a of type lx_real is written to the output channel.
// The output has the form:  {2**(ex),lr}
{
	lx_interval a(b);
	real p;
	l_interval m;
	l_real x;
	
	Bin2Dec(a,p,m);
	x = mid(m);
	
	s << "{ " 
	  << "10**("
	  << SaveOpt << SetPrecision(0,0) << Fixed << p << RestoreOpt
	  << ")" 
	  << "*" 
	  << x
	  << " }";
	return s;
}

} // end namespace cxsc
