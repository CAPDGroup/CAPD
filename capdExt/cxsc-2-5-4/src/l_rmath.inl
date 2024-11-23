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

/* CVS $Id: l_rmath.inl,v 1.30 2014/01/30 17:23:46 cxsc Exp $ */

#include "l_interval.hpp"
#include "l_imath.hpp"

namespace cxsc {

inline l_real sqr(const l_real &x) throw() 
   { return (x*x); }
inline l_real sqrt(const l_real & x,int n) 
   { return mid(sqrt(_l_interval(x),n)); }
inline l_real sqrtp1m1(const l_real & x) throw()
   { return mid(sqrtp1m1(_l_interval(x))); }
inline l_real sqrtx2m1(const l_real & x) throw() 
   { return mid(sqrtx2m1(_l_interval(x))); } // sqrt(x^2-1);
inline l_real sqrt1mx2(const l_real& x) throw() 
   { return mid(sqrt1mx2(_l_interval(x)));}
inline l_real ln_sqrtx2y2(const l_real& x, const l_real& y) throw()
   { return mid(ln_sqrtx2y2(_l_interval(x),_l_interval(y)));}
inline l_real acoshp1(const l_real& x)
   { return mid(acoshp1(l_interval(x)));}

inline l_real sin(const l_real & x) throw() 
   { return mid(sin(_l_interval(x))); }
inline l_real cos(const l_real & x) throw() 
   { return mid(cos(_l_interval(x))); }
inline l_real tan(const l_real & x) throw() 
   { return mid(tan(_l_interval(x))); }
inline l_real cot(const l_real & x) throw() 
   { return mid(cot(_l_interval(x))); }

inline l_real asin(const l_real & x) 
   { return mid(asin(_l_interval(x))); }
inline l_real acos(const l_real & x)
   { return mid(acos(_l_interval(x))); }
inline l_real atan(const l_real & x)
   { return mid(atan(_l_interval(x))); }
inline l_real acot(const l_real & x)
   { return mid(acot(_l_interval(x))); }

inline l_real exp(const l_real & x) throw() 
   { return mid(exp(_l_interval(x))); }
inline l_real exp2(const l_real & x) throw() 
{ return mid(exp2(_l_interval(x))); }
inline l_real exp10(const l_real & x) throw() 
{ return mid(exp10(_l_interval(x))); }
inline l_real expm1(const l_real &x) throw()
   { return mid(expm1(l_interval(x))); }
inline l_real expmx2(const l_real& x) throw() 
   { return mid(expmx2(l_interval(x)));}

inline l_real ln(const l_real & x)
   { return mid(ln(_l_interval(x))); }
inline l_real log2(const l_real & x)
{ return mid(log2(_l_interval(x))); }
inline l_real log10(const l_real & x)
{ return mid(log10(_l_interval(x))); }

inline l_real sinh(const l_real & x) throw() 
   { return mid(sinh(_l_interval(x))); }
inline l_real cosh(const l_real & x) throw() 
   { return mid(cosh(_l_interval(x))); }
inline l_real tanh(const l_real & x) throw() 
   { return mid(tanh(_l_interval(x))); }
inline l_real coth(const l_real & x) throw() 
   { return mid(coth(_l_interval(x))); }

inline l_real asinh(const l_real & x)
   { return mid(asinh(_l_interval(x))); }
inline l_real acosh(const l_real & x)
   { return mid(acosh(_l_interval(x))); }
inline l_real atanh(const l_real & x)
   { return mid(atanh(_l_interval(x))); }
inline l_real acoth(const l_real & x)
   { return mid(acoth(_l_interval(x))); }

inline l_real pow(const l_real & x,const l_real & expo) 
   { return mid(pow(_l_interval(x),_l_interval(expo))); }

} // namespace cxsc

