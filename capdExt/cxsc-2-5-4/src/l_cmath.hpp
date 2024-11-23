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

/* CVS $Id: l_cmath.hpp,v 1.13 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_L_CMATH_HPP_INCLUDED
#define _CXSC_L_CMATH_HPP_INCLUDED

#include <l_cimath.hpp>

namespace cxsc {
	
//! Calculates an approximation of \f$ z^2 \f$
inline l_complex sqr (const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt(z) \f$
l_complex sqrt(const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt(1+z)-1 \f$
l_complex sqrtp1m1(const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt(1+z^2) \f$
l_complex sqrt1px2(const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt(z^2-1) \f$
l_complex sqrtx2m1(const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt(1-z^2) \f$
l_complex sqrt1mx2(const l_complex&) throw();
//! Calculates an approximation of \f$ \exp(z) \f$
l_complex exp(const l_complex&) throw();
//! Calculates an approximation of \f$ \exp(z)-1 \f$
l_complex expm1(const l_complex&) throw();
//! Calculates an approximation of \f$ 2^z \f$
l_complex exp2(const l_complex&) throw();
//! Calculates an approximation of \f$ 10^z \f$
l_complex exp10(const l_complex&) throw();
//! Calculates an approximation of \f$ \sin(z) \f$
l_complex sin(const l_complex&) throw();
//! Calculates an approximation of \f$ \cos(z) \f$
l_complex cos(const l_complex&) throw();
//! Calculates an approximation of \f$ \tan(z) \f$
l_complex tan(const l_complex&) throw();
//! Calculates an approximation of \f$ \cot(z) \f$
l_complex cot(const l_complex&) throw();
//! Calculates an approximation of \f$ \arcsin(z) \f$
l_complex asin(const l_complex&) throw();
//! Calculates an approximation of \f$ \arccos(z) \f$
l_complex acos(const l_complex&) throw();
//! Calculates an approximation of \f$ \arctan(z) \f$
l_complex atan(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{arccot}(z) \f$
l_complex acot(const l_complex&) throw();
//! Calculates an approximation of \f$ \sinh(z) \f$
l_complex sinh(const l_complex&) throw();
//! Calculates an approximation of \f$ \cosh(z) \f$
l_complex cosh(const l_complex&) throw();
//! Calculates an approximation of \f$ \tanh(z) \f$
l_complex tanh(const l_complex&) throw();
//! Calculates an approximation of \f$ \coth(z) \f$
l_complex coth(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{arcsinh}(z) \f$
l_complex asinh(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{arccosh}(z) \f$
l_complex acosh(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{arctanh}(z) \f$
l_complex atanh(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{arccoth}(z) \f$
l_complex acoth(const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt{z}  \f$ and returns all possible solutions
std::list<l_complex>sqrt_all(const l_complex&);
//! Calculates an approximation of \f$ \sqrt[n]{z} \f$
l_complex sqrt(const l_complex&, int) throw();
//! Calculates an approximation of \f$ \mbox{arg}(z) \f$
l_real arg(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{arg}(z) \f$
l_real Arg(const l_complex&) throw();
//! Calculates an approximation of \f$ \sqrt[n]{z} \f$ and returns all possible solutions
std::list<l_complex>sqrt_all(const l_complex&, int);
//! Calculates an approximation of \f$ \ln(z) \f$
l_complex ln(const l_complex&) throw();
//! Calculates an approximation of \f$ \ln(1+z) \f$
l_complex lnp1(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{log2}(z) \f$
l_complex log2(const l_complex&) throw();
//! Calculates an approximation of \f$ \mbox{log10}(z) \f$
l_complex log10(const l_complex&) throw();
//! Calculates an approximation of \f$ z^n \f$
l_complex power(const l_complex&, int) throw();
//! Calculates an approximation of \f$ z^n \f$
l_complex power_fast(const l_complex&, int) throw();
//! Calculates an approximation of \f$ z^y \f$
l_complex pow(const l_complex&, const l_real&) throw();
//! Calculates an approximation of \f$ z_1^{z_2} \f$
l_complex pow(const l_complex&, const l_complex&) throw();

} // namespace cxsc 

#include "l_cmath.inl"
#endif // _CXSC_L_CMATH_HPP_INCLUDED
