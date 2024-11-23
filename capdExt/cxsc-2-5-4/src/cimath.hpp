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

/* CVS $Id: cimath.hpp,v 1.24 2014/01/30 17:23:43 cxsc Exp $ */

/*
**
**  COmplex interval STandard functions LibrarY, CoStLy Version 1.0.3
**
**  Copyright (C) Markus Neher,        markus.neher@math.uni-karlsruhe.de
**                Ingo Eble,           ingoeble@web.de
**                Frithjof Blomquist,  Blomquist@math.uni-wuppertal.de
**
**  The complex interval elementary functions in C-XSC are based on the 
**  CoStLy library written by Markus Neher. Additional improvements have
**  been done by Frithjof Blomquist.
**  
**  References:
**  - Neher, M: "Complex Standard Functions and their Implementation in
**    the CoStLy Library", Preprint Nr. 04/18, Fakultaet fuer Mathematik,
**    Universitaet Karlsruhe, 2004.
**  - Blomquist, F.; Hofschuster, W.; Kraemer, W.: "Complex Interval Functions 
**    in C-XSC", Preprint BUW-WRSWT 2005/2, Universitaet Wuppertal, 2005.
**  - Neher, M: "Complex Standard Functions and their Implementation in
**    the CoStLy Library", to appear in ACM Transactions on Mathematical 
**    Software (TOMS).
**
*/


#ifndef _CXSC_CIMATH_HPP_INCLUDED
#define _CXSC_CIMATH_HPP_INCLUDED

#include "cinterval.hpp"
#include <list>
#include <string>
#include <cmath>

namespace cxsc{

//! Calculates \f$ \exp([z]) \f$
cinterval exp(const cinterval&) throw();
//! Calculates \f$ \exp([z])-1 \f$
cinterval expm1(const cinterval&) throw();
//! Calculates \f$ 2^{[z]} \f$
cinterval exp2(const cinterval&) throw();
//! Calculates \f$ 10^{[z]} \f$
cinterval exp10(const cinterval&) throw();
//! Calculates \f$ \cos([z]) \f$
cinterval cos(const cinterval&) throw();
//! Calculates \f$ \sin([z]) \f$
cinterval sin(const cinterval&) throw();
//! Calculates \f$ \cosh([z]) \f$
cinterval cosh(const cinterval&) throw();
//! Calculates \f$ \sinh([z]) \f$
cinterval sinh(const cinterval&) throw();

//! Calculates \f$ \tan([z]) \f$
cinterval tan(const cinterval&) throw();
//! Calculates \f$ \cot([z]) \f$
cinterval cot(const cinterval&) throw();
//! Calculates \f$ \tanh([z]) \f$
cinterval tanh(const cinterval&) throw();
//! Calculates \f$ \coth([z]) \f$
cinterval coth(const cinterval&) throw();

//! Calculates \f$ \mbox{arg}([z]) \f$
interval arg(const cinterval&) throw();
interval arg_inclmon(const cinterval&) throw();
//! Calculates \f$ \mbox{arg}([z]) \f$
interval Arg(const cinterval&) throw();

//! Calculates \f$ \ln([z]) \f$
cinterval ln(const cinterval&) throw();
//! Calculates \f$ \ln([z]) \f$
cinterval Ln(const cinterval&) throw();
//! Calculates \f$ \ln(1+[z]) \f$
cinterval lnp1(const cinterval&) throw();

//! Calculates \f$ \mbox{log2}([z]) \f$
cinterval log2(const cinterval&) throw();
//! Calculates \f$ \mbox{log2}([z]) \f$
cinterval log10(const cinterval&) throw();

//! Calculates \f$ [z]^2  \f$
cinterval sqr(const cinterval&) throw();

//! Calculates \f$ \sqrt{[z]}  \f$
cinterval sqrt(const cinterval&) throw();
//! Calculates \f$ \sqrt{1+[z]}-1  \f$
cinterval sqrtp1m1(const cinterval&) throw();
//! Calculates \f$ \sqrt{1+[z]^2}  \f$
cinterval sqrt1px2(const cinterval&) throw();
//! Calculates \f$ \sqrt{[z]^2-1}  \f$
cinterval sqrtx2m1(const cinterval&) throw();
//! Calculates \f$ \sqrt{1-[z]^2}  \f$
cinterval sqrt1mx2(const cinterval&) throw();

//! Calculates \f$ \sqrt{[z]}  \f$ and returns all possible solutions
std::list<cinterval>sqrt_all(const cinterval&);
//! Calculates \f$ \sqrt[n]{[z]} \f$
cinterval sqrt(const cinterval&, int) throw();
//! Calculates \f$ \sqrt[n]{[z]} \f$ and returns all possible solutions
std::list<cinterval>sqrt_all(const cinterval&, int);

//! Calculates \f$ [z]^n \f$
cinterval power_fast(const cinterval&,int) throw();
//! Calculates \f$ [z]^n \f$
cinterval power(const cinterval&,int) throw();
//! Calculates \f$ [z]^{[y]} \f$
cinterval pow(const cinterval&, const interval&) throw();
//! Calculates \f$ [z_1]^{[z_2]} \f$
cinterval pow(const cinterval&, const cinterval&) throw();
//! Calculates \f$ [z]^{[y]} \f$ and returns all possible solutions
std::list<cinterval>pow_all(const cinterval&, const interval&) throw();
//! Fast multiplication of reference parameter [z] with \f$ 2^n \f$
void times2pown(cinterval& x, int n) throw();

//! Calculates \f$ \arcsin([z]) \f$
cinterval asin(const cinterval&) throw();
//! Calculates \f$ \arccos([z]) \f$
cinterval acos(const cinterval&) throw();
//! Calculates \f$ \mbox{arcsinh}([z]) \f$
cinterval asinh(const cinterval&) throw();
//! Calculates \f$ \mbox{arccosh}([z]) \f$
cinterval acosh(const cinterval&) throw();
//! Calculates \f$ \arctan([z]) \f$
cinterval atan(const cinterval&) throw();
//! Calculates \f$ \mbox{arccot}([z]) \f$
cinterval acot(const cinterval&) throw();
//! Calculates \f$ \mbox{arctanh}([z]) \f$
cinterval atanh(const cinterval&) throw();
//! Calculates \f$ \mbox{arccoth}([z]) \f$
cinterval acoth(const cinterval&) throw();

} // namespace cxsc

#endif // _CXSC_CIMATH_HPP_INCLUDED

/*

  End of File: cimath.hpp

*/
