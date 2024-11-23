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

/* CVS $Id: l_cimath.hpp,v 1.18 2014/01/30 17:23:46 cxsc Exp $ */

/*
**
**  File: l_cimath.hpp, 2007/03/04
**
**  Complex l_interval Standard functions LibrarY, Version 1.0
**
**  Copyright (C) Markus Neher, markus.neher@math.uni-karlsruhe.de
**                Ingo Eble,    ingoeble@web.de
**                Frithjof Blomquist,  Blomquist@math.uni-wuppertal.de
*/


#ifndef _CXSC_L_CIMATH_HPP_INCLUDED
#define _CXSC_L_CIMATH_HPP_INCLUDED

#include "l_cinterval.hpp"
#include <list>
#include <string>

namespace cxsc{

   //! Calculates \f$ \exp([z]) \f$
   l_cinterval exp(const l_cinterval&) throw();
   //! Calculates \f$ 2^{[z]} \f$
   l_cinterval exp2(const l_cinterval&) throw();
   //! Calculates \f$ 10^{[z]} \f$
   l_cinterval exp10(const l_cinterval&) throw();
   //! Calculates \f$ \exp([z])-1 \f$
   l_cinterval expm1(const l_cinterval&) throw();
   //! Calculates \f$ \cos([z]) \f$
   l_cinterval cos(const l_cinterval&) throw();
   //! Calculates \f$ \sin([z]) \f$
   l_cinterval sin(const l_cinterval&) throw();
   //! Calculates \f$ \cosh([z]) \f$
   l_cinterval cosh(const l_cinterval&) throw();
   //! Calculates \f$ \sinh([z]) \f$
   l_cinterval sinh(const l_cinterval&) throw();

   //! Calculates \f$ [z]^2  \f$
   l_cinterval sqr(const l_cinterval&) throw();

   //! Calculates \f$ \tan([z]) \f$
   l_cinterval tan(const l_cinterval&) throw();
   //! Calculates \f$ \cot([z]) \f$
   l_cinterval cot(const l_cinterval&) throw();
   //! Calculates \f$ \tanh([z]) \f$
   l_cinterval tanh(const l_cinterval&) throw();
   //! Calculates \f$ \coth([z]) \f$
   l_cinterval coth(const l_cinterval&) throw();

   // l_interval Atan(const l_interval& y, const l_interval& x) throw(); 

   //! Calculates \f$ \mbox{arg}([z]) \f$
   l_interval arg(const l_cinterval&) throw();
   l_interval arg_inclmon(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{arg}([z]) \f$
   l_interval Arg(const l_cinterval&) throw();

   //! Calculates \f$ \ln([z]) \f$
   l_cinterval Ln(const l_cinterval&) throw();
   //! Calculates \f$ \ln([z]) \f$
   l_cinterval ln(const l_cinterval&) throw();
   //! Calculates \f$ \ln(1+[z]) \f$
   l_cinterval lnp1(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{log2}([z]) \f$
   l_cinterval log2(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{log10}([z]) \f$
   l_cinterval log10(const l_cinterval&) throw();

   //! Calculates \f$ \sqrt{[z]}  \f$
   l_cinterval sqrt(const l_cinterval&) throw();
   //! Calculates \f$ \sqrt{1+[z]}-1  \f$
   l_cinterval sqrtp1m1(const l_cinterval&) throw();
   //! Calculates \f$ \sqrt{1+[z]^2}  \f$
   l_cinterval sqrt1px2(const l_cinterval&) throw();
   //! Calculates \f$ \sqrt{[z]^2-1}  \f$
   l_cinterval sqrtx2m1(const l_cinterval&) throw();
   //! Calculates \f$ \sqrt{1-[z]^2}  \f$
   l_cinterval sqrt1mx2(const l_cinterval&) throw();
   //! Calculates \f$ \sqrt{[z]}  \f$ and returns all possible solutions
   std::list<l_cinterval>sqrt_all(const l_cinterval&);
   //! Calculates \f$ \sqrt[n]{[z]} \f$
   l_cinterval sqrt(const l_cinterval&, int) throw();
   //! Calculates \f$ \sqrt[n]{[z]} \f$ and returns all possible solutions
   std::list<l_cinterval>sqrt_all(const l_cinterval&, int);
   //! Calculates \f$ [z]^n \f$
   l_cinterval power_fast(const l_cinterval&, int) throw();
   //! Calculates \f$ [z]^n \f$
   l_cinterval power(const l_cinterval&,int) throw();
   //! Calculates \f$ [z]^[y] \f$
   l_cinterval pow(const l_cinterval&, const l_interval&) throw();
   //! Calculates \f$ [z_1]^{[z_2]} \f$
   l_cinterval pow(const l_cinterval&, const l_cinterval&) throw();
   //! Calculates \f$ [z]^{[y]} \f$ and returns all possible solutions
   std::list<l_cinterval>pow_all(const l_cinterval&, 
                                 const l_interval&) throw();
   //! Calculates \f$ \arcsin([z]) \f$
   l_cinterval asin(const l_cinterval&) throw();
   //! Calculates \f$ \arccos([z]) \f$
   l_cinterval acos(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{arcsinh}([z]) \f$
   l_cinterval asinh(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{arccosh}([z]) \f$
   l_cinterval acosh(const l_cinterval&) throw();

   //! Calculates \f$ \arctan([z]) \f$
   l_cinterval atan(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{arccot}([z]) \f$
   l_cinterval acot(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{arctanh}([z]) \f$
   l_cinterval atanh(const l_cinterval&) throw();
   //! Calculates \f$ \mbox{arccoth}([z]) \f$
   l_cinterval acoth(const l_cinterval&) throw();

} // namespace cxsc

#endif // _CXSC_L_CIMATH_HPP_INCLUDED

/*

  End of File: l_cimath.hpp

*/
