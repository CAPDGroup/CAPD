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

/* CVS $Id: imath.hpp,v 1.38 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_IMATH_HPP_INCLUDED
#define _CXSC_IMATH_HPP_INCLUDED

#include "except.hpp"
#include "interval.hpp"

namespace cxsc {

//! Calculates \f$ [x]^2  \f$
       interval sqr     (const interval&) throw(); // sqr(x)
//! Calculates \f$ \sqrt{[x]}  \f$
inline interval sqrt    (const interval&);         // sqrt(x)
//! Calculates \f$ \sqrt[n]{[x]} \f$
       interval sqrt    (const interval &, int) throw(STD_FKT_OUT_OF_DEF);  // sqrt(x, n)
//! Calculates \f$ \sqrt{1+[x]^2} \f$
       interval sqrt1px2(const interval&) throw(); // sqrt(1+x^2)
//! Calculates \f$ \sqrt{([x]+1)-1} \f$
       interval sqrtp1m1(const interval&) throw(); // sqrt(x+1)-1
//! Calculates \f$ \sqrt{[x]^2-1} \f$
       interval sqrtx2m1(const interval&);         // sqrt(x^2-1)
//! Calculates \f$ \sqrt{1-[x]^2} \f$
       interval sqrt1mx2(const interval&);         // sqrt(1-x^2)

//! Calculates \f$ \sin([x]) \f$
inline interval sin     (const interval&) throw(); // sin(x)
//! Calculates \f$ \mbox{sin}(\pi*x)/\pi \f$;
       interval sinpix_pi(const interval& x);      // sin(pi*x)/pi
//! Calculates \f$ \cos([x]) \f$
inline interval cos     (const interval&) throw(); // cos(x)
//! Calculates \f$ \tan([x]) \f$
inline interval tan     (const interval&) throw(); // tan(x)
//! Calculates \f$ \cot([x]) \f$
inline interval cot     (const interval&) throw(); // cot(x)

//! Calculates \f$ \arcsin([x]) \f$
inline interval asin    (const interval&);         // ASin(x)
//! Calculates \f$ \arccos([x]) \f$
inline interval acos    (const interval&);         // ACos(x)
//! Calculates \f$ \arctan([x]) \f$
inline interval atan    (const interval&);         // ATan(x) 
//! Calculates \f$ \mbox{arccot}([x]) \f$
inline interval acot    (const interval&);         // ACot(x)

//! Calculates \f$ \exp([x]) \f$
inline interval exp     (const interval&) throw(); // exp(x)
//! Calculates \f$ \exp(-[x]^2) \f$
       interval expmx2  (const interval&);         // exp(-x^2)
//! Calculates \f$ \exp([x])-1 \f$
       interval expm1   (const interval&);         // exp(x)-1
//! Calculates \f$ \exp([x]^2) \f$
       interval expx2   (const interval& x);       // e^{+x^2}
//! Calculates \f$ \exp([x]^2)-1 \f$
       interval expx2m1 (const interval& x);       // e^{+x^2}-1

//! Calculates \f$ \ln([x]) \f$
inline interval ln      (const interval&);         // ln(x)
//! Calculates \f$ \ln(1+[x]) \f$
       interval lnp1    (const interval&) throw(); // ln(1+x)
//! Calculates \f$ \mbox{log}_2([x]) \f$
inline interval log2    (const interval&);         // log2(x)
//! Calculates \f$ \mbox{log}_{10}([x]) \f$
inline interval log10   (const interval&);         // log10(x)

//! Calculates \f$ \sinh([x]) \f$
inline interval sinh    (const interval&) throw(); // Sinh(x)
//! Calculates \f$ \cosh([x]) \f$
inline interval cosh    (const interval&) throw(); // Cosh(x)
//! Calculates \f$ \tanh([x]) \f$
inline interval tanh    (const interval&) throw(); // Tanh(x) 
//! Calculates \f$ \coth([x]) \f$
inline interval coth    (const interval&) throw(); // Coth(x)

//! Calculates \f$ \mbox{arcsinh}([x]) \f$
inline interval asinh   (const interval&);         // ASinh(x)
//! Calculates \f$ \mbox{arccosh}([x]) \f$
inline interval acosh   (const interval&);         // ACosh(x)
//! Calculates \f$ \arccos(1+[x]) \f$
       interval acoshp1 (const interval&);         // acosh(1+x)
//! Calculates \f$ \mbox{arctanh}([x]) \f$
inline interval atanh   (const interval&);         // ATanh(x)
//! Calculates \f$ \mbox{arccoth}([x]) \f$
inline interval acoth   (const interval&);         // ACoth(x)

//! The Gauss error function \f$ \mbox{erf}([x]) = \frac{2}{\sqrt{\pi}} \int \limits_0^{[x]} e^{-t^2} dt \f$ 
interval erf     (const interval&);         // error function
//! The complementary Gauss error function \f$ \mbox{erfc}([x]) = 1 - \mbox{erf}([x]) = \frac{2}{\sqrt{\pi}} \int \limits_{[x]}^\infty e^{-t^2} dt \f$ 
interval erfc    (const interval&);         // complementary error function
//! The Gamma function 
interval gamma (const interval& x);     // Gamma(x)
//! The inverse Gamma function: 1/Gamma(x)
interval gammar(const interval& x);     // 1/Gamma(x)

//! Calculates \f$ [x]^{[y]} \f$
interval pow     (const interval&, const interval&) throw();  // Pow(x,y)
//! Calculates \f$ [x]^n \f$
interval power   (const interval&, int);    // Power(x,n)
//! Calculates \f$ [x]^n \f$
interval Power   (const interval&, int);    // Power(x,n)

//! Calculates \f$ \sqrt{[x]^2+[y]^2} \f$
interval sqrtx2y2(const interval&, const interval&) throw(); // sqrt(x^2+y^2)
//! Calculates \f$ \ln{\sqrt{[x]^2+[y]^2}} \f$
interval ln_sqrtx2y2(const interval&, const interval&) throw();
                                            // ln( sqrt(x^2+y^2) ) == 0.5*ln(x^2+y^2);

//! Enclosure-Interval for \f$ \pi \f$
interval Pi         ( );
} // namespace cxsc 

#include "imath.inl"
#endif // _CXSC_IMATH_HPP_INCLUDED
