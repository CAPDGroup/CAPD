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

/* CVS $Id: rmath.hpp,v 1.36 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_RMATH_HPP_INCLUDED
#define _CXSC_RMATH_HPP_INCLUDED

#include "real.hpp"
#include <dot.hpp> // Blomquist

namespace cxsc {

//! Calculates \f$ x^2  \f$
inline real sqr    (const real&) throw(); // Sqr(x)
//! Calculates \f$ \sqrt{x}  \f$
inline real sqrt   (const real&);         // Sqrt(x)
//! Calculates \f$ \sqrt[n]{x} \f$
inline real sqrt   (const real &, int);   // Sqrt(x, n)
//! Calculates \f$ \sqrt{1+x^2} \f$
       real sqrt1px2(const real&) throw();  // Sqrt(1+x^2); Blomquist 13.12.02
//! Calculates \f$ \sqrt{(x+1)-1} \f$
inline real sqrtm1  (const real&); // Sqrt(x+1)-1; ohne Fehlerabsch�zung!
//! Calculates \f$ \sqrt{(x+1)-1} \f$
       real sqrtp1m1(const real&) throw(); // Blomquist 05.08.03
//! Calculates \f$ \sqrt{x^2-1} \f$
       real sqrtx2m1(const real&);
//! Calculates \f$ \sqrt{1-x^2} \f$
       real sqrt1mx2(const real&) throw(STD_FKT_OUT_OF_DEF);
              

//! Calculates \f$ \sin(x) \f$
inline real sin    (const real&) throw();        // Sin(x)
//! Calculates \f$ \mbox{sin}(\pi*x)/\pi \f$;
       real sinpix_pi(const real& x);            // sin(pi*x)/pi;
//! Calculates \f$ \cos(x) \f$
inline real cos    (const real&) throw();        // Cos(x)
//! Calculates \f$ \tan(x) \f$
inline real tan    (const real&) throw();        // Tan(x)
//! Calculates \f$ \cot(x) \f$
inline real cot    (const real&) throw();        // Cot(x)

//! Calculates \f$ \arcsin(x) \f$
inline real asin   (const real&);                // ASin(x)
//! Calculates \f$ \arccos(x) \f$
inline real acos   (const real&);                // ACos(x)
//! Calculates \f$ \arctan(x) \f$
inline real atan   (const real&);                // ATan(x) 
//! Calculates \f$ \mbox{arccot}(x) \f$
inline real acot   (const real&);                // ACot(x)

//! Calculates \f$ \exp(x) \f$
inline real exp    (const real&) throw();        // Exp(x)
//! Calculates \f$ \exp(x)-1 \f$
inline real expm1  (const real&) throw();        // Exp(x)-1
//! Calculates \f$ \exp(-x^2) \f$
       real expmx2 (const real&) throw();        // Exp(-x^2)
//! Calculates \f$ \exp(x^2) \f$
       real expx2  (const real& x);              // e^{+x^2}
//! Calculates \f$ \exp(x^2)-1 \f$
       real expx2m1(const real& x);              // e^{+x^2}-1

//! Calculates \f$ \ln(1+x) \f$
inline real lnp1   (const real&);                // Ln(1+x)
//! Calculates \f$ \ln(x) \f$
inline real ln     (const real&);                // Ln(x)
//! Calculates \f$ \mbox{log}_2(x) \f$
inline real log2     (const real&);              // Log2(x)
//! Calculates \f$ \mbox{log}_{10}(x) \f$
inline real log10     (const real&);             // Log10(x)

//! Calculates \f$ \sinh(x) \f$
inline real sinh   (const real&) throw();        // Sinh(x)
//! Calculates \f$ \cosh(x) \f$
inline real cosh   (const real&) throw();        // Cosh(x)
//! Calculates \f$ \arccos(1+x) \f$
       real acoshp1(const real& x) throw();      // acosh(1+x)
//! Calculates \f$ \tanh(x) \f$
inline real tanh   (const real&) throw();        // Tanh(x) 
//! Calculates \f$ \coth(x) \f$
inline real coth   (const real&) throw();        // Coth(x)

//! Calculates \f$ \mbox{arcsinh}(x) \f$
inline real asinh  (const real&);                // ASinh(x)
//! Calculates \f$ \mbox{arccosh}(x) \f$
inline real acosh  (const real&);                // ACosh(x)
//! Calculates \f$ \mbox{arctanh}(x) \f$
inline real atanh  (const real&);                // ATanh(x)
//! Calculates \f$ \mbox{arccoth}(x) \f$
inline real acoth  (const real&);                // ACoth(x)
//! The Gauss error function \f$ \mbox{erf}(x) = \frac{2}{\sqrt{\pi}} \int \limits_0^x e^{-t^2} dt \f$ 
inline real erf    (const real&);                // error function
//! The complementary Gauss error function \f$ \mbox{erfc}(x) = 1 - \mbox{erf}(x) = \frac{2}{\sqrt{\pi}} \int \limits_x^\infty e^{-t^2} dt \f$ 
inline real erfc   (const real&);                // complementary error function
//! The Gamma function 
       real gamma (const real& x);       // Gamma(x)
//! The inverse Gamma function: 1/Gamma(x)
       real gammar(const real& x);       // 1/Gamma(x)
  
//! Calculates \f$ x^y \f$
inline real pow    (const real&, const real&);  // Pow(x,y)
//! Calculates \f$ x^n \f$
inline real power  (const real&, const int);    // Power(x,n)

//! Calculates \f$ \sqrt{x^2+y^2} \f$
       real sqrtx2y2(const real&, const real&) throw(); // Sqrt(x^2+y^2)
//! Calculates \f$ \ln{\sqrt{x^2+y^2}} \f$
       real ln_sqrtx2y2(const real&, const real&) throw(STD_FKT_OUT_OF_DEF);  
                                                 // ln( sqrt(x^2+y^2) )

//! Returns a real value, which corresponds with the first 24 mantissa bits of x
real Cut24(const real&);
//! Returns a real value, which corresponds with the first 25 mantissa bits of x
real Cut25(const real&);
//! Returns a real value, which corresponds with the first 26 mantissa bits of x
real Cut26(const real&);
//! Rouding to the next integer; |x| < 2147483647.5
int Round(const real& x) throw(); 
//! Rounding to the smallest integer greater or equal x; -2147483649 < x <= 2147483647.0;
int ceil(const real& x) throw();
//! Rounding to the greates integer smaller or equal x; -2147483649 < x <= 2147483647.0;
int ifloor(const real& x) throw();

extern "C" {
   void r_lfsr(void); // Siehe real.hpp in real_ari...?!?!
}

} // namespace cxsc 

#include "rmath.inl"
#endif // _CXSC_RMATH_HPP_INCLUDED

