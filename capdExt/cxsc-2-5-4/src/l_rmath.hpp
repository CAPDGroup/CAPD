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

/* CVS $Id: l_rmath.hpp,v 1.32 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_L_RMATH_HPP_INCLUDED
#define _CXSC_L_RMATH_HPP_INCLUDED

#include "l_real.hpp"

namespace cxsc {

//! Calculates \f$ x^2  \f$
inline l_real sqr     (const l_real&) throw(); // Sqr(x)
//! Calculates \f$ \sqrt{x}  \f$
       l_real sqrt    (const l_real&) throw(ERROR_LREAL_STD_FKT_OUT_OF_DEF);
                                               // Sqrt(x)
//! Calculates \f$ \sqrt{x^2+y^2} \f$
       l_real sqrtx2y2(const l_real&, const l_real&) throw(); 
                                               // Sqrt(x^2+y^2)
//! Calculates \f$ \sqrt{1+x^2} \f$
       l_real sqrt1px2(const l_real&) throw(); // Sqrt(1+x^2)
//! Calculates \f$ \sqrt{(x+1)-1} \f$
inline l_real sqrtp1m1(const l_real &) throw(); 
//! Calculates \f$ \sqrt{x^2-1} \f$
inline l_real sqrtx2m1(const l_real &) throw();
//! Calculates \f$ \sqrt{1-x^2} \f$
inline l_real sqrt1mx2(const l_real &) throw();
//! Calculates \f$ \exp(x)-1 \f$
inline l_real expm1   (const l_real &x) throw();
//! Calculates \f$ \exp(-x^2) \f$
inline l_real expmx2  (const l_real&) throw();
//! Calculates \f$ \ln{\sqrt{x^2+y^2}} \f$
inline l_real ln_sqrtx2y2(const l_real& x, const l_real& y) throw();
//! Calculates \f$ \arccos(1+x) \f$
inline l_real acoshp1 (const l_real& x);

// inline l_real sqrt   (const l_real &, int);  // Sqrt(x, n)
// inline l_real sin    (const l_real&) throw();        // Sin(x)
// inline l_real cos    (const l_real&) throw();        // Cos(x)
// inline l_real tan    (const l_real&) throw();        // Tan(x)
// inline l_real cot    (const l_real&) throw();        // Cot(x)
// inline l_real asin   (const l_real&);        // ASin(x)
// inline l_real acos   (const l_real&);        // ACos(x)
// inline l_real atan   (const l_real&);        // ATan(x) 
// inline l_real acot   (const l_real&);        // ACot(x)
// inline l_real exp    (const l_real&) throw();        // Exp(x)
// inline l_real ln     (const l_real&);        // Ln(x)
// inline l_real sinh   (const l_real&) throw();        // Sinh(x)
// inline l_real cosh   (const l_real&) throw();        // Cosh(x)
// inline l_real tanh   (const l_real&) throw();        // Tanh(x) 
// inline l_real coth   (const l_real&) throw();        // Coth(x)
// inline l_real asinh  (const l_real&);        // ASinh(x)
// inline l_real acosh  (const l_real&);        // ACosh(x)
// inline l_real atanh  (const l_real&);        // ATanh(x)
// inline l_real acoth  (const l_real&);        // ACoth(x)

//! Calculates \f$ x^y \f$
inline l_real pow    (const l_real&, const l_real&); // Pow(x,y)
//! Calculates \f$ x^n \f$
l_real power         (const l_real&, const int);     // Power(x,n)

// real staggered constants (the same as in l_interval.hpp):
l_real Ln2_l_real()   throw();   // ln(2) 
l_real Ln10_l_real()  throw();   // ln(10)
l_real Ln10r_l_real() throw();   // 1/ln(10)
l_real Pid4_l_real()  throw();   // Pi/4
l_real Sqrt2_l_real() throw();   // sqrt(2)
l_real Sqrt5_l_real() throw();   // sqrt(5)
l_real Sqrt7_l_real() throw();   // sqrt(7)
l_real Ln2r_l_real() throw();     // 1/ln(2)
l_real Pi_l_real() throw();       // Pi
l_real Pid2_l_real() throw();     // Pi/2
l_real Pi2_l_real() throw();      // 2*Pi
l_real Pid3_l_real() throw();     // Pi/3
l_real Pir_l_real() throw();      // 1/Pi
l_real Pi2r_l_real() throw();     // 1/(2*Pi)
l_real SqrtPi_l_real() throw();   // sqrt(Pi)
l_real Sqrt2Pi_l_real() throw();  // sqrt(2*Pi)
l_real SqrtPir_l_real() throw();  // 1/sqrt(Pi)
l_real Sqrt2Pir_l_real() throw(); // 1/sqrt(2*Pi)
l_real Pip2_l_real() throw();     // Pi^2
l_real Sqrt2r_l_real() throw();   // 1/sqrt(2)
l_real Sqrt3_l_real() throw();    // sqrt(3)
l_real Sqrt3d2_l_real() throw();  // sqrt(3)/2
l_real Sqrt3r_l_real() throw();   // 1/sqrt(3)
l_real LnPi_l_real() throw();     // ln(Pi)
l_real Ln2Pi_l_real() throw();    // ln(2*Pi)
l_real E_l_real() throw();        // e = exp(1)
l_real Er_l_real() throw();       // 1/e
l_real Ep2_l_real() throw();      // e^2
l_real Ep2r_l_real() throw();     // 1/e^2
l_real EpPi_l_real() throw();     // e^Pi
l_real Ep2Pi_l_real() throw();    // e^(2*Pi)
l_real EpPid2_l_real() throw();   // e^(Pi/2)
l_real EpPid4_l_real() throw();   // e^(Pi/4)
l_real EulerGa_l_real() throw();  // EulerGamma
l_real Catalan_l_real() throw();  // Catalan

} // namespace cxsc 

#include "l_rmath.inl"
#endif // _CXSC_L_RMATH_HPP_INCLUDED
