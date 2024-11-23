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

/* CVS $Id: l_imath.hpp,v 1.33 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_L_IMATH_HPP_INCLUDED
#define _CXSC_L_IMATH_HPP_INCLUDED

#include "l_interval.hpp"

namespace cxsc {

//! Calculates \f$ [x]^{[y]} \f$
l_interval pow     (const l_interval&, const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF); // Pow(x,y)
//! Calculates \f$ [x]^n \f$
l_interval power   (const l_interval&, int);    // Power(x,n)
//! Calculates \f$ [x]^2  \f$
l_interval sqr     (const l_interval&);         // Sqr(x)

//! Calculates \f$ \sqrt{[x]}  \f$
l_interval sqrt    (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF); 
                                                // Sqrt(x)
//! Calculates \f$ \sqrt[n]{[x]} \f$
l_interval sqrt    (const l_interval&, int) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF); 
                                                // NSqrt(n,x)

//! Calculates \f$ \sin([x]) \f$
l_interval sin     (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW);    
                                                // Sin(x)
//! Calculates \f$ \cos([x]) \f$
l_interval cos     (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW);    
                                                // Cos(x)
//! Calculates \f$ \tan([x]) \f$
l_interval tan     (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW,ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF);
                                                // Tan(x)
//! Calculates \f$ \cot([x]) \f$
l_interval cot     (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW,ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF);
                                                // Cot(x)

//! Calculates \f$ \arcsin([x]) \f$
l_interval asin    (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF);   
                                                // ASin(x)
//! Calculates \f$ \arccos([x]) \f$
l_interval acos    (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF);   
                                                // ACos(x)
//! Calculates \f$ \arctan([x]) \f$
l_interval atan    (const l_interval&) throw(); // ATan(x)
//! Calculates \f$ \mbox{arccot}([x]) \f$
l_interval acot    (const l_interval&) throw(); // ACot(x)

//! Calculates \f$ \exp([x]) \f$
l_interval exp     (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW);    
                                                // exp(x)
//! Calculates \f$ \exp2([x]) \f$
l_interval exp2(const l_interval &); // 2^x

//! Calculates \f$ \exp10([x]) \f$
l_interval exp10(const l_interval &); // 10^x

//! Calculates \f$ \exp([x])-1 \f$
l_interval expm1(const l_interval & x) throw(); // exp(x)-1;
//! Calculates \f$ \exp(-[x]^2) \f$
l_interval expmx2  (const l_interval&);         // e^(-x^2);
//! Calculates \f$ \ln([x]) \f$
l_interval ln      (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF);
                                                // Ln(x)
//! Calculates \f$ \log2([x]) \f$
l_interval log2(const l_interval &);
//! Calculates \f$ \log10([x]) \f$
l_interval log10(const l_interval &);
//! Calculates \f$ \ln(1+[x]) \f$
l_interval lnp1    (const l_interval&) throw(); // ln(1+x)
//! Calculates \f$ \sinh([x]) \f$
l_interval sinh    (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW);   
                                                // Sinh(x)
//! Calculates \f$ \cosh([x]) \f$
l_interval cosh    (const l_interval&) throw(ERROR_LINTERVAL_FAK_OVERFLOW);   
                                                // Cosh(x)
//! Calculates \f$ \tanh([x]) \f$
l_interval tanh    (const l_interval&) throw(); // Tanh(x)
//! Calculates \f$ \coth([x]) \f$
l_interval coth    (const l_interval&) throw(); // Coth(x)           
 
//! Calculates \f$ \mbox{arcsinh}([x]) \f$
l_interval asinh   (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF,ERROR_LINTERVAL_FAK_OVERFLOW);  
                                                // ASinh(x)
//! Calculates \f$ \mbox{arccosh}([x]) \f$
l_interval acosh   (const l_interval&) throw(); // ACosh(x)
//! Calculates \f$ \arccos(1+[x]) \f$
l_interval acoshp1 (const l_interval& x);
//! Calculates \f$ \mbox{arctanh}([x]) \f$
l_interval atanh   (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF,ERROR_LINTERVAL_FAK_OVERFLOW);
                                                // ATanh(x)
//! Calculates \f$ \mbox{arccoth}([x]) \f$
l_interval acoth   (const l_interval&) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF,ERROR_LINTERVAL_FAK_OVERFLOW);  
                                                // ACoth(x)

//! Calculates \f$ \sqrt{1+[x]^2} \f$
l_interval sqrt1px2(const l_interval&) throw(); // Sqrt(1+x^2); 
//! Calculates \f$ \sqrt{[x]^2+[y]^2} \f$
l_interval sqrtx2y2(const l_interval&, const l_interval&) throw(); 
                                                // Sqrt(x^2+y^2); 
//! Calculates \f$ \sqrt{([x]+1)-1} \f$
l_interval sqrtp1m1(const l_interval&) throw(STD_FKT_OUT_OF_DEF);
                                                // sqrtp1m1(x) calculates an inclusion of sqrt(x+1)-1
//! Calculates \f$ \sqrt{[x]^2-1} \f$
l_interval sqrtx2m1(const l_interval&);         // sqrt(x^2-1)
//! Calculates \f$ \sqrt{1-[x]^2} \f$
l_interval sqrt1mx2(const l_interval&);         // sqrt(1-x^2) 

//! Calculates \f$ \ln{\sqrt{[x]^2+[y]^2}} \f$
l_interval ln_sqrtx2y2(const l_interval&, const l_interval&) throw();

// some constants as functions for l_interval
// l_interval li_ln2() throw();                    // ln(2)
// l_interval li_ln10() throw();                   // ln(10)
// l_interval li_Rln10() throw();                  // 1/ln(10)
// l_interval li_pi4() throw();                    // Pi/4
// l_interval li_sqrt2() throw();                  // sqrt(2)
//! Enclosure-Interval for \f$ \ln 2 \f$
l_interval Ln2_l_interval()   throw();   // ln(2) 
//! Enclosure-Interval for \f$ \ln 10 \f$
l_interval Ln10_l_interval()  throw();   // ln(10)
//! Enclosure-Interval for \f$ \frac{1}{\ln 10} \f$
l_interval Ln10r_l_interval() throw();   // 1/ln(10)
//! Enclosure-Interval for \f$ \frac{\pi}{4} \f$
l_interval Pid4_l_interval()  throw();   // Pi/4
//! Enclosure-Interval for \f$ \sqrt{2} \f$
l_interval Sqrt2_l_interval() throw();   // sqrt(2)
//! Enclosure-Interval for \f$ \sqrt{5} \f$
l_interval Sqrt5_l_interval() throw();   // sqrt(5)
//! Enclosure-Interval for \f$ \sqrt{7} \f$
l_interval Sqrt7_l_interval() throw();   // sqrt(7)

// obsolete functions, for compability
//! Enclosure-Interval for \f$ \ln 2 \f$
inline l_interval li_ln2()  {return Ln2_l_interval();}    // ln(2)
      //! Enclosure-Interval for \f$ \ln 10 \f$
inline l_interval li_ln10() {return Ln10_l_interval();}   // ln(10)
      //! Enclosure-Interval for \f$ \frac{1}{\ln 10} \f$
inline l_interval li_Rln10(){return Ln10r_l_interval();}  // 1/ln(10)
//! Enclosure-Interval for \f$ \frac{\pi}{4} \f$
inline l_interval li_pi4()  {return Pid4_l_interval();}   // Pi/4
//! Enclosure-Interval for \f$ \sqrt{2} \f$
inline l_interval li_sqrt2(){return Sqrt2_l_interval();}  // sqrt(2)

      //! Enclosure-Interval for \f$ \frac{1}{\ln 2} \f$
l_interval Ln2r_l_interval() throw();     // 1/ln(2)
      //! Enclosure-Interval for \f$ \pi \f$
l_interval Pi_l_interval() throw();       // Pi
      //! Enclosure-Interval for \f$ \frac{\pi}{2} \f$
l_interval Pid2_l_interval() throw();     // Pi/2
      //! Enclosure-Interval for \f$ 2\pi \f$
l_interval Pi2_l_interval() throw();      // 2*Pi
      //! Enclosure-Interval for \f$ \frac{\pi}{3} \f$
l_interval Pid3_l_interval() throw();     // Pi/3
      //! Enclosure-Interval for \f$ \frac{1}{\pi} \f$
l_interval Pir_l_interval() throw();      // 1/Pi
      //! Enclosure-Interval for \f$ \frac{1}{2\pi} \f$
l_interval Pi2r_l_interval() throw();     // 1/(2*Pi)
      //! Enclosure-Interval for \f$ \sqrt{\pi} \f$
l_interval SqrtPi_l_interval() throw();   // sqrt(Pi)
      //! Enclosure-Interval for \f$ \sqrt{2\pi} \f$
l_interval Sqrt2Pi_l_interval() throw();  // sqrt(2*Pi)
      //! Enclosure-Interval for \f$ \frac{1}{\sqrt{\pi}} \f$
l_interval SqrtPir_l_interval() throw();  // 1/sqrt(Pi)
      //! Enclosure-Interval for \f$ \frac{1}{\sqrt{2\pi}} \f$
l_interval Sqrt2Pir_l_interval() throw(); // 1/sqrt(2*Pi)
      //! Enclosure-Interval for \f$ 2^\pi \f$
l_interval Pip2_l_interval() throw();     // Pi^2
      //! Enclosure-Interval for \f$ \frac{1}{\sqrt{2}} \f$
l_interval Sqrt2r_l_interval() throw();   // 1/sqrt(2)
      //! Enclosure-Interval for \f$ \sqrt{3} \f$
l_interval Sqrt3_l_interval() throw();    // sqrt(3)
      //! Enclosure-Interval for \f$ \frac{\sqrt{3}}{2} \f$
l_interval Sqrt3d2_l_interval() throw();  // sqrt(3)/2
      //! Enclosure-Interval for \f$ \frac{1}{\sqrt{3}} \f$
l_interval Sqrt3r_l_interval() throw();   // 1/sqrt(3)
      //! Enclosure-Interval for \f$ \ln \pi \f$
l_interval LnPi_l_interval() throw();     // ln(Pi)
      //! Enclosure-Interval for \f$ \ln 2\pi \f$
l_interval Ln2Pi_l_interval() throw();    // ln(2*Pi)
      //! Enclosure-Interval for \f$ e \f$
l_interval E_l_interval() throw();        // e = exp(1)
      //! Enclosure-Interval for \f$ \frac{1}{e} \f$
l_interval Er_l_interval() throw();       // 1/e
      //! Enclosure-Interval for \f$ e^2 \f$
l_interval Ep2_l_interval() throw();      // e^2
      //! Enclosure-Interval for \f$ \frac{1}{e^2} \f$
l_interval Ep2r_l_interval() throw();     // 1/e^2
      //! Enclosure-Interval for \f$ e^\pi \f$
l_interval EpPi_l_interval() throw();     // e^Pi
      //! Enclosure-Interval for \f$ e^{2\pi} \f$
l_interval Ep2Pi_l_interval() throw();    // e^(2*Pi)
      //! Enclosure-Interval for \f$ e^{\frac{\pi}{2}} \f$
l_interval EpPid2_l_interval() throw();   // e^(Pi/2)
      //! Enclosure-Interval for \f$ e^{\frac{\pi}{4}} \f$
l_interval EpPid4_l_interval() throw();   // e^(Pi/4)
      //! Enclosure-Interval for Euler Gamma
l_interval EulerGa_l_interval() throw();  // EulerGamma
      //! Enclosure-Interval for Catalan Numbers
l_interval Catalan_l_interval() throw();  // Catalan
} // namespace cxsc 

#endif // _CXSC_L_IMATH_HPP_INCLUDED
