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

/* CVS $Id: l_rmath.cpp,v 1.28 2014/01/30 17:23:46 cxsc Exp $ */

#include <l_rmath.hpp>
#include <rmath.hpp>
#include <l_interval.hpp>

namespace cxsc {

l_real sqrt(const l_real &x) throw(ERROR_LREAL_STD_FKT_OUT_OF_DEF)
// Blomquist, additional scaling, 10.12.02
{
   int stagsave = stagprec, stagmax = 19, stagcalc;
   l_real y;
   if (sign(x[1])<0) // Blomquist: is faster than  x < 0.0;
   cxscthrow(ERROR_LREAL_STD_FKT_OUT_OF_DEF("l_real sqrt(const l_real &x)"));
   else if (zero_(x) || x == 1.0 ) // Blomquist, faster than x == 0;
      y = x;
   else 
   {
       l_real x1 = x;
       int ex = expo(x1[1]);
       ex = 1021-ex;        // ex: optimal scaling factor
       if (ex%2) ex--;      // ex is now even;
       times2pown(x1,ex);   // scaling x1 with 2^ex

       // Einsatz des Newton-Verfahrens y = y-f(y)/f'(y)
       if (stagprec < stagmax) stagcalc = stagprec+1;  
       else stagcalc = stagmax+1;  // Blomquist with +1; 30.11.02;
       y = sqrt(_real(x1));
       stagprec = 1; 
       while (stagprec < stagcalc) 
       {
	   stagprec += stagprec; // now calculation in double precision:
           if (stagprec > stagcalc) stagprec = stagcalc;
           // y = 0.5*((x1/y)+y);
	   y += x1/y;  // Blomquist, 30.11.02;
	   times2pown(y,-1); // Blomquist, this multiplication is faster!!
       }
       times2pown(y,-ex/2); // backscaling
       stagprec = stagsave; // restore old stagprec
       adjust(y);
   }
   return y;
}

l_real sqrtx2y2(const l_real& x, 
                       const l_real& y) throw() // Blomquist 6.12.02
// Calculation of an approximation of sqrt(x^2+y^2).
// In general the maximum precision is stagprec=19, predifined by the used
// sqrt-function declared in l_rmath.hpp.
// If the difference |exa-exb| of the exponents (base 2) is sufficient high,
// precision and accuracy can be choosen greater 19.
{ 
    int stagsave = stagprec, stagcalc, stagmax = 19;
    l_real a,b,r,r1;
    int exa,exb,ex;
    a = x;  b = y;
    exa = expo(a[1]);
    exb = expo(b[1]);
    if (exb > exa)
    {  // Permutation of a,b:
	r = a;  a = b;  b = r;
	ex = exa;  exa = exb;  exb = ex;
    } // |a| >= |b|
    if (sign(a[1]) < 0) a = -a;  // a == abs(a);
    if (sign(b[1]==0)) return a;

    // |a| = a >= |b| > 0:
    ex = 0;  // initialization
    if (6*exb < 5*exa-1071)
    {   // approximation with a + 0.5*b^2/a - b^4/(8*a^3), if the next
        // Taylor term b^6/(16*a^5) is less than 2^-1074==minreal;
        // minreal is the smallest positive number in IEEE-format.
	r1 = (b/a);
	r = r1*r1;  
	times2pown(r,-2);   // r = r/4;
	r = 1 - r;
	r1 *= b;
	times2pown(r1,-1);  // r = r/2;
	r *= r1;
	r += a; // r = a + 0.5*b^2/a - b^4/(8*a^3);
    } else
    {   // Scaling ubwards or downwards to achiev optimal accuracy and 
        // to avoid overflow by accu-reading:
        stagcalc = stagprec;  
	if (stagcalc > stagmax) stagcalc = stagmax; // using sqrt(...) a 
	stagprec = stagcalc; // higher precision than stagmax makes no sense!

	ex = 511 - exa; // scaling factor to avoide overflow by accu-reading
	if (ex < 0)
	{   // sqrt(a^2+b^2) = a*sqrt( 1+(b/a)^2 )
	    r = b/a;
	    r = a*sqrt(1+r*r);  // sqrt(...) declared in l_rmath.hpp
	} else
	{   
            times2pown(a,ex);   // exact scaling with ex >= 0
	    times2pown(b,ex);   // exact scaling eith ex >= 0
	    dotprecision dot(0.0);
	    accumulate(dot,a,a);
	    accumulate(dot,b,b);
	    r = dot; // r with no overflow!
	    r = sqrt(r);  // sqrt(...) declared in l_rmath.hpp
	    times2pown(r,-ex);  // back-scaling 
	}
	stagprec = stagsave;  // restore old stagprec
    }
    return r;
} // sqrtx2y2(...)

l_real sqrt1px2(const l_real& x) throw()
// Inclusion of sqrt(1+x^2); Blomquist, 13.12.02;
// With stagmax=19 we get about 16*19=304 exact decimal digits.
{
    l_real y,t;
    int stagsave, stagmax=19;
    stagsave = stagprec;
    if (stagprec > stagmax) stagprec = stagmax;
    if (expo(x[1]) > 260)
    {  // sqrt(1+x^2) = |x| + 1/(2*|x|)
	y = abs(x);
	t = 1/y;
	times2pown(t,-1);
	y += t;
    } else y = sqrt(1+x*x);
    stagprec = stagsave;
    y = adjust(y);
    return y;
}

l_real power(const l_real& x, int n)
{
   int       stagsave = stagprec, stagmax = 19;
   long int  zhi = 2;
   l_real    y, neu;
   bool      neg=false;

   if (x == 1.0)
      y = x;
   else if (n == 0)
      y = adjust(_l_real(1.0));
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else
         stagprec = stagmax;

      if (n == 1)
         y = x;
      else if (n == 2)
         y = x*x;
      else 
      {
         if (n < 0) 
         {
            neg = true;
            n = -n;
         }

         // Initialisierung
         if (n%2)   
            y = x;
         else      
            y = 1.0;  // Praezision wird bei 1 Mult. auf
                      // aktuellen Wert gesetzt;
                      // Berechnung durch binaere Darstellung der n
         neu = x*x;
         do {
            if ((n/zhi)%2)
               y *= neu;
            zhi += zhi;
            if (zhi <= n)  // letzte Mult. entfaellt --> schneller
               neu *= neu;
         } while (zhi <= n);
         if (neg) 
            y = 1/(y);
      }
      stagprec = stagsave;
      y = adjust(y);
   }
   return y;
}

// real staggered constants (the same as in l_interval.hpp):
l_real Ln2_l_real()   throw()   // ln(2) 
{ return mid( Ln2_l_interval() ); }
l_real Ln10_l_real()  throw()   // ln(10)
{ return mid( Ln10_l_interval() ); }
l_real Ln10r_l_real() throw()   // 1/ln(10)
{ return mid( Ln10r_l_interval()); }
l_real Pid4_l_real()  throw()   // Pi/4
{ return mid( Pid4_l_interval() ); }
l_real Sqrt2_l_real() throw()   // sqrt(2)
{ return mid( Sqrt2_l_interval() ); }
l_real Sqrt5_l_real() throw()   // sqrt(5)
{ return mid( Sqrt5_l_interval() ); }
l_real Sqrt7_l_real() throw()   // sqrt(7)
{ return mid( Sqrt7_l_interval() ); }
l_real Ln2r_l_real() throw()     // 1/ln(2)
{ return mid( Ln2r_l_interval() ); }
l_real Pi_l_real() throw()       // Pi
{ return mid( Pi_l_interval() ); }
l_real Pid2_l_real() throw()     // Pi/2
{ return mid( Pid2_l_interval() ); }
l_real Pi2_l_real() throw()      // 2*Pi
{ return mid( Pi2_l_interval() ); }
l_real Pid3_l_real() throw()     // Pi/3
{ return mid( Pid3_l_interval() ); }
l_real Pir_l_real() throw()      // 1/Pi
{ return mid( Pir_l_interval() ); }
l_real Pi2r_l_real() throw()     // 1/(2*Pi)
{ return mid( Pi2r_l_interval() ); }
l_real SqrtPi_l_real() throw()   // sqrt(Pi)
{ return mid( SqrtPi_l_interval() ); }
l_real Sqrt2Pi_l_real() throw()  // sqrt(2*Pi)
{ return mid( Sqrt2Pi_l_interval() ); }
l_real SqrtPir_l_real() throw()  // 1/sqrt(Pi)
{ return mid( SqrtPir_l_interval() ); }
l_real Sqrt2Pir_l_real() throw() // 1/sqrt(2*Pi)
{ return mid( Sqrt2Pir_l_interval() ); }
l_real Pip2_l_real() throw()     // Pi^2
{ return mid( Pip2_l_interval() ); }
l_real Sqrt2r_l_real() throw()   // 1/sqrt(2)
{ return mid( Sqrt2r_l_interval() ); }
l_real Sqrt3_l_real() throw()    // sqrt(3)
{ return mid( Sqrt3_l_interval() ); }
l_real Sqrt3d2_l_real() throw()  // sqrt(3)/2
{ return mid( Sqrt3d2_l_interval() ); }
l_real Sqrt3r_l_real() throw()   // 1/sqrt(3)
{ return mid( Sqrt3r_l_interval() ); }
l_real LnPi_l_real() throw()     // ln(Pi)
{ return mid( LnPi_l_interval() ); }
l_real Ln2Pi_l_real() throw()    // ln(2*Pi)
{ return mid( Ln2Pi_l_interval() ); }
l_real E_l_real() throw()        // e = exp(1)
{ return mid( E_l_interval() ); }
l_real Er_l_real() throw()       // 1/e
{ return mid( Er_l_interval() ); }
l_real Ep2_l_real() throw()      // e^2
{ return mid( Ep2_l_interval() ); }
l_real Ep2r_l_real() throw()     // 1/e^2
{ return mid( Ep2r_l_interval() ); }
l_real EpPi_l_real() throw()     // e^Pi
{ return mid( EpPi_l_interval() ); }
l_real Ep2Pi_l_real() throw()    // e^(2*Pi)
{ return mid( Ep2Pi_l_interval() ); }
l_real EpPid2_l_real() throw()   // e^(Pi/2)
{ return mid( EpPid2_l_interval() ); }
l_real EpPid4_l_real() throw()   // e^(Pi/4)
{ return mid( EpPid4_l_interval() ); }
l_real EulerGa_l_real() throw()  // EulerGamma
{ return mid(EulerGa_l_interval()  ); }
l_real Catalan_l_real() throw()  // Catalan
{ return mid( Catalan_l_interval() ); }

} // namespace cxsc

