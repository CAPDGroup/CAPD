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

/* CVS $Id: l_imath.cpp,v 1.38 2014/01/30 17:23:46 cxsc Exp $ */

#include <math.h>
#include "l_imath.hpp"
#include "imath.hpp"
#include "rmath.hpp"

namespace cxsc {

#define  CXSC_One       1.
#define  CXSC_Zero      0.
#define  CXSC_MinusOne  -1.

l_interval pow(const l_interval & x, const l_interval & e) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF)
{
   int         stagsave = stagprec,
               stagmax = 19,
               intexp;
               
   bool        fertig;

   l_interval  y;
   interval    dx = interval(x),
               de = interval(e),
               einfachgenau;
   real        supabsde = Sup(abs(de));

   einfachgenau = pow(dx,de);

   fertig = false;
   if (Inf(de) == Sup(de)) // &&
      if (supabsde < 32768.0) 
      {
         intexp = int(_double(real(Sup(e))));
         if (real(intexp) == Sup(e)) 
         {
            y = power(x,intexp);   // Integer-Potenz wesentlich schneller
            fertig = true;
         }
      }
      
   if (!fertig) 
   {
      if (Inf(x) < 0.0) 
         cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval pow(const l_interval & x, const l_interval & e)"));
      else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_One)
         y = x;
      else if (Inf(de) == Sup(de) && Sup(de) == CXSC_One)
         y = x;
      else if (Inf(de) == Sup(de) && Sup(de) == CXSC_Zero)
         y = 1.0;
      else 
      {
         if (stagprec < stagmax) 
            stagprec++;
         else                    
            stagprec = stagmax;
         y = exp(e*ln(x));
         stagprec = stagsave;
         y = adjust(y);
         y = y & einfachgenau;
      }
   }

   return y;
}

l_interval power(const l_interval &x, int n)       // Power(x,n)
{
   int         stagsave = stagprec,
               stagmax = 19;
              
   bool        neg = false;
              
   long int    zhi = 2;
   interval    dx = interval(x),
               einfachgenau;
   l_interval  y, neu;

   einfachgenau = Power(dx,n);

   if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_One)
      y = x;
   else if (n == 0)
      y = adjust(l_interval(1.0));
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;
      
      if (n == 1)
         y = x;
      else if (n == 2)
         y = sqr(x);
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
            y = l_interval(1.0);  // Praezision wird bei 1 Mult. auf
                                        // aktuellen Wert gesetzt;
         // Berechnugn durch binaere Darstellung der n
         neu = sqr(x);   // neu = x*x;
         do {
            if ((n/zhi)%2)  y *= neu;
            zhi += zhi;
            if (zhi <= n)  // letzte Mult. entfaellt --> schneller
               neu *= neu;
         } while (zhi <= n);

         if (neg) 
            y = 1.0/(y);
      }
      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;
   }

   return y;
}

l_interval sqr(const l_interval & x)              // Sqr(x)
{
   l_interval	y;

   if (Inf(x) >= 0.0)       /* result = [x.inf*x.inf, x.sup*x.sup] */
      y = x * x;
   else if (Sup(x)<= 0.0) 
   { 
      /* result = [x.sup*x.sup, x.inf*x.inf] */
      y = l_interval (-Sup(x) , -Inf(x));
      y = (y) * (y);
   } else 
   {                   
      /* result = [0.0, max(x.sup*x.sup,x.inf*x.inf)] */
      if (abs(Inf(x)) >= abs(Sup(x)))
         y = l_interval(0.0, abs(Inf(x)));
      else
         y = l_interval(0.0, abs(Sup(x)));
      y = (y) * (y);
   }
	return y;
}

l_interval sqrt(const l_interval & x) 
                         throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF) // Sqrt(x)
{  // Blomquist: scaling with 2^ex is necessary if expo(Sup(dx)) is too small!
   int         stagsave = stagprec,
               stagmax = 30,
               stagsave2;
   interval    dx = interval(x),
               einfachgenau;
   l_interval  a1,y,t,mt;
   bool Inf_Zero;  

   einfachgenau = sqrt(dx);

   if (Inf(x) < 0.0) 
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval sqrt(const l_interval & x)"));
   else if (Inf(dx) == Sup(dx) && 
                               (Sup(dx) == CXSC_Zero || Sup(dx) == CXSC_One))
      y = x;
   else 
   {   // scaling necessary if exponent ex < 0
       l_interval x1 = x;  
       Inf_Zero = (Inf(dx)==0);    
       if (Inf_Zero) x1 = Sup(x1); 
       int ex = expo(Sup(dx));
       if (ex>0) ex = 0;  // scaling unnecessary if ex>0
       else
       {
	   ex = -ex; // ex >= 0
	   if (ex > 1023) ex = 1023; // ex==1023 is sufficient
	   if (ex%2) ex--;  // ex>=0 is even
       }
       if (ex) times2pown(x1,ex);  // ex >= 0  --->  exact scaling!
      // Interval-Newton-methode: y = m(y)-f(m(y))/f'(y)
      t = sqrt(interval(x1));
      if (stagprec < stagmax) 
         stagsave2 = stagprec+1;
      else                    
         stagsave2 = stagmax;
      stagprec = 1;
      while (stagprec < stagsave2) 
      { 
         stagprec += stagprec;
         if (stagprec > stagmax) 
            stagprec = stagmax;
         mt = mid(t);  times2pown(t,1);
         t = mt-((mt*mt-x1)/t);
      }
      if (ex) times2pown(t,-ex/2); // ex!=0 --> backscaling with 2^(-ex/2)
      stagprec = stagsave;   // restore the previous precision
      y = adjust(t);         // matching to previous stagprec.
      if (Inf_Zero) SetInf(y,0.0); 
      y = y & einfachgenau;  // seeking optimal inclusion with intersection
  }
  return y;
}

l_interval sqrt(const l_interval &x, int n) 
                throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF) // Sqrt(x,n)
// -1073741822 <= n <= +1073741823, sonst autom. Programm-Abbruch
// sqrt(x,n) jetzt mit Skalierung --> Hohe Genauigkeit in allen Bereichen!
// Blomquist, 28.12.03;
{  
    int stagsave = stagprec,
	stagshort, staglong, i, ex, N,
	stagmax = 19,
	max = 2;
    l_interval  my, fy, corr, y, xx;
    interval dx = interval(x), einfachgenau;

    einfachgenau = sqrt(dx,n);

    if (Inf(x) < 0.0) 
	cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
                        ("l_interval sqrt(const l_interval &x, int n)"));
    else if (stagprec == 1)
	y = pow(dx, interval(1.0)/n);
    else if (Inf(dx) == Sup(dx) && (Sup(dx) == CXSC_Zero || 
                                                   Sup(dx) == CXSC_One))
	y = x;
    else 
    {
	if (stagprec < stagmax) 
	    stagprec++;
	else stagprec = stagmax;
         
	while (max < stagprec) 
	    max += max;  // quadratische Konvergenz kann
                         // nicht eingehalten werden,
	max += max;      // deshalb eine Schleife mehr

	xx = x;
	ex = expo(Sup(dx));
	N = -ex;
        // Skalierung mit 2^N, so dass dx,xx etwa = 1; ----------------
	times2pown(dx,N);
	if (N>1023) { 
	    times2pown(xx,1023);
	    times2pown(xx,N-1023);
	}
	else times2pown(xx,N);
        // Skalierung beendet, Blomquist 28.12.03 --------------------- 

	y = pow(dx, interval(1.0/n));
	stagprec = 1;
	// Intervall-Newton-Verfahren
	for (i = 2; i <= max; i += i) 
	{
	    // Verdoppelung der Genauigkeit:
	    stagshort = stagprec;
	    stagprec += stagprec; // Verdoppelung hier
	    if (stagprec > stagmax) 
		stagprec = stagmax;
	    staglong  = stagprec;
	    my = l_interval(mid(y));
	    fy = power(my, n)-xx;
	    // Fehlerauswertung nur in halber Genauigkeit notwendig!
	    stagprec  = stagshort;
	    corr = fy/(real(n)*power(y, n-1));
	    stagprec  = staglong;
	    // Fehlerkorrektur in normaler Genauigkeit:
	    y = my-corr;
	}
        // Rueckskalierung mit dem Faktor 2^(-N/n):
	fy = l_interval(-N)/n;
	y *= exp(fy*li_ln2());  // li_ln2() = ln(2)
	    
	stagprec = stagsave;
	y = adjust(y);  // Anpassung an Ausgangs-stagprec = stagsave.
	y = y & einfachgenau;  // Falls y breiter ist als einfachgenau
    }

    return y;
} // sqrt(x,n)

l_interval sqrt1px2(const l_interval& x) throw()
// Calculation of an optimal inclusion of sqrt(1+x^2); Blomquist, 13.12.02;
// With stagmax=19 we get about 16*19=304 exact decimal digits.
{   // First step: Calculation of sqrt(1+x*x) in simple type interval:
    interval einfach = sqrt1px2( interval(x) );  // Only interval types
    int stagsave, stagmax=19;
    stagsave = stagprec;
    if (stagprec > stagmax) stagprec = stagmax;
    // Second step: Inclusion of sqrt(1+x^2) with stagprec <= 19
    const int exmax=512;
    l_interval y = abs(x);
    l_real lr = Sup( 1 + l_interval(Sup(y)) );
    interval z = interval(Sup(x));
    int ex = expo(Sup(z));
    if (ex > exmax) 
    {   // scaling to avoid overflow by sqr(y):
	ex = exmax - ex;    // ex = 512 - ex;
	times2pown(y,ex);   // scaling to avoid overflow by sqr(y)
	y = sqrt( comp(0.5,2*ex+1) + sqr(y) ); // sqr(y) without overflow!
	times2pown(y,-ex);  // backscaling of the result with 2^(-ex)
    } else
    { // no scaling necessary:
        y = sqrt(1+sqr(x));
    }
    if (Inf(y)<1.0) SetInf(y,1.0); // improvement, Blomquist, 25.02.07;
    if (Sup(y)>lr) SetSup(y,lr);   // improvement, Blomquist, 26.02.07;
    stagprec = stagsave;  // restore the old stagprec value
    y = adjust(y);       // y gets the previous precision
    y = einfach & y; // This intersection delivers for intervals x with large
                     // diameters the optimal result inclusion
    return y;
} // sqrt1px2(...)

l_interval sqrtx2y2(const l_interval& x, const l_interval& y) throw()
// Inclusion of sqrt(x^2+y^2); Blomquist, 14.12.02;
{
    interval ia = abs(interval(x)), ib = abs(interval(y)), einfach;
    einfach = sqrtx2y2(ia,ib); // Inclusion only with type interval.
    if (!einfach) return l_interval(0.0);

    int stagsave=stagprec, stagmax=19;
    if (stagprec>stagmax) stagprec = stagmax;

    l_interval a=abs(x), b=abs(y), r;
    int exa=expo(Sup(ia)), exb=expo(Sup(ib)), ex;
    if (exb > exa)
    {  // Permutation of a,b:
	r = a;  a = b;  b = r;
	ex = exa;  exa = exb;  exb = ex;
    }
    ex = 511 - exa;   exa = 0;
    if (ex>1022)
    { exa = ex - 1022;  ex = 1022; }  // ex > 1022 --> scaling in two steps
    times2pown(a,ex);                 // is necessary!
    if (exa) times2pown(a,exa);
    times2pown(b,ex);            // First step:  scaling b with 2^ex;
    if (exa) times2pown(b,exa);  // Second step: scaling b with 2^exa;
    r = sqrt(a*a + b*b);
    times2pown(r,-ex);           // Backscaling, first step
    if (exa) times2pown(r,-exa); // Backscaling, second step
    stagprec = stagsave;
    r = adjust(r);
    r = einfach & r;
    return r;
} // sqrtx2y2

l_interval sqrtp1m1(const l_interval& x) throw(STD_FKT_OUT_OF_DEF)
// sqrtp1m1(x) calculates an inclusion of sqrt(x+1)-1;
// Blomquist, 05.08.03;
{
    int stagsave=stagprec, stagmax=19;
    stagprec++;
    if (stagprec>stagmax) stagprec = stagmax;
    l_interval y,tmp;
    interval z = interval(x); // z is an inclusion of x
    real r = Inf(z);
    if (r < -1)
    cxscthrow(STD_FKT_OUT_OF_DEF("l_interval sqrtp1m1(const l_interval&)"));
    const real c = 1e-10;
    tmp = x+1;
    y = x<=interval(-c,c) ? x / (sqrt(tmp)+1) : sqrt(tmp)-1;
    stagprec = stagsave;
    y = adjust(y);
    return y; 
} // sqrtp1m1


l_interval sin(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW)    // Sin(x)
{
   int         stagsave = stagprec,
               stagmax = 19;
   l_interval  pihalbe,
               y;
   interval    dx = interval(x),
               einfachgenau;

   einfachgenau = sin(dx);

   if (stagprec == 1) 
      y = sin(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = adjust(l_interval(0.0));
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;
      
//      pihalbe = 2.0*atan(l_interval(1.0));
      pihalbe = li_pi4();
      times2pown(pihalbe,1); // Blomquist, 05.12.03;
      y = x-pihalbe;
      
      try {
         y = cos(y);
      }
      catch(const ERROR_LINTERVAL_FAK_OVERFLOW &) // Damit Funktionsname stimmt
      {
         cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval sin(const l_interval & x)")); 
      }
      
      
      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;
   }

   return y;
}

l_interval cos(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW)   // Cos(x)
{
   long int    mm = 6;
   int         stagsave = stagprec,
               stagmax = 19,
               n = 0,
               digits = 53,
               sign = 0,
               degree, k,
               m2;
               
   bool        xinf = false,
               xsup = false,
               fertig = false;

   real        abst, m, eps,
//               bas, tn, t4,
               lneps, lnt,
               zk,
               zhn = 1.0,
               lnb = 0.69314718,
               fak = 720.0;    // 6!
   interval    dx = interval(x),
               einfachgenau,
               dt, extr, error;
   l_interval  zwopi, t, t2, p, ph, test,
               y;

   einfachgenau = cos(dx);
   if (stagprec == 1) 
      y = cos(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = adjust(l_interval(1.0));
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;

      // stagprec++;
      // zwopi = 8.0*atan(l_interval(1.0));
      zwopi = li_pi4();
      times2pown(zwopi,3); // Blomquist, 05.12.03;
      // stagprec--;

      // erste Argumentreduktion ( cos(x) = cos(x+2k*pi) )
      if (Sup(zwopi) < Sup(abs(x))) 
      {
         m = floor(_double(real(Sup((x/zwopi)))));  // floor rundet zur naechsten
         t = x-m*zwopi;                              // ganzen Zahl kleiner x
         if (Sup(zwopi) < Sup(abs(t))) 
         { // das ganze nochmal, falls m uebergelaufen ist!
                                      // ohne Sup wird Inclusion geprueft!
            m = floor(_double(real(Sup((t/zwopi)))));// rundet zur naechsten Zahl < x
            t = t-m*zwopi;
         }
      } else 
         t = x;

      // ueberpruefen, ob Maximum oder Minimum im Inneren des Intervalls
      extr = interval(2.0/zwopi*t);
      m2 = int(double(floor(_double(Sup(extr)))));
      if (interval(real(m2)) <= extr) 
      {
         if (interval(real(m2-1),real(m2)) <= extr) 
         {
            y = l_interval(-1.0,1.0);
            fertig = true;
         } else 
            if (m2%2)  
               xinf = TRUE;
            else            
               xsup = TRUE;
      }

      if (!(fertig)) 
      {
         // zweite Argumentreduktion
         dt = interval(t);
         // eps = 0.01/(stagprec*stagprec*stagprec);
         eps = 0.01;
         while (Sup(abs(dt))/zhn >= eps) 
         {
            n++;
            zhn += zhn;
         }
         t /= zhn;

         // Abschaetzung der Genauigkeit

         t2 = t*t;
         abst = real(Sup(abs(t)));
         if (abst < MinReal) 
            abst = MinReal; // nur zur Sicherheit
         lnt = ln(abst);
         lneps = (1.0-digits*stagprec)*lnb;
         while (lneps-(real(mm)*lnt+ln(2.0/fak)) <= 0.0) 
         {
            mm  += 4;
            if (mm > 170) 
            {   
               // 170! = 7.26*10^306
               cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval cos(const l_interval & x)"));
               mm = 170;
               break;
            }
            fak *= mm*(mm-1.0)*(mm-2.0)*(mm-3.0);
         }
         /*
            bas = real(Inf(power(l_interval(2.0),(1-digits*stagprec))));
            tn = abs(real(Sup(power(t,6))));
            t4 = real(mid(power(t,4)));
            while (bas-2.0*tn/fak <= 0.0) 
            {
               mm += 4;
               if (mm > 170) 
               {    // 170! = 7.26*10^306
                  errmon (ERROR_LINTERVAL(FAKOVERFLOW));
                  mm = 170;
                  break;
            }
            tn *= t4;
            fak *= mm*(mm-1.0)*(mm-2.0)*(mm-3.0);
         }
         */

         degree = mm-2;     // Achtung mm := 2n+2 !

         // Polynomauswertung

         sign = (degree/2)%2;
         zk  = real(degree)*real(degree-1);
         if (sign) 
            p = -t2/zk;
         else      
            p = t2/zk;
         for (k = degree-2; k >= 2; k -= 2) 
         {
            sign = 1-sign;
            if (sign)   
               p -= 1.0;
            else       
               p += 1.0;
            zk = real(k)*real(k-1);
            p *= t2/zk;
         }

         error = pow(interval(2.0), interval(1.0-digits*stagprec))
                 * interval(-0.5,0.5);
         p = l_interval(1.0)+error+p;

         // Rueckgaengigmachung der zweiten Argumentreduktion
         for (int i = 0; i < n; i++) 
            p = 2.0*p*p-1.0;

         stagprec = stagsave;
         y = adjust(p);
         if (Inf(y) < -1.0) 
            SetInf(y,-1.0);
         if (Sup(y) > 1.0) 
            SetSup(y,1.0);
         if (xinf) 
            SetInf(y,-1.0);
         if (xsup) 
            SetSup(y,1.0);
      }
      y = y & einfachgenau;
   }
   return y;
}

l_interval tan(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW,ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF)   // Tan(x)
{
   interval    dx = interval(x),
               einfachgenau;
   l_interval  s, c, y;

   einfachgenau = tan(dx);

   if (stagprec == 1) 
      y = tan(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
         y = adjust(l_interval(0.0));
   else 
   {
      try 
      {
         c = cos(x);
      }
      catch(const ERROR_LINTERVAL_FAK_OVERFLOW &) // Damit Funktionsname stimmt
      {
         cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval tan(const l_interval &x)"));
      }

      if (interval(0.0) <= c) 
      {
         cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval tan(const l_interval &x)"));
      }
      try
      {
         s = sin(x);
      }
      catch(const ERROR_LINTERVAL_FAK_OVERFLOW &)
      {
         cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval tan(const l_interval &x)"));
      }
      stagprec++;
      y = s/c;
      stagprec--;
      y = adjust(y);
      y = y & einfachgenau;
   }
   return y;
}

l_interval cot(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW,ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF)   // Cot(x)
{
   interval    dx = interval(x),
               einfachgenau;
   l_interval  s, c, y;

   einfachgenau = cot(dx);

   if (stagprec == 1) 
      y = tan(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = adjust(l_interval(0.0));
   else 
   {
      try
      {
         s = sin(x);
      }
      catch(const ERROR_LINTERVAL_FAK_OVERFLOW &) // Damit Funktionsname stimmt
      {
         cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval cot(const l_interval &x)"));
      }

      if (interval(0.0) <= s) 
      {
         cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval cot(const l_interval &x)"));
      }

      try 
      {
         c = cos(x);
      }
      catch(const ERROR_LINTERVAL_FAK_OVERFLOW &) // Damit Funktionsname stimmt
      {
         cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval cot(const l_interval &x)"));
      }
   
      c = cos(x);
      stagprec++;
      y = c/s;
      stagprec--;
      y = adjust(y);
      y = y & einfachgenau;
   }

   return y;
}

l_interval asin(const l_interval & x) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF)  // ASin(x)
{
   l_interval  t, ta, u, pihalbe,
               y;
   interval    dx = interval(x),
               einfachgenau;
   real        supabsdx = Sup(abs(dx)),
               infdx = Inf(dx),
               supdx = Sup(dx);

   einfachgenau = asin(dx);

   stagprec++;
//   pihalbe = 2.0*atan(l_interval(1.0));
   pihalbe = li_pi4();
   times2pown(pihalbe,1); // Blomquist, 05.12.03;
   stagprec--;
   
   if (Inf(x) < CXSC_MinusOne || Sup(x) > CXSC_One) 
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval asin(const l_interval & x)"));
   else if (stagprec == 1) 
      y = asin(dx);
   else if (infdx == supdx && supdx == CXSC_Zero)
      y = adjust(l_interval(0.0));
   else if (infdx == supdx && supabsdx == CXSC_One) 
   {
      if (supdx == 1.0) 
         y = pihalbe;
      else              
         y = -pihalbe;
   } else 
   {
      stagprec++;
      try 
      {
         if (supabsdx <= 0.75) 
            u = x;
         else                  
            u = 2.0*x*sqrt((1.0-x)*(1.0+x));
         t = u/sqrt((1.0-u)*(1.0+u));
         stagprec--;
         ta = atan(t);
      }
      catch(const ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF &)
      {
         cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval asin(const l_interval & x)"));
      }
      stagprec++;
      if (supabsdx <= 0.75) 
         y = ta;
      else if (Sup(t) >= 0.0)  
         y = pihalbe-0.5*ta;
      else 
         y = -pihalbe-0.5*ta;
      stagprec--;
      y = adjust(y);
      y = y & einfachgenau;
   }

   return y;
}

l_interval acos(const l_interval & x) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF)   // ACos(x)
{
   bool        neg=false;
   l_interval  pi, s, y;
   interval    dx = interval(x),
               einfachgenau;
   real        supabsdx = Sup(abs(dx)),
               infdx = Inf(dx),
               supdx = Sup(dx);

   try {

      einfachgenau = acos(dx);

      stagprec++;
//      pi = 4.0*atan(l_interval(1.0));
      pi = li_pi4();
      times2pown(pi,2); // Blomquist, 05.12.03;
      stagprec--;
      if (Inf(x) < CXSC_MinusOne || Sup(x) > CXSC_One) 
         cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF());
      else if (stagprec == 1) 
         y = acos(dx);
      else if (infdx == supdx && supabsdx == CXSC_One) 
      {
         if (supdx == 1.0) 
            y = adjust(l_interval(0.0));
         else              
            y = adjust(pi);
      } else 
      {
         if (supdx < 0.0) 
         {
            y = -x;
            neg = true;
         } else 
         {
            y = x;
            neg = false;
         }
         stagprec++;
         if (supabsdx > 0.75) 
         {
            y = (1.0-(y))*(1.0+(y));
            //      stagprec--;
            y = asin(sqrt(y));
            //      stagprec++;
            if (neg)  
               y = pi-(y);
         } else 
         {
            //      stagprec--;
            s = asin(y);
            //      stagprec++;
            if (neg)  
               y = 0.5*pi+s;
            else      
               y = 0.5*pi-s;
         }

         // Fehler in der Systemumgebung, deshalb die Aufblaehung des Intervalls
         real      err1, err2;
         interval  error;
         err1 = 5.0*real(diam(y));
         err2 = 5.0*power(10.0,-16*(stagprec-1)-1);
         if (err1 > err2)  
            error = interval(-err1, err1);
         else              
            error = interval(-err2, err2);
         y += error;

         stagprec--;
         y = adjust(y);
         y = y & einfachgenau;
      }
   }
   catch(const ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF &)
   {
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval acos(const l_interval & x)"));
   }
   return y;
}

static real CXSC_ln2[21]; // CXSC_ln2[0], ... CXSC_ln2[20] 
static bool CXSC_ln2_initialized = false;

// l_interval li_ln2() throw()
l_interval Ln2_l_interval() throw()
// Inclusion of ln(2), Blomquist, 04.12.03;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_ln2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+162E42FEFA39EFe3FE";
       str >> CXSC_ln2[0];
       str = "+1ABC9E3B39803Fe3C7";
       str >> CXSC_ln2[1];
       str = "+17B57A079A1934e390";
       str >> CXSC_ln2[2];
       str = "-1ACE93A4EBE5D1e35A";
       str >> CXSC_ln2[3];
       str = "-123A2A82EA0C24e324";
       str >> CXSC_ln2[4];
       str = "+1D881B7AEB2615e2ED";
       str >> CXSC_ln2[5];
       str = "+19552FB4AFA1B1e2B7";
       str >> CXSC_ln2[6];
       str = "+1DA5D5C6B82704e27E";
       str >> CXSC_ln2[7];
       str = "+14427573B29117e247";
       str >> CXSC_ln2[8];
       str = "-191F6B05A4D7A7e211";
       str >> CXSC_ln2[9];
       str = "-1DB5173AE53426e1DB";
       str >> CXSC_ln2[10];
       str = "+11317C387EB9EBe1A3";
       str >> CXSC_ln2[11];
       str = "-190F13B267F137e16D";
       str >> CXSC_ln2[12];
       str = "+16FA0EC7657F75e137";
       str >> CXSC_ln2[13];
       str = "-1234C5E1398A6Be101";
       str >> CXSC_ln2[14];
       str = "+1195EBBF4D7A70e0CA";
       str >> CXSC_ln2[15];
       str = "+18192432AFD0C4e094";
       str >> CXSC_ln2[16];
       str = "-1A1BE38BA4BA4De05E";
       str >> CXSC_ln2[17];
       str = "-1D7860151CFC06e024";
       str >> CXSC_ln2[18];
       str = "+1000032847ED6Fe000";
       str >> CXSC_ln2[19];
       str = "+1000032847ED70e000";
       str >> CXSC_ln2[20];

       CXSC_ln2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_ln2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

static real CXSC_ln10[21]; // CXSC_ln10[0], ... CXSC_ln10[20] 
static bool CXSC_ln10_initialized = false;

// l_interval li_ln10() throw()
l_interval Ln10_l_interval() throw()
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_ln10_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+126BB1BBB55516e400";
       str >> CXSC_ln10[0];
       str = "-1F48AD494EA3E9e3CA";
       str >> CXSC_ln10[1];
       str = "-19EBAE3AE0260Ce394";
       str >> CXSC_ln10[2];
       str = "-12D10378BE1CF1e35E";
       str >> CXSC_ln10[3];
       str = "+10403E05AE52C6e328";
       str >> CXSC_ln10[4];
       str = "-1FA509CAFDF466e2F0";
       str >> CXSC_ln10[5];
       str = "-1C79A1FE9D0795e2BA";
       str >> CXSC_ln10[6];
       str = "+1058C448308218e284";
       str >> CXSC_ln10[7];
       str = "-1D250470877BFDe24D";
       str >> CXSC_ln10[8];
       str = "-1AE92987D3075De215";
       str >> CXSC_ln10[9];
       str = "-1D5CDBB8626956e1DF";
       str >> CXSC_ln10[10];
       str = "-13C4F27CE0410Ae1A9";
       str >> CXSC_ln10[11];
       str = "+1B3AC12ACF1BE9e173";
       str >> CXSC_ln10[12];
       str = "+1161BB49D219C8e13D";
       str >> CXSC_ln10[13];
       str = "-110D6613293728e107";
       str >> CXSC_ln10[14];
       str = "+142163A4CDA351e0CF";
       str >> CXSC_ln10[15];
       str = "+1E2713D6C22C16e097";
       str >> CXSC_ln10[16];
       str = "-15090EF85CB0ADe05E";
       str >> CXSC_ln10[17];
       str = "-1C5B3E859F876Ee027";
       str >> CXSC_ln10[18];
       str = "-10000703552C52e000";
       str >> CXSC_ln10[19];
       str = "-10000703552C51e000";
       str >> CXSC_ln10[20];

       CXSC_ln10_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_ln10[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

static real CXSC_Rln10[21]; // CXSC_Rln10[0], ... CXSC_Rln10[20] 
static bool CXSC_Rln10_initialized = false;

// l_interval li_Rln10() throw()
l_interval Ln10r_l_interval() throw()
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Rln10_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1BCB7B1526E50Ee3FD";
       str >> CXSC_Rln10[0];
       str = "+195355BAAAFAD3e3C6";
       str >> CXSC_Rln10[1];
       str = "+1EE191F71A3012e38F";
       str >> CXSC_Rln10[2];
       str = "+17268808E8FCB5e358";
       str >> CXSC_Rln10[3];
       str = "+13DE3A94F1D509e320";
       str >> CXSC_Rln10[4];
       str = "+1DF42805E7E524e2E9";
       str >> CXSC_Rln10[5];
       str = "+11AAC96323250Be2B3";
       str >> CXSC_Rln10[6];
       str = "-1CE63884C058E4e27D";
       str >> CXSC_Rln10[7];
       str = "-1A1C82EA3969BAe247";
       str >> CXSC_Rln10[8];
       str = "+1B4F6686AD7A33e211";
       str >> CXSC_Rln10[9];
       str = "-1B97C8035FFC70e1DB";
       str >> CXSC_Rln10[10];
       str = "+1630771369962Ee1A0";
       str >> CXSC_Rln10[11];
       str = "-1E15BD37B295AFe16A";
       str >> CXSC_Rln10[12];
       str = "-132484B432318Be134";
       str >> CXSC_Rln10[13];
       str = "+15430212AE68C0e0FE";
       str >> CXSC_Rln10[14];
       str = "+1351923B322731e0C8";
       str >> CXSC_Rln10[15];
       str = "+11F934D794D64Fe092";
       str >> CXSC_Rln10[16];
       str = "+13E4B475D9FF20e05B";
       str >> CXSC_Rln10[17];
       str = "+185D9B63ED9A24e025";
       str >> CXSC_Rln10[18];
       str = "+1000035B8CA18Ce000";
       str >> CXSC_Rln10[19];
       str = "+1000035B8CA18De000";
       str >> CXSC_Rln10[20];

       CXSC_Rln10_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Rln10[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

static real    CXSC_pi4[21];
static bool    CXSC_pi4_initialized = false;

// l_interval li_pi4() throw()
l_interval Pid4_l_interval() throw()
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_pi4_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1921FB54442D18e3FE";
       str >> CXSC_pi4[0];
       str = "+11A62633145C06e3C8";
       str >> CXSC_pi4[1];
       str = "+1C1CD129024E08e393";
       str >> CXSC_pi4[2];
       str = "+114CF98E804178e35E";
       str >> CXSC_pi4[3];
       str = "-159C4EC64DDAECe327";
       str >> CXSC_pi4[4];
       str = "+1410F31C6809BCe2F2";
       str >> CXSC_pi4[5];
       str = "-106AE64C32C5BCe2BB";
       str >> CXSC_pi4[6];
       str = "-1C99FA9EB241B4e286";
       str >> CXSC_pi4[7];
       str = "-1D791603D95252e24E";
       str >> CXSC_pi4[8];
       str = "-1571EDD0DBD254e218";
       str >> CXSC_pi4[9];
       str = "-133B4302721768e1E2";
       str >> CXSC_pi4[10];
       str = "+10BA698DFB5AC2e1AD";
       str >> CXSC_pi4[11];
       str = "+1FFAE5B7A035C0e178";
       str >> CXSC_pi4[12];
       str = "-1211C79404A576e143";
       str >> CXSC_pi4[13];
       str = "-1816945836FBA0e10D";
       str >> CXSC_pi4[14];
       str = "-1DA700CDB6BCCEe0D8";
       str >> CXSC_pi4[15];
       str = "+11ECE45B3DC200e0A3";
       str >> CXSC_pi4[16];
       str = "+1F2E2858EFC166e06D";
       str >> CXSC_pi4[17];
       str = "+1B4906C38ABA73e036";
       str >> CXSC_pi4[18];
       str = "+19A458FEA3F493e000";
       str >> CXSC_pi4[19];
       str = "+19A458FEA3F494e000";
       str >> CXSC_pi4[20];

       CXSC_pi4_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_pi4[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

static real    CXSC_sqrt2[21];
static bool    CXSC_sqrt2_initialized = false;

// l_interval li_sqrt2() throw()
l_interval Sqrt2_l_interval() throw()
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_sqrt2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+16A09E667F3BCDe3FF";
       str >> CXSC_sqrt2[0];
       str = "-1BDD3413B26456e3C9";
       str >> CXSC_sqrt2[1];
       str = "+157D3E3ADEC175e393";
       str >> CXSC_sqrt2[2];
       str = "+12775099DA2F59e35B";
       str >> CXSC_sqrt2[3];
       str = "+160CCE64552BF2e322";
       str >> CXSC_sqrt2[4];
       str = "+1821D5C5161D46e2E9";
       str >> CXSC_sqrt2[5];
       str = "-1C032046F8498Ee2B3";
       str >> CXSC_sqrt2[6];
       str = "+1EE950BC8738F7e27B";
       str >> CXSC_sqrt2[7];
       str = "-1AC3FDBC64E103e245";
       str >> CXSC_sqrt2[8];
       str = "+13B469101743A1e20D";
       str >> CXSC_sqrt2[9];
       str = "+15E3E9CA60B38Ce1D7";
       str >> CXSC_sqrt2[10];
       str = "+11BC337BCAB1BDe19C";
       str >> CXSC_sqrt2[11];
       str = "-1BBA5DEE9D6E7De166";
       str >> CXSC_sqrt2[12];
       str = "-1438DD083B1CC4e130";
       str >> CXSC_sqrt2[13];
       str = "+1B56A28E2EDFA7e0FA";
       str >> CXSC_sqrt2[14];
       str = "+1CCB2A634331F4e0C4";
       str >> CXSC_sqrt2[15];
       str = "-1BD9056876F83Ee08D";
       str >> CXSC_sqrt2[16];
       str = "-1234FA22AB6BEFe057";
       str >> CXSC_sqrt2[17];
       str = "+19040CA4A81395e020";
       str >> CXSC_sqrt2[18];
       str = "-1000002A493818e000";
       str >> CXSC_sqrt2[19];
       str = "-1000002A493817e000";
       str >> CXSC_sqrt2[20];

       CXSC_sqrt2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_sqrt2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

static real CXSC_sqrt5[21]; // CXSC_sqrt5[0], ... CXSC_sqrt5[20] 
static bool CXSC_sqrt5_initialized = false;

l_interval Sqrt5_l_interval() throw()
// Inclusion of sqrt(5), Blomquist, 30.11.2008;
{
	l_interval y;
	int stagsave = stagprec,
 stagmax = 20;
   
 if (!CXSC_sqrt5_initialized) 
 {
	 std::string str;
	 std::cout << SaveOpt;
	 std::cout << Hex;
	 str = "+11E3779B97F4A8e400";
	 str >> CXSC_sqrt5[0];
	 str = "-1F506319FCFD19e3C9";
	 str >> CXSC_sqrt5[1];
	 str = "+1B906821044ED8e393";
	 str >> CXSC_sqrt5[2];
	 str = "-18BB1B5C0F272Ce35B";
	 str >> CXSC_sqrt5[3];
	 str = "+11D0C18E952768e324";
	 str >> CXSC_sqrt5[4];
	 str = "-1E9D585B0901F9e2EB";
	 str >> CXSC_sqrt5[5];
	 str = "-1C7DD252073EC0e2B5";
	 str >> CXSC_sqrt5[6];
	 str = "-1FCEF21EDAF7FAe27F";
	 str >> CXSC_sqrt5[7];
	 str = "+160EB25D20799Be241";
	 str >> CXSC_sqrt5[8];
	 str = "-1C90F95285168Fe208";
	 str >> CXSC_sqrt5[9];
	 str = "+1E1DFA160E75BCe1D2";
	 str >> CXSC_sqrt5[10];
	 str = "-10A08E66CB368Ce196";
	 str >> CXSC_sqrt5[11];
	 str = "+1C5371682CADD1e160";
	 str >> CXSC_sqrt5[12];
	 str = "-1998100220F4EDe129";
	 str >> CXSC_sqrt5[13];
	 str = "+1C6771A0968663e0F3";
	 str >> CXSC_sqrt5[14];
	 str = "+1DFB9E3C86CA7Ce0BD";
	 str >> CXSC_sqrt5[15];
	 str = "-18AE38ED5304B1e086";
	 str >> CXSC_sqrt5[16];
	 str = "+182A5FEC507706e050";
	 str >> CXSC_sqrt5[17];
	 str = "-1B5191A18C5647e018";
	 str >> CXSC_sqrt5[18];
	 str = "+100000000F9D52e000";
	 str >> CXSC_sqrt5[19];
	 str = "+100000000F9D53e000";
	 str >> CXSC_sqrt5[20];

	 CXSC_sqrt5_initialized = true;
	 std::cout << RestoreOpt;
 }
 stagprec = stagmax;
 y = adjust(l_interval(0.0)); 
 for (int i=0; i <= stagmax; i++)  
    y.data[i] = CXSC_sqrt5[i];
 stagprec = stagsave;
 y = adjust(y);
 return y;             
}

// ************************************************************************

static real CXSC_sqrt7[21]; // CXSC_sqrt7[0], ... CXSC_sqrt7[20] 
static bool CXSC_sqrt7_initialized = false;

l_interval Sqrt7_l_interval() throw()
// Inclusion of sqrt(7), Blomquist, 30.11.2008;
{
	l_interval y;
	int stagsave = stagprec,
 stagmax = 20;
   
 if (!CXSC_sqrt7_initialized) 
 {
	 std::string str;
	 std::cout << SaveOpt;
	 std::cout << Hex;
	 str = "+152A7FA9D2F8EAe400";
	 str >> CXSC_sqrt7[0];
	 str = "-121C62B033C079e3CA";
	 str >> CXSC_sqrt7[1];
	 str = "-177CAAD6200612e391";
	 str >> CXSC_sqrt7[2];
	 str = "-1EFA880DC72D64e359";
	 str >> CXSC_sqrt7[3];
	 str = "-171D206D5B1A4Ce31F";
	 str >> CXSC_sqrt7[4];
	 str = "+119392FA9B0494e2E6";
	 str >> CXSC_sqrt7[5];
	 str = "+17BB8A64890057e2AD";
	 str >> CXSC_sqrt7[6];
	 str = "-17E89300383DDEe277";
	 str >> CXSC_sqrt7[7];
	 str = "+130FB7AF68A6FBe241";
	 str >> CXSC_sqrt7[8];
	 str = "+1322281D303D36e209";
	 str >> CXSC_sqrt7[9];
	 str = "+1996109A16D3B1e1D3";
	 str >> CXSC_sqrt7[10];
	 str = "+1F239C301DFBB4e19C";
	 str >> CXSC_sqrt7[11];
	 str = "-1B5CA40AB771A2e163";
	 str >> CXSC_sqrt7[12];
	 str = "-1675711487FEAAe12A";
	 str >> CXSC_sqrt7[13];
	 str = "+122CB7FA26ABA5e0F4";
	 str >> CXSC_sqrt7[14];
	 str = "+1059211B7D5398e0BD";
	 str >> CXSC_sqrt7[15];
	 str = "-10F15BFA46EB7Fe087";
	 str >> CXSC_sqrt7[16];
	 str = "+15AB71566CE72Be051";
	 str >> CXSC_sqrt7[17];
	 str = "-1386BDCA3845C7e01A";
	 str >> CXSC_sqrt7[18];
	 str = "+10000000AC4BC7e000";
	 str >> CXSC_sqrt7[19];
	 str = "+10000000AC4BC8e000";
	 str >> CXSC_sqrt7[20];

	 CXSC_sqrt7_initialized = true;
	 std::cout << RestoreOpt;
 }
 stagprec = stagmax;
 y = adjust(l_interval(0.0)); 
 for (int i=0; i <= stagmax; i++)  
	 y.data[i] = CXSC_sqrt7[i];
 stagprec = stagsave;
 y = adjust(y);
 return y;
}


// ************************************************************************

static real CXSC_ln2r[21]; // CXSC_ln2r[0], ... CXSC_ln2r[20] 
static bool CXSC_ln2r_initialized = false;

l_interval Ln2r_l_interval() throw()
// Inclusion of 1/ln(2), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_ln2r_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+171547652B82FEe3FF";
       str >> CXSC_ln2r[0];
       str = "+1777D0FFDA0D24e3C7";
       str >> CXSC_ln2r[1];
       str = "-160BB8A5442AB9e391";
       str >> CXSC_ln2r[2];
       str = "-14B52D3BA6D74De359";
       str >> CXSC_ln2r[3];
       str = "+19A342648FBC39e323";
       str >> CXSC_ln2r[4];
       str = "-1E0455744994EEe2ED";
       str >> CXSC_ln2r[5];
       str = "+1B25EEB82D7C16e2B7";
       str >> CXSC_ln2r[6];
       str = "+1F5485CF306255e281";
       str >> CXSC_ln2r[7];
       str = "-1EC07680A1F958e24B";
       str >> CXSC_ln2r[8];
       str = "-106326680EB5B6e215";
       str >> CXSC_ln2r[9];
       str = "-1B3D04C549BC98e1DF";
       str >> CXSC_ln2r[10];
       str = "+1EABCEAD10305Be1A9";
       str >> CXSC_ln2r[11];
       str = "-14440C57D7AB97e170";
       str >> CXSC_ln2r[12];
       str = "-17185D42A4E6D6e139";
       str >> CXSC_ln2r[13];
       str = "-1F332B5BE48526e101";
       str >> CXSC_ln2r[14];
       str = "+12CE4F199E108De0CB";
       str >> CXSC_ln2r[15];
       str = "-18DAFCC6077F2Ae092";
       str >> CXSC_ln2r[16];
       str = "+19ABB71EC25E12e05B";
       str >> CXSC_ln2r[17];
       str = "-11473D7A3366BDe022";
       str >> CXSC_ln2r[18];
       str = "-1000004977D38Be000";
       str >> CXSC_ln2r[19];
       str = "-1000004977D38Ae000";
       str >> CXSC_ln2r[20];

       CXSC_ln2r_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_ln2r[i];
//       y[i+1] = CXSC_ln2r[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

       
// **********************************************************************

static real CXSC_Pi[21]; // CXSC_Pi[0], ... CXSC_Pi[20] 
static bool CXSC_Pi_initialized = false;

l_interval Pi_l_interval() throw()
// Inclusion of Pi, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1921FB54442D18e400";
       str >> CXSC_Pi[0];
       str = "+11A62633145C07e3CA";
       str >> CXSC_Pi[1];
       str = "-1F1976B7ED8FBCe392";
       str >> CXSC_Pi[2];
       str = "+14CF98E804177De35C";
       str >> CXSC_Pi[3];
       str = "+131D89CD9128A5e326";
       str >> CXSC_Pi[4];
       str = "+10F31C6809BBDFe2EC";
       str >> CXSC_Pi[5];
       str = "+1519B3CD3A431Be2B5";
       str >> CXSC_Pi[6];
       str = "+18158536F92F8Ae27E";
       str >> CXSC_Pi[7];
       str = "+1BA7F09AB6B6A9e246";
       str >> CXSC_Pi[8];
       str = "-1EDD0DBD2544CFe20E";
       str >> CXSC_Pi[9];
       str = "+179FB1BD1310BAe1D7";
       str >> CXSC_Pi[10];
       str = "+1A637ED6B0BFF6e1A1";
       str >> CXSC_Pi[11];
       str = "-1A485FCA40908Ee16A";
       str >> CXSC_Pi[12];
       str = "-1E501295D98169e133";
       str >> CXSC_Pi[13];
       str = "-1160DBEE83B4E0e0FD";
       str >> CXSC_Pi[14];
       str = "-19B6D799AE131Ce0C5";
       str >> CXSC_Pi[15];
       str = "+16CF70801F2E28e08F";
       str >> CXSC_Pi[16];
       str = "+163BF0598DA483e059";
       str >> CXSC_Pi[17];
       str = "+1871574E69A459e023";
       str >> CXSC_Pi[18];
       str = "-10000005702DB4e000";
       str >> CXSC_Pi[19];
       str = "-10000005702DB3e000";
       str >> CXSC_Pi[20];

       CXSC_Pi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pi[i]; 
//       y[i+1] = CXSC_Pi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

// **********************************************************************

static real CXSC_Pid2[21]; // CXSC_Pid2[0], ... CXSC_Pid2[20] 
static bool CXSC_Pid2_initialized = false;

l_interval Pid2_l_interval() throw()
// Inclusion of Pi/2, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pid2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1921FB54442D18e3FF";
       str >> CXSC_Pid2[0];
       str = "+11A62633145C07e3C9";
       str >> CXSC_Pid2[1];
       str = "-1F1976B7ED8FBCe391";
       str >> CXSC_Pid2[2];
       str = "+14CF98E804177De35B";
       str >> CXSC_Pid2[3];
       str = "+131D89CD9128A5e325";
       str >> CXSC_Pid2[4];
       str = "+10F31C6809BBDFe2EB";
       str >> CXSC_Pid2[5];
       str = "+1519B3CD3A431Be2B4";
       str >> CXSC_Pid2[6];
       str = "+18158536F92F8Ae27D";
       str >> CXSC_Pid2[7];
       str = "+1BA7F09AB6B6A9e245";
       str >> CXSC_Pid2[8];
       str = "-1EDD0DBD2544CFe20D";
       str >> CXSC_Pid2[9];
       str = "+179FB1BD1310BAe1D6";
       str >> CXSC_Pid2[10];
       str = "+1A637ED6B0BFF6e1A0";
       str >> CXSC_Pid2[11];
       str = "-1A485FCA40908Ee169";
       str >> CXSC_Pid2[12];
       str = "-1E501295D98169e132";
       str >> CXSC_Pid2[13];
       str = "-1160DBEE83B4E0e0FC";
       str >> CXSC_Pid2[14];
       str = "-19B6D799AE131Ce0C4";
       str >> CXSC_Pid2[15];
       str = "+16CF70801F2E28e08E";
       str >> CXSC_Pid2[16];
       str = "+163BF0598DA483e058";
       str >> CXSC_Pid2[17];
       str = "+1871574E69A459e022";
       str >> CXSC_Pid2[18];
       str = "-10000002B816DAe000";
       str >> CXSC_Pid2[19];
       str = "-10000002B816D0e000";
       str >> CXSC_Pid2[20];

       CXSC_Pid2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pid2[i]; 
//       y[i+1] = CXSC_Pid2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

// **********************************************************************

static real CXSC_Pi2[21];  // CXSC_Pi2[0], ... CXSC_Pi2[20] 
static bool CXSC_Pi2_initialized = false;

l_interval Pi2_l_interval() throw()
// Inclusion of 2*Pi, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pi2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1921FB54442D18e401";
       str >> CXSC_Pi2[0];
       str = "+11A62633145C07e3CB";
       str >> CXSC_Pi2[1];
       str = "-1F1976B7ED8FBCe393";
       str >> CXSC_Pi2[2];
       str = "+14CF98E804177De35D";
       str >> CXSC_Pi2[3];
       str = "+131D89CD9128A5e327";
       str >> CXSC_Pi2[4];
       str = "+10F31C6809BBDFe2ED";
       str >> CXSC_Pi2[5];
       str = "+1519B3CD3A431Be2B6";
       str >> CXSC_Pi2[6];
       str = "+18158536F92F8Ae27F";
       str >> CXSC_Pi2[7];
       str = "+1BA7F09AB6B6A9e247";
       str >> CXSC_Pi2[8];
       str = "-1EDD0DBD2544CFe20F";
       str >> CXSC_Pi2[9];
       str = "+179FB1BD1310BAe1D8";
       str >> CXSC_Pi2[10];
       str = "+1A637ED6B0BFF6e1A2";
       str >> CXSC_Pi2[11];
       str = "-1A485FCA40908Ee16B";
       str >> CXSC_Pi2[12];
       str = "-1E501295D98169e134";
       str >> CXSC_Pi2[13];
       str = "-1160DBEE83B4E0e0FE";
       str >> CXSC_Pi2[14];
       str = "-19B6D799AE131Ce0C6";
       str >> CXSC_Pi2[15];
       str = "+16CF70801F2E28e090";
       str >> CXSC_Pi2[16];
       str = "+163BF0598DA483e05A";
       str >> CXSC_Pi2[17];
       str = "+1871574E69A459e024";
       str >> CXSC_Pi2[18];
       str = "-1000000AE05B67e000";
       str >> CXSC_Pi2[19];
       str = "-1000000AE05B66e000";
       str >> CXSC_Pi2[20];

       CXSC_Pi2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pi2[i];
//       y[i+1] = CXSC_Pi2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

// **********************************************************************

static real CXSC_Pid3[21]; // CXSC_Pid3[0], ... CXSC_Pid3[20] 
static bool CXSC_Pid3_initialized = false;

l_interval Pid3_l_interval() throw()
// Inclusion of Pi/3, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pid3_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+10C152382D7366e3FF";
       str >> CXSC_Pid3[0];
       str = "-1EE6913347C2A6e3C9";
       str >> CXSC_Pid3[1];
       str = "-14BBA47A9E5FD2e391";
       str >> CXSC_Pid3[2];
       str = "-1CCAEF65529B02e35B";
       str >> CXSC_Pid3[3];
       str = "+197CB7BCC18B87e324";
       str >> CXSC_Pid3[4];
       str = "-13EBBDA1FF3058e2EE";
       str >> CXSC_Pid3[5];
       str = "-11D10CB320F4D1e2B6";
       str >> CXSC_Pid3[6];
       str = "+1958EB892987ECe27F";
       str >> CXSC_Pid3[7];
       str = "+167C54B11CF247e249";
       str >> CXSC_Pid3[8];
       str = "+12C2E985923A44e210";
       str >> CXSC_Pid3[9];
       str = "+1945484A2DD81Fe1D8";
       str >> CXSC_Pid3[10];
       str = "+1197A9E475D54Fe1A0";
       str >> CXSC_Pid3[11];
       str = "-1E181FEE158585e16A";
       str >> CXSC_Pid3[12];
       str = "+1047FCE7066A6Ee134";
       str >> CXSC_Pid3[13];
       str = "+1D1A8602EA0C85e0FE";
       str >> CXSC_Pid3[14];
       str = "+14430C5998BF34e0C8";
       str >> CXSC_Pid3[15];
       str = "+173BF40AAD43D9e091";
       str >> CXSC_Pid3[16];
       str = "-137B014DDEDCF5e05B";
       str >> CXSC_Pid3[17];
       str = "-1A5F1B210EE7C5e022";
       str >> CXSC_Pid3[18];
       str = "+100000A8DA9B6Ee000";
       str >> CXSC_Pid3[19];
       str = "+100000A8DA9B6Fe000";
       str >> CXSC_Pid3[20];

       CXSC_Pid3_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pid3[i];
//       y[i+1] = CXSC_Pid3[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Pir[21]; // CXSC_Pir[0], ... CXSC_Pir[20] 
static bool CXSC_Pir_initialized = false;

l_interval Pir_l_interval() throw()
// Inclusion of 1/Pi, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pir_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+145F306DC9C883e3FD";
       str >> CXSC_Pir[0];
       str = "-16B01EC5417056e3C7";
       str >> CXSC_Pir[1];
       str = "-16447E493AD4CEe391";
       str >> CXSC_Pir[2];
       str = "+1E21C820FF28B2e35B";
       str >> CXSC_Pir[3];
       str = "-1508510EA79237e324";
       str >> CXSC_Pir[4];
       str = "+1B8E909374B802e2EC";
       str >> CXSC_Pir[5];
       str = "-1B6D115F62E6DEe2B6";
       str >> CXSC_Pir[6];
       str = "-180F10A71A76B3e27F";
       str >> CXSC_Pir[7];
       str = "+1CFBA208D7D4BBe248";
       str >> CXSC_Pir[8];
       str = "-12EDEC598E3F65e210";
       str >> CXSC_Pir[9];
       str = "-1741037D8CDC54e1D9";
       str >> CXSC_Pir[10];
       str = "+1CC1A99CFA4E42e1A3";
       str >> CXSC_Pir[11];
       str = "+17E2EF7E4A0EC8e16C";
       str >> CXSC_Pir[12];
       str = "-1DA00087E99FC0e130";
       str >> CXSC_Pir[13];
       str = "-10D0EE74A5F593e0FA";
       str >> CXSC_Pir[14];
       str = "+1F6D367ECF27CBe0C2";
       str >> CXSC_Pir[15];
       str = "+136E9E8C7ECD3De089";
       str >> CXSC_Pir[16];
       str = "-100AE9456C229Ce053";
       str >> CXSC_Pir[17];
       str = "-141A0E84C2F8C6e01A";
       str >> CXSC_Pir[18];
       str = "-1000000010EB5Be000";
       str >> CXSC_Pir[19];
       str = "-1000000010EB5Ae000";
       str >> CXSC_Pir[20];

       CXSC_Pir_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pir[i];
//       y[i+1] = CXSC_Pir[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}

// **********************************************************************

static real CXSC_Pi2r[21]; // CXSC_Pi2r[0], ... CXSC_Pi2r[20] 
static bool CXSC_Pi2r_initialized = false;

l_interval Pi2r_l_interval() throw()
// Inclusion of 1/(2*Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pi2r_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+145F306DC9C883e3FC";
       str >> CXSC_Pi2r[0];
       str = "-16B01EC5417056e3C6";
       str >> CXSC_Pi2r[1];
       str = "-16447E493AD4CEe390";
       str >> CXSC_Pi2r[2];
       str = "+1E21C820FF28B2e35A";
       str >> CXSC_Pi2r[3];
       str = "-1508510EA79237e323";
       str >> CXSC_Pi2r[4];
       str = "+1B8E909374B802e2EB";
       str >> CXSC_Pi2r[5];
       str = "-1B6D115F62E6DEe2B5";
       str >> CXSC_Pi2r[6];
       str = "-180F10A71A76B3e27E";
       str >> CXSC_Pi2r[7];
       str = "+1CFBA208D7D4BBe247";
       str >> CXSC_Pi2r[8];
       str = "-12EDEC598E3F65e20F";
       str >> CXSC_Pi2r[9];
       str = "-1741037D8CDC54e1D8";
       str >> CXSC_Pi2r[10];
       str = "+1CC1A99CFA4E42e1A2";
       str >> CXSC_Pi2r[11];
       str = "+17E2EF7E4A0EC8e16B";
       str >> CXSC_Pi2r[12];
       str = "-1DA00087E99FC0e12F";
       str >> CXSC_Pi2r[13];
       str = "-10D0EE74A5F593e0F9";
       str >> CXSC_Pi2r[14];
       str = "+1F6D367ECF27CBe0C1";
       str >> CXSC_Pi2r[15];
       str = "+136E9E8C7ECD3De088";
       str >> CXSC_Pi2r[16];
       str = "-100AE9456C229Ce052";
       str >> CXSC_Pi2r[17];
       str = "-141A0E84C2F8C6e019";
       str >> CXSC_Pi2r[18];
       str = "-100000000875AEe000";
       str >> CXSC_Pi2r[19];
       str = "-100000000875ADe000";
       str >> CXSC_Pi2r[20];

       CXSC_Pi2r_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pi2r[i];
//       y[i+1] = CXSC_Pi2r[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_SqrtPi[21]; // CXSC_SqrtPi[0], ... CXSC_SqrtPi[20] 
static bool CXSC_SqrtPi_initialized = false;

l_interval SqrtPi_l_interval() throw()
// Inclusion of Sqrt(Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_SqrtPi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1C5BF891B4EF6Be3FF";
       str >> CXSC_SqrtPi[0];
       str = "-1618F13EB7CA89e3C9";
       str >> CXSC_SqrtPi[1];
       str = "-1B1F0071B7AAE4e391";
       str >> CXSC_SqrtPi[2];
       str = "-1389B5A46BDFE8e35A";
       str >> CXSC_SqrtPi[3];
       str = "-160AF5C5C89448e324";
       str >> CXSC_SqrtPi[4];
       str = "-14835F07122994e2E8";
       str >> CXSC_SqrtPi[5];
       str = "+1CEC283C18EE8Fe2B2";
       str >> CXSC_SqrtPi[6];
       str = "-13ADEBB9223CA8e27B";
       str >> CXSC_SqrtPi[7];
       str = "+1454912430D291e245";
       str >> CXSC_SqrtPi[8];
       str = "-1E8B2345020EF6e20F";
       str >> CXSC_SqrtPi[9];
       str = "-17262982556291e1D8";
       str >> CXSC_SqrtPi[10];
       str = "+1196FA9B140CABe1A1";
       str >> CXSC_SqrtPi[11];
       str = "-175EEE59D91D39e16B";
       str >> CXSC_SqrtPi[12];
       str = "+1789268B7D9D48e130";
       str >> CXSC_SqrtPi[13];
       str = "+17162E2F06B89Ce0FA";
       str >> CXSC_SqrtPi[14];
       str = "+1EC9C08F40A3DBe0C3";
       str >> CXSC_SqrtPi[15];
       str = "+1B6048DD0729E2e08D";
       str >> CXSC_SqrtPi[16];
       str = "+1471CF4C33FF6Be056";
       str >> CXSC_SqrtPi[17];
       str = "+1D75FBD8B36F94e020";
       str >> CXSC_SqrtPi[18];
       str = "+1000002D74B3A2e000";
       str >> CXSC_SqrtPi[19];
       str = "+1000002D74B3A3e000";
       str >> CXSC_SqrtPi[20];

       CXSC_SqrtPi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_SqrtPi[i];
//       y[i+1] = CXSC_SqrtPi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Sqrt2Pi[21]; // CXSC_Sqrt2Pi[0], ... CXSC_Sqrt2Pi[20] 
static bool CXSC_Sqrt2Pi_initialized = false;

l_interval Sqrt2Pi_l_interval() throw()
// Inclusion of sqrt(2*Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Sqrt2Pi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+140D931FF62706e400";
       str >> CXSC_Sqrt2Pi[0];
       str = "-1A6A0D6F814637e3CA";
       str >> CXSC_Sqrt2Pi[1];
       str = "-1311D073060ACEe394";
       str >> CXSC_Sqrt2Pi[2];
       str = "+16000B50DC2F41e35B";
       str >> CXSC_Sqrt2Pi[3];
       str = "+16EF75CA45A834e324";
       str >> CXSC_Sqrt2Pi[4];
       str = "+19BDB2B4C39342e2EC";
       str >> CXSC_Sqrt2Pi[5];
       str = "+1F5582E2063EE6e2B5";
       str >> CXSC_Sqrt2Pi[6];
       str = "+183F879BEA150Ce27C";
       str >> CXSC_Sqrt2Pi[7];
       str = "-1F1EA3CA289B00e244";
       str >> CXSC_Sqrt2Pi[8];
       str = "-1699CDA77736F9e20D";
       str >> CXSC_Sqrt2Pi[9];
       str = "-11A379D298B55Ee1D4";
       str >> CXSC_Sqrt2Pi[10];
       str = "-1A6DDB0152BA94e19E";
       str >> CXSC_Sqrt2Pi[11];
       str = "-1957E2E58A02FEe167";
       str >> CXSC_Sqrt2Pi[12];
       str = "-1D6160F18E604De131";
       str >> CXSC_Sqrt2Pi[13];
       str = "+1311860CDF7215e0F8";
       str >> CXSC_Sqrt2Pi[14];
       str = "+12271F44C50274e0C1";
       str >> CXSC_Sqrt2Pi[15];
       str = "-100BF5C5497A21e08A";
       str >> CXSC_Sqrt2Pi[16];
       str = "+1E94B6E6AD51E2e052";
       str >> CXSC_Sqrt2Pi[17];
       str = "-1C910B5F3D27CEe019";
       str >> CXSC_Sqrt2Pi[18];
       str = "+100000007C99B0e000";
       str >> CXSC_Sqrt2Pi[19];
       str = "+100000007C99B1e000";
       str >> CXSC_Sqrt2Pi[20];

       CXSC_Sqrt2Pi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Sqrt2Pi[i];
//       y[i+1] = CXSC_Sqrt2Pi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_SqrtPir[21]; // CXSC_SqrtPir[0], ... CXSC_SqrtPir[20] 
static bool CXSC_SqrtPir_initialized = false;

l_interval SqrtPir_l_interval() throw()
// Inclusion of 1/sqrt(Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_SqrtPir_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+120DD750429B6De3FE";
       str >> CXSC_SqrtPir[0];
       str = "+11AE3A914FED80e3C6";
       str >> CXSC_SqrtPir[1];
       str = "-13CBBEBF65F145e38F";
       str >> CXSC_SqrtPir[2];
       str = "-1E0C574632F53Ee358";
       str >> CXSC_SqrtPir[3];
       str = "-1E6633BE9E7F15e322";
       str >> CXSC_SqrtPir[4];
       str = "+1CF859270F1141e2EB";
       str >> CXSC_SqrtPir[5];
       str = "-1FE4FB499C328Ae2B4";
       str >> CXSC_SqrtPir[6];
       str = "-10B82C446DC78De27D";
       str >> CXSC_SqrtPir[7];
       str = "-1878B089078800e247";
       str >> CXSC_SqrtPir[8];
       str = "-13DAEADA9E233Ee20F";
       str >> CXSC_SqrtPir[9];
       str = "+1137197A708BD2e1D9";
       str >> CXSC_SqrtPir[10];
       str = "-109009506D5BA2e19E";
       str >> CXSC_SqrtPir[11];
       str = "+17C9F0B5951E94e168";
       str >> CXSC_SqrtPir[12];
       str = "-1735F4949633A4e131";
       str >> CXSC_SqrtPir[13];
       str = "-146014DBC90D0Ee0FB";
       str >> CXSC_SqrtPir[14];
       str = "+1CAB0B222EEEA0e0C5";
       str >> CXSC_SqrtPir[15];
       str = "+1B1C750754B40Ae08F";
       str >> CXSC_SqrtPir[16];
       str = "-16B2CD2E72C16Ee057";
       str >> CXSC_SqrtPir[17];
       str = "-148C024FF194B2e021";
       str >> CXSC_SqrtPir[18];
       str = "+10000073E19B74e000";
       str >> CXSC_SqrtPir[19];
       str = "+10000073E19B75e000";
       str >> CXSC_SqrtPir[20];

       CXSC_SqrtPir_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_SqrtPir[i];
//       y[i+1] = CXSC_SqrtPir[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Sqrt2Pir[21]; // CXSC_Sqrt2Pir[0], ... CXSC_Sqrt2Pir[20] 
static bool CXSC_Sqrt2Pir_initialized = false;

l_interval Sqrt2Pir_l_interval() throw()
// Inclusion of 1/sqrt(2*Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Sqrt2Pir_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+19884533D43651e3FD";
       str >> CXSC_Sqrt2Pir[0];
       str = "-1CBC0D30EBFD15e3C7";
       str >> CXSC_Sqrt2Pir[1];
       str = "-1C7402C7D60CFBe38F";
       str >> CXSC_Sqrt2Pir[2];
       str = "+12706D8C0471B5e357";
       str >> CXSC_Sqrt2Pir[3];
       str = "-1FF6718B45881De321";
       str >> CXSC_Sqrt2Pir[4];
       str = "-13AABB82C248DCe2EB";
       str >> CXSC_Sqrt2Pir[5];
       str = "-1458A899162EE4e2B2";
       str >> CXSC_Sqrt2Pir[6];
       str = "-14EBD8868F41EBe27B";
       str >> CXSC_Sqrt2Pir[7];
       str = "+13278E993445F1e243";
       str >> CXSC_Sqrt2Pir[8];
       str = "-1CC019F5F4780Ae20D";
       str >> CXSC_Sqrt2Pir[9];
       str = "+147CE4B4ECDBD7e1D7";
       str >> CXSC_Sqrt2Pir[10];
       str = "-19A3DCC6A3534Be19F";
       str >> CXSC_Sqrt2Pir[11];
       str = "+11379A7BA8CB0Ae169";
       str >> CXSC_Sqrt2Pir[12];
       str = "-12D909C875312Ee132";
       str >> CXSC_Sqrt2Pir[13];
       str = "+1C1CEC4882C77Be0FB";
       str >> CXSC_Sqrt2Pir[14];
       str = "-14C4078263DF36e0C5";
       str >> CXSC_Sqrt2Pir[15];
       str = "+1AB3FC8D2AB243e08F";
       str >> CXSC_Sqrt2Pir[16];
       str = "+17B9172454310Ae059";
       str >> CXSC_Sqrt2Pir[17];
       str = "-1444B6B781B7F2e023";
       str >> CXSC_Sqrt2Pir[18];
       str = "-100001DB5C6774e000";
       str >> CXSC_Sqrt2Pir[19];
       str = "-100001DB5C6773e000";
       str >> CXSC_Sqrt2Pir[20];

       CXSC_Sqrt2Pir_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Sqrt2Pir[i];
//       y[i+1] = CXSC_Sqrt2Pir[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Pip2[21]; // CXSC_Pip2[0], ... CXSC_Pip2[20] 
static bool CXSC_Pip2_initialized = false;

l_interval Pip2_l_interval() throw()
// Inclusion of Pi^2, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Pip2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+13BD3CC9BE45DEe402";
       str >> CXSC_Pip2[0];
       str = "+1692B71366CC04e3CC";
       str >> CXSC_Pip2[1];
       str = "+18358E10ACD480e396";
       str >> CXSC_Pip2[2];
       str = "-1F2F5DD7997DDFe35F";
       str >> CXSC_Pip2[3];
       str = "+129E39B47B884Ee324";
       str >> CXSC_Pip2[4];
       str = "-12CF7459DD5DAFe2EE";
       str >> CXSC_Pip2[5];
       str = "-11842F87B5FE0Fe2B8";
       str >> CXSC_Pip2[6];
       str = "+1FFD8A79616A21e282";
       str >> CXSC_Pip2[7];
       str = "+12492A6663E899e24C";
       str >> CXSC_Pip2[8];
       str = "-1A15F4352CC511e215";
       str >> CXSC_Pip2[9];
       str = "-1301AA1792FF3Ce1DE";
       str >> CXSC_Pip2[10];
       str = "+122B6F31626EFEe1A8";
       str >> CXSC_Pip2[11];
       str = "+1B317FA13BDD8Fe172";
       str >> CXSC_Pip2[12];
       str = "+16F83B49040075e13C";
       str >> CXSC_Pip2[13];
       str = "-1B1890A945FE17e106";
       str >> CXSC_Pip2[14];
       str = "+12DCD389B96CDBe0D0";
       str >> CXSC_Pip2[15];
       str = "-1743F5DDE2F157e097";
       str >> CXSC_Pip2[16];
       str = "-153F96FFD4AEB5e060";
       str >> CXSC_Pip2[17];
       str = "+13CD6F5847D569e028";
       str >> CXSC_Pip2[18];
       str = "+10001471E79A7Be000";
       str >> CXSC_Pip2[19];
       str = "+10001471E79A8Be000";
       str >> CXSC_Pip2[20];

       CXSC_Pip2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Pip2[i];
//       y[i+1] = CXSC_Pip2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Sqrt2r[21]; // CXSC_Sqrt2r[0], ... CXSC_Sqrt2r[20] 
static bool CXSC_Sqrt2r_initialized = false;

l_interval Sqrt2r_l_interval() throw()
// Inclusion of 1/sqrt(2), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Sqrt2r_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+16A09E667F3BCDe3FE";
       str >> CXSC_Sqrt2r[0];
       str = "-1BDD3413B26456e3C8";
       str >> CXSC_Sqrt2r[1];
       str = "+157D3E3ADEC175e392";
       str >> CXSC_Sqrt2r[2];
       str = "+12775099DA2F59e35A";
       str >> CXSC_Sqrt2r[3];
       str = "+160CCE64552BF2e321";
       str >> CXSC_Sqrt2r[4];
       str = "+1821D5C5161D46e2E8";
       str >> CXSC_Sqrt2r[5];
       str = "-1C032046F8498Ee2B2";
       str >> CXSC_Sqrt2r[6];
       str = "+1EE950BC8738F7e27A";
       str >> CXSC_Sqrt2r[7];
       str = "-1AC3FDBC64E103e244";
       str >> CXSC_Sqrt2r[8];
       str = "+13B469101743A1e20C";
       str >> CXSC_Sqrt2r[9];
       str = "+15E3E9CA60B38Ce1D6";
       str >> CXSC_Sqrt2r[10];
       str = "+11BC337BCAB1BDe19B";
       str >> CXSC_Sqrt2r[11];
       str = "-1BBA5DEE9D6E7De165";
       str >> CXSC_Sqrt2r[12];
       str = "-1438DD083B1CC4e12F";
       str >> CXSC_Sqrt2r[13];
       str = "+1B56A28E2EDFA7e0F9";
       str >> CXSC_Sqrt2r[14];
       str = "+1CCB2A634331F4e0C3";
       str >> CXSC_Sqrt2r[15];
       str = "-1BD9056876F83Ee08C";
       str >> CXSC_Sqrt2r[16];
       str = "-1234FA22AB6BEFe056";
       str >> CXSC_Sqrt2r[17];
       str = "+19040CA4A81395e01F";
       str >> CXSC_Sqrt2r[18];
       str = "-10000015249C0Ce000";
       str >> CXSC_Sqrt2r[19];
       str = "-10000015249C0Be000";
       str >> CXSC_Sqrt2r[20];

       CXSC_Sqrt2r_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Sqrt2r[i]; 
//       y[i+1] = CXSC_Sqrt2r[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Sqrt3[21]; // CXSC_Sqrt3[0], ... CXSC_Sqrt3[20] 
static bool CXSC_Sqrt3_initialized = false;

l_interval Sqrt3_l_interval() throw()
// Inclusion of sqrt(3), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Sqrt3_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1BB67AE8584CAAe3FF";
       str >> CXSC_Sqrt3[0];
       str = "+1CEC95D0B5C1E3e3C9";
       str >> CXSC_Sqrt3[1];
       str = "-1F11DB689F2CCFe391";
       str >> CXSC_Sqrt3[2];
       str = "+13DA4798C720A6e35B";
       str >> CXSC_Sqrt3[3];
       str = "+121B9169B89243e325";
       str >> CXSC_Sqrt3[4];
       str = "-1813508751212Be2EC";
       str >> CXSC_Sqrt3[5];
       str = "-1B3D547B775C1Ee2B5";
       str >> CXSC_Sqrt3[6];
       str = "-19D986D92E2F0Ae27C";
       str >> CXSC_Sqrt3[7];
       str = "+1A34334CE806B6e245";
       str >> CXSC_Sqrt3[8];
       str = "+1A383B9E122E61e20F";
       str >> CXSC_Sqrt3[9];
       str = "+1C61D736F2F6F2e1D8";
       str >> CXSC_Sqrt3[10];
       str = "-10AF49233F9250e1A1";
       str >> CXSC_Sqrt3[11];
       str = "-1558A109EC0523e16A";
       str >> CXSC_Sqrt3[12];
       str = "+1F799D4D4FF2BCe134";
       str >> CXSC_Sqrt3[13];
       str = "-1AD7B219E34EDBe0FE";
       str >> CXSC_Sqrt3[14];
       str = "+15AB940B6677E3e0C8";
       str >> CXSC_Sqrt3[15];
       str = "-1D9B2A8203B8F0e091";
       str >> CXSC_Sqrt3[16];
       str = "-1DB0C8975A3834e05B";
       str >> CXSC_Sqrt3[17];
       str = "-1BCAAB3F6BE884e025";
       str >> CXSC_Sqrt3[18];
       str = "+100000531C2B6Ce000";
       str >> CXSC_Sqrt3[19];
       str = "+100000531C2B6De000";
       str >> CXSC_Sqrt3[20];

       CXSC_Sqrt3_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Sqrt3[i];
//       y[i+1] = CXSC_Sqrt3[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Sqrt3d2[21]; // CXSC_Sqrt3d2[0], ... CXSC_Sqrt3d2[20] 
static bool CXSC_Sqrt3d2_initialized = false;

l_interval Sqrt3d2_l_interval() throw()
// Inclusion of Sqrt(3)/2, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Sqrt3d2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1BB67AE8584CAAe3FE";
       str >> CXSC_Sqrt3d2[0];
       str = "+1CEC95D0B5C1E3e3C8";
       str >> CXSC_Sqrt3d2[1];
       str = "-1F11DB689F2CCFe390";
       str >> CXSC_Sqrt3d2[2];
       str = " +13DA4798C720A6e35A";
       str >> CXSC_Sqrt3d2[3];
       str = "+121B9169B89243e324";
       str >> CXSC_Sqrt3d2[4];
       str = " -1813508751212Be2EB";
       str >> CXSC_Sqrt3d2[5];
       str = "-1B3D547B775C1Ee2B4";
       str >> CXSC_Sqrt3d2[6];
       str = "-19D986D92E2F0Ae27B";
       str >> CXSC_Sqrt3d2[7];
       str = "+1A34334CE806B6e244";
       str >> CXSC_Sqrt3d2[8];
       str = "+1A383B9E122E61e20E";
       str >> CXSC_Sqrt3d2[9];
       str = "+1C61D736F2F6F2e1D7";
       str >> CXSC_Sqrt3d2[10];
       str = "-10AF49233F9250e1A0";
       str >> CXSC_Sqrt3d2[11];
       str = " -1558A109EC0523e169";
       str >> CXSC_Sqrt3d2[12];
       str = "+1F799D4D4FF2BCe133";
       str >> CXSC_Sqrt3d2[13];
       str = "-1AD7B219E34EDBe0FD";
       str >> CXSC_Sqrt3d2[14];
       str = "+15AB940B6677E3e0C7";
       str >> CXSC_Sqrt3d2[15];
       str = "-1D9B2A8203B8F0e090";
       str >> CXSC_Sqrt3d2[16];
       str = "-1DB0C8975A3834e05A";
       str >> CXSC_Sqrt3d2[17];
       str = "-1BCAAB3F6BE884e024";
       str >> CXSC_Sqrt3d2[18];
       str = "+100000298E15B6e000";
       str >> CXSC_Sqrt3d2[19];
       str = "+100000298E15B7e000";
       str >> CXSC_Sqrt3d2[20];

       CXSC_Sqrt3d2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Sqrt3d2[i];
//       y[i+1] = CXSC_Sqrt3d2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Sqrt3r[21]; // CXSC_Sqrt3r[0], ... CXSC_Sqrt3r[20] 
static bool CXSC_Sqrt3r_initialized = false;

l_interval Sqrt3r_l_interval() throw()
// Inclusion of 1/Sqrt(3), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Sqrt3r_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1279A74590331Ce3FE";
       str >> CXSC_Sqrt3r[0];
       str = "+134863E0792BEDe3C8";
       str >> CXSC_Sqrt3r[1];
       str = "-1A82F9E6C53222e392";
       str >> CXSC_Sqrt3r[2];
       str = "-1CB0F41134253Ae35C";
       str >> CXSC_Sqrt3r[3];
       str = "+1859ED919EC30Be326";
       str >> CXSC_Sqrt3r[4];
       str = "+1454874FB1F3F4e2EF";
       str >> CXSC_Sqrt3r[5];
       str = "-1DE69C6D3D2741e2B9";
       str >> CXSC_Sqrt3r[6];
       str = "+17EEC450C48BE1e283";
       str >> CXSC_Sqrt3r[7];
       str = "-16F743EEE65D53e24D";
       str >> CXSC_Sqrt3r[8];
       str = "-1887B505D7E7C2e215";
       str >> CXSC_Sqrt3r[9];
       str = "-1484D2E10C1161e1DE";
       str >> CXSC_Sqrt3r[10];
       str = "-1A0B1F86177FB7e1A8";
       str >> CXSC_Sqrt3r[11];
       str = "+1FE389D3F2C54Ee170";
       str >> CXSC_Sqrt3r[12];
       str = "+1F29F77C671544e13A";
       str >> CXSC_Sqrt3r[13];
       str = "-16CE74ED77D9BEe104";
       str >> CXSC_Sqrt3r[14];
       str = "-1E38708FF0CCB5e0CE";
       str >> CXSC_Sqrt3r[15];
       str = "-1F13BCC70157D1e098";
       str >> CXSC_Sqrt3r[16];
       str = "+17EC34CF9B1930e062";
       str >> CXSC_Sqrt3r[17];
       str = "-117A638EFF3A8Be02B";
       str >> CXSC_Sqrt3r[18];
       str = "-10016A8EF69C32e000";
       str >> CXSC_Sqrt3r[19];
       str = "-10016A8EF69C31e000";
       str >> CXSC_Sqrt3r[20];

       CXSC_Sqrt3r_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Sqrt3r[i]; 
//       y[i+1] = CXSC_Sqrt3r[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_LnPi[21]; // CXSC_LnPi[0], ... CXSC_LnPi[20] 
static bool CXSC_LnPi_initialized = false;

l_interval LnPi_l_interval() throw()
// Inclusion of ln(Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_LnPi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1250D048E7A1BDe3FF";
       str >> CXSC_LnPi[0];
       str = "+17ABF2AD8D5088e3C6";
       str >> CXSC_LnPi[1];
       str = "-16CCF43244818Ae38E";
       str >> CXSC_LnPi[2];
       str = "+1F9303719C0176e358";
       str >> CXSC_LnPi[3];
       str = "+15DF52611CB54Ee322";
       str >> CXSC_LnPi[4];
       str = "-1D9056E74F8C97e2EC";
       str >> CXSC_LnPi[5];
       str = "+100B095B6C2E1Ae2B5";
       str >> CXSC_LnPi[6];
       str = "-18C7557878A9E7e27F";
       str >> CXSC_LnPi[7];
       str = "+1B9BBBB4F4CEE7e248";
       str >> CXSC_LnPi[8];
       str = "+1B477FCC702F86e212";
       str >> CXSC_LnPi[9];
       str = "+141F1344A31799e1DC";
       str >> CXSC_LnPi[10];
       str = "+1B6740BE95CD58e1A6";
       str >> CXSC_LnPi[11];
       str = "-1F2C63904D27DBe16E";
       str >> CXSC_LnPi[12];
       str = "+1426F00B933976e136";
       str >> CXSC_LnPi[13];
       str = "+125703BE5FAA20e100";
       str >> CXSC_LnPi[14];
       str = "-1DADAE5397F95Be0C9";
       str >> CXSC_LnPi[15];
       str = "+17C9D110381543e091";
       str >> CXSC_LnPi[16];
       str = "-1259230E627FCAe05B";
       str >> CXSC_LnPi[17];
       str = "+191CEAB6B13A33e024";
       str >> CXSC_LnPi[18];
       str = "+10000109D49A14e000";
       str >> CXSC_LnPi[19];
       str = "+10000109D49A15e000";
       str >> CXSC_LnPi[20];

       CXSC_LnPi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_LnPi[i];
//       y[i+1] = CXSC_LnPi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Ln2Pi[21]; // CXSC_Ln2Pi[0], ... CXSC_Ln2Pi[20] 
static bool CXSC_Ln2Pi_initialized = false;

l_interval Ln2Pi_l_interval() throw()
// Inclusion of ln(2*Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Ln2Pi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1D67F1C864BEB5e3FF";
       str >> CXSC_Ln2Pi[0];
       str = "-165B5A1B7FF5DFe3C9";
       str >> CXSC_Ln2Pi[1];
       str = "-1B7F70C13DC1CCe392";
       str >> CXSC_Ln2Pi[2];
       str = "+13458B4DDEC6A3e35C";
       str >> CXSC_Ln2Pi[3];
       str = "+133DAA155D2130e324";
       str >> CXSC_Ln2Pi[4];
       str = "-18A007FC5E501Be2EE";
       str >> CXSC_Ln2Pi[5];
       str = "-15406FA3AA9644e2B4";
       str >> CXSC_Ln2Pi[6];
       str = "-13E8D52A392CC9e27E";
       str >> CXSC_Ln2Pi[7];
       str = "-1A43099131E88De248";
       str >> CXSC_Ln2Pi[8];
       str = "-114835B6623C4De212";
       str >> CXSC_Ln2Pi[9];
       str = "-1ABB7858CF827Ae1DC";
       str >> CXSC_Ln2Pi[10];
       str = "+1D8D7045A5A495e1A6";
       str >> CXSC_Ln2Pi[11];
       str = "+1A26094B3F6FC5e16F";
       str >> CXSC_Ln2Pi[12];
       str = "-1EF27932D0E3D0e137";
       str >> CXSC_Ln2Pi[13];
       str = "-12128804136AB6e100";
       str >> CXSC_Ln2Pi[14];
       str = "+15F8A4AC0BEE17e0C7";
       str >> CXSC_Ln2Pi[15];
       str = "+1892F2A5B69B5Fe091";
       str >> CXSC_Ln2Pi[16];
       str = "+1CC7C09477ADCEe05B";
       str >> CXSC_Ln2Pi[17];
       str = "-116DD579AF074Ae022";
       str >> CXSC_Ln2Pi[18];
       str = "+100000321C8783e000";
       str >> CXSC_Ln2Pi[19];
       str = "+100000321C8784e000";
       str >> CXSC_Ln2Pi[20];

       CXSC_Ln2Pi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Ln2Pi[i]; 
//       y[i+1] = CXSC_Ln2Pi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_E[21]; // CXSC_E[0], ... CXSC_E[20] 
static bool CXSC_E_initialized = false;

l_interval E_l_interval() throw()
// Inclusion of e=exp(1), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_E_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+15BF0A8B145769e400";
       str >> CXSC_E[0];
       str = "+14D57EE2B1013Ae3CA";
       str >> CXSC_E[1];
       str = "-1618713A31D3E2e392";
       str >> CXSC_E[2];
       str = "+1C5A6D2B53C26De35C";
       str >> CXSC_E[3];
       str = "-1F75CDE60219B6e326";
       str >> CXSC_E[4];
       str = "-188C76D93041A1e2EF";
       str >> CXSC_E[5];
       str = "+12FE363630C75Ee2B9";
       str >> CXSC_E[6];
       str = "-1C25F937F544EEe283";
       str >> CXSC_E[7];
       str = "-1E852C20E12A2Ae24D";
       str >> CXSC_E[8];
       str = "-14D4F6DE605705e212";
       str >> CXSC_E[9];
       str = "-1F3225EF539355e1D8";
       str >> CXSC_E[10];
       str = "-16109728625547e1A2";
       str >> CXSC_E[11];
       str = "-194301506D94CFe16C";
       str >> CXSC_E[12];
       str = "-1879C78F8CBA44e136";
       str >> CXSC_E[13];
       str = "-1D5976250C1018e0FD";
       str >> CXSC_E[14];
       str = "+1C877C56284DABe0C7";
       str >> CXSC_E[15];
       str = "+1E73530ACCA4F5e091";
       str >> CXSC_E[16];
       str = "-1F161A150FD53Ae05B";
       str >> CXSC_E[17];
       str = "+159927DB0E8845e022";
       str >> CXSC_E[18];
       str = "+10000094BB2C8Ee000";
       str >> CXSC_E[19];
       str = "+10000094BB2C8Fe000";
       str >> CXSC_E[20];

       CXSC_E_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_E[i]; 
//       y[i+1] = CXSC_E[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Er[21]; // CXSC_Er[0], ... CXSC_Er[20] 
static bool CXSC_Er_initialized = false;

l_interval Er_l_interval() throw()
// Inclusion of 1/e, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Er_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+178B56362CEF38e3FD";
       str >> CXSC_Er[0];
       str = "-1CA8A4270FADF5e3C6";
       str >> CXSC_Er[1];
       str = "-1837912B3FD2AAe390";
       str >> CXSC_Er[2];
       str = "-152711999FB68Ce35A";
       str >> CXSC_Er[3];
       str = "-17AD7C1289274Ee324";
       str >> CXSC_Er[4];
       str = "+17E8E56842B705e2E6";
       str >> CXSC_Er[5];
       str = "-1D24CB13796C2De2B0";
       str >> CXSC_Er[6];
       str = "-1456AABDA5C8F2e279";
       str >> CXSC_Er[7];
       str = "+1229F03C6276DDe243";
       str >> CXSC_Er[8];
       str = "-1569CFC4F53109e20D";
       str >> CXSC_Er[9];
       str = "-155B63C9B68091e1D5";
       str >> CXSC_Er[10];
       str = "+1580CF14DC087Ce19F";
       str >> CXSC_Er[11];
       str = "+1F9FF222313669e168";
       str >> CXSC_Er[12];
       str = "+15BC9CB1A22487e132";
       str >> CXSC_Er[13];
       str = "-1857E415C89B13e0FB";
       str >> CXSC_Er[14];
       str = "+13DF75706E3643e0C5";
       str >> CXSC_Er[15];
       str = "+13BDF5B7646234e08D";
       str >> CXSC_Er[16];
       str = "+1C956A5A3BE55De057";
       str >> CXSC_Er[17];
       str = "-167243FE9CD95Ee020";
       str >> CXSC_Er[18];
       str = "+1000002F30CCDBe000";
       str >> CXSC_Er[19];
       str = "+1000002F30CCDCe000";
       str >> CXSC_Er[20];

       CXSC_Er_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Er[i];
//       y[i+1] = CXSC_Er[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Ep2[21]; // CXSC_Ep2[0], ... CXSC_Ep2[20] 
static bool CXSC_Ep2_initialized = false;

l_interval Ep2_l_interval() throw()
// Inclusion of e^2, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Ep2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1D8E64B8D4DDAEe401";
       str >> CXSC_Ep2[0];
       str = "-19E62E22EFCA4Ce3CA";
       str >> CXSC_Ep2[1];
       str = "+1577508F5CF5EDe394";
       str >> CXSC_Ep2[2];
       str = "-186EF0294C2511e35E";
       str >> CXSC_Ep2[3];
       str = "+177D109F148782e327";
       str >> CXSC_Ep2[4];
       str = "+166BBC354AB700e2F0";
       str >> CXSC_Ep2[5];
       str = "-1273AEC0115969e2BA";
       str >> CXSC_Ep2[6];
       str = "-1C5AE00D3BEEF1e284";
       str >> CXSC_Ep2[7];
       str = "+15ACA3FDC9595Fe24C";
       str >> CXSC_Ep2[8];
       str = "-113FCDFE2B1F0Ce215";
       str >> CXSC_Ep2[9];
       str = "+10EEDFD1AE90C9e1DF";
       str >> CXSC_Ep2[10];
       str = "+1D2CB8EDC7078Be1A9";
       str >> CXSC_Ep2[11];
       str = "+11827A19F175F8e173";
       str >> CXSC_Ep2[12];
       str = "-10267512A9BFB2e13C";
       str >> CXSC_Ep2[13];
       str = "-19A1E2FC413AE3e105";
       str >> CXSC_Ep2[14];
       str = "+1170C7A5981ADBe0CF";
       str >> CXSC_Ep2[15];
       str = "-1FC991480067CFe099";
       str >> CXSC_Ep2[16];
       str = "-12E9A54CF5CFB5e062";
       str >> CXSC_Ep2[17];
       str = "-166FA6C468910Ae02A";
       str >> CXSC_Ep2[18];
       str = "+100043EA6DC142e000";
       str >> CXSC_Ep2[19];
       str = "+100043EA6DC143e000";
       str >> CXSC_Ep2[20];

       CXSC_Ep2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Ep2[i];
//       y[i+1] = CXSC_Ep2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Ep2r[21]; // CXSC_Ep2r[0], ... CXSC_Ep2r[20] 
static bool CXSC_Ep2r_initialized = false;

l_interval Ep2r_l_interval() throw()
// Inclusion of 1/e^2, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Ep2r_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1152AAA3BF81CCe3FC";
       str >> CXSC_Ep2r[0];
       str = "-1809224547B4BFe3C6";
       str >> CXSC_Ep2r[1];
       str = "-16A8E079134F13e390";
       str >> CXSC_Ep2r[2];
       str = "+14564CACF0994Ee358";
       str >> CXSC_Ep2r[3];
       str = "+1B796438129AF8e322";
       str >> CXSC_Ep2r[4];
       str = "-1ACFED57EF2AE5e2EC";
       str >> CXSC_Ep2r[5];
       str = "-1A968CBDBB5D9De2B5";
       str >> CXSC_Ep2r[6];
       str = "+1A7238CBD97B71e27C";
       str >> CXSC_Ep2r[7];
       str = "-146C53DB77BB01e245";
       str >> CXSC_Ep2r[8];
       str = "-1EEC161C3EBBD7e20C";
       str >> CXSC_Ep2r[9];
       str = "-12D084DC157ACEe1D5";
       str >> CXSC_Ep2r[10];
       str = "+12A61F46883347e19F";
       str >> CXSC_Ep2r[11];
       str = "+1993BAF10CAE0Be164";
       str >> CXSC_Ep2r[12];
       str = "+1F9224351178FFe12E";
       str >> CXSC_Ep2r[13];
       str = "-1C366D1C7BA64Ae0F7";
       str >> CXSC_Ep2r[14];
       str = "-17D9938EFA4657e0C0";
       str >> CXSC_Ep2r[15];
       str = "+1B6668DF0C1286e08A";
       str >> CXSC_Ep2r[16];
       str = "+1F7A4FFC9B48C6e050";
       str >> CXSC_Ep2r[17];
       str = "+1F3E3AF6F17591e01A";
       str >> CXSC_Ep2r[18];
       str = "+100000006C7831e000";
       str >> CXSC_Ep2r[19];
       str = "+100000006C7832e000";
       str >> CXSC_Ep2r[20];

       CXSC_Ep2r_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Ep2r[i];
//       y[i+1] = CXSC_Ep2r[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_EpPi[21]; // CXSC_EpPi[0], ... CXSC_EpPi[20] 
static bool CXSC_EpPi_initialized = false;

l_interval EpPi_l_interval() throw()
// Inclusion of e^(Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_EpPi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1724046EB0933Ae403";
       str >> CXSC_EpPi[0];
       str = "-184C962DD81952e3CD";
       str >> CXSC_EpPi[1];
       str = "-12D659C0BCD22Ee396";
       str >> CXSC_EpPi[2];
       str = "+117496B8A92F91e360";
       str >> CXSC_EpPi[3];
       str = "+16A8C4203E5FCDe32A";
       str >> CXSC_EpPi[4];
       str = "-166B11F99A663Be2F4";
       str >> CXSC_EpPi[5];
       str = "-118EC2076DABB1e2BE";
       str >> CXSC_EpPi[6];
       str = "+19776E5BEB18A5e288";
       str >> CXSC_EpPi[7];
       str = "+1AD4091E84B051e252";
       str >> CXSC_EpPi[8];
       str = "+1E89AA12909B40e21C";
       str >> CXSC_EpPi[9];
       str = "+1ACE3C0DDBB994e1E3";
       str >> CXSC_EpPi[10];
       str = "+141EC9379CBBFEe1AC";
       str >> CXSC_EpPi[11];
       str = "+1FC4E78D00A016e173";
       str >> CXSC_EpPi[12];
       str = "+1608BE35B9A409e13D";
       str >> CXSC_EpPi[13];
       str = "-1A0D8AA90EB6B9e103";
       str >> CXSC_EpPi[14];
       str = "+106FE8AFD21ACFe0CD";
       str >> CXSC_EpPi[15];
       str = "+1C072FEA1BFCAFe095";
       str >> CXSC_EpPi[16];
       str = "+1915B9F352EC68e05B";
       str >> CXSC_EpPi[17];
       str = "-13FA07C37897E9e024";
       str >> CXSC_EpPi[18];
       str = "-100003D8039138e000";
       str >> CXSC_EpPi[19];
       str = "-100003D8039137e000";
       str >> CXSC_EpPi[20];

       CXSC_EpPi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_EpPi[i]; 
//       y[i+1] = CXSC_EpPi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Ep2Pi[21]; // CXSC_Ep2Pi[0], ... CXSC_Ep2Pi[20] 
static bool CXSC_Ep2Pi_initialized = false;

l_interval Ep2Pi_l_interval() throw()
// Inclusion of e^(2*Pi), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Ep2Pi_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+10BBEEE9177E19e408";
       str >> CXSC_Ep2Pi[0];
       str = "+1C7DD9272526B1e3D0";
       str >> CXSC_Ep2Pi[1];
       str = "+15200F57AB89EDe39A";
       str >> CXSC_Ep2Pi[2];
       str = "-1FCCB6EDBE9C36e363";
       str >> CXSC_Ep2Pi[3];
       str = "+1BEA0BF179A589e32D";
       str >> CXSC_Ep2Pi[4];
       str = "-1F3AD5A6B77F9Ee2F7";
       str >> CXSC_Ep2Pi[5];
       str = "-1622F702B57637e2C0";
       str >> CXSC_Ep2Pi[6];
       str = "-100C09AE818734e287";
       str >> CXSC_Ep2Pi[7];
       str = "+10DA7ADA79EFE6e24D";
       str >> CXSC_Ep2Pi[8];
       str = "+1FF9BF48B72959e216";
       str >> CXSC_Ep2Pi[9];
       str = "-17AD7A3F6D2A14e1E0";
       str >> CXSC_Ep2Pi[10];
       str = "+1FCD4B0FA971E4e1A9";
       str >> CXSC_Ep2Pi[11];
       str = "+193A2CDC04526Be172";
       str >> CXSC_Ep2Pi[12];
       str = "-18CBE5FDFAF25Fe13C";
       str >> CXSC_Ep2Pi[13];
       str = "+1D47EEE171DA93e105";
       str >> CXSC_Ep2Pi[14];
       str = "-15B0F8DA29DB32e0CF";
       str >> CXSC_Ep2Pi[15];
       str = "-19207AD7E637D8e097";
       str >> CXSC_Ep2Pi[16];
       str = "+191CA743F265A6e061";
       str >> CXSC_Ep2Pi[17];
       str = "+1A15069182EF28e02A";
       str >> CXSC_Ep2Pi[18];
       str = "-1000EAC5FC05A9e000";
       str >> CXSC_Ep2Pi[19];
       str = "-1000EAC5FC05A8e000";
       str >> CXSC_Ep2Pi[20];

       CXSC_Ep2Pi_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Ep2Pi[i]; 
//       y[i+1] = CXSC_Ep2Pi[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_EpPid2[21]; // CXSC_EpPid2[0], ... CXSC_EpPid2[20] 
static bool CXSC_EpPid2_initialized = false;

l_interval EpPid2_l_interval() throw()
// Inclusion of e^(Pi/2), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_EpPid2_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+133DEDC855935Fe401";
       str >> CXSC_EpPid2[0];
       str = "+13E45A768FB73Ce3CB";
       str >> CXSC_EpPid2[1];
       str = "-1FB31CF300FF3Ce395";
       str >> CXSC_EpPid2[2];
       str = "-1E80D8BEB83F79e35F";
       str >> CXSC_EpPid2[3];
       str = "-14A3DE039142DDe326";
       str >> CXSC_EpPid2[4];
       str = "-18792D7A37282Be2EB";
       str >> CXSC_EpPid2[5];
       str = "-19DF43A5980C28e2B5";
       str >> CXSC_EpPid2[6];
       str = "-1C6F0F641C0D67e27F";
       str >> CXSC_EpPid2[7];
       str = "-1779C86C2DB5ACe249";
       str >> CXSC_EpPid2[8];
       str = "+168521EE91B16Fe211";
       str >> CXSC_EpPid2[9];
       str = "+12530F905D97BDe1DB";
       str >> CXSC_EpPid2[10];
       str = "+13498112CB7585e1A5";
       str >> CXSC_EpPid2[11];
       str = "+1BA4546B13A434e16F";
       str >> CXSC_EpPid2[12];
       str = "+14FF791C56421Ce138";
       str >> CXSC_EpPid2[13];
       str = "-1F375C223A2152e102";
       str >> CXSC_EpPid2[14];
       str = "-126AB0C8C77412e0CC";
       str >> CXSC_EpPid2[15];
       str = "-1B39C9C0B8C54Ae094";
       str >> CXSC_EpPid2[16];
       str = "-167741414E31E3e05D";
       str >> CXSC_EpPid2[17];
       str = "+1DEFB4462546C1e025";
       str >> CXSC_EpPid2[18];
       str = "-1000010F7B89CDe000";
       str >> CXSC_EpPid2[19];
       str = "-1000010F7B89CCe000";
       str >> CXSC_EpPid2[20];

       CXSC_EpPid2_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_EpPid2[i];
//       y[i+1] = CXSC_EpPid2[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_EpPid4[21]; // CXSC_EpPid4[0], ... CXSC_EpPid4[20] 
static bool CXSC_EpPid4_initialized = false;

l_interval EpPid4_l_interval() throw()
// Inclusion of e^(Pi/4), Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_EpPid4_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+118BD669471CAAe400";
       str >> CXSC_EpPid4[0];
       str = "+1F0ED609715756e3CA";
       str >> CXSC_EpPid4[1];
       str = "-1B9C7B871FE1DBe394";
       str >> CXSC_EpPid4[2];
       str = "+15C0FECE98F209e35D";
       str >> CXSC_EpPid4[3];
       str = "+18C9FACC5DF3CEe327";
       str >> CXSC_EpPid4[4];
       str = "+15EDE838B4A399e2EF";
       str >> CXSC_EpPid4[5];
       str = "-1C7EFACA363051e2B9";
       str >> CXSC_EpPid4[6];
       str = "-1A1EBEA1646411e283";
       str >> CXSC_EpPid4[7];
       str = "+1AEF54E68CE03Be24C";
       str >> CXSC_EpPid4[8];
       str = "-11250CB97FDDBFe212";
       str >> CXSC_EpPid4[9];
       str = "-169ADC0E65B8A7e1DB";
       str >> CXSC_EpPid4[10];
       str = "+198A501DB90EDDe1A5";
       str >> CXSC_EpPid4[11];
       str = "-1586909A3F6365e16E";
       str >> CXSC_EpPid4[12];
       str = "+1BE542410F8CE7e138";
       str >> CXSC_EpPid4[13];
       str = "+1E7EEC51889EECe102";
       str >> CXSC_EpPid4[14];
       str = "-1913C9FC19333Ce0CC";
       str >> CXSC_EpPid4[15];
       str = "+1112C71EA1E6F0e095";
       str >> CXSC_EpPid4[16];
       str = "-1C4CCF0F5D1E14e05E";
       str >> CXSC_EpPid4[17];
       str = "+1AC4A72310FA27e028";
       str >> CXSC_EpPid4[18];
       str = "-100013EC6A07AEe000";
       str >> CXSC_EpPid4[19];
       str = "-100013EC6A07ADe000";
       str >> CXSC_EpPid4[20];

       CXSC_EpPid4_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_EpPid4[i]; 
//       y[i+1] = CXSC_EpPid4[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_EulerGa[21]; // CXSC_EulerGa[0], ... CXSC_EulerGa[20] 
static bool CXSC_EulerGa_initialized = false;

l_interval EulerGa_l_interval() throw()
// Inclusion of EulerGamma, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_EulerGa_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+12788CFC6FB619e3FE";
       str >> CXSC_EulerGa[0];
       str = "-16CB90701FBFABe3C5";
       str >> CXSC_EulerGa[1];
       str = "-134A95E3133C51e38F";
       str >> CXSC_EulerGa[2];
       str = "+19730064300F7De359";
       str >> CXSC_EulerGa[3];
       str = "-171ECA0084E369e322";
       str >> CXSC_EulerGa[4];
       str = "-1302FE2B078898e2EC";
       str >> CXSC_EulerGa[5];
       str = "+192732D88415F4e2B5";
       str >> CXSC_EulerGa[6];
       str = "+11056AE9132136e27F";
       str >> CXSC_EulerGa[7];
       str = "-17DC6F12E630A3e249";
       str >> CXSC_EulerGa[8];
       str = "+175FD4B1BD70F2e212";
       str >> CXSC_EulerGa[9];
       str = "-19BC9466120C20e1DC";
       str >> CXSC_EulerGa[10];
       str = "-18FD5699260EADe1A6";
       str >> CXSC_EulerGa[11];
       str = "-12EA987665551Fe16F";
       str >> CXSC_EulerGa[12];
       str = "-1FB159BA4A423De138";
       str >> CXSC_EulerGa[13];
       str = "+1FA543D43BCC60e102";
       str >> CXSC_EulerGa[14];
       str = "-1E6F04E0F639F6e0C9";
       str >> CXSC_EulerGa[15];
       str = "-1A23768654F43De091";
       str >> CXSC_EulerGa[16];
       str = "-14F1C5CB4F55EBe058";
       str >> CXSC_EulerGa[17];
       str = "+1E71DF52EDAA7Fe020";
       str >> CXSC_EulerGa[18];
       str = "+1000001C398F9Be000";
       str >> CXSC_EulerGa[19];
       str = "+1000001C398F9Ce000";
       str >> CXSC_EulerGa[20];

       CXSC_EulerGa_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_EulerGa[i];
//       y[i+1] = CXSC_EulerGa[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


// **********************************************************************

static real CXSC_Catalan[21]; // CXSC_Catalan[0], ... CXSC_Catalan[20] 
static bool CXSC_Catalan_initialized = false;

l_interval Catalan_l_interval() throw()
// Inclusion of Catalan-Constant, Blomquist, 12.09.06;
{
   l_interval y;
   int stagsave = stagprec,
       stagmax = 20;
   
   if (!CXSC_Catalan_initialized) 
   {
       std::string str;
       std::cout << SaveOpt;
       std::cout << Hex;
       str = "+1D4F9713E8135De3FE";
       str >> CXSC_Catalan[0];
       str = "+11485608B8DF4De3C5";
       str >> CXSC_Catalan[1];
       str = "-12F39C13BC1EC8e38F";
       str >> CXSC_Catalan[2];
       str = "+1C2FF8094A263Ee357";
       str >> CXSC_Catalan[3];
       str = "+168F335DBE5370e321";
       str >> CXSC_Catalan[4];
       str = "+16291BBB16163Ee2E9";
       str >> CXSC_Catalan[5];
       str = "+124D663F739C43e2B3";
       str >> CXSC_Catalan[6];
       str = "-136A0725ED0E94e27B";
       str >> CXSC_Catalan[7];
       str = "-1D3A26F9C06FCEe240";
       str >> CXSC_Catalan[8];
       str = "-164E42486BFCD2e209";
       str >> CXSC_Catalan[9];
       str = "14F358CFDEC843e1D3";
       str >> CXSC_Catalan[10];
       str = "-11EB82210976ABe19D";
       str >> CXSC_Catalan[11];
       str = "-17D31F6DF5E801e167";
       str >> CXSC_Catalan[12];
       str = "+13FD19CE3E396Ae131";
       str >> CXSC_Catalan[13];
       str = "-1C8CBB3852FF3Fe0F8";
       str >> CXSC_Catalan[14];
       str = "+1A86EB34EAD01Ae0C2";
       str >> CXSC_Catalan[15];
       str = "+1C68C37800513Be087";
       str >> CXSC_Catalan[16];
       str = "+1D46EBB334D7C9e050";
       str >> CXSC_Catalan[17];
       str = "-1944C5E2711625e019";
       str >> CXSC_Catalan[18];
       str = "-100000005E2172e000";
       str >> CXSC_Catalan[19];
       str = "-100000005E2171e000";
       str >> CXSC_Catalan[20];

       CXSC_Catalan_initialized = true;
       std::cout << RestoreOpt;
   }
   stagprec = stagmax;
   y = adjust(l_interval(0.0)); 
   for (int i=0; i <= stagmax; i++)  
       y.data[i] = CXSC_Catalan[i];
//       y[i+1] = CXSC_Catalan[i];
   stagprec = stagsave;
   y = adjust(y);
   return y;             
}


l_interval atan(const l_interval & x) throw()  // aTan(x)
{
   int         stagsave = stagprec,
               stagmax = 19,
               digits = 53,
               m = 1,
               stagsave2,
               degree, sign, zk;
   real        lnt2, lneps, ln3m, eps,
   //            bas, tn, t3,
               zhn = 1.0,
               lnb = 0.69314718, // ln(2.0)
               ln3 = 1.098612289;  // ln(3.0), Blomquist 27.12.03
   l_interval  t, t2, p,
               y;
   interval    dx = interval(x),
               einfachgenau,
               err;

   einfachgenau = atan(dx);

   if (stagprec == 1) 
      y = atan(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = adjust(l_interval(0.0));
   else 
   {
       if (dx == interval(1.0)) y = li_pi4(); // y = Pi/4
       else 
       {
	   if (stagprec < stagmax) 
	       stagprec++;
	   else                    
	       stagprec = stagmax;

           // Argumentreduktion
	   t = x;
	   eps = 0.01/stagprec;    // vorher 0.1, aber zu gross
	   while (Sup(abs(interval(t))) > eps) 
	   {  // t = t/(1.0+sqrt(1.0+t*t));
	       t = t/(1.0+sqrt1px2(t)); // Blomquist, 13.12.02;
	       zhn += zhn;
	   }
	   t2 = t*t;

           // Bestimmung des Approximationsgrades
           // Beginn Blomquist, 27.12.03 ------------------------------
	   err = interval(abs(t2));
	   if (expo(Sup(err)) < -300) m = 4; // Vermeidet alte Fehlermeldg.
	   else 
	   {   
	       // lnt2 = Sup(ln(err)); // Erzeugt alten Fehler nicht mehr!! ALT!!
	       lnt2 = ln(Sup(err));    // Neu, 24.08.12; 
	       ln3m = ln(3.0-real(Sup(t2)));
	       lneps = (1.0-digits*stagprec)*lnb;
	       do {
		   m += 3;
	       } while (lneps-ln3-real(m+1)*lnt2+ln(2.0*m+3.0)+ln3m <= 0.0);
	   }
           // Ende Blomquist, 27.12.03 ---------------------------------

	   degree = m;

	   // Approximation
	   sign = (degree%2) ? -1 : 1;
	   zk = sign*(2*degree+1);
	   p = l_interval(1.0)/real(zk);
	   for (int k = degree-1; k >= 0; k--) 
	   {
	       sign = -sign;
	       zk = sign*(2*k+1);
	       p = p*t2+l_interval(1.0)/real(zk);
	   } // letztes t wird beim Fehler mit einmultipliziert

	   // Fehler bestimmen
	   stagsave2 = stagprec;
	   stagprec = 2;
	   l_interval relerr;
	   stagprec = stagsave2;
	   err = pow(interval(2.0), interval(1.0-digits*stagprec))
	       * interval(-1.0, 1.0);      // von 2^(2-d*s) ge\x{201E}ndert !
	   relerr[1] = 1.0;
	   relerr[2] = Inf(err);
	   relerr[3] = Sup(err);

           // Argumentreduktion rueckgaengig machen, mit t mult. 
           // und Fehler einbauen
	   y = zhn*t*p*relerr;

	   stagprec = stagsave;
	   y = adjust(y);
	   y = y & einfachgenau;
       }
   }

   return y;
}

l_interval acot(const l_interval &x) throw()  // ACot(x)
{
   l_interval  pihalbe, y;
   interval    dx = interval(x),
               einfachgenau;

   einfachgenau = acot(dx);
   stagprec++;
//   pihalbe = 2.0*atan(l_interval(1.0));
   pihalbe = li_pi4();
   times2pown(pihalbe,1); // Blomquist, 05.12.03;
   stagprec--;

   if (stagprec == 1)
      y = acot(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = adjust(pihalbe);
   else 
   {
      y = pihalbe-atan(x);
      // y = atan(1.0/x);  andere M�glichkeit der Berechnung
      y = adjust(y);
      y = y & einfachgenau;
   }

   return y;
}

l_interval exp(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW) // exp(x)
{
   long int    n = 2;
   int         stagsave = stagprec,
               stagmax = 19,
               degree,
               rednum = 0,
               digits = 53;
   real        fak = 2.0,
               zhn = 1.0,       // kann groesser 32768 werden
               tmp, lny, lneps, eps,
               lnb = 0.69314718;
   l_interval  p, t,
               y;
   interval    dx = interval(x),
               einfachgenau,
               dt, error;
   //            ilnb =interval(0.6931471, 0.6931572),
   //            ilny, ilneps;

   einfachgenau = exp(dx);

   // gueltigen Bereich pruefen
   if (stagprec == 1) 
      y = exp(dx);
   else if (Inf(dx) == Sup(dx) && Inf(dx) == CXSC_Zero)
      y = adjust(l_interval(1.0));
   else if (Inf(dx)<-708.396418532264) y = einfachgenau; // Blomquist,12.07.04
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;

      if (Sup(dx) <= 0.0)  
         t = -x;  // Werte unter -1e308 nicht auswertbar,
      else                 
         t =  x;  // dafuer andere Ergebnisse exakter
    
      dt = interval(t);

      // Argumentreduktion
      //    eps = 0.000001*real(Sup(t))/stagprec;
      //    eps = 1.0/pow10(stagprec/4);
      eps = 0.1;
      while(Sup(dt)/zhn > eps) 
      {
         rednum++;
         zhn += zhn;
      }
      t /= l_interval(zhn);

      // Taylor Approximation
      dt  = interval(t);

      // Anzahl der Durchlaeufe bei der Approximation
      tmp = Sup(abs(dt));
      if (MinReal > tmp) 
         tmp = MinReal;
      lny = ln(tmp);
      lneps = (1.0-digits*stagprec)*lnb;
      
      while(lneps-lnb+ln(fak)-real(n)*lny <= 0.0) 
      {
         n += 3;
         if (n > 170) 
         {    // 170! = 7.26*10306
            cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval exp(const l_interval & x)"));
            n = 170;
            break;
         }
         fak = fak*real(n)*(n-1.0)*(n-2.0);
       }
      degree = n;

      p = t/real(degree);
      int i;
      for (i=degree-1; i >= 1; i--)
         p = (p+1.0)*t/real(i);

      // Fehlerabschaetzung
      error = interval(-2.0, 2.0)*pow(dt, interval(real(n)))/fak;

      // Fehler einbauen
      p += 1.0+l_interval(error);

      // Argumentreduktion wird rueckgaengig gemacht
      for (i = 1; i <= rednum; i++)
         p *= p;

      if (Sup(dx) <= 0.0)  
         p = 1.0/p;

      stagprec = stagsave;
      y = adjust(p);
      y = y & einfachgenau;
   }

   return y;
}

l_interval exp2(const l_interval & x) // 2^x
{
	int stagsave = stagprec,
   stagmax = 19;
	if (stagprec>stagmax)
		stagprec = stagmax;
	l_interval y;
	
	y = exp(x*Ln2_l_interval());
	
	stagprec = stagsave;
	y = adjust(y);
	
	return y;
}

l_interval exp10(const l_interval & x) // 10^x
{
	int stagsave = stagprec,
   stagmax = 19;
	if (stagprec>stagmax)
		stagprec = stagmax;
   l_interval y;
	
   y = exp(x*Ln10_l_interval());
	
   stagprec = stagsave;
   y = adjust(y);
	
   return y;
}

l_interval expm1(const l_interval & x) throw()
// exp(x)-1;            
{
    l_interval y(0.0);
    real B,S_absdx,Sbru,abserr;
    interval dx = interval(x),
             einfachgenau,sx,bru;
    int stagsave = stagprec,ex,N, 
        stagmax = 19; // from exponential function 
    einfachgenau = expm1(dx); // produces all necessary error messages!

    if (stagprec == 1) y = einfachgenau;
    else 
    {
        if (stagprec>stagmax) stagprec = stagmax;
        S_absdx = Sup( abs(dx) );
        ex = expo(S_absdx);
        if (ex<-49) { // Taylor series
            // P_N(x) := x + x^2/2! + ... + x^N/N! 
	    sx = interval(S_absdx); // sx: point interval
	    N = ex-53*stagprec; // This N is here not the polynomial degree!
	    if (N<-1022) N = -1022;
	    if (N > 0) goto Fertig; // because S_absdx = 0;
	    B = comp(0.5,N); // B: rough bound of the absolute error
            // The absolute error should be smaller than B.
	    N=0;  bru = sx; // N is now the polynomial degree! 
            // Calculation of the polynomia degree N:
            // Calculation of the absolute error:
            // |absolute error| := |(exp(x)-1) - P_N(x)| <= AE;
            // AE = S_absdx^(N+1)/[(N+1)!*(1-S_absdx)]
	    do 
	    {
		N++;
		bru = (bru * S_absdx)/(N+1);
		Sbru = Sup(bru);
	    }
	    while(Sbru > B);
	    bru = bru * 1.00000001; // for  1/(1-S_absdx)
	    abserr = Sup(bru); 
            // |absolute error| <= abserr;
	    // Caculating an inclusion y of P_N(x):
	    y = x/N;
	    for (int i=N-1; i>=1; i--)
		y = (y+1)*x/i; 
            // x + x^2/2! + ... + x^N/N! = P_N(x) included by y
            // Conserning the absolute error:
	    y = y + interval(-abserr,abserr);
            
        } else {
	    stagprec++;
	    if (stagprec>stagmax) stagprec = stagmax;
	    y = exp(x) - 1;
	}
    Fertig:
	stagprec = stagsave; // restore old stagprec value
        y = adjust(y);  // adjust precision of y to old stagprec value
        y = y & einfachgenau;
    }
    return y;
} // expm1

l_interval expmx2(const l_interval& x)
// e^(-x^2);  Blomquist, 13.04.04;
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval y,z = abs(x);
    l_real r1,r2;
    interval dx = interval(z),
	     einfachgenau;
    einfachgenau = expmx2(dx);
    if (stagprec>stagmax) stagprec = stagmax;

    if (stagprec == 1) y = expmx2(dx); // simple precision 
    else if (Inf(z)>30) y = einfachgenau; // verhindert Overflow bei z*z!
    else y = exp(-z*z); 

    stagprec = stagsave;
    y = adjust(y);
    y = y & einfachgenau;
    return y;
} // expmx2




int cxsc_zweihoch(int) throw();

l_interval ln(const l_interval & x) throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF)    // Ln(x)
{ 
   int         stagsave = stagprec,
               stagmax = 19,
               k,
               m = 0,
               n = 0,
               digits = 53,   // Anzahl der Mantissenstellen
               cmp,
               zhn1 = 1,
               zhn2 = 1;
   bool        zwei=false;
   real        mx, tmp, 
   //            zkp1,
               ln2 = 0.69314718,
               lnb = 0.69314718,
               lny, lneps;
   l_interval  y, t, t2, p;
   interval    dx = interval(x),
               einfachgenau,
               dy, error;

   einfachgenau = ln(dx);
   // ---------------------------------------------------------------------
   if (Sup(dx) > succ(succ(Inf(dx)))) {  // Dieser Teil vermeidet eine
       y = einfachgenau; y = adjust(y);  // Fehlermeldung der pow(..)-Fkt,
       goto fertig;                      // wenn x zu breit ist und die 1
   }                                     // enthaelt. Bei zu breiten Inter-
   // vallen wird y daher mit einfachgenau berechnet.  Blomquist,05.12.03;
   // ---------------------------------------------------------------------
   if (Inf(x) <= 0.0) 
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval ln(const l_interval & x)"));
   else if (stagprec == 1) 
      y = ln(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_One)
      y = adjust(l_interval(0.0));
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;

      // erste Argumentreduktion
      /*
         if (Inf(dx) > 1.0) 
         {
            y = 1.0/x;
            dx = interval(x);
            neg = true;
         } else 
            y = x;
      */
      y = x;

      // zweite Argumentreduktion
      mx = Sup(dx);
      tmp = 1.0/Sup(dx);
      if (tmp > mx) 
         mx=tmp;
      tmp = 1.0/Inf(dx);
      if (tmp > mx) 
         mx = tmp;
      if (mx > 1.1) 
      {
         cmp = int(_double(ln(ln(mx)/ln(1.1))/ln2));
         if (cmp > 0) 
         {
            if (cmp > 14) 
               n = 14;
            else 
               n = cmp;
            zhn1 = cxsc_zweihoch(n);
         }
      }
//      y = sqrt(y, zhn1); // Old slow version
      if (zhn1 != 1)
	  for (int k=1; k<=n; k++)
	      y = sqrt(y); // Blomquist, 26.06.2007;
      mx = abs(real(Sup(y)));
      if (mx > 1.1) 
      {
         cmp = int(_double(ln(ln(mx)/ln(1.1))/ln2));
         if (cmp > 0) 
         {
            if (cmp > 14) 
               n = 14;
            else 
               n = cmp;
            zhn2 = cxsc_zweihoch(n);
         }
//         y = sqrt(y, zhn2); // Old slow version
	 if (zhn2 != 1)
	     for (int k=1; k<=n; k++)
		 y = sqrt(y); // Blomquist, 26.06.2007;
         zwei = TRUE;
      }

      // dritte Argumentreduktion
      t = (y-1.0)/(y+1.0);

      // Abschaetzung der Genauigkeit
      t2 = t*t;
      tmp = Sup(interval(t2));
      if (tmp == 0.0) 
         tmp = MinReal;
      dy = tmp;
      lny = Sup(ln(dy));
      lneps = (1.0-digits*stagprec)*lnb;
      
      do {
         m += 2;
      } while (lneps-lnb+ln(2.0*m+3.0)-real(m+1)*lny <= 0.0);

      // Polynomauswertung
      p = 0.0;
      for (k = 2*m+1; k >= 1; k -= 2)
         p = p*t2+l_interval(2.0)/real(k);
      p *= t;

      // Fehlereinschliessung
      error = interval(-4.0, 4.0)*pow(interval(t2), interval(real(m+1)))/(2.0*m+3.0);

      // Fehler einbauen
      p = p+error;

      // Argumentreduktion 2 wird rueckgaengig gemacht
      y = real(zhn1)*p;
      if (zwei) 
         y *= real(zhn2);

      // Argumentreduktion 1 wird rueckgaengig gemacht
      /*    if (neg) y = -(y);*/

      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;

      // Zusaetzliche Fehlerbehandlung noetig, da Underflow-Fehler 
      // nicht abgefangen werden kann
      /*    
         error = interval(-4.0*pow10(-stagprec*16-1), 4.0*pow10(-stagprec*16-1));
         y += error;
      */

   }
 fertig: // Blomquist,05.12.03;
   return y;
}

int cxsc_zweihoch(int n) throw()
{ 
   // liefert 2^n

   int res = 1;
   switch (n) 
   {
      case 0:
         return 1;
      case 1:
         return 2;
      case 2:
         return 4;
      default:
         {
            int y = 1,
            x = 2,
            zhi = 2;

            if (n%2)  
               res = 2;
            y = x*x;
            do {
               if ((n/zhi)%2)  
                  res *= y;
               zhi += zhi;
               if (zhi <= n)  
                  y *= y;
            } while (zhi <= n);
         }
   }
   return res;
}

l_interval log2(const l_interval & x) // log2(x)
{
	int stagsave = stagprec,
   stagmax = 19;
	if (stagprec>stagmax)
		stagprec = stagmax;
   l_interval y;
	
   y = ln(x)/Ln2_l_interval();
	
   stagprec = stagsave;
   y = adjust(y);
	
 return y;
}

l_interval log10(const l_interval & x) // log10(x)
{
	int stagsave = stagprec,
   stagmax = 19;
	if (stagprec>stagmax)
		stagprec = stagmax;
   l_interval y;
	
   y = ln(x)/Ln10_l_interval();
	
	stagprec = stagsave;
	y = adjust(y);
	
	return y;
}

l_interval sinh(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW)   // Sinh(x)
{
   long int    n = 1;
   int         stagsave = stagprec,
               stagmax = 19,
               sign = 1,
               digits = 53,
               degree, stagsave2;
   real        tmp,
               lnt, lneps,
               fak = 1.0,
               lnb = 0.69314718;
   l_interval  y,
               t, t2, p, pot;
   interval    dy,
               dx = interval(x),
               einfachgenau,
               ibase = interval(2.0),
   //            ilnb = interval(0.6931471,0.6931472),
   //            ilnx, ilneps, in,
               err;
   //            ifak = interval(fak);

   einfachgenau = sinh(dx);

   if (stagprec == 1) 
      y = sinh(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = x;
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;

      // Ausnuztung der Punktsymmetrie
      if (Sup(dx) < 0.0) 
      {
         y = -x;
         sign = -1;
      } else  
         y =  x;
    
      dy = interval(y);

      // Bei Werten �ber 0.5 -- Auswertung �ber e-Funktion
      if (Sup(dy) > 0.5) 
      {
         try {
            t = exp(y);
         }
         catch(const ERROR_LINTERVAL_FAK_OVERFLOW &)
         {
            cxscthrow( ERROR_LINTERVAL_FAK_OVERFLOW("l_interval sinh(const l_interval & x)")); 
         }
         y = sign*0.5*(t-1.0/t);
      }
      // Auswertung �ber Potenzreihenentwicklung
      else 
      {
         t = y;
         t2 = t*t;

         // Absch�tzung der Genauigkeit
         tmp = real(Sup(abs(t)));       // bei Andreas und Christian "t2"
         if (tmp < MinReal)  
            tmp = MinReal;
         lnt = ln(tmp);
         lneps = (1.0-digits*stagprec)*lnb;
         do {
            n += 3;
            if (n > 170) 
            {    // 170! = 7.26*10^306
               cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval sinh(const l_interval & x)"));
               n = 170;
               fak *= n*(n-1.0);
               break;
            }
            fak *= n*(n-1.0)*(n-2.0);
         } while(lneps+ln(fak)-real(n)*lnt-lnb <= 0.0);
         /*
            real bas, tn, t2d;
            bas = double(real(power(l_real(2.0),1-digits*stagprec)));
            tn = abs(double(real(mid(t))));
            t2d = double(real(mid(t2)));
            do {
               n += 2.0;
               if (n > 170) 
               {    // 170! = 7.26*10^306
                  errmon(ERROR_LINTERVAL(FAKOVERFLOW));
                  n = 170;
                  fak *= n;
                  break;
               }
               tn *= t2d;
               fak *= n*(n-1.0);
            } while(bas-tn/fak <= 0.0);
         */
         degree = 2*(n/2+1);  // degree ist gerade
                              // Achtung degree = 2n+2!

         // Auswertung mit Horner-Schema
         p = 1.0;
         for (int i = degree; i >= 2; i -= 2) 
         {
            pot = interval((i+1.0)*i);
            p = (p*(t2))/pot+1.0;
         }
         p *= t;

         // Negative Werte zur�ckwandeln
         if (sign == -1)  
            p = -p;

         // Fehlerauswertung
         stagsave2 = stagprec;
         stagprec = 2;
         l_interval  relerr;
         stagprec = stagsave2;
         err = pow(ibase , interval(1.0-digits*stagprec)) * interval(-1.0, 1.0);
         relerr[1] = 1.0;
         relerr[2] = Inf(err);
         relerr[3] = Sup(err);

         // Fehler einbauen
         y = p*relerr;
      }
      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;
   }
   return y;
}

l_interval cosh(const l_interval & x) throw(ERROR_LINTERVAL_FAK_OVERFLOW)   // Cosh(x)
{
   int         stagsave = stagprec,
               stagmax = 19;
   l_interval  y, t;
   interval    dx = interval(x),
               einfachgenau;

   einfachgenau = cosh(dx);

   if (stagprec == 1) 
      y = cosh(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = adjust(l_interval(1.0));
   else 
   {
      if (stagprec < stagmax) 
         stagprec++;
      else                    
         stagprec = stagmax;
      if (Sup(dx) <= 0.0)  
         y = -x;
      else                 
         y =  x;
      
      try 
      {
         t = exp(y);
      } 
      catch(const  ERROR_LINTERVAL_FAK_OVERFLOW &)
      {
         cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval cosh(const l_interval & x)"));
      }
      
      y = 0.5*(t+1.0/t);
      if (interval(0.0) <= dx) 
         SetInf(y,1.0);
      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;
   }

   return y;
}

// Fields for tanh-function
int  ex_tanh[8] = { 1,-1,-2,-4,-5,-6,-8,-9 };  
int  c_tanhN[8] = { 1,-1,2,-17,62,-1382,21844,-929569 };
int  c_tanhD[8] = { 1,3,15,315,2835,155925,6081075,638512875 }; 

l_interval tanh(const l_interval & x) throw()
// tanh(x) with Taylor-series and x --> MaxReal; Blomquist 24.05.04;           
{
    l_interval s,t,y;
    interval dx = interval(x),
             einfachgenau,r,r2;
    int stagsave = stagprec,ex,m,N,k,
	stagmax = 19;
    bool neg;
    real Sdx;
    einfachgenau = tanh(dx); 
    if (stagprec>stagmax) stagprec = stagmax;
    if (stagprec == 1) y = tanh(dx);
    else if (dx==0) y = 0;
    else if (0<=dx) y = einfachgenau;
    else // 0 not in x and not in dx:
    {   
	t = x;  
	neg = Sup(dx)<0;  Sdx = Sup(dx);
	if (neg) {t = -t; dx = -dx; Sdx = Sup(dx);} // Inf(dx) > 0:
	ex = expo(Sdx);
	if (ex < -70) {
	    m = ex - 53*stagprec; // m >= -1074 necessary !!
	    if (m < -1074) m = -1074;
	    r = interval(Sdx); 
	    r2 = r*r;
	    N = 0;
	    do {
		N++;
		r = r*r2;
		k = expo(Sup(r)) + ex_tanh[N];
	    } while (k>m);
	    r2 = interval(c_tanhN[N])/c_tanhD[N];
	    r2 = r * abs(r2); // r2: inclusion of the absolute error
	    N--;  // N: Polynomial degree
	    y = l_interval(c_tanhN[N])/c_tanhD[N];
	    s = t*t;
	    for (k=N-1; k>=0; k--)  // Horner-scheme
		y = y*s + l_interval(c_tanhN[k])/c_tanhD[k];
	    y = y * t;
	    y = y + interval(-Sup(r2),Sup(r2)); // Approxim. error
	} else if (ex<-4) {
	    s = sinh(t);  y = cosh(t);
	    if (stagprec<stagmax) stagprec++;
	    y = s/y;
	} else if (Sdx<352) {
	    times2pown(t,1);
	    y = 1 - 2/(exp(t)+1); // --> Sup(y) <= 1 !!
	} else {
	    if (Inf(dx)<352) {
		t = Inf(t);  times2pown(t,1);
		y = 1 - 2/(exp(t)+1);
	    } else y = l_interval(1) - comp(0.5,-1013);
	    SetSup(y,1); 
	}
	if (neg) y = -y;
    } // else

    stagprec = stagsave; // restore old stagprec value
    y = adjust(y);  // adjust precision of y to old stagprec value
    y = y & einfachgenau; // optimal inclusion in case of too
                          // large relative diameters of y;
    return y;
} // tanh

// ************************************************************************
// Fields for coth-function
int  ex_coth[8] = { -1,-5,-8,-12,-15,-18,-22,-25 };  
int  c_cothN[8] = { 1,-1,2,-1,2,-1382,4,-3617 };
real c_cothD[8] = { 1,15,315,1575,31185,212837625,
                    6081075,54273594375.0 }; // components exactly stored!
// ************************************************************************

l_interval coth(const l_interval & x) throw()
// coth(x); Blomquist 17.04.04;           
{
    l_interval t,s,c,y;
    interval dx = interval(x),
             einfachgenau,r,r2;
    int stagsave = stagprec,ex,m,N,k,
	stagmax = 19;
    bool neg;
    einfachgenau = coth(dx); // produces all necessary error messages!

    if (stagprec == 1) y = coth(dx);
    else 
    {
	if (stagprec>stagmax) stagprec = stagmax;
	neg = Sup(dx)<0.0;
	t = x;
	if (neg) { t = -t; dx = -dx; } // Inf(t),Inf(dx) > 0;
	ex = expo(Inf(dx));
	if (ex<-66) { // Laurent series
	    m = -ex - 53*stagprec;
	    r = interval(Sup(dx)); r2 = r*r;
	    N = 0;
	    do { 
		N++;
		if (N>7) { y = einfachgenau; goto Fertig; }
		r = r*r2;
		k = expo(Sup(r)) + ex_coth[N];
	    } while (k>m);
	    r2 = interval(c_cothN[N])/c_cothD[N]; 
	    r2 = r*abs(r2)/3;  // r2: inclusion of the absolute error
	    N--; // Polynomial degree
	    y = l_interval(c_cothN[N])/c_cothD[N]; 
	    s = t*t;
	    for (k=N-1; k>=0; k--)
		y = y*s + l_interval(c_cothN[k])/c_cothD[k];
	    y = (y*t)/3 + 1/t;
	    y = y + interval(-Sup(r2),Sup(r2)); // with approx. error
	} else if (ex < 2) {
	    if (stagprec<stagmax) stagprec++; // stagprec <= 19
	    c = cosh(t);
	    s = sinh(t);
	    y = c/s;
	} else if (Inf(dx)<353.0) {
	    if (stagprec<stagmax) stagprec++; // stagprec <= 19
	    times2pown(t,1);
	    y = 1.0 + 2.0/(exp(t)-1.0);
	} else { // Inf(dx) >= 353.0
	    y = l_interval(1.0) + comp(0.5,-1016);
	    SetInf(y,1.0);
	}
	if (interval(1.0) <= y)  SetInf(y,1.0);
	if (neg) y = -y;
    Fertig:	stagprec = stagsave; // restore old stagprec value
	y = adjust(y);  // adjust precision of y to old stagprec value
	y = y & einfachgenau;
    }
    return y;
} // coth


l_interval acosh(const l_interval & x) throw()  
// acosh(x); Blomquist 14.04.04;
{
    int         stagsave = stagprec,ex1,ex2,
                stagmax = 19;
    interval    dx = interval(x),
                einfachgenau;
    l_interval  y,t;

    einfachgenau = acosh(dx); // false definition range stops program
      
    if (stagprec == 1) y = acosh(dx);
    else if (Inf(dx) == Sup(dx) && Sup(dx) == 1.0)
	y = adjust(l_interval(0.0));
    else  
    {
	if (stagprec < stagmax) 
            stagprec++;
	else stagprec = stagmax;
	ex1 = expo(Inf(dx));  ex2 = expo(Sup(dx));
	if (ex1>500) { // acosh(x) approx ln(2) + ln(x):
	    y = li_ln2() + ln(x);
	    // Absolute approximation error:
	    t = 1.0/Inf(x);
	    y += l_interval(-Sup(t*t),0); // approxim. error realized here
	} else if (ex2<2) {
	    t = x-1;
	    y = lnp1( t+sqrt(t*(2+t)) );
	} else y = ln(x+sqrtx2m1(x));

	stagprec = stagsave;
	y = adjust(y);
	y = y & einfachgenau;
    } 
   return y;
}

l_interval asinh(const l_interval & x) 
    throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF,ERROR_LINTERVAL_FAK_OVERFLOW) 
// ASinh(x) hier ohne zeitaufwendiges Newton-Verfahren, da die Funktion 
// sqrt1px2(x) = sqrt(1+x*x) zur Verfuegung steht, die Overflow vermeidet!
// Blomquist, 28.12.03;
{
   int         stagsave = stagprec,
               stagmax = 19;
   l_interval  y, my;
   interval    dx = interval(x), einfachgenau, 
               error = abs(dx);
   real absSupdx = Sup(error), r;

   try 
   {
      einfachgenau = asinh(dx);

      if (stagprec == 1) 
            y = asinh(dx);
      else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
            y = x;
      else // if (absSupdx < MaxReal) 
      {
         if (stagprec < stagmax) stagprec++;
         else stagprec = stagmax; // Erhoete Praezision:

	 if (absSupdx < 2E-108) { // Taylorreihe von ln(x+sqrt(1+x2))
         // ln(x+sqrt(1+x2)) = x - x3/6 +- .....
	     y = x;
	     dx = interval(absSupdx);
	     error = dx*dx*dx/6;
	     r = Sup(error);  // r: Oberschranke des absoluten Fehlers
         // Jetzt folgt die Addition einer garantierten Einschliessung 
         // des abs. Fehlers:
	     y += l_interval(-r,r); 
	 } else
	     if (Sup(x) < 0.0) 
		 if (absSupdx < 1E+10) y = lnp1( -x+sqrtp1m1(x*x) );  
		 else y = -ln(-x+sqrt1px2(x));
	     else 
		 if (absSupdx < 1E+10) y = lnp1( x+sqrtp1m1(x*x) );
		 else y = ln(x+sqrt1px2(x));

         stagprec = stagsave;
         y = adjust(y);
	 y = y & einfachgenau;
      } 
   } // try 
   catch(const ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF &)
   {
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval asinh(const l_interval & x)"));
   }
   catch(const ERROR_LINTERVAL_FAK_OVERFLOW &)
   {
      cxscthrow(ERROR_LINTERVAL_FAK_OVERFLOW("l_interval asinh(const l_interval & x)"));
   }
   return y;
} // asinh(x)

l_interval atanh(const l_interval & x) 
   throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF,ERROR_LINTERVAL_FAK_OVERFLOW) 
// ATanh(x), Blomquist 29.12.03;
{
   int         stagsave = stagprec,
               stagmax = 19;
   l_interval  y, my;
   interval    dx = interval(x), 
               error = abs(dx),
               einfachgenau;
   real absSupdx = Sup(error);

   einfachgenau = atanh(dx);

   if (Inf(x) <= CXSC_MinusOne || Sup(x) >= CXSC_One) 
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("l_interval atanh(const l_interval & x)"));
   else if (stagprec == 1)
      y = atanh(dx);
   else if (Inf(dx) == Sup(dx) && Sup(dx) == CXSC_Zero)
      y = x;
   else 
   {
      if (stagprec < stagmax) stagprec++;
      else stagprec = stagmax; // Ab jetzt hoehere Praezision
      if (absSupdx < 0.125) {
	  y = x / (1-x);
	  times2pown(y,1);  // Schnelle Multiplikation mit 2
	  y = lnp1(y);
	  times2pown(y,-1); // Schnelle Division durch 2
      } else { 
	  y = ln((1.0+x)/(1.0-x));
	  times2pown(y,-1);
      }
      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;
   } 

   return y;
} // atanh

l_interval acoth(const l_interval & x) 
    throw(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF,ERROR_LINTERVAL_FAK_OVERFLOW)
// Acoth(x) hier ohne das langsame Newtonverfahren; Fue grosse x gilt:
// Acoth = 0.5*ln(1+2/(x-1));  Blomquist, 28.12.03;
{
   int         stagsave = stagprec,
               stagmax = 19;
   l_interval  y, my;
   interval    dx = interval(x),
               einfachgenau,
               error = abs(dx);
   real absSupdx = Sup(error);

   einfachgenau = acoth(dx);

   if ((l_interval(Inf(x))) < l_interval(CXSC_MinusOne,CXSC_One)
    || (l_interval(Sup(x))) < l_interval(CXSC_MinusOne,CXSC_One))
      cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
                       ("l_interval acoth(const l_interval & x)"));
   else if (stagprec == 1)
      y = acoth(dx);
   else  
   {
      if (stagprec < stagmax) stagprec++;
      else stagprec = stagmax; // Rechg. in erhoeter Praezision:

      if (absSupdx > 1E10) {
	  y = l_interval(2)/(x-1);
	  y = lnp1(y);
	  times2pown(y,-1);
      } else { 
	  y = ln((x+1.0)/(x-1.0));
	  times2pown(y,-1);
      }

      stagprec = stagsave;
      y = adjust(y);
      y = y & einfachgenau;
   } 

   return y;
} // acoth

l_interval lnp1(const l_interval& x) throw()
// ln(1+x) = zeta * P(zeta);  zeta := x/(2+x);
// P(zeta) := 2 + (2/3)zeta^2 + (2/5)zeta^4 + ...
// Approximation of P(zeta) with P_N(zeta), defined by
// P_N(zeta) := 2 + (2/3)zeta^(2*1) + 
//                + (2/5)zeta^(2*2) + ... + (2/2N+1)zeta^(2*N);
// Absolute approximation error Delta := P(zeta) - P_N(zeta);
// With z:=zeta^2  it holds:    Delta <= (2*z^(N+1))/{(2N+3)*(1-z)};
{
    const real c1(0.693147181);
    int stagsave(stagprec),
	stagmax(19);  // For ln-function
    if (stagprec>stagmax) stagprec = stagmax; // stagprec now <= 19;
    interval dx(x),    // dx of type interval inclusing x.
	einfachgenau(lnp1(dx)), // einfachgenau inclusing lnp1(x).
	ax( abs(dx) ); // ax inclusing the absolute values of dx.
    real gr(Sup(ax)),  // gr: maximum of the absolute values from x.
	lngr,u;
    l_interval t(x);
    if (gr < 1E-8)  // Series expansion, gr = 1E-8 is Ok!
	if (sign(gr)==0) t=0;
	else { // Calculation of N for the approximation of P(zeta)
               // using the polynomial P_N(zeta) of degree N: 
	    int N(0),TwoNp3(3),k;
	    k = expo(gr);
	    if (k>-1019) {
	    k = k - 1 - 53*stagprec;
	    lngr = ln(gr/2);
	    if (k < -1074) k = -1074;
	    k--; 
	    u = k*c1;  // u = (k-1)*ln(2);
	    do {
		N++;
		TwoNp3 += 2;
	    } while(TwoNp3*lngr - ln(TwoNp3) > u);
	    }  	
            // Polynomial degree 0<= N <=17 is now calculated!
            // Calculation of the N+1 polynomial coefficients:
	    l_interval a_lnp1[18], // Declaration of the coeffits.
		zeta,z2;
	    a_lnp1[0] = 2;
	    for (int i=1; i<=N; ++i) {
		a_lnp1[i] = a_lnp1[0]/(2*i + 1); 
	    }
            // Polyn.-coeff. a_lnp1[i] are now calculated; 0<=i<=N;
            // The calculation of P_N(zeta) follows now:
	    zeta = t/(2+t);
	    z2 = sqr(zeta);  
	    dx = z2; // is used for the approximation error!
	    t = a_lnp1[N]; // Horner-scheme:
	    for (int i=N-1; i>=0; --i) {
		t = t*z2 + a_lnp1[i];
	    } // t = P_N(zeta) is now calculated
            // Calculating now the approximation error 
            // with dx of type interval:
	    dx = interval(Sup(dx));
	    ax = Power(dx,N+1); 
	    times2pown(ax,1); // Multiplication with 2;
	    dx = ax / ((2*N+3)*(1-dx));
	    t = t + l_interval(0.0,Sup(dx)); 
            // Approximation error implemented now 
	    t = zeta * t;
	}
    else 
	if (gr<1) { // Calculation in higher accuracy:
	    stagprec++;
	    if (stagprec>stagmax) stagprec = stagmax;
	    t = ln(1+t);
	} else t = ln(1+t); // Normal accuracy:
    stagprec = stagsave;
    t = adjust(t);
    t = t & einfachgenau; // intersection of t and einfachgenau;
    // If x is too wide and therefore t because of the inevitable 
    // interval overestimations, t is the best possible inclusion.  
    return t;
}

l_interval sqrtx2m1(const l_interval& x)
// sqrt(x^2-1);  Blomquist, 13.04.04;
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval y,z = abs(x);
    l_real r1,r2;
    interval dx = interval(z),
	     einfachgenau;
    einfachgenau = sqrtx2m1(dx);
    if (stagprec>stagmax) stagprec = stagmax;

    if (stagprec == 1) y = sqrtx2m1(dx); // simple interval 
    else if (Inf(z)==1) {
	r1=0;  z = Sup(z);  z = sqrt(z*z-1);
	r2 = Sup(z);
	y = l_interval(r1,r2);
    } else if (expo(Sup(dx))<500) y = sqrt(z*z-1.0); // no overflow!
    else { // avoiding overflow using:  x-1/(2x) < sqrt(x^2-1) < x
	r2 = Sup(z);
	z = Inf(z);  // z: Point interval
	y = 1.0/z;
	times2pown(y,-1);
	y = z - y;
	r1 = Inf(y);
	y = l_interval(r1,r2);
    } 
    stagprec = stagsave;
    y = adjust(y);
    y = y & einfachgenau;
    return y;
} // sqrtx2m1

l_interval sqrt1mx2(const l_interval& x)
// sqrt(1-x^2);  Blomquist, 13.04.04;
{
    int stagsave = stagprec,
	stagmax = 19;
    l_interval y,z = abs(x);
    l_real r1,r2;
    interval dx = interval(z),
	     einfachgenau;
    einfachgenau = sqrt1mx2(dx);
    if (stagprec>stagmax) stagprec = stagmax;

    if (stagprec == 1) y = sqrt1mx2(dx); // simple interval 
    else {
	y = comp(0.5,1023); // y = 2^(+1022)
	times2pown(z,511);  
	y = sqrt(y-z*z);
	times2pown(y,-511);
    } 
    stagprec = stagsave;
    y = adjust(y);
    y = y & einfachgenau;
    return y;
} // sqrt1mx2


l_interval ln_sqrtx2y2(const l_interval& x, const l_interval& y) throw()
{ // Inclusion of ln(sqrt{x^2+y^2}) in staggered arithmetic
    int stagsave = stagprec;
//	stagmax = 19;  
    interval dx = interval(x),
	dy = interval(y),
	einfachgenau = ln_sqrtx2y2(dx,dy);
    dx = abs(dx);  dy = abs(dy);
    l_interval ax(abs(x)),ay(abs(y)),ar; 
    int ex_x(expo(Sup(dx))), ex_y(expo(Sup(dy))),
	N,N1;
    if (ex_y>ex_x) ex_x = ex_y;
    if (ex_x>508) { // To avoid overflow:
	N = ex_x-500;
	times2pown(ax,-N); times2pown(ay,-N);
	ar = ax*ax + ay*ay;
	ar = ln(ar);
	times2pown(ar,-1);
	ar += N*li_ln2();
    } else 
	if (ex_x<-20) { // To avoid underflow:
	    N = 500 - ex_x;
	    if (N>1023) { // Avoiding an range error with times2pown(...)  
		N1 = N-1023;
		times2pown(ax,1023);  times2pown(ax,N1); 
		times2pown(ay,1023);  times2pown(ay,N1);
	    } else {
		times2pown(ax,N); 
		times2pown(ay,N);
	    }
	    ar = ax*ax + ay*ay;
	    ar = ln(ar);
	    times2pown(ar,-1);
	    ar -= N*li_ln2();  // ar = ar - N*ln(2)
	} else { // Calculation with function lnp1:
	    ar = sqr(ax) + sqr(ay) - 1;
	    ar = lnp1(ar);
	    times2pown(ar,-1); 
	}
    stagprec = stagsave;
    ar = adjust(ar);
    ar = ar & einfachgenau;

    return ar;
}

l_interval acoshp1(const l_interval& x)
// acoshp1(x) = acosh(1+x);  Blomquist, 20.04.05;
{
    int stagsave = stagprec,ex,
	stagmax = 19;
    l_interval y,t;
    l_real lr;
    interval dx = interval(x),
	     einfachgenau;
    einfachgenau = acoshp1(dx);
    if (stagprec>stagmax) stagprec = stagmax;

    if (stagprec == 1) y = acoshp1(dx); // simple interval 
    else {
	ex = expo(Sup(dx));
	if (ex<=-1016) {
	    t = x;
	    times2pown(t,1);  // fast multiplication with 2 = 2^1;
	    t = sqrt(t);      // t = sqrt(2x);
	    lr = Inf(t);      // Calculating the infimum (Leibniz-series)
	    y = l_interval(lr)*(1.0 - l_interval(lr)/12.0);
	    y = l_interval(Inf(y),Sup(t)); // using alternating Leibniz-series
	}
	else 
	    if (ex<-400) y = lnp1(x*(1.0+sqrt(1+2.0/x)));
	    else y = acosh(1+x); // x = MaxReal not allowed
    } 
    stagprec = stagsave;
    y = adjust(y);
    y = y & einfachgenau;
    return y;
} // acoshp1


} // namespace cxsc

