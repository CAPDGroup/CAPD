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

/* CVS $Id: cinterval.cpp,v 1.20 2014/01/30 17:23:44 cxsc Exp $ */

#include "cinterval.hpp"
#include "cidot.hpp"
#include "dot.hpp"
#include "rmath.hpp"
#include "imath.hpp"

namespace cxsc {

#define CXSC_Zero 0.0

cinterval::cinterval(const dotprecision &a) throw() : re(a),im(0) {}
cinterval::cinterval(const idotprecision &a) throw() : re(a),im(0) {}
cinterval::cinterval(const cdotprecision &a) throw() : re(Re(a)),im(Im(a)) {}
cinterval::cinterval(const cidotprecision &a) throw() : 
   re(rnd(InfRe(a),RND_DOWN),rnd(SupRe(a),RND_UP)),
   im(rnd(InfIm(a),RND_DOWN),rnd(SupIm(a),RND_UP))  
{
}


cinterval mult_operator(const cinterval & a,const cinterval & b) throw()
{
   cidotprecision akku;
   akku=0.0;
   accumulate(akku,a,b);
   return rnd(akku);
} 
 
// In complex.cpp
void product  (real, real, real, real, int&, real&, interval&);
real quotient (real, interval, real, interval, int, int, int);

// optimale komplexe Intervalldivision

bool cxsc_complex_division_p[5];

real cxsc_complex_division_f(real a, real b, real c, real d, int round)
{
   int      zoverfl,noverfl;
   real     z1,n1;
   interval z2,n2;

   //  f:=(a*c+b*d)/(SQR(c)+SQR(d))

   product( a, c, b, d, zoverfl, z1,z2 );
   product( c, c, d, d, noverfl, n1,n2 );

   return quotient( z1,z2, n1,n2, round, zoverfl, noverfl );
}

static real minmax(int minimum, real a, real b, real y0,
                   interval x, int i, int j)
// Calculates the inner minimum or maximum of f = (ac+bd)/(cc+dd)
// on the interval x = [c.inf,c.sup] ( a,b,d=y0 fixated ).
// If minimum = true the minimum will be calculated, otherwise the maximum.

{
   real q1,q2,w1,w2,t1,t2,x1,x2,ay0, minmax;
   dotprecision akku;
   bool scaling = false;  // scaling = true <==> scaling is necessary

   if (minimum) 
      minmax = MaxReal;
   else         
      minmax = -MaxReal;

   if (Inf(x) == Sup(x)) 
   {
      if (cxsc_complex_division_p[i] && cxsc_complex_division_p[j])
         minmax = cxsc_complex_division_f( a, b, Inf(x), y0, 1-2*minimum );
      
      cxsc_complex_division_p[i] = false;
      cxsc_complex_division_p[j] = false;
  } else 
  if (a == 0.0) 
  {
      if ( b == CXSC_Zero  ||  y0 == CXSC_Zero ) 
      {
         minmax = 0.0;
         cxsc_complex_division_p[i]   = false;
         cxsc_complex_division_p[j]   = false;
      } else 
      {
         if (0.0 < x) {
            if (minimum  && sign(b) != sign(y0) ) 
            {
               minmax = divd(b, y0);
               cxsc_complex_division_p[i]   = false;
               cxsc_complex_division_p[j]   = false;
            } else 
            if (!minimum  &&  sign(b) == sign(y0) ) 
            {
               minmax =  divu(b, y0);
               cxsc_complex_division_p[i]   = false;
               cxsc_complex_division_p[j]   = false;
            }
	 }   
      }
   } else 
   { 
      //  a != 0.0
      if (y0 == 0.0) 
      {
         if (minimum) 
         {
            if (a > 0.0) 
               minmax = divd(a, Sup(x));
            else         
               minmax = divd(a, Inf(x));
         } else 
         {
            if (a > 0.0) 
               minmax = divu(a, Inf(x));
            else         
               minmax = divu(a, Sup(x));
         }
         cxsc_complex_division_p[i] = false;
         cxsc_complex_division_p[j] = false;
      } else 
      { 
         // y0 != 0.0, Calculation of extrema points and minimum|maximum
         // values: 
         //  IF NOTBEKANNT THEN
         //  Calculation of: t = sign(a) * ( |b/a| + sqrt( 1+|b/a|^2 ) )
         //  in staggered presentation:  t ~ t1 + t2

         // Exponent over-/undeflow in |b/a| is now considered from 
         // Blomquist/Hofschuster 06.11.02.

         real invf2=1, a_skal;
	 int exf1=0,  exf2=0,  exinf1=0,  exinf2=0;

         if (sign(b)==0) { t1 = 1;  t2 = 0; }
         else
	 {  // b != 0, ---> expo(b) != -2147483647;
         if (a>0.0) 
            q2 =  abs(b);
         else       
            q2 = -abs(b); 

         // Skaling to avoid overflow by division b/a: 
         int expo_diff = expo(b)-expo(a), ex;
         if (expo_diff >= 512) 
	 {
	     int exdiff5 = expo_diff - 500;
	     scaling = true;
	     if (exdiff5 > 1024) // Two scaling factors necessary!
	     {   
		 ex = exdiff5 / 2;
                 exf1 = ex-1;
		 exf2 = exdiff5 - ex;
		 exinf1 = -exf1;
		 invf2 = comp(0.5,1-exf2);
		 exinf2 = -exf2;
		 a_skal = a;
		 times2pown(a_skal,exf1);
                 times2pown(a_skal,exf2);
             }
	     else // Scaling with only one factor!
	     {   
		 exf2 = exdiff5-1;    // exf1 = 0;
		 invf2 = comp(0.5,2-exdiff5);  // invf1 = 1;
		 exinf2 = -exf2;  // exinf1 = 0;
		 a_skal = a;  
		 times2pown(a_skal,exf2);
             }
	 }
	 else // Scaling not necessary!
	 { a_skal = a; }

         q1   = q2/a_skal;
         akku = q2;
         accumulate(akku, -q1, a_skal);
         q2   = rnd(akku) / a_skal;

         akku = 0.0;
	 if (exinf1 == 0) accumulate(akku, invf2, invf2);	 
         accumulate(akku, q1, q1);
         accumulate(akku, q1, q2);
         accumulate(akku, q1, q2);
         accumulate(akku, q2, q2);

         w1   = sqrt(rnd(akku));

         accumulate(akku, -w1, w1);
         w2 = rnd(akku) / (2.0*w1);

         akku  = q1;
         akku += q2;
         akku += w1;
         akku += w2;

         t1 = rnd(akku);

         akku -= t1;
         t2 = rnd(akku);
	 }
	 if (a<0.0)  // if (a_skale<0.0)
         {
            t1 = -t1;
            t2 = -t2;
         }

         // Fall differentiation for min-,max- calculation:
         ay0 = abs(y0);
         if (( sign(b) == sign(y0) ) == minimum) 
         {
            //   Calculation of:  x1 + x2  =  |y0| * ( t1 + t2 )
            akku  = 0.0;
            accumulate(akku,ay0,t1);
            accumulate(akku,ay0,t2);
            x1    = rnd(akku);
	    if (expo(x1) == 2147483647) goto Ende;
            akku -= x1;
            x2    = rnd(akku);
	    if (scaling) 
            { 
		if (expo(x1)+exf1 > 1023) goto Ende;
		times2pown(x1,exf1);
		if (expo(x1)+exf2 > 1023) goto Ende;
                times2pown(x1,exf2); 
	        times2pown(x2,exf1);
                times2pown(x2,exf2);
            }
         } else 
         {
            //  Calculation of:  x1 + x2  =  |y0| / ( t1 + t2 )
            x1 = ay0 / t1;
            akku = ay0;
            accumulate(akku, -t1, x1);
            accumulate(akku, -t2, x1);
            x2 = rnd(akku) / t1;
	    if (scaling) 
            {  
		if (expo(x1)+exinf1 > 1023) goto Ende;
	        times2pown(x1,exinf1);
		if (expo(x1)+exinf2 > 1023) goto Ende;
                times2pown(x1,exinf2);
	        times2pown(x2,exinf1);
                times2pown(x2,exinf2);
            }
         }
         if (minimum) 
         {
            x1 = -x1;
            x2 = -x2;
         }

         if (x1 < x) 
         {
            //  Calculation of:  a / ( 2*(x1+x2) )
            q1   = a/(2*x1);
            akku = 0.0;
            accumulate(akku, -x1, q1);  // vorher: accumulate(akku, x1, q1);
            accumulate(akku, -x2, q1);

            // exact calculation of (a + akku + akku) in new variable akku:
	    akku += akku;  
	    akku += a; 
	    q2 = rnd(akku) / (2.0*x1);
            if (minimum) 
            {
		if (sign(q2)==0 && sign(akku)!=0) minmax = pred(q1);
		else minmax = addd(q1, q2);
	    }
            else 
	    {   
                if (sign(q2)==0 && sign(akku)!=0) minmax = succ(q1);
                else minmax = addu(q1, q2);
	    }
        
            cxsc_complex_division_p[i] = false;
            cxsc_complex_division_p[j] = false;
         }
      Ende:;
      }  // y0 != 0.0
   }
   return minmax;
} // *** minmax ***

cinterval cidiv(const cinterval& A, const cinterval& B)
{
   real     realteilINF, realteilSUP,
            imagteilINF, imagteilSUP;
   // da sonst eventuell zwischendurch illegale Intervalle entstehen
   real     a0,b0;
   bool     a_repeat,b_repeat; 
   int      i, rep, j;
   real     AREINF, ARESUP, AIMINF, AIMSUP,
            BREINF, BRESUP, BIMINF, BIMSUP;
   interval ARE, AIM, BRE, BIM;

   // keine Fehlerabfrage -> ist schon in CINTVAL.CPP
   //  IF ( 0.0 IN B.RE ) AND ( 0.0 IN B.IM ) THEN
   //   CIDIVISION:= COMPL( 1.0 / INTVAL(-1.0,1.0), INTVAL(0.0) );
   //   Fehlerabbruch erzwingen: Intervall enthaelt 0

   //  *** Berechnung des Realteils ***

   AREINF = Inf(Re(A));
   ARESUP = Sup(Re(A));
   AIMINF = Inf(Im(A));
   AIMSUP = Sup(Im(A));
   BREINF = Inf(Re(B));
   BRESUP = Sup(Re(B));
   BIMINF = Inf(Im(B));
   BIMSUP = Sup(Im(B));
   ARE    = Re(A);
   AIM    = Im(A);
   BRE    = Re(B);
   BIM    = Im(B);

   a_repeat = ( BREINF < CXSC_Zero ) && ( CXSC_Zero < BRESUP );
   b_repeat = ( BIMINF < CXSC_Zero ) && ( CXSC_Zero < BIMSUP );

   if (a_repeat || b_repeat) 
      rep = 2;
   else                      
      rep = 1;

   if (BREINF >= 0.0) 
      a0 = ARESUP;
   else               
      a0 = AREINF;
  
   if (BIMINF >= 0.0) 
      b0 = AIMSUP;
   else               
      b0 = AIMINF;

   realteilSUP = -MaxReal;

   for (j=1; j<=rep; j++) 
   {
      for (i=1; i<=4; cxsc_complex_division_p[i++] = true);
    
      realteilSUP =
             max( realteilSUP,
                   max( max( minmax( false, a0, b0, BIMINF, BRE, 1,2 ),
                             minmax( false, a0, b0, BIMSUP, BRE, 3,4 ) ),
                        max( minmax( false, b0, a0, BREINF, BIM, 1,3 ),
                             minmax( false, b0, a0, BRESUP, BIM, 2,4 ) ) )

                );
                
      if (cxsc_complex_division_p[1])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BREINF, BIMINF, +1 ) );
      if (cxsc_complex_division_p[2])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BRESUP, BIMINF, +1 ) );
      if (cxsc_complex_division_p[3])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BREINF, BIMSUP, +1 ) );
      if (cxsc_complex_division_p[4])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BRESUP, BIMSUP, +1 ) );

      if (a_repeat) 
         a0 = ARESUP; 
      else if (b_repeat) 
         b0 = AIMSUP;
   }

   if (BREINF >= 0.0) 
      a0 = AREINF;
   else               
      a0 = ARESUP;
   if (BIMINF >= 0.0) 
      b0 = AIMINF;
   else               
      b0 = AIMSUP;

   realteilINF = MaxReal;

   for (j=1; j<=rep; j++) 
   {
      for (i=1; i<=4; cxsc_complex_division_p[i++] = true);
      
      realteilINF =
              min( realteilINF,
                   min( min( minmax( true, a0, b0, BIMINF, BRE, 1,2 ),
                             minmax( true, a0, b0, BIMSUP, BRE, 3,4 ) ),
                        min( minmax( true, b0, a0, BREINF, BIM, 1,3 ),
                             minmax( true, b0, a0, BRESUP, BIM, 2,4 ) ) )
                 );
      if (cxsc_complex_division_p[1])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BREINF, BIMINF, -1 ) );
      if (cxsc_complex_division_p[2])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BRESUP, BIMINF, -1 ) );
      if (cxsc_complex_division_p[3])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BREINF, BIMSUP, -1 ) );
      if (cxsc_complex_division_p[4])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BRESUP, BIMSUP, -1 ) );

      if (a_repeat) 
         a0 = AREINF; 
      else if (b_repeat) 
         b0 = AIMINF;
   }


   //  Calculation of the img. part: 
   //  g(a, b, c, d) = cxsc_complex_division_f(b, -a, c, d) 

   a_repeat = ( BIMINF < CXSC_Zero ) && ( CXSC_Zero < BIMSUP );
   b_repeat = ( BREINF < CXSC_Zero ) && ( CXSC_Zero < BRESUP );

   //  IF a_repeat OR b_repeat THEN rep:= 2 ELSE rep:= 1;  

   if (BREINF >= 0.0) 
      b0 = AIMSUP;
   else 
      b0 = AIMINF;
   
   if (BIMINF >= 0.0) 
      a0 = AREINF;
   else 
      a0 = ARESUP;

   imagteilSUP = -MaxReal;

   for (j=1; j<=rep; j++) 
   {
      for (i=1; i<=4; cxsc_complex_division_p[i++] = true) ;
    
      imagteilSUP =
              max( imagteilSUP,
                   max( max( minmax( false,  b0, -a0, BIMINF, BRE, 1,2 ),
                             minmax( false,  b0, -a0, BIMSUP, BRE, 3,4 ) ),
                        max( minmax( false, -a0,  b0, BREINF, BIM, 1,3 ),
                             minmax( false, -a0,  b0, BRESUP, BIM, 2,4 ) ) )
                 );
    
      if (cxsc_complex_division_p[1])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BREINF, BIMINF, +1 ) );
      if (cxsc_complex_division_p[2])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BRESUP, BIMINF, +1 ) );
      if (cxsc_complex_division_p[3])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BREINF, BIMSUP, +1 ) );
      if (cxsc_complex_division_p[4])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BRESUP, BIMSUP, +1 ) );

      if (b_repeat) 
         b0 = AIMSUP;  
      else if (a_repeat) 
         a0 = AREINF  ;
   }

   if (BREINF >= 0.0) 
      b0 = AIMINF;
   else 
      b0 = AIMSUP;
  
   if (BIMINF >= 0.0) 
      a0 = ARESUP;
   else 
      a0 = AREINF;

   imagteilINF = MaxReal;

   for (j=1; j<=rep; j++) 
   {
      for (i=1; i<=4; cxsc_complex_division_p[i++] = true) ;
    
      imagteilINF =
              min( imagteilINF,
                   min( min( minmax( true,  b0, -a0, BIMINF, BRE, 1,2 ),
                             minmax( true,  b0, -a0, BIMSUP, BRE, 3,4 ) ),
                        min( minmax( true, -a0,  b0, BREINF, BIM, 1,3 ),
                             minmax( true, -a0,  b0, BRESUP, BIM, 2,4 ) ) )
                 );
    
      if (cxsc_complex_division_p[1])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BREINF, BIMINF, -1 ) );
      if (cxsc_complex_division_p[2])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BRESUP, BIMINF, -1 ) );
      if (cxsc_complex_division_p[3])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BREINF, BIMSUP, -1 ) );
      if (cxsc_complex_division_p[4])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BRESUP, BIMSUP, -1 ) );

      if (b_repeat) 
         b0 = AIMINF;  
      else if (a_repeat) 
         a0 = ARESUP;
   }

   return cinterval(interval(realteilINF, realteilSUP),
                     interval(imagteilINF, imagteilSUP));
}  //    CIDIVISION

cinterval C_point_div(const cinterval& z, const cinterval& n)
// Division of complex point intervals; 
// z,n must be point intervals!!  Blomquist, 07,11.02
// This function only for internal use!
{
    complex a,b,q1,q2;
    a = complex(InfRe(z),InfIm(z));
    b = complex(InfRe(n),InfIm(n));
    q1 = divd(a,b);
    q2 = divu(a,b);

    interval re, im;
    re = interval( Re(q1),Re(q2) );
    im = interval( Im(q1),Im(q2) );

    return cinterval(re,im);
}  // C_point_div

cinterval div_operator (const cinterval & a, const cinterval & b) throw(DIV_BY_ZERO)
{
    bool a_point, b_point;
    a_point = InfRe(a)==SupRe(a) && InfIm(a)==SupIm(a);
    b_point = InfRe(b)==SupRe(b) && InfIm(b)==SupIm(b);
    if(a_point && b_point) return C_point_div(a,b); // a,b are point intervals
    else return cidiv(a,b);
}

interval abs(const cinterval &a) throw()
{
//    idotakku[2]=0;
//    accumulate(idotakku[2],a.re,a.re);
//    accumulate(idotakku[2],a.im,a.im);
//    return sqrt(rnd(idotakku[2]));
    return sqrtx2y2(a.re,a.im);
}


// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const cinterval& a) throw()
{
   s << '('          
     << a.re << ','  
     << a.im       
     << ')';
   return s;
}
std::string & operator << (std::string &s, const cinterval &a) throw()
{
   s+='(';
   s << a.re;
   s+=',';
   s << a.im; 
   s+=')';
   return s;
}

std::istream & operator >> (std::istream &s, cinterval &a) throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{ // New version for cinterval input; Blomquist, 27.10.02;
    char c;
    skipeolnflag = inpdotflag = true;
    c = skipwhitespacessinglechar (s, '(');
    if (inpdotflag) s.putback(c);
    c = skipwhitespacessinglechar (s, '[');
    if (inpdotflag) s.putback(c);
    s >> SaveOpt >> RndDown >> Inf(a.re);
    skipeolnflag = inpdotflag = true; 
    c = skipwhitespacessinglechar (s, ','); 
    if (inpdotflag) s.putback(c);
    s  >> RndUp >> Sup(a.re);
    c = skipwhitespacessinglechar (s, ']');
    if (inpdotflag) s.putback(c);
    c = skipwhitespacessinglechar (s, ',');
    if (inpdotflag) s.putback(c);
    
    c = skipwhitespacessinglechar (s, '[');
    if (inpdotflag) s.putback(c);
    s >> RndDown >> Inf(a.im);
    skipeolnflag = inpdotflag = true; 
    c = skipwhitespacessinglechar (s, ','); 
    if (inpdotflag) s.putback(c);
    s  >> RndUp >> Sup(a.im) >> RestoreOpt;

   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ']') 
         s.putback(c);
   }
   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ')') 
         s.putback(c);
	 }
   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
   cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("std::istream & operator >> (std::istream &s, cinterval &a)"));
      
   return s;
}

std::string & operator >> (std::string &s, cinterval &a) throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   s = skipwhitespacessinglechar (s, '(');
   s = skipwhitespacessinglechar (s, '[');
   s = s >> SaveOpt >> RndDown >> Inf(a.re);
   s = skipwhitespacessinglechar (s, ',');
   s = s >> RndUp >> Sup(a.re);
   s = skipwhitespacessinglechar (s, ']');
   s = skipwhitespacessinglechar (s, ',');
   s = skipwhitespacessinglechar (s, '[');
   s = s >> RndDown >> Inf(a.im);
   s = skipwhitespacessinglechar (s, ',');
   s = s >> RndUp >> Sup(a.im) >> RestoreOpt;
   s = skipwhitespaces (s);
   if (s[0] == ']') 
      s.erase(0,1);
   s = skipwhitespaces (s);
   if (s[0] == ')') 
      s.erase(0,1);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("std::string & operator >> (std::string &s, cinterval &a)"));

   return s;
}

void operator >>(const std::string &s,cinterval &a) throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   std::string r(s);
   r>>a;
}
void operator >>(const char *s,cinterval &a) throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   std::string r(s);
   r>>a;
}

int in ( const cinterval& x, const cinterval& y )    // Contained-in-the-interior relation
{                                        //-----------------------------------
  return ( in(Re(x),Re(y)) && in(Im(x),Im(y)) );
}
/*!
\param x The complex interval for which the epsilon inflation should be computed
\param eps The real value of epsilon
\return The inflated complex interval

\sa cxsc::Blow(const interval& x, const real& eps )
*/
cinterval Blow ( cinterval x, const real& eps )           // Epsilon inflation
{                                                         //------------------
  return cinterval(Blow(Re(x),eps),Blow(Im(x),eps));
}
} // namespace cxsc















