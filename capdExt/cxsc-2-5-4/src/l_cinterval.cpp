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

/* CVS $Id: l_cinterval.cpp,v 1.18 2014/01/30 17:23:46 cxsc Exp $ */

#include "l_cinterval.hpp"
#include "cidot.hpp"
#include "dot.hpp"
#include "l_rmath.hpp"
#include "l_imath.hpp"

namespace cxsc {


#define CXSC_Zero 0.0

cinterval::cinterval(const l_cinterval & a) throw()
{
    interval u,v;
    u = a.re;
    v = a.im;
    *this = cinterval(u,v);
}

cinterval & cinterval::operator = (const l_cinterval & a) throw()
{
    interval u,v;
    u = a.re;
    v = a.im;
return *this = cinterval(u,v);
}

l_cinterval::l_cinterval(const dotprecision  &a) throw() : re(a),im(0) {}
l_cinterval::l_cinterval(const idotprecision &a) throw() : re(a),im(0) {}
l_cinterval::l_cinterval(const cdotprecision &a) 
                                         throw() : re(Re(a)),im(Im(a)) {}
l_cinterval::l_cinterval(const cidotprecision &a) throw() : 
                         re( l_interval(Re(a))),im(l_interval(Im(a)) ) {}

l_cinterval operator * (const l_cinterval & a, const l_cinterval & b) throw()
{
    idotprecision akku;
    l_cinterval res;
    l_interval u,v;
    akku = 0.0;
    accumulate(akku,a.re,b.re);
    accumulate(akku,-a.im,b.im);
    u = akku;
    if (Inf(u)>Sup(u)) 
    {   std::cerr << "Error in l_cinterval * l_cinterval" << std::endl;
	exit(1);
    }
    akku = 0.0;
    accumulate(akku,a.im,b.re);
    accumulate(akku,a.re,b.im);
    v = akku; // v: Imaginaerteil
    if (Inf(v)>Sup(v)) 
    {   std::cerr << "Error in l_cinterval * l_cinterval" << std::endl;
	exit(1);
    }
    res = l_cinterval(u,v);
    return res;
}


// *********************************************************************
// In l_complex.cpp implemented: (In l_complex.hpp not declared!)
void product(const l_real& a, const l_real& b, const l_real& c, 
             const l_real& d, int& ex, l_interval& res);
void product(const l_real& c, const l_real& d, int& ex, l_interval& res);
l_real quotient(const l_interval& z, const l_interval& n, int round, 
		int ex_z, int ex_n);
//void Times2pown(l_interval& a, int n) throw();

// *********************************************************************

static const int max_expo  = 1020, max_expo1 = 1023;

// optimale komplexe Intervalldivision

bool cxsc_l_complex_division_p[5];

l_real cxsc_complex_division_f(l_real a, l_real b, l_real c, l_real d, 
                               int round)
{
   int ex1, ex2;
   l_interval z,n;

   //  f:=(a*c+b*d)/(SQR(c)+SQR(d))

   product(a, c, b, d, ex1, z);
   product(c, d, ex2, n);
   return quotient(z, n, round, ex1, ex2);
}

// *************************************************************************
// *************************************************************************

static l_real minmax(int minimum, l_real a, l_real b, l_real y0,
                     l_interval x, int i, int j)
// Calculates the inner minimum or maximum of f = (ac+bd)/(cc+dd)
// on the interval x = [c.inf,c.sup] ( a,b,d=y0 fixated ).
// If minimum = true the minimum will be calculated, otherwise the maximum.

{
   l_real ay0, minmax;
   l_interval t,q,x0,two_Da;

   int Da(0);

   a += 0.0;   b += 0.0;  y0 += 0.0; 

   if (minimum) 
      minmax = MaxReal;
   else         
      minmax = -MaxReal;

   if (Inf(x) == Sup(x)) 
   {
      if (cxsc_l_complex_division_p[i] && cxsc_l_complex_division_p[j])
         minmax = cxsc_complex_division_f( a, b, Inf(x), y0, 1-2*minimum );
      
      cxsc_l_complex_division_p[i] = false;
      cxsc_l_complex_division_p[j] = false;
  } else 
  if (a == 0.0) 
  {
      if ( b == CXSC_Zero  ||  y0 == CXSC_Zero ) 
      {
         minmax = 0.0;
         cxsc_l_complex_division_p[i]   = false;
         cxsc_l_complex_division_p[j]   = false;
      } else 
      { // b*y0 <> 0:
         if (0.0 < x) {
            if (minimum  && sign(b) != sign(y0) ) 
            {
		// minmax = divd(b, y0);
		minmax = Inf(l_interval(b)/y0);
		cxsc_l_complex_division_p[i]   = false;
		cxsc_l_complex_division_p[j]   = false;
            } else 
            if (!minimum  &&  sign(b) == sign(y0) ) 
            {
		// minmax =  divu(b, y0);
		minmax = Sup(l_interval(b)/y0);
		cxsc_l_complex_division_p[i]   = false;
		cxsc_l_complex_division_p[j]   = false;
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
		// minmax = divd(a, Sup(x));
		minmax = Inf(l_interval(a)/l_interval(Sup(x)));
            else         
		// minmax = divd(a, Inf(x));
		minmax = Inf(l_interval(a)/l_interval(Inf(x)));
         } else 
         {
            if (a > 0.0) 
		// minmax = divu(a, Inf(x));
		minmax = Sup(l_interval(a)/l_interval(Inf(x)));
            else         
		// minmax = divu(a, Sup(x));
		minmax = Sup(l_interval(a)/l_interval(Sup(x)));
         }
         cxsc_l_complex_division_p[i] = false;
         cxsc_l_complex_division_p[j] = false;
      } else 
      {  // a  != 0.0 and
         // y0 != 0.0, Calculation of extrema points and minimum|maximum
         // values: 
         //  Calculation of: t = sign(a) * ( |b/a| + sqrt( 1+|b/a|^2 ) )

         l_real invf2(1.0), a_skal;
	 // int exf1=0,  exf2=0,  exinf1=0,  exinf2=0; // unused variable

         // We first calculate:   t = |b/a| + sqrt( 1+|b/a|^2 );

         if (sign(b)==0) t = 1.0;
         else
	 {   // a != 0.0  and  b != 0;
             // To avoid overflow by calculating |b/a| + sqrt( 1+|b/a|^2 )
             // we must multiply a with 2^Da:
	     int expo_diff = expo(b[1]) - expo(a[1]), ex;
	     a_skal = a;
	     if (expo_diff > max_expo) 
	     {
		 Da = expo_diff-max_expo; // Da > 0;
                 // a must be multiplied with 2^Da to avoid overflow 
                 // by calculating |b/a| + sqrt( 1+|b/a|^2 ) :
		 if (Da>max_expo1)
		 {
		     times2pown(a_skal,max_expo1);
		     ex = Da - max_expo1;
		     times2pown(a_skal,ex);
		 }
		 else times2pown(a_skal,Da);

                 // Now calculating an inclusion t of 2^(-Da):
		 if (Da>1022)
		 {
		     two_Da = l_interval( comp(0.5,-1021) ); 
		     times2pown(two_Da,1022-Da);     
		 }
		 else two_Da = l_interval( comp(0.5,1-Da) ); 
                 // Now two_Da is am inclusion of 2^(-Da);
	     }
	     q = l_interval(b)/a_skal;
	     if (sign(q[1])<0) q = -q;
	     // q: Inclusion of |b/(a*2^Da)|;

	     t = (Da > 0)? q + sqrtx2y2(two_Da,q) : q + sqrt1px2(q);
	 }

	 if (a<0.0)  t = -t;

// if (Da > 0) the value t from the last line must additionally be 
// multiplied with 2^Da, to get an inclusion of the expression:    
//               sign(a) * ( |b/a| + sqrt( 1+|b/a|^2 ) );

// Now to a fall differentiation for min-,max- calculation:
// First we will calculate an inclusion x0 of the point of the
// relative minimum or maximum:

         ay0 = abs(y0);

         if ( (sign(b) == sign(y0)) == minimum ) 
         {   // Calculation of x0 = |y0| * t :
	     if (expo(ay0[1]) + expo(t[1]) + Da > max_expo1) goto Ende;
	     else x0 = ay0 * t;
	     if (Da>0) Times2pown(x0,Da);
         } 
	 else  //  Calculation of x0 = |y0| / t :
         {   
	     if (expo(ay0[1]) - expo(t[1]) - Da > max_expo1) goto Ende;
	     else x0 = ay0 / t;             
	     if (Da>0) Times2pown(x0,-Da); 
         }

         if (minimum) x0 = -x0;

         if (x0 < x) // The minimum or maximum point lies in 
         {           // the interior of x.
	     //  Calculation of:  a / ( 2*x0 )
	     q = a/x0;
	     times2pown(q,-1); // q: inclusion of a / ( 2*x0 );

	     if (minimum) minmax = Inf(q);
	     else minmax = Sup(q);

	     cxsc_l_complex_division_p[i] = false;
	     cxsc_l_complex_division_p[j] = false;
         }
      Ende:;
      }  // y0 != 0.0
  }
   return minmax;
} // *** minmax ***

l_real max(const l_real& u, const l_real& v)
{
    l_real res(u);
    if (v>u) res = v;
    return res;
}

l_real min(const l_real& u, const l_real& v)
{
    l_real res(u);
    if (v<u) res = v;
    return res;
}

l_cinterval cidiv(const l_cinterval& A, const l_cinterval& B)
{
   l_real   realteilINF, realteilSUP,
            imagteilINF, imagteilSUP;
   // da sonst eventuell zwischendurch illegale Intervalle entstehen
   l_real     a0,b0;
   bool     a_repeat,b_repeat; 
   int      i, rep, j;
   l_real   AREINF, ARESUP, AIMINF, AIMSUP,
            BREINF, BRESUP, BIMINF, BIMSUP;
   l_interval ARE, AIM, BRE, BIM;

   // keine Fehlerabfrage -> ist schon in CINTVAL.CPP
   //  IF ( 0.0 IN B.RE ) AND ( 0.0 IN B.IM ) THEN
   //   CIDIVISION:= COMPL( 1.0 / INTVAL(-1.0,1.0), INTVAL(0.0) );
   //   Fehlerabbruch erzwingen: Intervall enthaelt 0

   //  *** Berechnung des Realteils ***

   AREINF = Inf(Re(A));    AREINF += 0.0;
   ARESUP = Sup(Re(A));    ARESUP += 0.0;
   AIMINF = Inf(Im(A));    AIMINF += 0.0;
   AIMSUP = Sup(Im(A));    AIMSUP += 0.0;
   BREINF = Inf(Re(B));    BREINF += 0.0;
   BRESUP = Sup(Re(B));    BRESUP += 0.0;
   BIMINF = Inf(Im(B));    BIMINF += 0.0;
   BIMSUP = Sup(Im(B));    BIMSUP += 0.0;
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
      for (i=1; i<=4; cxsc_l_complex_division_p[i++] = true);
    
      realteilSUP =
             max( realteilSUP,
                   max( max( minmax( false, a0, b0, BIMINF, BRE, 1,2 ),
                             minmax( false, a0, b0, BIMSUP, BRE, 3,4 ) ),
                        max( minmax( false, b0, a0, BREINF, BIM, 1,3 ),
                             minmax( false, b0, a0, BRESUP, BIM, 2,4 ) ) )

                );

      if (cxsc_l_complex_division_p[1])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BREINF, BIMINF, +1 ) );
      if (cxsc_l_complex_division_p[2])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BRESUP, BIMINF, +1 ) );
      if (cxsc_l_complex_division_p[3])
         realteilSUP = max( realteilSUP, cxsc_complex_division_f( a0, b0, BREINF, BIMSUP, +1 ) );
      if (cxsc_l_complex_division_p[4])
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
      for (i=1; i<=4; cxsc_l_complex_division_p[i++] = true);
      
      realteilINF =
              min( realteilINF,
                   min( min( minmax( true, a0, b0, BIMINF, BRE, 1,2 ),
                             minmax( true, a0, b0, BIMSUP, BRE, 3,4 ) ),
                        min( minmax( true, b0, a0, BREINF, BIM, 1,3 ),
                             minmax( true, b0, a0, BRESUP, BIM, 2,4 ) ) )
                 );

      if (cxsc_l_complex_division_p[1])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BREINF, BIMINF, -1 ) );
      if (cxsc_l_complex_division_p[2])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BRESUP, BIMINF, -1 ) );
      if (cxsc_l_complex_division_p[3])
         realteilINF = min( realteilINF, cxsc_complex_division_f( a0, b0, BREINF, BIMSUP, -1 ) );
      if (cxsc_l_complex_division_p[4])
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
      for (i=1; i<=4; cxsc_l_complex_division_p[i++] = true) ;
    
      imagteilSUP =
              max( imagteilSUP,
                   max( max( minmax( false,  b0, -a0, BIMINF, BRE, 1,2 ),
                             minmax( false,  b0, -a0, BIMSUP, BRE, 3,4 ) ),
                        max( minmax( false, -a0,  b0, BREINF, BIM, 1,3 ),
                             minmax( false, -a0,  b0, BRESUP, BIM, 2,4 ) ) )
                 );
    
      if (cxsc_l_complex_division_p[1])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BREINF, BIMINF, +1 ) );
      if (cxsc_l_complex_division_p[2])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BRESUP, BIMINF, +1 ) );
      if (cxsc_l_complex_division_p[3])
         imagteilSUP = max( imagteilSUP, cxsc_complex_division_f( b0, -a0, BREINF, BIMSUP, +1 ) );
      if (cxsc_l_complex_division_p[4])
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
      for (i=1; i<=4; cxsc_l_complex_division_p[i++] = true) ;
    
      imagteilINF =
              min( imagteilINF,
                   min( min( minmax( true,  b0, -a0, BIMINF, BRE, 1,2 ),
                             minmax( true,  b0, -a0, BIMSUP, BRE, 3,4 ) ),
                        min( minmax( true, -a0,  b0, BREINF, BIM, 1,3 ),
                             minmax( true, -a0,  b0, BRESUP, BIM, 2,4 ) ) )
                 );
    
      if (cxsc_l_complex_division_p[1])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BREINF, BIMINF, -1 ) );
      if (cxsc_l_complex_division_p[2])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BRESUP, BIMINF, -1 ) );
      if (cxsc_l_complex_division_p[3])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BREINF, BIMSUP, -1 ) );
      if (cxsc_l_complex_division_p[4])
         imagteilINF = min( imagteilINF, cxsc_complex_division_f( b0, -a0, BRESUP, BIMSUP, -1 ) );

      if (b_repeat) 
         b0 = AIMINF;  
      else if (a_repeat) 
         a0 = ARESUP;
   }

   return l_cinterval(l_interval(realteilINF, realteilSUP),
                      l_interval(imagteilINF, imagteilSUP));
}  // cidiv

l_cinterval C_point_div(const l_cinterval& z, const l_cinterval& n)
// Division of complex point intervals; 
// z,n must be point intervals!!  Blomquist, 07,11.02
// This function only for internal use!
{
    l_complex a,b,q1,q2;
    a = l_complex(InfRe(z),InfIm(z));
    b = l_complex(InfRe(n),InfIm(n));
    q1 = divd(a,b);
    q2 = divu(a,b);

    l_interval re, im;
    re = l_interval( Re(q1),Re(q2) );
    im = l_interval( Im(q1),Im(q2) );

    return l_cinterval(re,im);
}  // C_point_div


l_cinterval operator / (const l_cinterval & a, const l_cinterval & b) 
                                                      throw(DIV_BY_ZERO)
{
    if (0.0 <= b.re && 0.0 <= b.im ) {
//    if (0.0 <= (sqr(b.re) + sqr(b.im))) {
      cxscthrow(DIV_BY_ZERO("l_cinterval operator / (const l_cinterval&, const l_cinterval&)"));
      return a; // dummy result
    }      
    bool a_point, b_point;
    a_point = InfRe(a)==SupRe(a) && InfIm(a)==SupIm(a);
    b_point = InfRe(b)==SupRe(b) && InfIm(b)==SupIm(b);
    if(a_point && b_point) return C_point_div(a,b); // a,b are point intervals
    else return cidiv(a,b);
}

l_interval abs(const l_cinterval &a) throw()
{
    return sqrtx2y2(a.re,a.im);
}


// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const l_cinterval& a) throw()
{
    s << '('          
      << a.re << ','  
      << a.im       
      << ')';
    return s;
}

std::string & operator << (std::string &s, const l_cinterval& a) throw()
{
// string s; l_cinterval a;
// s << a; s delivers the string of the value a in the form:
// ([Inf(real-part(a)),Sup(real-part(a))],[Inf(img-part(a)),Sup(img-part(a))])
    s+='(';
    s << a.re;
    s+=',';
    s << a.im; 
    s+=')';
    return s;
}

std::string & operator >> (std::string &s, l_cinterval &a) 
                                     throw(EMPTY_INTERVAL)
// With: 
//       l_cinterval a;
//       string("([1.234,1.234],[2.567,2.567])") >> a;
// the value a will be an inclusion of the above string.
// The actual precisions of the staggered intervals a.re and a.im 
// will not be affected by the operator >> !
// The above braces, brackets and commas must not be used in the 
// string, however the four numbers must then be seperated by spaces!
// Thus, the following string will produce the same inclusion a:
//       string("1.234 1.234 2.567 2.567 ") >> a;
// Blomquist, 15.11.2006;
{
    l_real Iar,Sar,Iai,Sai;
    l_interval lr,li;
    int stagprec_old(stagprec);
    dotprecision dot;

    s = skipwhitespacessinglechar (s, '(');
    s = skipwhitespacessinglechar (s, '[');
    s = s >> dot;
    stagprec = StagPrec(a.re);
    lr = l_interval(dot);
    Iar = Inf(lr);
    s = skipwhitespacessinglechar (s, ',');
    s = s >> dot; 
    lr = l_interval(dot);
    Sar = Sup(lr);
    lr = l_interval(Iar,Sar);

    stagprec = StagPrec(a.im);
    s = skipwhitespacessinglechar (s, ']');
    s = skipwhitespacessinglechar (s, ',');
    s = skipwhitespacessinglechar (s, '[');
    s = s >> dot; 

    li = l_interval(dot);
    Iai = Inf(li);
    s = skipwhitespacessinglechar (s, ',');
    s = s >> dot; 
    li = l_interval(dot);
    Sai = Sup(li);
    li = l_interval(Iai,Sai);

    a = l_cinterval(lr,li );
    s = skipwhitespaces (s);
    if (s[0] == ']') 
        s.erase(0,1);
    s = skipwhitespaces (s);
    if (s[0] == ')') 
        s.erase(0,1);
    stagprec = stagprec_old;
    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(EMPTY_INTERVAL
        ("std::string & operator >> (std::string &s, cinterval &a)"));

   return s;
}

std::istream & operator >> (std::istream & s, l_cinterval& a) 
                                                 throw(EMPTY_INTERVAL)
// With: 
//       l_cinterval lc;
//       cout << "([a,b],[c,d]) = ?" << endl;
//       cin >> lc;
// the input string ([1.23,1.23],[3.45,3.45]) will be included by lc.
// The actual precisions of the staggered intervals  lc.re  and  lc.im 
// will not be affected by the operator >> !
// The above braces, brackets and commas must not be used in the 
// string, however the four numbers a,b,c,d must then be seperated by 
// spaces! Thus, the following input string   1.23 1.23 3.45 3.45
// will produce the same inclusion lc:
// Blomquist, 15.11.2006;
{ 
    l_real Iar,Sar,Iai,Sai;
    l_interval lr,li;
    dotprecision dot;
    // int stagprec_old(stagprec); // unused variable

    char c;
    skipeolnflag = inpdotflag = true;
    stagprec = StagPrec(a.re);
    c = skipwhitespacessinglechar (s, '(');
    if (inpdotflag) s.putback(c);
    c = skipwhitespacessinglechar (s, '[');
    if (inpdotflag) s.putback(c);
    s >> dot;
    lr = l_interval(dot);
    Iar = Inf(lr);
    skipeolnflag = inpdotflag = true; 
    c = skipwhitespacessinglechar (s, ','); 
    if (inpdotflag) s.putback(c);
    s >> dot;
    lr = l_interval(dot);
    Sar = Sup(lr);
    lr = l_interval(Iar,Sar);
    c = skipwhitespacessinglechar (s, ']');
    if (inpdotflag) s.putback(c);
    c = skipwhitespacessinglechar (s, ',');
    if (inpdotflag) s.putback(c);
    
    c = skipwhitespacessinglechar (s, '[');
    if (inpdotflag) s.putback(c);
//    s >> RndDown >> Inf(a.im);
    stagprec = StagPrec(a.im);
    s >> dot;
    li = l_interval(dot);
    Iai = Inf(li);
    skipeolnflag = inpdotflag = true; 
    c = skipwhitespacessinglechar (s, ','); 
    if (inpdotflag) s.putback(c);
    s >> dot;
    li = l_interval(dot);
    Sai = Sup(li);
    li = l_interval(Iai,Sai);

    a = l_cinterval(lr,li);

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
   cxscthrow(EMPTY_INTERVAL
        ("std::istream & operator >> (std::istream &s, cinterval &a)"));
      
   return s;
}

void operator >> (const std::string &s, l_cinterval &a) throw(EMPTY_INTERVAL)
{
   std::string r(s);
   r >> a;
}

void operator >> (const char *s, l_cinterval &a) throw(EMPTY_INTERVAL)
{
   std::string r(s);
   r >> a;
}

} // namespace cxsc


























