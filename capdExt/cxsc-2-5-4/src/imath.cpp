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

/* CVS $Id: imath.cpp,v 1.43 2014/01/30 17:23:45 cxsc Exp $ */

#include "imath.hpp"
#include "rmath.hpp"

// Auch fi_lib hat ein eigenes interval (wie RTS sein dotprecision)
//typedef struct fi_interval { double INF, SUP;} fi_interval;
#undef LINT_ARGS
#define CXSC_INCLUDE
#include <fi_lib.hpp>

extern "C" {
#ifndef rfcth_included
#define rfcth_included
#include "r_fcth.h"
#endif
}

namespace cxsc {


interval sqr    (const interval &a) throw() 
{ 
  interval res;
  res= a*a;
  if (Inf(res)<0) Inf(res)=0;
  return res;
}

interval sqrt   (const interval &a, int n)  throw(STD_FKT_OUT_OF_DEF)
{ 
   if ( ((n>0) && (Inf(a)>=0.0)) || ((n<0) && (Inf(a)>0.0)) ) 
      return pow(a,interval(1.0,1.0)/n); 
   else { 
      cxscthrow(STD_FKT_OUT_OF_DEF("interval sqrt (const interval &a, int n)"));
      return interval(-1.0); // dummy result
   }
}

interval sqrt1px2(const interval& x) throw()
// Inclusion of sqrt(1+x^2); Blomquist, 13.12.02;
{
    interval t = abs(x),y;
    if (expo(Inf(t)) > 33)
    {
	y = t;
	Sup(y) = succ(Sup(y));
    } else if (expo(Sup(t)) > 33)
    {
	y = interval(Inf(t));  // interval(Inf(t)) is a point interval!
	y = sqrt(1+y*y);       // --->  y*y == sqr(y);
	y = interval(Inf(y),succ(Sup(t)));
     } else y = sqrt(1+sqr(t));
    return y;
}

interval sqrtx2y2(const interval& x, const interval& y) throw()
// Inclusion of sqrt(x^2+y^2); Blomquist, 13.12.02;
{
    interval a=abs(x), b=abs(y), r;
    int exa=expo(Sup(a)), exb=expo(Sup(b)), ex;
    if (exb > exa)
    {  // Permutation of a,b:
	r = a;  a = b;  b = r;
	ex = exa;  exa = exb;  exb = ex;
    }
    ex = 511 - exa;
    times2pown(a,ex);
    times2pown(b,ex);
    r = sqrt(a*a + b*b);
    times2pown(r,-ex);
    return r;
} // sqrtx2y2

//********************************************************************
//  Constants for:   interval sqrtp1m1(const interval& x) throw()
//  Blomquist 05.08.03
const real Delta_f = 2*minreal;
const real q_sqrtp1m1m = 9007199254740984.0 / 9007199254740992.0;
const real q_sqrtp1m1p = 4503599627370502.0 / 4503599627370496.0;
//********************************************************************
interval sqrtp1m1(const interval& x) throw()
// interval(a,b) is an inclusion of sqrt(x+1)-1;
// Exported by imath.hpp;       Blomquist, 05.08.03;
{
    real a=0,b=0,ix=Inf(x),sx=Sup(x);
    int ex_ix,ex_sx,sgn_ix,sgn_sx;

    // Calculating the lower bound:
    ex_ix = expo(ix);  sgn_ix = sign(ix);
    if (ex_ix<=-1021) //  ==> |ix| < 2^(-1021)
        { if (sgn_ix) a = sqrtp1m1(ix) - Delta_f; }
    else // |ix| >= 2^(-1021)
	if (ex_ix<=-51) 
	{
	    times2pown(ix,-1); // exact division by 2
	    a = pred(ix);
	}
	else 
	    if (sgn_ix>0) a = (ix>0.67) ? 
	       Inf( sqrt(interval(ix)+1)-1 ) : sqrtp1m1(ix)*q_sqrtp1m1m;
	    else a = (ix<-0.25) ? 
	       Inf( sqrt(interval(ix)+1)-1 ) : sqrtp1m1(ix)*q_sqrtp1m1p;
    // Calculating the upper bound:
    if (ix == sx) { ex_sx = ex_ix;  sgn_sx = sgn_ix; }
    else { ex_sx = expo(sx);  sgn_sx = sign(sx); }
    if (ex_sx<=-1021) //  ==> |sx| < 2^(-1021)
        { if (sgn_sx) b = sqrtp1m1(sx) + Delta_f; }
    else // |sx| >= 2^(-1021)
	if (ex_sx<=-47) { b = sx; times2pown(b,-1); }
	else 
	    if (sgn_sx>0) b = (sx>0.58) ? 
	        Sup( sqrt(interval(sx)+1)-1 ) : sqrtp1m1(sx)*q_sqrtp1m1p;
	    else b = (sx<-0.32) ? 
	         Sup( sqrt(interval(sx)+1)-1 ) : sqrtp1m1(sx)*q_sqrtp1m1m;
    return interval(a,b);
} // sqrtp1m1()

//********************************************************************
// Konstanten fr die Intervallfunktion sqrt(x^2-1):
real q_sqrtx2m1p(4503599627370501.0/4503599627370496.0);  // (1+e(f))
real q_sqrtx2m1m(9007199254740986.0/9007199254740992.0);  // (1-e(f))
//********************************************************************

interval sqrtx2m1(const interval& x)
// sqrt(x^2-1);  Blomquist, 13.04.04;
{
    interval z = abs(x);
    real r1,r2;
    r1 = sqrtx2m1(Inf(z)) * q_sqrtx2m1m;
    r2 = Sup(z); 
    if (expo(r2)<26)  
        r2 = sqrtx2m1(r2) * q_sqrtx2m1p; 
    // expo(r2) >= 26 --> r2 = Sup(z) ist optimale Oberschranke!
    return interval(r1,r2);
} // sqrtx2m1 (Intervallargumente)


//********************************************************************
// Konstanten fr die Intervallfunktion sqrt(1-x2):
real q_sqrt1mx2p(4503599627370501.0/4503599627370496.0);  // (1+e(f))
real q_sqrt1mx2m(9007199254740985.0/9007199254740992.0);  // (1-e(f))
//********************************************************************

interval sqrt1mx2(const interval& x)
// sqrt(1-x2);  Blomquist, 13.04.04;
{
    interval z = abs(x);  // sqrt(1-t2) monoton fallend in t aus z;
    real r1,r2,sz,iz;
    sz = Sup(z);  iz = Inf(z);
    // Berechnung der Unterschranke r1:
    r2 = sqrt1mx2(sz);
    r1 = (sz==0)? 1 : r2 * q_sqrt1mx2m;
    // Berechnung der Oberrschranke r2:
    if (iz<4.81e-8) r2 = 1; // r2=1 ist immer korrekte Oberschranke!
    else r2 = (sz==iz)? r2*q_sqrt1mx2p : sqrt1mx2(iz)*q_sqrt1mx2p;
    return interval(r1,r2);

} // sqrtx2m1 (Intervallargumente)


// Konstanten fuer die Intervallfunktion exp(-x^2):
const real q_exp_x2p = 4503599627370502.0 / 4503599627370496.0;  // (1+e(f))
const real q_exp_x2m = 9007199254740984.0 / 9007199254740992.0;  // (1-e(f))
  // Oberschranke expmx2_UB fuer |x| > expmx2_x0:  
const real expmx2_UB = 2.225073858507447856659E-308;
const real expmx2_x0 = 7491658466053896.0 / 281474976710656.0;

interval expmx2(const interval& x)
// e^(-x^2);  Blomquist, 05.07.04;
{
    real y,r1,r2,Sz,Iz;
    interval z = abs(x); 
    // Berechnung einer Unterschranke:
    Sz = Sup(z);  Iz = Inf(z);
    y = expmx2(Sz);
    r1 = y;
    if (Sz<=expmx2_x0) r1 = r1 * q_exp_x2m; // Abrunden
    if (Sz==0) r1 = y;
    // Berechnung einer Oberschranke:
    if (Iz>expmx2_x0) r2 = expmx2_UB;
    else r2 = (Sz==Iz)? y * q_exp_x2p : expmx2(Iz) * q_exp_x2p; 
    if (r2>1) r2 = 1.0;
    return interval(r1,r2);

} // expmx2 (Intervallargumente)


interval expm1(const interval& x)
// Interval function of exp(x)-1
{
    real y,r1,r2,Sx,Ix;
    Sx = Sup(x);  Ix = Inf(x);
    y = expm1(Ix);
    // calculation of a lower bound r1:
    r1 = (y>0)? y*q_exmm : y*q_exmp;  // rounding downwards
    if (r1<-1) r1 = -1;
    // calculation of an upper bound r2:
    if (Sx!=Ix) y = expm1(Sx);
    r2 = (y>0)? y*q_exmp : y*q_exmm;

    return interval(r1,r2);
} // expm1(...) 

// Die folgenden Konstanten q_expx2_p, q_expx2_m beziehen sich 
// auf eps(expx2) = 4.618958e-16;    f(x) = e^{+x^2};

const real q_expx2_p = 4503599627370500.0 / 4503599627370496.0;  // (1+e(f))
const real q_expx2_m = 9007199254740985.0 / 9007199254740992.0;  // (1-e(f))

interval expx2(const interval& x)
// e^(+x^2);  Blomquist, 25.07.06;
{
    real y,r1,r2,Sz,Iz;
    interval z = abs(x);
    Sz = Sup(z);  Iz = Inf(z);
    // Berechnung einer Unterschranke:
    y = expx2(Iz);
    r1 = y;
    r1 = r1 * q_expx2_m; // Abrunden
    if (r1<1.0) r1 = 1.0;
    // Berechnung einer Oberschranke:
    r2 = (Sz==Iz)? y * q_expx2_p : expx2(Sz) * q_expx2_p;
    if (Sz==0) r2 = 1.0;

    return interval(r1,r2);
} // expx2 (Intervallfunktion)

// Die folgenden Konstanten q_expx2m1_p, q_expx2m1_m beziehen sich 
// auf eps(expx2m1) = 4.813220e-16;    f(x) = e^{+x^2}-1;

const real q_expx2m1_p = 4503599627370500.0 / 4503599627370496.0; // (1+e(f))
const real q_expx2m1_m = 9007199254740985.0 / 9007199254740992.0; // (1-e(f))
const real expx2m1_0   = comp(0.5,-510);  // 2^(-511)

void sqr2uv(const real&, real&, real&);

real expx2m1_intv(const real& x)
// Zur Implementierung der Intervallfunktion;
// e^(+x^2)-1; rel. Fehlerschranke: eps = 4.813220E-16 = e(f) gilt
// fuer alle x, mit:  2^(-511) <= x <= x0 = 26.64174755704632....
// x0 = 7498985273150791.0 / 281474976710656.0;
// Fuer x > x0 --> Programmabbruch wegen Overflow;
// Fuer 0 <= x < 2^(-511) werden die Funktionswerte auf Null gesetzt.
// Ausfuehrlich getestet;  Blomquist, 10.08.2006;
{
    real t(x),u,v,y,res(0);
    int ex;
    if (t<0) t = -t;  // t >= 0;
    
    if (t>=6.5) res = expx2(t);
    else 
    {
	ex = expo(t);
	sqr2uv(x,u,v);  // u := x*x und v aus S(2,53); 
	if (ex>=2) // g4(x)
	{
	    y = exp(u); 
	    res = 1 - v*y;
	    res = y - res;
	}
	else 
	    if (ex>=-8) res = expm1(u) + v*exp(u); // g3(x)
	    else 
		if (ex>=-25) { // g2(x)
		    y = u*u;
                    times2pown(y,-1);
		    res = (1+u/3)*y + u;
		}
		else 
		    if(ex>=-510) res = u;  // g1(x)
    }

    return res;
} // expx2m1_intv

interval expx2m1(const interval& x)
// e^(+x^2)-1;  Blomquist, 10.08.06;
{
    real y,r1,r2,Sz,Iz;
    interval z = abs(x);
    Sz = Sup(z);  Iz = Inf(z);
    // Berechnung einer Unterschranke:
    y = expx2m1_intv(Iz);
    r1 = y;
    r1 = r1 * q_expx2m1_m; // Abrunden
    // Berechnung einer Oberschranke:
    if (Sz < expx2m1_0) 
    { 
	r2 = MinReal;
	if (Sz==0) r2 = 0.0;
    }
    else r2 = (Sz==Iz)? y * q_expx2m1_p : expx2m1_intv(Sz) * q_expx2m1_p;

    return interval(r1,r2);
} // expx2m1 (Intervallfunktion)

// ------  1-eps and 1+eps for function lnp1, Blomquist 28,07,03;  -----------
// ------------------------  eps = 2.5082e-16  -------------------------------
static real q_lnp1m = 9007199254740986.0 / 9007199254740992.0; // 1-eps
static real q_lnp1p = 4503599627370501.0 / 4503599627370496.0; // 1+eps
// ---------------------------------------------------------------------------

interval lnp1(const interval& x) throw()
// returns an inclusion of ln(1+t), with t in x; Blomquist 28.07.03;
{ 
    real ix=Inf(x), sx=Sup(x),a,b;   // ln(1+x) <= [a,b]
    // Calculating the lower bound a:
    int sgn_ix = sign(ix), ex_ix = expo(ix), sgn_sx,ex_sx;
    if (!sgn_ix) a = 0;  // optimal lower bound!
    else if (sgn_ix>0)
	     a = (ex_ix<=-53) ? pred(ix) : lnp1(ix) * q_lnp1m;
         else 
	     a = (ex_ix<=-54) ? pred(ix) : lnp1(ix) * q_lnp1p;
    // Calculating the upper bound b:
    if (ix == sx) { sgn_sx = sgn_ix; ex_sx = ex_ix; }
    else { sgn_sx = sign(sx); ex_sx = expo(sx); }
    if (sgn_sx>=0)
	b = (ex_sx<=-49) ? sx : lnp1(sx) * q_lnp1p;
    else // sx < 0:
	b = (ex_sx<=-50) ? sx : lnp1(sx) * q_lnp1m;
    return interval(a,b); // ln(1+x) in [a,b]
} // lnp1

/*!
\param a The value for which to compute the value of the error function
\return The computed result of the error function

\sa erf(const real & arg)
*/
interval erf    (const interval &a)         { return j_erf(a);  }
/*!
\param a The value for which to compute the value of the complementary error function
\return The computed result of the complementary error function

\sa erf(const real & arg)
*/
interval erfc   (const interval &a)         { return j_erfc(a); }

//interval pow    (const interval &a, const interval &b) throw(ERROR_INTERVAL_STD_FKT_OUT_OF_DEF)
//{
//	if(Inf(a)>0)
//		return j_exp(b*ln(a));
//	else if(Inf(a)==0 && Inf(b)>=0)
//	{
//		if(Sup(a)>=1)
//			return interval(0,pow(Sup(a),Sup(b)));
//		else
//			return interval(0,pow(Sup(a),Inf(b)));
//	}
//	else
//	{
//		cxscthrow(ERROR_INTERVAL_STD_FKT_OUT_OF_DEF("interval pow(const interval &,const interval &)"));
//		return interval(0,0);
//	}
//}

inline a_intv _a_intv(const interval &x)
{       
       return *((const a_intv *)(&x));
}
inline interval _interval(const a_intv &x)
{       
       return *((const interval *)(&x));
}

interval pow    (const interval &a, const interval &b) throw()
       { 
         interval res; 
         if(Inf(a)==0 && Inf(b)>=0)
         {
           if(Sup(a)>=1)
             res=interval(0,pow(Sup(a),Sup(b)));
           else
             res=interval(0,pow(Sup(a),Inf(b)));
         }
	 else res = _interval(i_pow(_a_intv(a),_a_intv(b))); 

         if (Inf(res) <= Sup(res))
           return res;
         else
           return interval(Sup(res),Inf(res));
       }
       
//----------------------------------------------------------------------------
// Purpose: The local function 'Power()' is used to compute a lower or an
//    upper bound for the power function with real argument and integer
//    exponent, respectively.
// Parameters:
//    In: 'x'      : real argument
//        'n'      : integer exponent
//        'RndMode': rounding mode,
//                   (-1 = downwardly directed, +1 = upwardly directed)
// Description:
//    This function is used to speed up the interval power function defined
//    below. The exponentiation is reduced to multiplications using the
//    binary shift method. Depending on 'n', this function is up to 40 times
//    as fast as the standard power function for real argument and real
//    exponent. However, its accuracy is less than one ulp (unit in the last
//    place of the mantissa) since about log2(n) multiplications are executed
//    during computation. Since directed roundings are antisymmetric, one
//    gets
//
//       down(x^n) = -up((-x)^n)   and   up(x^n) = -down((-x)^n)
//
//    for x < 0 and odd n, where 'down' and 'up' denote the downwardly and
//    upwardly directed roundings, respectively.
//----------------------------------------------------------------------------
static real Power (const real & x, int n, int RndMode )
{                         // Signals change of the rounding mode
  int  ChangeRndMode;     // for x < 0 and odd n
  real p, z;

  ChangeRndMode = ( (x < 0.0) && (n % 2 == 1) );

  if (ChangeRndMode) {
    z = -x;
    RndMode = -RndMode;
  }
  else
    z = x;

  p = 1.0;
  switch (RndMode) {                             // Separate while-loops used
    case -1 : while (n > 0) {                    // to gain speed at runtime
                if (n % 2 == 1) p = muld(p,z);   //--------------------------
                n = n / 2;
                if (n > 0) z = muld(z,z);
              }
              break;
    case +1 : while (n > 0) {
                if (n % 2 == 1) p = mulu(p,z);
                n = n / 2;
                if (n > 0) z = mulu(z,z);
              }
              break;
  }

  if (ChangeRndMode)
    return -p;
  else
    return p;
}

interval power(const interval& a, int n)
// Calculating a^n; 
// Examples: [-1,4]^2 = [0,16];  [-1,4]^3 = [-16,+64];   
{
    bool neg(n<0);
    int N(n);          //,k(-1),r; 
    interval res,h;
    real Lower, Upper;
    if (neg) N = -N;
    if (N==0) res = 1;
    else if (N==1) res = a;
    else // N > 1:
	if (Inf(a)>=MinReal) res = exp(N*ln(a));
	else if (Sup(a)<=-MinReal) {
	    h = interval(-Sup(a),-Inf(a));
	    res = exp(N*ln(h));
	    if (N%2 == 1) res = -res;
	} 
	else 
	{
//	    h = a;
//	    while(N>0) 
//	    {
//		k++;
//		r = N % 2;
//		if (k==0) 
//		    if (r==1) res=a; else res=1;
//		else {
//		    h = sqr(h);
//		    if (r==1) res *= h;
//		}
//		N = N / 2;
//	    }
          if ( (0.0 < Inf(a)) || (N % 2 == 1) ) {
            Lower = Power(Inf(a),N,-1);
            Upper = Power(Sup(a),N,+1);
          }
          else if (0.0 > Sup(a)) {
            Lower = Power(Sup(a),N,-1);
            Upper = Power(Inf(a),N,+1);
          }
          else {
          Lower = 0.0;
          Upper = Power(AbsMax(a),N,+1);
          }
          res = interval(Lower,Upper);
	}
    if (neg) res = 1/res;
    return res;
} // power


//----------------------------------------------------------------------------
// Purpose: This version of the function 'Power()' is used to compute an
//    enclosure for the power function with interval argument and integer
//    exponent.
// Parameters:
//    In: 'x': interval argument
//        'n': integer exponent
// Description:
//    In general, this implementation does not deliver a result of maximum
//    accuracy, but it is about 30-40 times faster than the standard power
//    function for interval arguments and interval exponents. The resulting
//    interval has a width of approximately 2*log2(n) ulps. Since x^n is
//    considered as a monomial, we define x^0 := 1. For negative exponents
//    and 0 in 'x', the division at the end of the function will cause a
//    runtime error (division by zero).
//----------------------------------------------------------------------------
interval Power (const interval & x, int n )
{
  int  m;
  real Lower, Upper;

  if (n == 0) return(interval(1.0,1.0));

  if (n > 0)  m = n;  else  m = -n;

  if ( (0.0 < Inf(x)) || (m % 2 == 1) ) {
    Lower = Power(Inf(x),m,-1);
    Upper = Power(Sup(x),m,+1);
  }
  else if (0.0 > Sup(x)) {
    Lower = Power(Sup(x),m,-1);
    Upper = Power(Inf(x),m,+1);
  }
  else {
    Lower = 0.0;
    Upper = Power(AbsMax(x),m,+1);
  }

  if (n > 0)
    return(interval(Lower,Upper));
  else                                    // Causes a runtime error
    return(1.0/interval(Lower,Upper));    // if 0 in 'x'.
}

//----------------------------------------------------------------------------

interval Pi ( )                                    // Enclosure of constant pi
  { return(4.0*atan(_interval(1.0,1.0))); }        //-------------------------

// Error bounds for the interval function ln_sqrtx2y2:
real ln_x2y2_abs(2.225076E-308); // Absolute error bond
real q_lnx2y2p(4503599627370502.0 / 4503599627370496.0); // 1+e(f)
real q_lnx2y2m(9007199254740984.0 / 9007199254740992.0); // 1-e(f)
// With the following b0 
// real b0 = 6369051672525773.0 / 30191699398572330817932436647906151127335369763331523427009650401964993299137190816689013801421270140331747000246110759198164677039398341060491474011461568349195162615808.0;
real b0 = MakeHexReal(0,1022-510,0x0016A09E,0x667F3BCD);
// it holds:
// 1. b < bo  ==> g(b) := (0.5*b)*b < MinReal with rounding downwards
// 2. b >= b0 ==> g(b) := (0.5*b)*b >= MinReal with arbitrary rounding
//                                   modus by the two multiplications.

interval ln_sqrtx2y2(const interval& x, const interval& y) throw()
// ln( sqrt(x^2+y^2) ) == 0.5*ln(x^2+y^2);   Blomquist, 22.11.03;
{
    interval ax=abs(x), ay=abs(y);
    real Ix=Inf(ax), Sx=Sup(ax), Iy=Inf(ay), Sy=Sup(ay),f,u1,u2;
    // Calculating the lower bound u1:
    f = ln_sqrtx2y2(Ix,Iy);
    if ((Ix==1 && Iy<b0) || (Iy==1 && Ix<b0)) {
    // f in the denormalized range!
	u1 = f - ln_x2y2_abs;   // directed rounding not necessary!
	if (sign(u1)<0) u1 = 0;
    } else  u1 = (sign(f)<0) ? f*q_lnx2y2p : f*q_lnx2y2m;
    // Calculating the upper bound u2:
    if (Ix==Sx && Iy==Sy) // x and y are point-intervals
	if ((Sx==1 && Sy<b0) || (Sy==1 && Sx<b0)) {
        // f in the denormalized range!
	    u2 = (Sy==0 || Sx==0) ? f : f+ln_x2y2_abs; 
	} else  u2 = (sign(f)<0) ? f*q_lnx2y2m : f*q_lnx2y2p;
    else // x or y is no point-interval:
    {
	f = ln_sqrtx2y2(Sx,Sy);
	if ((Sx==1 && Sy<b0) || (Sy==1 && Sx<b0))
        // f in the denormalized range!
	    u2 = (sign(Sy)==0 || sign(Sx)==0) ? f : f+ln_x2y2_abs; 
	else  u2 = (sign(f)<0) ? f*q_lnx2y2m : f*q_lnx2y2p;
    }
    return interval(u1,u2);
} // ln_sqrtx2y2


// Constants for the interval function acoshp1(x) = acosh(1+x):
// (1+e(f)):
static const real q_acoshp1p(4503599627370503.0/4503599627370496.0);  
// (1-e(f))
static const real q_acoshp1m(9007199254740981.0/9007199254740992.0);  

interval acoshp1(const interval& x)
// acoshp1;  Blomquist, 28.03.2005;
{
    real r1,r2,sx,ix;
    sx = Sup(x);  ix = Inf(x);
    // Calculating of the lower bound r1:
    r2 = acoshp1(ix);
    r1 = r2 * q_acoshp1m;
    // Calculating of the upper bound r2:
    r2 = (sx==ix)? r2*q_acoshp1p : acoshp1(sx)*q_acoshp1p;
    return interval(r1,r2);

} // acoshp1 (interval arguments)

// *****************************************************************************
//                               sin(Pi*x)/Pi                                  *
// ********************************************************************************

// Die folgenden Konstanten q_sinpix_p, q_sinpix_m beziehen sich 
// auf eps(sinpix_pi) = 3.401895e-16;    f(x) = sin(pi*x)/pi;

const real q_sinpix_p = 4503599627370499.0 / 4503599627370496.0;  // (1+e(f))
const real q_sinpix_m = 9007199254740986.0 / 9007199254740992.0;  // (1-e(f))

real rounded_up(const real& x)
{
	real y;
	y = (x>=0)? x*q_sinpix_p : x*q_sinpix_m;
	return y;
}

real rounded_down(const real& x)
{
	real y;
	y = (x>=0)? x*q_sinpix_m : x*q_sinpix_p;
	return y;
}

interval sinpix_pi(const interval& x)
{
	const real Pir = Sup(Pir_interval); // 1/pi upwards rounded 
	interval y;
	int ma,mb;
	real y1,y2,a(Inf(x)),b(Sup(x)),ya,yb;
	bool bl, ya_klg_yb;
		
	ma = Round(a);  mb = Round(b);
	bl = (ma%2 != 0);
	if (mb-ma>1)
	{
		y1 = -Pir;  y2 = Pir;
	}
	else // 0 <= mb-ma <=1
		if (mb==ma)
			if (a==b)  // x: Point interval
			{
				ya = sinpix_pi(a);
            y1 = rounded_down(ya);
				y2 = rounded_up(ya);
			}
			else // (mb==ma)  and  (a!=b)
			{
				ya = sinpix_pi(a);
				yb = sinpix_pi(b);
				if (!bl)
				{
					y1 = rounded_down(ya);
					y2 = rounded_up(yb);
				}
				else
				{
					y1 = rounded_down(yb);
					y2 = rounded_up(ya);
				}
			}
		else // mb-ma=1;
		{
			ya = sinpix_pi(a);
			yb = sinpix_pi(b);
			ya_klg_yb = (ya <= yb);
			
			if (bl)
			{
				if (!ya_klg_yb)
					yb = ya;
				ya = -Pir;
			}
			else
			{
				if (!ya_klg_yb)
					ya = yb;
				yb = Pir;
			}
			y1 = rounded_down(ya);
			if (y1<-Pir) 
				y1 = -Pir;
			y2 = rounded_up(yb);
			if (y2>Pir)
				y2 = Pir;
		}
		
	y = interval(y1,y2);
	return y;
}

// ********************************************************************************
// ***************  Konstanten bez. Gamma(x) und 1/Gamma(x) ********************
// ********************************************************************************

/* worst case relative error bound for 1/Gamma(x)            */    
/* eps(gammar)   = 2.866906e-15;                             */
/*     q_gammarm  = 1 - eps(gammar)                          */
const real   q_gammarm  = 9007199254740964.0 / 9007199254740992.0;
/*     q_gammarp  = 1 + eps(gammar)                       */
const	real   q_gammarp  = 4503599627370510.0 / 4503599627370496.0;

// Inclusion of the extremes of 1/gamma(x):
interval pow2(const interval& x, int ex)
{
	interval y(x);
   times2pown(y,ex);
	return y;
}
									 
const real Ne =  9007199254740992.0;
const real Ne1 = 1125899906842624.0;
const real Ne2 = 562949953421312.0;
const real Ne3 = 281474976710656.0;  
const real Ne4 = 140737488355328.0;  
const real Ne5 = 70368744177664.0;   
const real Ne6 = 35184372088832.0;
		    
const interval gam_rxi[171] =
{
	interval( 6582605572834349.0 / 4503599627370496.0,6582606400588712.0 / 
		       4503599627370496.0 ),
	interval( -4540376432147063.0 / 9007199254740992.0,-4540375772996112.0 /
			    9007199254740992.0 ),
	interval( -7086407292338520.0 / 4503599627370496.0,-7086406981597106.0 /
			    4503599627370496.0 ),
	interval( -5878820838740338.0 / 2251799813685248.0,-5878820690102701.0 /
			    2251799813685248.0 ),
	interval( -8185952996852629.0 / 2251799813685248.0,-8185952850644519.0 /
			    2251799813685248.0 ),
	interval( -5239079997162568.0 / Ne1,-5239079928185648.0 /
			    Ne1 ),
	interval( -6380657697812205.0 / Ne1,-6380657632438250.0 /
			    Ne1 ),
	interval( -7519230477777525.0 / Ne1,-7519230410301402.0 /
			    Ne1 ),
   interval( -8655680190901081.0 / Ne1,-8655680125714323.0 /
			    Ne1 ),
	interval( -4895280046470312.0 / Ne2,-4895280015191181.0 /
			    Ne2 ),
	interval( -5462119069950045.0 / Ne2,-5462119039144308.0 /
			    Ne2 ),
	interval( -6028485171921533.0 / Ne2,-6028485140574985.0 /
			    Ne2 ),
	interval( -6594470676196825.0 / Ne2,-6594470646005838.0 /
			    Ne2 ),
	interval( -7160144161412306.0 / Ne2,-7160144131821972.0 /
			    Ne2 ),
	interval( -7725557826948019.0 / Ne2,-7725557796923813.0 /
			    Ne2 ),
	interval( -8290752238810453.0 / Ne2,-8290752208368199.0 /
			    Ne2 ),
	interval( -8855759486553113.0 / Ne2,-8855759457969009.0 /
			    Ne2 ),
	interval( -4710302676530551.0 / Ne3,-4710302661730969.0 /
			    Ne3 ),
	interval( -4992655414739558.0 / Ne3,-4992655400196254.0 /
			    Ne3 ),
	interval( -5274946608174960.0 / Ne3,-5274946594248303.0 /
			    Ne3 ),
	interval( -5557183461054268.0 / Ne3,-5557183446978818.0 /
			    Ne3 ),
	interval( -5839372030862353.0 / Ne3,-5839372016359958.0 /
			    Ne3 ),
	interval( -6121517453741464.0 / Ne3,-6121517439898965.0 /
			    Ne3 ),
	interval( -6403624121061720.0 / Ne3,-6403624107245033.0 /
			    Ne3 ),
	interval( -6685695811452843.0 / Ne3,-6685695797850064.0 /
			    Ne3 ),
	interval( -6967735798799276.0 / Ne3,-6967735784842353.0 /
			    Ne3 ),
	interval( -7249746935136834.0 / Ne3,-7249746921647546.0 /
			    Ne3 ),
	interval( -7531731721183113.0 / Ne3,-7531731707547301.0 /
			    Ne3 ),
	interval( -7813692357395941.0 / Ne3,-7813692343778691.0 /
			    Ne3 ),
	interval( -8095630791428232.0 / Ne3,-8095630778040390.0 /
			    Ne3 ),
	interval( -8377548753853165.0 / Ne3,-8377548740538224.0 /
			    Ne3 ),
	interval( -8659447788741678.0 / Ne3,-8659447775645579.0 /
			    Ne3 ),
	interval( -8941329278863280.0 / Ne3,-8941329265761967.0 /
			    Ne3 ),
	interval( -4611597233397515.0 / Ne4,-4611597226897947.0 /
			    Ne4 ),
	interval( -4752522236545681.0 / Ne4,-4752522230097182.0 /
			    Ne4 ),
	interval( -4893440155674552.0 / Ne4,-4893440149271021.0 /
			    Ne4 ),
	interval( -5034351450592848.0 / Ne4,-5034351444042898.0 /
			    Ne4 ),
	interval( -5175256539394880.0 / Ne4,-5175256533033078.0 /
			    Ne4 ),
	interval( -5316155803584886.0 / Ne4,-5316155797241740.0 /
			    Ne4 ),
	interval( -5457049591969558.0 / Ne4,-5457049585698186.0 /
			    Ne4 ),
	interval( -5597938224227692.0 / Ne4,-5597938218047102.0 /
			    Ne4 ),
	interval( -5738821993949958.0 / Ne4,-5738821987610287.0 /
			    Ne4 ),
	interval( -5879701171316002.0 / Ne4,-5879701164978161.0 /
			    Ne4 ),
	interval( -6020576005634620.0 / Ne4,-6020575999374578.0 /
			    Ne4 ),
	interval( -6161446727231484.0 / Ne4,-6161446721061003.0 /
			    Ne4 ),
	interval( -6302313549465809.0 / Ne4,-6302313543147939.0 /
			    Ne4 ),
	interval( -6443176669927728.0 / Ne4,-6443176663720294.0 /
			    Ne4 ),
	interval( -6584036272403005.0 / Ne4,-6584036266213061.0 /
			    Ne4 ),
	interval( -6724892527807298.0 / Ne4,-6724892521515390.0 /
			    Ne4 ),
	interval( -6865745595463323.0 / Ne4,-6865745589218114.0 /
			    Ne4 ),
	interval( -7006595624078029.0 / Ne4,-7006595617954995.0 /
			    Ne4 ),
	interval( -7147442752315627.0 / Ne4,-7147442746313355.0 /
			    Ne4 ),
	interval( -7288287110549528.0 / Ne4,-7288287104399241.0 /
			    Ne4 ),
	interval( -7429128820262002.0 / Ne4,-7429128814290427.0 /
			    Ne4 ),
	interval( -7569967996009183.0 / Ne4,-7569967989919521.0 /
			    Ne4 ),
	interval( -7710804744971319.0 / Ne4,-7710804738911590.0 /
			    Ne4 ),
	interval( -7851639168099204.0 / Ne4,-7851639162048862.0 /
			    Ne4 ),
	interval( -7992471360380410.0 / Ne4,-7992471354478088.0 /
			    Ne4 ),
	interval( -8133301411685542.0 / Ne4,-8133301405639704.0 /
			    Ne4 ),
	interval( -8274129406251900.0 / Ne4,-8274129400330745.0 /
			    Ne4 ),
	interval( -8414955423947592.0 / Ne4,-8414955418025784.0 /
			    Ne4 ),
	interval( -8555779540237542.0 / Ne4,-8555779534343191.0 /
			    Ne4 ),
	interval( -8696601826519818.0 / Ne4,-8696601820560983.0 /
			    Ne4 ),
	interval( -8837422350374443.0 / Ne4,-8837422344547779.0 /
			    Ne4 ),
	interval( -8978241175812537.0 / Ne4,-8978241170031881.0 /
			    Ne4 ),
	interval( -4559529181955286.0 / Ne5,-4559529178973419.0 /
			    Ne5 ),
	interval( -4629936985962592.0 / Ne5,-4629936983062621.0 /
			    Ne5 ),
	interval( -4700344027557972.0 / Ne5,-4700344024658010.0 /
			    Ne5 ),
	interval( -4770750332748164.0 / Ne5,-4770750329806175.0 /
			    Ne5 ),
	interval( -4841155926312136.0 / Ne5,-4841155923514622.0 /
			    Ne5 ),
	interval( -4911560831959402.0 / Ne5,-4911560829205074.0 /
			    Ne5 ),
	interval( -4981965072279533.0 / Ne5,-4981965069314643.0 /
			    Ne5 ),
	interval( -5052368668591817.0 / Ne5,-5052368665644419.0 /
			    Ne5 ),
	interval( -5122771641448337.0 / Ne5,-5122771638497776.0 /
			    Ne5 ),
	interval( -5193174010445884.0 / Ne5,-5193174007591888.0 /
			    Ne5 ),
	interval( -5263575794318673.0 / Ne5,-5263575791488264.0 /
			    Ne5 ),
	interval( -5333977010931459.0 / Ne5,-5333977008101358.0 /
			    Ne5 ),
	interval( -5404377677505176.0 / Ne5,-5404377674566742.0 /
			    Ne5 ),
	interval( -5474777810213221.0 / Ne5,-5474777807391522.0 /
			    Ne5 ),
	interval( -5545177424917448.0 / Ne5,-5545177422065289.0 /
			    Ne5 ),
	interval( -5615576536591204.0 / Ne5,-5615576533728944.0 /
			    Ne5 ),
	interval( -5685975159660620.0 / Ne5,-5685975156829327.0 /
			    Ne5 ),
	interval( -5756373307993982.0 / Ne5,-5756373305121310.0 /
			    Ne5 ),
	interval( -5826770994840416.0 / Ne5,-5826770992034928.0 /
			    Ne5 ),
	interval( -5897168232948637.0 / Ne5,-5897168230130823.0 /
			    Ne5 ),
	interval( -5967565034702908.0 / Ne5,-5967565031698652.0 /
			    Ne5 ),
	interval( -6037961411567024.0 / Ne5,-6037961408739972.0 /
			    Ne5 ),
	interval( -6108357375160195.0 / Ne5,-6108357372312820.0 /
			    Ne5 ),
	interval( -6178752936267084.0 / Ne5,-6178752933424618.0 /
			    Ne5 ),
	interval( -6249148105390399.0 / Ne5,-6249148102629437.0 /
			    Ne5 ),
	interval( -6319542892697349.0 / Ne5,-6319542889886519.0 /
			    Ne5 ),
	interval( -6389937307844225.0 / Ne5,-6389937305053873.0 /
			    Ne5 ),
	interval( -6460331360217357.0 / Ne5,-6460331357417915.0 /
			    Ne5 ),
	interval( -6530725058889375.0 / Ne5,-6530725056107292.0 /
			    Ne5 ),
	interval( -6601118412624317.0 / Ne5,-6601118409842826.0 /
			    Ne5 ),
	interval( -6671511429748937.0 / Ne5,-6671511426971253.0 /
			    Ne5 ),
	interval( -6741904118401788.0 / Ne5,-6741904115722630.0 /
			    Ne5 ),
	interval( -6812296486576567.0 / Ne5,-6812296483762819.0 /
			    Ne5 ),
	interval( -6882688541678557.0 / Ne5,-6882688538949795.0 /
			    Ne5 ),
	interval( -6953080291125250.0 / Ne5,-6953080288362282.0 /
			    Ne5 ),
	interval( -7023471742013192.0 / Ne5,-7023471739228496.0 /
			    Ne5 ),
	interval( -7093862901097948.0 / Ne5,-7093862898353711.0 /
			    Ne5 ),
	interval( -7164253775147266.0 / Ne5,-7164253772348037.0 /
			    Ne5 ),
	interval( -7234644370421069.0 / Ne5,-7234644367656947.0 /
			    Ne5 ),
	interval( -7305034693228153.0 / Ne5,-7305034690493625.0 /
			    Ne5 ),
	interval( -7375424749560875.0 / Ne5,-7375424746816433.0 /
			    Ne5 ),
	interval( -7445814545246703.0 / Ne5,-7445814542529259.0 /
			    Ne5 ),
	interval( -7516204085881418.0 / Ne5,-7516204083182192.0 /
			    Ne5 ),
	interval( -7586593377039902.0 / Ne5,-7586593374291972.0 /
			    Ne5 ),
	interval( -7656982423878545.0 / Ne5,-7656982421139464.0 /
			    Ne5 ),
	interval( -7727371231629585.0 / Ne5,-7727371228950751.0 /
			    Ne5 ),
	interval( -7797759805243431.0 / Ne5,-7797759802601748.0 /
			    Ne5 ),
	interval( -7868148149631789.0 / Ne5,-7868148146919441.0 /
			    Ne5 ),
	interval( -7938536269389156.0 / Ne5,-7938536266702540.0 /
			    Ne5 ),
	interval( -8008924169145439.0 / Ne5,-8008924166482991.0 /
			    Ne5 ),
	interval( -8079311853288492.0 / Ne5,-8079311850638146.0 /
			    Ne5 ),
	interval( -8149699326208468.0 / Ne5,-8149699323496057.0 /
			    Ne5 ),
	interval( -8220086591946405.0 / Ne5,-8220086589180599.0 /
			    Ne5 ),
	interval( -8290473654625668.0 / Ne5,-8290473651911160.0 /
			    Ne5 ),
	interval( -8360860518207539.0 / Ne5,-8360860515490326.0 /
			    Ne5 ),
	interval( -8431247186492146.0 / Ne5,-8431247183838485.0 /
			    Ne5 ),
	interval( -8501633663256230.0 / Ne5,-8501633660616442.0 /
			    Ne5 ),
	interval( -8572019952114999.0 / Ne5,-8572019949515073.0 /
			    Ne5 ),
	interval( -8642406056614155.0 / Ne5,-8642406053936017.0 /
			    Ne5 ),
	interval( -8712791980171444.0 / Ne5,-8712791977516241.0 /
			    Ne5 ),
	interval( -8783177726114270.0 / Ne5,-8783177723474143.0 /
			    Ne5 ),
	interval( -8853563297747803.0 / Ne5,-8853563295082565.0 /
			    Ne5 ),
	interval( -8923948698211540.0 / Ne5,-8923948695603642.0 /
			    Ne5 ),
	interval( -8994333930651793.0 / Ne5,-8994333928013544.0 /
			    Ne5 ),
	interval( -4532359499005981.0 / Ne6,-4532359497669548.0 /
			    Ne6 ),
	interval( -4567551951590474.0 / Ne6,-4567551950297444.0 /
			    Ne6 ),
	interval( -4602744324609335.0 / Ne6,-4602744323253914.0 /
			    Ne6 ),
	interval( -4637936619312858.0 / Ne6,-4637936618035482.0 /
			    Ne6 ),
	interval( -4673128837201202.0 / Ne6,-4673128835890259.0 /
			    Ne6 ),
	interval( -4708320979539414.0 / Ne6,-4708320978220509.0 /
			    Ne6 ),
	interval( -4743513047586994.0 / Ne6,-4743513046289186.0 /
			    Ne6 ),
	interval( -4778705042649671.0 / Ne6,-4778705041357978.0 /
			    Ne6 ),
	interval( -4813896965939740.0 / Ne6,-4813896964638916.0 /
			    Ne6 ),
	interval( -4849088818685547.0 / Ne6,-4849088817383583.0 /
			    Ne6 ),
	interval( -4884280602023663.0 / Ne6,-4884280600723731.0 /
			    Ne6 ),
	interval( -4919472317096448.0 / Ne6,-4919472315799988.0 /
			    Ne6 ),
	interval( -4954663965077128.0 / Ne6,-4954663963749440.0 /
			    Ne6 ),
	interval( -4989855546965370.0 / Ne6,-4989855545671197.0 /
			    Ne6 ),
	interval( -5025047063935622.0 / Ne6,-5025047062637366.0 /
			    Ne6 ),
	interval( -5060238516963350.0 / Ne6,-5060238515663821.0 /
			    Ne6 ),
	interval( -5095429907060327.0 / Ne6,-5095429905761378.0 /
			    Ne6 ),
	interval( -5130621235251224.0 / Ne6,-5130621233930590.0 /
			    Ne6 ),
	interval( -5165812502482099.0 / Ne6,-5165812501185488.0 /
			    Ne6 ),
	interval( -5201003709724055.0 / Ne6,-5201003708430144.0 /
			    Ne6 ),
	interval( -5236194857910350.0 / Ne6,-5236194856593808.0 /
			    Ne6 ),
	interval( -5271385947936424.0 / Ne6,-5271385946609048.0 /
			    Ne6 ),
	interval( -5306576980681870.0 / Ne6,-5306576979376506.0 /
			    Ne6 ),
	interval( -5341767957039505.0 / Ne6,-5341767955734187.0 /
			    Ne6 ),
	interval( -5376958877810379.0 / Ne6,-5376958876528096.0 /
			    Ne6 ),
	interval( -5412149743940434.0 / Ne6,-5412149742614602.0 /
			    Ne6 ),
	interval( -5447340556059975.0 / Ne6,-5447340554811246.0 /
			    Ne6 ),
	interval( -5482531315163654.0 / Ne6,-5482531313857339.0 /
			    Ne6 ),
	interval( -5517722021905912.0 / Ne6,-5517722020635122.0 /
			    Ne6 ),
	interval( -5552912677100886.0 / Ne6,-5552912675823002.0 /
			    Ne6 ),
	interval( -5588103281480912.0 / Ne6,-5588103280171229.0 /
			    Ne6 ),
	interval( -5623293835745456.0 / Ne6,-5623293834493516.0 /
			    Ne6 ),
	interval( -5658484340702630.0 / Ne6,-5658484339399074.0 /
			    Ne6 ),
	interval( -5693674796958315.0 / Ne6,-5693674795683212.0 /
			    Ne6 ),
	interval( -5728865205222665.0 / Ne6,-5728865203973973.0 /
			    Ne6 ),
	interval( -5764055566229013.0 / Ne6,-5764055564966520.0 /
			    Ne6 ),
	interval( -5799245880597279.0 / Ne6,-5799245879305275.0 /
			    Ne6 ),
	interval( -5834436148937784.0 / Ne6,-5834436147684294.0 /
			    Ne6 ),
	interval( -5869626371953711.0 / Ne6,-5869626370690622.0 /
			    Ne6 ),
	interval( -5904816550239413.0 / Ne6,-5904816548965896.0 /
			    Ne6 ),
	interval( -5940006684383290.0 / Ne6,-5940006683093915.0 /
			    Ne6 ),
	interval( -5975196775000579.0 / Ne6,-5975196773734016.0 /
			    Ne6 ) };
// Inclusion of the extremes of 1/gamma(x): 
const interval gam_ryi[171] = { 
pow2( interval( 5085347089749720.0 / Ne,5085347089749823.0 / Ne ) , 1 ) ,
pow2( interval( -5082146609264467.0 / Ne,-5082146609264314.0 / Ne ) , -1 ) ,
pow2( interval( 7824158147621733.0 / Ne,7824158147621966.0 / Ne ) , -1 ) ,
pow2( interval( -5070842539852372.0 / Ne,-5070842539852221.0 / Ne ) , 1 ) ,
pow2( interval( 4593118780547419.0 / Ne,4593118780547576.0 / Ne ) , 3 ) ,
pow2( interval( -5333021955274733.0 / Ne,-5333021955274575.0 / Ne ) , 5 ) ,
pow2( interval( 7546574203185105.0 / Ne,7546574203185319.0 / Ne ) , 7 ) ,
pow2( interval( -6294628859031764.0 / Ne,-6294628859031469.0 / Ne ) , 10 ) ,
pow2( interval( 6045310252810166.0 / Ne,6045310252811273.0 / Ne ) , 13 ) ,
pow2( interval( -6568078652156336.0 / Ne,-6568078652156148.0 / Ne ) , 16 ) ,
pow2( interval( 7963169065060572.0 / Ne,7963169065060801.0 / Ne ) , 19 ) ,
pow2( interval( -5328217018030122.0 / Ne,-5328217018029960.0 / Ne ) , 23 ) ,
pow2( interval( 7800142897041864.0 / Ne,7800142897042089.0 / Ne ) , 26 ) ,
pow2( interval( -6199437664213474.0 / Ne,-6199437664213297.0 / Ne ) , 30 ) ,
pow2( interval( 5316470282961123.0 / Ne,5316470282961284.0 / Ne ) , 34 ) ,
pow2( interval( -4892929765135337.0 / Ne,-4892929765135165.0 / Ne ) , 38 ) ,
pow2( interval( 4810107119289947.0 / Ne,4810107119290088.0 / Ne ) , 42 ) ,
pow2( interval( -5030373421375086.0 / Ne,-5030373421374834.0 / Ne ) , 46 ) ,
pow2( interval( 5576144001185310.0 / Ne,5576144001185479.0 / Ne ) , 50 ) ,
pow2( interval( -6530685487420963.0 / Ne,-6530685487420774.0 / Ne ) , 54 ) ,
pow2( interval( 8057940169576582.0 / Ne,8057940169576818.0 / Ne ) , 58 ) ,
pow2( interval( -5223648494045513.0 / Ne,-5223648494045349.0 / Ne ) , 63 ) ,
pow2( interval( 7099855957135674.0 / Ne,7099855957135885.0 / Ne ) , 67 ) ,
pow2( interval( -5047359382236272.0 / Ne,-5047359382236084.0 / Ne ) , 72 ) ,
pow2( interval( 7492585872478835.0 / Ne,7492585872479188.0 / Ne ) , 76 ) ,
pow2( interval( -5795835662380422.0 / Ne,-5795835662380242.0 / Ne ) , 81 ) ,
pow2( interval( 4664800910382651.0 / Ne,4664800910382790.0 / Ne ) , 86 ) ,
pow2( interval( -7801058080117709.0 / Ne,-7801058080117472.0 / Ne ) , 90 ) ,
pow2( interval( 6767162072327001.0 / Ne,6767162072327282.0 / Ne ) , 95 ) ,
pow2( interval( -6082121514218736.0 / Ne,-6082121514218554.0 / Ne ) , 100 ) ,
pow2( interval( 5656800000052189.0 / Ne,5656800000052359.0 / Ne ) , 105 ) ,
pow2( interval( -5438268378952110.0 / Ne,-5438268378951951.0 / Ne ) , 110 ) ,
pow2( interval( 5398375606367166.0 / Ne,5398375606367329.0 / Ne ) , 115 ) ,
pow2( interval( -5527713447587841.0 / Ne,-5527713447587674.0 / Ne ) , 120 ) ,
pow2( interval( 5833125895912623.0 / Ne,5833125895912799.0 / Ne ) , 125 ) ,
pow2( interval( -6337936184674347.0 / Ne,-6337936184674153.0 / Ne ) , 130 ) ,
pow2( interval( 7084743510515278.0 / Ne,7084743510515501.0 / Ne ) , 135 ) ,
pow2( interval( -8141214882701327.0 / Ne,-8141214882701088.0 / Ne ) , 140 ) ,
pow2( interval( 4804968547193877.0 / Ne,4804968547194018.0 / Ne ) , 146 ) ,
pow2( interval( -5822137580509526.0 / Ne,-5822137580509355.0 / Ne ) , 151 ) ,
pow2( interval( 7236772755227956.0 / Ne,7236772755228162.0 / Ne ) , 156 ) ,
pow2( interval( -4610758665056508.0 / Ne,-4610758665056369.0 / Ne ) , 162 ) ,
pow2( interval( 6019530845699084.0 / Ne,6019530845699266.0 / Ne ) , 167 ) ,
pow2( interval( -8047036389398365.0 / Ne,-8047036389398123.0 / Ne ) , 172 ) ,
pow2( interval( 5504580189086749.0 / Ne,5504580189086968.0 / Ne ) , 178 ) ,
pow2( interval( -7703001513324420.0 / Ne,-7703001513324183.0 / Ne ) , 183 ) ,
pow2( interval( 5510183009440391.0 / Ne,5510183009440581.0 / Ne ) , 189 ) ,
pow2( interval( -8055535954952413.0 / Ne,-8055535954952173.0 / Ne ) , 194 ) ,
pow2( interval( 6014315232803007.0 / Ne,6014315232803294.0 / Ne ) , 200 ) ,
pow2( interval( -4584378555360492.0 / Ne,-4584378555360260.0 / Ne ) , 206 ) ,
pow2( interval( 7132212380084113.0 / Ne,7132212380084326.0 / Ne ) , 211 ) ,
pow2( interval( -5659549393054692.0 / Ne,-5659549393054526.0 / Ne ) , 217 ) ,
pow2( interval( 4579461117155838.0 / Ne,4579461117155977.0 / Ne ) , 223 ) ,
pow2( interval( -7554216840666713.0 / Ne,-7554216840666493.0 / Ne ) , 228 ) ,
pow2( interval( 6348787715758027.0 / Ne,6348787715758222.0 / Ne ) , 234 ) ,
pow2( interval( -5434979980476367.0 / Ne,-5434979980476204.0 / Ne ) , 240 ) ,
pow2( interval( 4737681191908824.0 / Ne,4737681191908967.0 / Ne ) , 246 ) ,
pow2( interval( -8407842664867513.0 / Ne,-8407842664867267.0 / Ne ) , 251 ) ,
pow2( interval( 7592052521188700.0 / Ne,7592052521188935.0 / Ne ) , 257 ) ,
pow2( interval( -6974119252551297.0 / Ne,-6974119252551090.0 / Ne ) , 263 ) ,
pow2( interval( 6515520808385677.0 / Ne,6515520808385874.0 / Ne ) , 269 ) ,
pow2( interval( -6188946869743481.0 / Ne,-6188946869743300.0 / Ne ) , 275 ) ,
pow2( interval( 5975502808844840.0 / Ne,5975502808845020.0 / Ne ) , 281 ) ,
pow2( interval( -5862842897072874.0 / Ne,-5862842897072704.0 / Ne ) , 287 ) ,
pow2( interval( 5843967448508660.0 / Ne,5843967448508828.0 / Ne ) , 293 ) ,
pow2( interval( -5916517001341501.0 / Ne,-5916517001341321.0 / Ne ) , 299 ) ,
pow2( interval( 6082464626325325.0 / Ne,6082464626325503.0 / Ne ) , 305 ) ,
pow2( interval( -6348157530347044.0 / Ne,-6348157530346858.0 / Ne ) , 311 ) ,
pow2( interval( 6724699799057619.0 / Ne,6724699799057843.0 / Ne ) , 317 ) ,
pow2( interval( -7228705737680202.0 / Ne,-7228705737679999.0 / Ne ) , 323 ) ,
pow2( interval( 7883493269720206.0 / Ne,7883493269720561.0 / Ne ) , 329 ) ,
pow2( interval( -8720834785364833.0 / Ne,-8720834785364561.0 / Ne ) , 335 ) ,
pow2( interval( 4891722644546351.0 / Ne,4891722644546502.0 / Ne ) , 342 ) ,
pow2( interval( -5564236710028970.0 / Ne,-5564236710028799.0 / Ne ) , 348 ) ,
pow2( interval( 6416191129172903.0 / Ne,6416191129173091.0 / Ne ) , 354 ) ,
pow2( interval( -7498890927628704.0 / Ne,-7498890927628487.0 / Ne ) , 360 ) ,
pow2( interval( 8881515552460572.0 / Ne,8881515552460999.0 / Ne ) , 366 ) ,
pow2( interval( -5328950915550370.0 / Ne,-5328950915550206.0 / Ne ) , 373 ) ,
pow2( interval( 6478093314396794.0 / Ne,6478093314397089.0 / Ne ) , 379 ) ,
pow2( interval( -7976303366065662.0 / Ne,-7976303366065426.0 / Ne ) , 385 ) ,
pow2( interval( 4972846688449830.0 / Ne,4972846688450017.0 / Ne ) , 392 ) ,
pow2( interval( -6278401907481090.0 / Ne,-6278401907480879.0 / Ne ) , 398 ) ,
pow2( interval( 8024854758356088.0 / Ne,8024854758356345.0 / Ne ) , 404 ) ,
pow2( interval( -5191277948909595.0 / Ne,-5191277948909444.0 / Ne ) , 411 ) ,
pow2( interval( 6797621462551740.0 / Ne,6797621462551941.0 / Ne ) , 417 ) ,
pow2( interval( -4503636668393666.0 / Ne,-4503636668393518.0 / Ne ) , 424 ) ,
pow2( interval( 6037997262493341.0 / Ne,6037997262493523.0 / Ne ) , 430 ) ,
pow2( interval( -8189485306115383.0 / Ne,-8189485306115130.0 / Ne ) , 436 ) ,
pow2( interval( 5617805845426844.0 / Ne,5617805845427124.0 / Ne ) , 443 ) ,
pow2( interval( -7795192616785187.0 / Ne,-7795192616784477.0 / Ne ) , 449 ) ,
pow2( interval( 5469175405734180.0 / Ne,5469175405734422.0 / Ne ) , 456 ) ,
pow2( interval( -7759929987383324.0 / Ne,-7759929987383086.0 / Ne ) , 462 ) ,
pow2( interval( 5565727978288701.0 / Ne,5565727978288876.0 / Ne ) , 469 ) ,
pow2( interval( -8070914994857895.0 / Ne,-8070914994857635.0 / Ne ) , 475 ) ,
pow2( interval( 5914931467943193.0 / Ne,5914931467943373.0 / Ne ) , 482 ) ,
pow2( interval( -8762204548045716.0 / Ne,-8762204548045455.0 / Ne ) , 488 ) ,
pow2( interval( 6558513517606168.0 / Ne,6558513517606353.0 / Ne ) , 495 ) ,
pow2( interval( -4960305627886271.0 / Ne,-4960305627886120.0 / Ne ) , 502 ) ,
pow2( interval( 7580642983583672.0 / Ne,7580642983583897.0 / Ne ) , 508 ) ,
pow2( interval( -5851844804194595.0 / Ne,-5851844804194367.0 / Ne ) , 515 ) ,
pow2( interval( 4563038858728436.0 / Ne,4563038858728577.0 / Ne ) , 522 ) ,
pow2( interval( -7187477492053316.0 / Ne,-7187477492052964.0 / Ne ) , 528 ) ,
pow2( interval( 5716852908386950.0 / Ne,5716852908387214.0 / Ne ) , 535 ) ,
pow2( interval( -4591808630269563.0 / Ne,-4591808630269411.0 / Ne ) , 542 ) ,
pow2( interval( 7448102539955649.0 / Ne,7448102539955986.0 / Ne ) , 548 ) ,
pow2( interval( -6098770429791387.0 / Ne,-6098770429791204.0 / Ne ) , 555 ) ,
pow2( interval( 5041550443966798.0 / Ne,5041550443966946.0 / Ne ) , 562 ) ,
pow2( interval( -8413996086583072.0 / Ne,-8413996086582821.0 / Ne ) , 568 ) ,
pow2( interval( 7086939987269423.0 / Ne,7086939987269731.0 / Ne ) , 575 ) ,
pow2( interval( -6024570065319942.0 / Ne,-6024570065319682.0 / Ne ) , 582 ) ,
pow2( interval( 5168535487082451.0 / Ne,5168535487082609.0 / Ne ) , 589 ) ,
pow2( interval( -8949051953781375.0 / Ne,-8949051953781115.0 / Ne ) , 595 ) ,
pow2( interval( 7817344426895164.0 / Ne,7817344426895996.0 / Ne ) , 602 ) ,
pow2( interval( -6889843867972878.0 / Ne,-6889843867972674.0 / Ne ) , 609 ) ,
pow2( interval( 6126229646423302.0 / Ne,6126229646423484.0 / Ne ) , 616 ) ,
pow2( interval( -5495122334906381.0 / Ne,-5495122334906222.0 / Ne ) , 623 ) ,
pow2( interval( 4971972094727164.0 / Ne,4971972094727314.0 / Ne ) , 630 ) ,
pow2( interval( -4537480959802395.0 / Ne,-4537480959802254.0 / Ne ) , 637 ) ,
pow2( interval( 8352835047353300.0 / Ne,8352835047353555.0 / Ne ) , 643 ) ,
pow2( interval( -7753443787904532.0 / Ne,-7753443787904298.0 / Ne ) , 650 ) ,
pow2( interval( 7257653550749169.0 / Ne,7257653550749382.0 / Ne ) , 657 ) ,
pow2( interval( -6850281165773769.0 / Ne,-6850281165773570.0 / Ne ) , 664 ) ,
pow2( interval( 6519305845448896.0 / Ne,6519305845449168.0 / Ne ) , 671 ) ,
pow2( interval( -6255266499085062.0 / Ne,-6255266499084872.0 / Ne ) , 678 ) ,
pow2( interval( 6050802311308162.0 / Ne,6050802311308350.0 / Ne ) , 685 ) ,
pow2( interval( -5900304762620398.0 / Ne,-5900304762620223.0 / Ne ) , 692 ) ,
pow2( interval( 5799657649647993.0 / Ne,5799657649648165.0 / Ne ) , 699 ) ,
pow2( interval( -5746047975553302.0 / Ne,-5746047975553134.0 / Ne ) , 706 ) ,
pow2( interval( 5737835419331524.0 / Ne,5737835419331693.0 / Ne ) , 713 ) ,
pow2( interval( -5774471890994117.0 / Ne,-5774471890993944.0 / Ne ) , 720 ) ,
pow2( interval( 5856465763387432.0 / Ne,5856465763387600.0 / Ne ) , 727 ) ,
pow2( interval( -5985387992102590.0 / Ne,-5985387992102406.0 / Ne ) , 734 ) ,
pow2( interval( 6163919695584074.0 / Ne,6163919695584257.0 / Ne ) , 741 ) ,
pow2( interval( -6395943042753787.0 / Ne,-6395943042753502.0 / Ne ) , 748 ) ,
pow2( interval( 6686679647283150.0 / Ne,6686679647283350.0 / Ne ) , 755 ) ,
pow2( interval( -7042883260256940.0 / Ne,-7042883260256730.0 / Ne ) , 762 ) ,
pow2( interval( 7473096566380533.0 / Ne,7473096566380749.0 / Ne ) , 769 ) ,
pow2( interval( -7987985534527481.0 / Ne,-7987985534527243.0 / Ne ) , 776 ) ,
pow2( interval( 8600769311605383.0 / Ne,8600769311605633.0 / Ne ) , 783 ) ,
pow2( interval( -4663884705694464.0 / Ne,-4663884705694325.0 / Ne ) , 791 ) ,
pow2( interval( 5094554684614484.0 / Ne,5094554684614634.0 / Ne ) , 798 ) ,
pow2( interval( -5604802840349871.0 / Ne,-5604802840349701.0 / Ne ) , 805 ) ,
pow2( interval( 6209951739735886.0 / Ne,6209951739736072.0 / Ne ) , 812 ) ,
pow2( interval( -6928963530888061.0 / Ne,-6928963530887851.0 / Ne ) , 819 ) ,
pow2( interval( 7785368708274196.0 / Ne,7785368708274423.0 / Ne ) , 826 ) ,
pow2( interval( -8808459126256481.0 / Ne,-8808459126256060.0 / Ne ) , 833 ) ,
pow2( interval( 5017412797579486.0 / Ne,5017412797579638.0 / Ne ) , 841 ) ,
pow2( interval( -5755173329981532.0 / Ne,-5755173329981361.0 / Ne ) , 848 ) ,
pow2( interval( 6646385258439176.0 / Ne,6646385258439444.0 / Ne ) , 855 ) ,
pow2( interval( -7727539896552529.0 / Ne,-7727539896552294.0 / Ne ) , 862 ) ,
pow2( interval( 4522473425691912.0 / Ne,4522473425692052.0 / Ne ) , 870 ) ,
pow2( interval( -5328812572761788.0 / Ne,-5328812572761623.0 / Ne ) , 877 ) ,
pow2( interval( 6320558000502691.0 / Ne,6320558000502885.0 / Ne ) , 884 ) ,
pow2( interval( -7546265781200776.0 / Ne,-7546265781200489.0 / Ne ) , 891 ) ,
pow2( interval( 4534316912522546.0 / Ne,4534316912522688.0 / Ne ) , 899 ) ,
pow2( interval( -5484491485989575.0 / Ne,-5484491485989407.0 / Ne ) , 906 ) ,
pow2( interval( 6676632315202014.0 / Ne,6676632315202302.0 / Ne ) , 913 ) ,
pow2( interval( -8180074398476253.0 / Ne,-8180074398476014.0 / Ne ) , 920 ) ,
pow2( interval( 5042989707083422.0 / Ne,5042989707083666.0 / Ne ) , 928 ) ,
pow2( interval( -6257379418480333.0 / Ne,-6257379418480019.0 / Ne ) , 935 ) ,
pow2( interval( 7813097673618694.0 / Ne,7813097673619043.0 / Ne ) , 942 ) ,
pow2( interval( -4908325621754370.0 / Ne,-4908325621754217.0 / Ne ) , 950 ) ,
pow2( interval( 6205346227418363.0 / Ne,6205346227418597.0 / Ne ) , 957 ) ,
pow2( interval( -7893590972392525.0 / Ne,-7893590972392227.0 / Ne ) , 964 ) ,
pow2( interval( 5051411882876310.0 / Ne,5051411882876506.0 / Ne ) , 972 ) ,
pow2( interval( -6504655602059905.0 / Ne,-6504655602059583.0 / Ne ) , 979 ) ,
pow2( interval( 8426810051054742.0 / Ne,8426810051054986.0 / Ne ) , 986 ) ,
pow2( interval( -5491407534973626.0 / Ne,-5491407534973452.0 / Ne ) , 994 ) ,
pow2( interval( 7199960218142557.0 / Ne,7199960218142768.0 / Ne ) , 1001 ) ,
pow2( interval( -4748178637044143.0 / Ne,-4748178637044000.0 / Ne ) , 1009 ) ,
pow2( interval( 6299691458188962.0 / Ne,6299691458189149.0 / Ne ) , 1016 ) };

// ****************************************************************************
// ******************** Gamma(x), 1/Gamma(x) **********************************
// ****************************************************************************

inline int round_g(const real& x) throw() 
// Only for the internal use ( interval gammar(x) )
// For |x| < 2147483647.5 the assignment  y = round_g(x)  delivers:
// y = round_g(-0.1);            --->  y = 1;
// y = round_g(-0.5);            --->  y = 1;
// y = round_g(-1.0);            --->  y = 1;
// y = round_g(-1.4);            --->  y = 2;		
// y = round_g(-6.8);            --->  y = 7;
// y = round_g(+0.0);            --->  y = 0;
// y = round_g(+0.1);            --->  y = 0;
// y = round_g(+0.5);            --->  y = 0;
// y = round_g(+2.6);            --->  y = 0;
// Blomquist, 25.06.2009;
{
	int n = ifloor(_double(x));
	n = (n>=0)? 0:-n;
	
	return n;
}

real gamr_even_Ma(const real& x1, const real& x2, int n1)
{
	real y;
	
	if ( x2<Inf(gam_rxi[n1]) || Sup(gam_rxi[n1])<x1 ) 
	{  // [x1,x2] & gam_rxi[n1] is the empty set:
		y = (x1<Inf(gam_rxi[n1]))? gammar(x2) : gammar(x1);
		y *= q_gammarp;
	}
	else // [x1,x2] & gam_rxi[n1] is not the empty set:
		y = Sup(gam_ryi[n1]); 
	
	return y;
}

real gamr_even_Mi(const real& x1, const real& x2, int n1)
{
	real y,y1;
	
	if ( x2<Inf(gam_rxi[n1]) || Sup(gam_rxi[n1])<x1 ) 
	{  // [x1,x2] & gam_rxi[n1] is the empty set:
		std::cout << "Leere Menge:" << std::endl;
		y = (x1<Inf(gam_rxi[n1]))? gammar(x1) : gammar(x2);
		y *= q_gammarm;
	}
	else // [x1,x2] & gam_rxi[n1] is not the empty set:
	{
		y1 = gammar(x1)*q_gammarm;
		y  = gammar(x2)*q_gammarm;
		if (y1<y) y=y1;
	}
	
	return y;
}

real gamr_odd_Mi(const real& x1, const real& x2, int n1)
{
	real y;
	
	if ( x2<Inf(gam_rxi[n1]) || Sup(gam_rxi[n1])<x1 ) 
	{  // [x1,x2] & gam_rxi[n1] is the empty set:
		y = (x1<Inf(gam_rxi[n1]))? gammar(x2) : gammar(x1);
		y *= q_gammarp;
	}
	else // [x1,x2] & gam_rxi[n1] is not the empty set:
		y = Inf(gam_ryi[n1]); 
	
	return y;
}

real gamr_odd_Ma(const real& x1, const real& x2, int n1)
{
	real y,y1;
	
	if ( x2<Inf(gam_rxi[n1]) || Sup(gam_rxi[n1])<x1 ) 
	{  // [x1,x2] & gam_rxi[n1] is the empty set:
		std::cout << "Leere Menge:" << std::endl;
		y = (x1<Inf(gam_rxi[n1]))? gammar(x1) : gammar(x2);
		y *= q_gammarm;
	}
	else // [x1,x2] & gam_rxi[n1] is not the empty set:
	{
		y1 = gammar(x1)*q_gammarm;
		y  = gammar(x2)*q_gammarm;
		if (y1>y) y=y1;
	}
	
	return y;
}

interval gammar(const interval& x)
// Calculating inclusions of 1/Gamma(x);
// Blomquist, 01.07.09
{
	interval y;
	real x0, x1(Inf(x)), x2(Sup(x)), y0,y1(0),y2;
	int n1,n2;
	
	n1 = round_g(x1); 
	n2 = round_g(x2);
	if (x1==x2) // x is point interval
		if (x1==-n1) y2 = y1; // y = [0,0];
	else 
	{
		y1 = gammar(x1);
		y2 = y1;
			// Lower bound y1, Upper bound y2:
		if (y1<0) 
		{
			y1 = y1*q_gammarp;
			y2 = y2*q_gammarm;
		}
		else 
		{
			y1 = y1*q_gammarm;
			y2 = y2*q_gammarp;
		}
	}
	else // x2>x1:
	{
		if (n1%2==0) // n1 even, i.e. n1=0,2,4,6,...
		{
			if (n1==n2)
			{
				y1 = gamr_even_Mi(x1,x2,n1);
				y2 = gamr_even_Ma(x1,x2,n1);
			}
			else
				if (n2==n1-1)
			{
				y1 = gamr_odd_Mi(-n2,x2,n2);
				y2 = gamr_even_Ma(x1,-n2,n1);
			}
			else // n2 <= n1-2
			{	
				y1 = Inf(gam_ryi[n1-1]);  // Minimum OK
				y2 = gamr_even_Ma(x1,-n1+1,n1);
				x0 = x2;
				if (x2>n1-3 && x2<0) 
					x0 = n1-3;
				y0 = gamr_even_Ma(-n1+2,x0,n1-2);
				if (y0>y2) y2 = y0;
					
				if (n1==4 && n2==0)
				{
					y0 = gamr_even_Ma(0,x2,0);
					if (y0>y2) y2=y0;
				}
			}
		} // n1 even;
		else // n1 odd:
		{
			if (n1==n2)
			{
				y1 = gamr_odd_Mi(x1,x2,n1);
				y2 = gamr_odd_Ma(x1,x2,n1);
			}
			else
				if (n2==n1-1)
			{
				y1 = gamr_odd_Mi(x1,-n2,n1);
				y2 = gamr_even_Ma(-n2,x2,n2);
			}
			else
				if (n2==n1-2)
			{
				y1 = gamr_odd_Mi(x1,n1-1,n1);
				y2 = gamr_odd_Mi(-n1+2,x2,n1-2);
				if (y2<y1) y1 = y2; // Minimum calculated
				y2 = Sup( gam_ryi[n1-1] );
			}
			else // 0 <= n2 <= n1-3;
			{  // Calculating the minimum y1:
				y1 = gamr_odd_Mi(x1,n1-1,n1);
				y2 = Inf( gam_ryi[n1-2] );
				if (y2<y1) y1 = y2; // Minimum y1 calculated
						// Now calculating the maximum y2:
				if (n1==3)
				{
					y2 = Sup( gam_ryi[n1-1]);
					y0 = gamr_even_Ma(0,x2,0);
					if (y0>y2) y2=y0;
				}
				else // n1 = 5,7,9,....
					y2 = Sup( gam_ryi[n1-1] );
			}
						
		} // n1 odd
	}
	
	y = interval(y1,y2);	
	return y;
}

interval gamma(const interval& x)
// Calculating inclusions of 1/Gamma(x);
// Blomquist, 01.07.09
{
	return 1/gammar(x);
}


} // namespace cxsc

