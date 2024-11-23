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

/* CVS $Id: rmath.cpp,v 1.30 2014/01/30 17:23:48 cxsc Exp $ */

#include "rmath.hpp"

#include "rtsrmath.h"

namespace cxsc {
// Constants and values for the function ln_sqrtx2y2:
real ln2_1067(6505485212531678.0/8796093022208.0); // 1067*ln(2)
// Exponents of the interval bounds:
int B_lnx2y2_1[22] = {21, 31, 51, 101, 151, 201, 251, 301, 351, 401,
                      451, 501, 551, 601, 651, 701, 751, 801, 851, 
                      901, 951, 1025}; 
int B_lnx2y2_2[22] = {-1021,-949,-899,-849,-799,-749,-699,-649,-599,
                      -549,-499,-449,-399,-349,-299,-249,-199,-149,
                      -99,-49,-29,-19}; 
// Optimal values N for N*ln(2):
int B_lnx2y2_N1[21] ={20, 40, 61, 122, 160, 229, 259, 320, 366, 427,
                      488, 518, 549, 610, 671, 732, 763, 825, 885, 
                      945, 976};
real B_lnx2y2_c1[21] = 
{
    7804143460206699.0 / 562949953421312.0, // N*ln(2) with N = 20;
    7804143460206699.0 / 281474976710656.0, // N*ln(2) with N = 40;
    5950659388407608.0 / 140737488355328.0, // N*ln(2) with N = 61;
    5950659388407608.0 / 70368744177664.0,  // N*ln(2) with N = 122;
    7804143460206699.0 / 70368744177664.0,  // N*ln(2) with N = 160;
    5584840163710419.0 / 35184372088832.0,  // N*ln(2) with N = 229;
    6316478613104797.0 / 35184372088832.0,  // N*ln(2) with N = 259;
    7804143460206699.0 / 35184372088832.0,  // N*ln(2) with N = 320;
    8925989082611412.0 / 35184372088832.0,  // N*ln(2) with N = 366;
    5206826964856657.0 / 17592186044416.0,  // N*ln(2) with N = 427;
    5950659388407608.0 / 17592186044416.0,  // N*ln(2) with N = 488;
    6316478613104797.0 / 17592186044416.0,  // N*ln(2) with N = 518;
    6694491811958559.0 / 17592186044416.0,  // N*ln(2) with N = 549;
    7438324235509510.0 / 17592186044416.0,  // N*ln(2) with N = 610;
    8182156659060461.0 / 17592186044416.0,  // N*ln(2) with N = 671;
    8925989082611412.0 / 17592186044416.0,  // N*ln(2) with N = 732;
    4652001140732587.0 / 8796093022208.0,   // N*ln(2) with N = 763;
    5030014339586349.0 / 8796093022208.0,   // N*ln(2) with N = 825;
    5395833564283538.0 / 8796093022208.0,   // N*ln(2) with N = 885;
    5761652788980727.0 / 8796093022208.0,   // N*ln(2) with N = 945;
    5950659388407608.0 / 8796093022208.0    // N*ln(2) with N = 976;
};

int uint_trail(const unsigned int& n)
//  Liefert die Anzahl aufeinanderfolgender Null-Bits am Ende von n;
//  n=0 --> 32;   n=1 --> 0;   n=2 --> 1;   n=4 --> 2;   n=1024 --> 10;
{
    int m1,p=0;
    unsigned int m=n;
    if (m==0) p = 32;
    else { // m != 0, d.h. die folgende while-Schleife bricht ab!
       do
       {
	   m1 = m & 1;  // m1=0(1), wenn letztes Bit von m gleich 0(1) ist.
	   if (m1==0) p++;
	   m = m >> 1;  // Bit-Muster um 1 nach rechts schieben.
       } while(m1==0);
    }
    return p;
} // uint_trail

void sqr2uv(const real& x, real& u, real& v)
// Liefert u,v für: x2 = u + v; EXAKTE Darstellung, falls kein overflow 
// auftritt und v im normalisierten Bereich liegt. u > |v|
// Vorsicht: Funktioniert zunächst nur auf INTEL-Systemen!!!
{
    real a,b,t,y1,y2;
    a = Cut26(x);
    b = x-a;  // x = a+b;
    u = x*x;
    t = u - a*a; // exakte Auswertung!
    y2 = a*b;
    times2pown(y2,1); // y2 = 2*a*b, exakt!
    t -= y2;
// Jetzt fehlt noch:  t2 - b*b, aber b*b wird nicht immer korrekt berechnet,
// daher nochmalige Aufspaltung von b in y1+y2!!
    y1 = Cut25(b);
    y2 = b - y1;   // b = y1+y2, exakt;
    t -= y1*y1;
    if (sign(y2)!=0)
    {
	a = y1*y2;
	times2pown(a,1); // a = 2*y1*y2, exakt!
	t -= a;
	t -= y2*y2;
    }
    v = -t;
} // sqr2uv


//------------------------------------------------------------------

// real-Konstante  expmx2_x0  fuer die Funktion e^(-x2);
const real expmx2_x0 = 7491658466053896.0 / 281474976710656.0;
// Fuer x > expmx2_x0 werden die Funktionswerte auf Null gesetzt.
// Die relative Fehlerschranke e(f) := 4.618919E-16 gilt fuer
// alle |x| <= expmx2_x0 = 26.61571750925.... 

real expmx2(const real& x) throw()
// e^(-x^2);  rel. Fehlerschranke:  eps = 4.618919E-16 = e(f) gilt
// fuer alle |x| <= expmx2_x0 = 26.61571750925....
// Fuer |x| > expmx2_x0 --> expmx2(x) = 0;
// Blomquist, 05.07.04;
{ 
    real t=x,u,v,res=0;
    int ex;
    if (t<0) t = -t;  // t >= 0;
    ex = expo(t);
    if (ex<=-26) res = 1;  // t < 2^(-26)
    else if (ex<=-6)       // t < 2^(-6)
    {
	u = t*t;  v = u; 
	times2pown(v,-1);  // v: 0.5*x2
	res = 1-u*( 1-v*(1-u/3) );
    } else if (t<=expmx2_x0) {
	sqr2uv(x,u,v);  // u:= x*x,v aus S(2,53); x2 = u+v (exakt!)
	res = exp(-u); 
	if (v!=0) { 
	    times2pown(res,500);  // Die Skalierung verhindert, dass 
	    v *= res; // v*exp(-u) in den denormalisierten Bereich faellt
	    res -= v;
	    times2pown(res,-500); // Rueckskalierung
	}
    }
    return res;
} // expmx2

real expx2(const real& x)
// e^(+x^2);  rel. Fehlerschranke:  eps = 4.618958E-16 = e(f) gilt
// fuer alle |x| <= x0 = 26.64174755704632....
// x0 = 7498985273150791.0 / 281474976710656.0;
// Fuer |x| > x0 --> Programmabbruch
// Ausfuehrlich getestet;  Blomquist, 26.07.06;
{
    real t=x,u,v,res;
    int ex;
    if (t<0) t = -t;  // t >= 0;
    ex = expo(t);
    if (ex<=-26) res = 1;  // t < 2^(-26)
    else 
	if (ex<=-6)       // t < 2^(-6)
	{
	    u = t*t;  v = u;
	    times2pown(v,-1);  // v: 0.5*x^2
	    res = 1+u*( 1+v*(1+u/3) );
	} 
	else 
	{
	    sqr2uv(x,u,v);  // u := x*x und v aus S(2,53); 
                            // x^2 = u+v ist exakt!
	    res = exp(u);
	    v *= res;       // v*exp(+u)
	    res += v;
	}
    return res;
} // expx2

real expx2m1(const real& x)
// e^(+x^2)-1; rel. Fehlerschranke: eps = 4.813220E-16 = e(f) gilt
// fuer alle x, mit:  2^(-511) <= x <= x0 = 26.64174755704632....
// x0 = 7498985273150791.0 / 281474976710656.0;
// Fuer x > x0 --> Programmabbruch wegen Overflow;
// Fuer 0 < x < 2^(-511) --> Programmabbruch, denorm. Bereich!!
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
		    else if (ex>=-1073) 
		    {
			std::cerr << "expx2m1: denormalized range!" 
                             << std::endl; exit(1);
		    }
    }

    return res;
} // expx2m1



//------------------------------------------------------------------

real sqrt1px2(const real& x) throw()
// sqrt(1+x^2); Blomquist 13.12.02;
{
    if (expo(x) > 33) return abs(x);
    else return sqrt(1+x*x);
}

// Folgende Konstante sqrtp1m1_s wird gebraucht in 
// real sqrtp1m1(const real& x) throw();  Blomquist, 05.08.03
const real sqrtp1m1_s = 9007199254740984.0 / 1125899906842624.0;

real sqrtp1m1(const real& x) throw()
// sqrtp1m1(x) = sqrt(x+1)-1;
// Blomquist, 05.08.03; 
{
    real y = x;
    int ex = expo(x);
    if (ex<=-50) times2pown(y,-1);  // |x|<2^(-50); fast division by 2 
    else if (ex>=105) y = sqrt(x);  // x >= 2^(+104) = 2.02824...e+31
    else if (ex>=53) y = sqrt(x)-1; // x >= 2^(+52)  = 4.50359...e+15
    else if (x>-0.5234375 && x<=sqrtp1m1_s ) y = x / (sqrt(x+1) + 1);
    else y = sqrt(x+1)-1;
    return y;
}


//------------------------------------------------------------------

real sqrtx2y2(const real& x, const real& y) throw()
// calculating sqrt(x^2 + y^2) in high accuracy. Blomquist 01.12.02
{
    real a,b,r;
    dotprecision dot;
    int exa,exb,ex;
    a = x;   b = y;
    exa = expo(a);  exb = expo(b);
    if (exb > exa)
    { // Permutation of a,b:
	r = a;  a = b;  b = r;
	ex = exa;  exa = exb;  exb = ex;
    }  // |a| >= |b|
    // a = |a| >= |b| > 0:
    ex = 511 - exa;   // scaling a with 2^ex --> expo(a) == 511, --> a*a and
    times2pown(a,ex); // b*b without overflow. An underflow of b*b will not
    times2pown(b,ex); // affect the accuracy of the result!
    dot = 0;
    accumulate(dot,a,a);
    accumulate(dot,b,b);
    r = rnd(dot);
    r = sqrt(r); // sqrt(...) declared in rmath.hpp
    times2pown(r,-ex); // back-scaling
    return r;
} // sqrtx2y2

//------------------------------------------------------------------

real sqrtx2m1(const real& x) 
// sqrt(x^2-1);  rel. Fehlerschranke:  eps = 2.221305E-16 = e(f)
// Blomquist, 13.04.04;
{    const real c1 = 1.000732421875, 
	       c2 = 44000.0,
               c3 = 1024.0;  // c1,c2,c3 werden exakt gespeichert!
    real res,ep,ep2,s1,s2,x1,x2,arg=x;
    if (sign(arg)<0) arg = -arg; // arg = |x| >= 0
    if (arg <= c1) { // x = 1+ep; x^2-1 = 2*ep + ep^2;
	ep = x - 1; // Differenz rundungsfehlerfrei!
	ep2 = ep*ep;  // ep2 i.a. fehlerbehaftet!
	times2pown(ep,1);  // ep = 2*ep;
	res = sqrt(ep+ep2); // res=y0: Startwert;
        // x - y0^2 = (2*eps - s1^2) + [eps^2 - s2*(y0 + s1)]
	s1 = Cut26(res);  s2 = res - s1; // Startwert y0 = s1 + s2;
	arg = ep - s1*s1;  // arg = 2*eps - s1^2;
	arg += (ep2 - s2*(res+s1));  // arg = x - y0^2
	if (sign(arg)>0) {
	arg = arg / res;
	times2pown(arg,-1);
	res += arg;  // 1. Newton-Schritt beendet; eps = 2.221261E-16
	}
    } else if (arg<c2) { 
        // x-y0^2 = [(x1^2-1)-s1^2] + [x2*(x+x1)-s2*(y0+s1)]
	x1 = Cut26(arg);  x2 = arg - x1;  // arg = x = x1 + x2;
	ep2 = x2*(arg+x1);  // ep2 ist fehlerbehaftet
	x2 = x1*x1;  ep = x2-1;
	res = sqrt(ep+ep2); // res ist Startwert für Newton-Verfahren
	s1 = Cut26(res);  s2 = res - s1;  // Startwert = s1 + s2;
	ep2 = ep2 - s2 * (res+s1); // ep2 = [x2*(x+x1)-s2*(y0+s1)]
	if (arg<c3) ep -= s1*s1;   // ep = (x1^2-1) - s1^2;
	else {
	    x2 -= s1*s1;  // x2 = x1^2-s1^2
	    ep = x2 - 1; }         // ep = (x1^2-s1^2) - 1;
	ep += ep2;        // ep = x - y0^2;
	ep /= res;
	times2pown(ep,-1);
	res = res + ep;  // 1. Newton-Schritt in hoher Genauigkeit
                         // beendet;  eps = 2.221305E-16
    } else { // arg = |x| >= 44000;
	res = -1/arg;
	times2pown(res,-1); // Multiplikation mit 0.5;
	res += arg;  // res = x - 1/(2*x);  eps = 2.221114E-16
    }
    return res;
} // sqrtx2m1 (Punktargumente)

//------------------------------------------------------------------

real sqrt1mx2(const real& x) throw(STD_FKT_OUT_OF_DEF)
// sqrt(1-x2);  rel. Fehlerschranke:  eps = 3.700747E-16 = e(f)
// Blomquist, 19.06.04;
{ 
    real t=x,res;
    int ex;
    if (sign(t)<0) t = -t; // Argument t >=0;
    if (t>1) cxscthrow(STD_FKT_OUT_OF_DEF("real sqrt1mx2(const real&)"));
    // For argument t now it holds: 0 <= t <=1;
    ex = expo(t);
    if (ex<=-26) res = 1; // t < 2^(-26) --> res = 1
    else if (ex<=-15) {   // t < 2^(-15) --> res = 1-x2/2
	res = t*t;  
	times2pown(res,-1);
	res = 1 - res;
    } else {
	if (ex>=0) 
	{  // ex>=0 --> t>=0.5; 
	    res = 1-t;   // res: delta = 1-t;
	    t = res * res;
	    times2pown(res,1); // res: 2*delta
	    res = res - t;     // res: 1-x2 = 2*delta - delta2
	} else res = 1-t*t;    // res: Maschinenwert von 1-x2
	res = sqrt(res);       // res: Nullte Naeherung 
    }
    return res;
} // sqrt1mx2


//------------------------------------------------------------------

int Interval_Nr(int* v, const int& n, const int& ex)
// n>0 subintervals:   |...\|...\|...\ ..... |...|
// subinterval Nr.:      0    1    2   .....  n-1  
{
    int i=0,j=n,k;  // n>0:  Number of subintervals
    do {
	k = (i+j)/2;
	if (ex < v[k]) j = k-1;
	else i = k+1;
    } while(i<=j);
    return j;  // x with ex=expo(x) lies in the subinterval number j  
}

//------------------------------------------------------------------

real ln_sqrtx2y2(const real& x, const real& y) 
                                throw(STD_FKT_OUT_OF_DEF)
// ln( sqrt(x^2+y^2) ) == 0.5*ln(x^2+y^2); Blomquist, 21.11.03;
// Relative error bound: 5.160563E-016;
// Absolute error bound: 2.225075E-308; if x=1 and 0<=y<=b0;
{
    int j,N;
    real a,b,r,r1;
    dotprecision dot;
   
    a = sign(x)<0 ? -x : x;  // a = |x| >= 0;
    b = sign(y)<0 ? -y : y;  // b = |y| >= 0;
    int exa=expo(a), exb=expo(b), ex;
    if (b > a)
    {
	r = a;   a = b;   b = r;
	ex = exa;   exa = exb;   exb = ex;
    }
    // It holds now:   0 <= b <= a 
    if (sign(a)==0)
        cxscthrow(STD_FKT_OUT_OF_DEF
                    ("real ln_sqrtx2y2(const real&, const real&)"));
    if (exa>20) // to avoid overflow by calculating a^2 + b^2
    {  // a>=2^(20):
	j = Interval_Nr(B_lnx2y2_1,21,exa); // j: No. of subinterval
	N = B_lnx2y2_N1[j];    // N: Optimal int value
	if (exb-exa > -25)
	{   // For (exb-exa>-25) we use the complete term:  
            // N*ln(2) + [ln(2^(-N)*a)+0.5*ln(1+(b/a)^2)]
	    b = b/a;  // a > 0
	    b = lnp1(b*b);
	    times2pown(b,-1);  // exact division by 2
	    times2pown(a,-N);
	    r = b + ln(a); // [ ... ] calculated!
	    r += B_lnx2y2_c1[j];
	}
	else { // For (exb-exa<=-25) only two summands!: 
	    times2pown(a,-N);
	    r = ln(a) + B_lnx2y2_c1[j];
	}
    }
    else  // exa<=20 or a<2^(20):
    {     // Now calculation of a^2+b^2 without overflow:
	if (exa<=-20) // to avoid underflow by calculating a^2+b^2
	    if (exa<=-1022) // a in the denormalized range
	    {
		r = b/a;
		r = lnp1(r*r);  times2pown(r,-1); // r: 0.5*ln(1+..)
		times2pown(a,1067);
		r += ln(a);     // [ .... ] ready
		r -= ln2_1067;  // rel. error = 2.459639e-16;
	    }
	    else // MinReal=2^(-1022) <= a < 2^(-20)
	    {    // Calculating the number j of the subinterval:
		j = 20 - Interval_Nr(B_lnx2y2_2,21,exa);
		r = a;  times2pown(r,B_lnx2y2_N1[j]);
		r = ln(r);  // r: ln(2^N*a);
		if (exb-exa > -25) { // calculating the complete term
		    b = b/a;
		    a = lnp1(b*b);
		    times2pown(a,-1);
		    r += a;  // [ ... ] ready now
		}
                // We now have: exb-exa<=-25,  ==>  b/a <= 2^(-24);
		r -= B_lnx2y2_c1[j]; // 0.5*ln(1+(b/a)^2) neglected!
                // relative error = 4.524090e-16 in both cases;
	    }
	else // calculation of a^2+b^2 without overflow or underflow:
	{   // exa>-20  respective  a>=2^(-20):
	    dot = 0;
	    accumulate(dot,a,a);
	    accumulate(dot,b,b);  // dot = a^2+b^2, exact!
	    real s = rnd(dot);    // s = a^2 + b^2, rounded!
	    if (s>=0.25 && s<=1.75)
		if (s>=0.828125 && s<=1.171875)
		{ // Series:
		    if (a==1 && exb<=-28)
		    {
			r = b;  times2pown(r,-1);
			r *= b;
		    }
		    else {
			dot -= 1;
			r = rnd(dot); // r = a^2+b^2-1 rounded!
			r = lnp1(r);
			times2pown(r,-1);
		    }
		}
		else { // Reading dot = a^2+b^2 twice:
		    r = rnd(dot);
		    dot -= r;
		    r1 = rnd(dot); // a^2+b^2 = r+r1, rounded!
		    r1 = lnp1(r1/r);
		    r = ln(r) + r1;
		    times2pown(r,-1); // exact division by 2
		}
	    else { // calculating straight from: 0.5*ln(x^2+y^2)
		r = ln(s);
		times2pown(r,-1);
	    } 
	}
    }
    return r;
} // ln_sqrtx2y2

typedef union { double f; char intern[8]; } help_real;

real Cut24(const real& x){
    // y = Cut24(x) liefert ein y, das mit den ersten 24 Mantissenbits
    // von x übereinstimmt, das hidden bit ist dabei mitgezählt!
    // Die restlichen 53-24=29 Mantissenbits werden auf Null gesetzt.
  help_real y;
  y.f = _double(x);
#if INTEL
  y.intern[3] &= 224;
  y.intern[0] = y.intern[1] = y.intern[2] = 0;
#else
  y.intern[4] &= 224;
  y.intern[7] = y.intern[6] = y.intern[5] = 0;
#endif
  return real(y.f);
}

real Cut25(const real& x){
    // y = Cut25(x) liefert ein y, das mit den ersten 25 Mantissenbits
    // von x übereinstimmt, das hidden bit ist dabei mitgezählt!
    // Die restlichen 53-25=28 Mantissenbits werden auf Null gesetzt.
  help_real y;
  y.f = _double(x);
#if INTEL
  y.intern[3] &= 240;
  y.intern[0] = y.intern[1] = y.intern[2] = 0;
#else
  y.intern[4] &= 240;
  y.intern[7] = y.intern[6] = y.intern[5] = 0;
#endif
  return real(y.f);
}

real Cut26(const real& x){
    // y = Cut26(x) liefert ein y, das mit den ersten 26 Mantissenbits
    // von x übereinstimmt, das hidden bit ist dabei mitgezählt!
    // Die restlichen 53-26=27 Mantissenbits werden auf Null gesetzt.
  help_real y;
  y.f = _double(x);
#if INTEL
  y.intern[3] &= 248;
  y.intern[0] = y.intern[1] = y.intern[2] = 0;
#else
  y.intern[4] &= 248;
  y.intern[7] = y.intern[6] = y.intern[5] = 0;
#endif
  return real(y.f);
}

int Round(const real& x) throw()
// y = Round(x) delivers the rounded value y of x.
// For |x| < 2147483647.5 the assignment  y = Round(x)  delivers
// the next integer value y.
// Examples:
// y = Round(-0.1);            --->  y = 0;
// y = Round(-123.5);          --->  y = -124.0;
// y = Round(+123.5);          --->  y = +124.0;
// y = Round(+123.499);        --->  y = +123.0;		
// y = Round(+2147483647.499); --->  y = +2147483647.0;
// y = Round(+2147483647.5);   --->  y = -2147483648.0;
// y = Round(-2147483647.499); --->  y = -2147483647.0;
// y = Round(-2147483647.5);   --->  y = -2147483648.0;
// y = Round(-2147483648.501); --->  y = -2147483648.0; wrong rounding!
// Blomquist, 29.10.2008;
{
	double dbl;
	
	dbl = _double(x);
	return dbl<0 ? int(dbl-0.5) : int(dbl+0.5);
}

int ceil(const real& x) throw()
{
	real y(x); 
	bool neg(y<0);
	if (neg) y = -y; // y >= 0;
	int res = int( _double(y) );
	y = y - res;
	if (!neg && y>0) res += 1;
	if (neg) res = -res;
	return res;
}

int ifloor(const real& x) throw()
{
	real y(x); 
	bool neg(y<0);
	if (neg) y = -y; // y >= 0;
	int res = int( _double(y) );
	y = y - res;
	if (neg) 
	{
		res = -res;
		if (y>0) res -= 1;
	}
	return res;
}

static real q_acoshp1[5] = // Polynomial coefficients of Q_4(x) 
           // roundet to nearest. acosh(1+x) = sqrt(2*x)*Q_4(x)
{ 1.0 / 1.0,                                    // q_acoshp1[0]
  -6004799503160661.0 / 72057594037927936.0,
  +5404319552844595.0 / 288230376151711744.0,
  -6433713753386423.0 / 1152921504606846976.0,
  +8756999275442631.0 / 4611686018427387904.0   // q_acoshp1[4]
};

static const real c_ln2_B = 6243314768165359.0 / 9007199254740992.0;
// c_ln2_B < ln(2) is the nearest machine number for ln(2) with an
// absolute error < 2.3190469E-17;

real acoshp1(const real& x) throw()
// acoshp1(x) = acosh(1+x);  rel. error: eps = 7.792706E-16 = e(f)
// Ausfuehrlich getestet;  Blomquist, 27.03.05;
{ 
    real res;
    int ex(expo(x));
    if (x<0) 
	cxscthrow(STD_FKT_OUT_OF_DEF("real acoshp1(const real&)"));
    // For argument x now it holds: 0 <= x <= MaxReal;
    if (ex<=-50) res = sqrt(2*x); // 0<=x<2^(-50): acoshp1(x)=sqrt(2x)
    else if (ex<=-9) // 2^(-50)<=x<2^{-9}: acoshp1(x)=sqrt(2x)*Q_4(x)
      res = sqrt(2*x)*((((q_acoshp1[4]*x+q_acoshp1[3])*x+q_acoshp1[2])
			    *x+q_acoshp1[1])*x + q_acoshp1[0]);
    else if (ex<=0) res = lnp1(x+sqrt(2*x+x*x));     // range A_3
    else if (ex<=50) res = lnp1(x*(1+sqrt(1+2/x)));  // range A_4
    else if (ex<=1022) res = ln(2*x);                // range A_5
    else res = ln(x) + c_ln2_B;                      // range A_6

    return res;
} // acoshp1


// *****************************************************************************
//                               sin(Pi*x)/Pi                                  *
// ********************************************************************************
//   Sinpix_pi:           Teilpunkte des Intervalls [0,0.5]:
// ********************************************************************************
real a_sinpix_pi[8] = {0,
                       7.450580596923828125e-9, // 2^(-27)
                       0.05078125,
                       0.1015625,
                       0.15234375,
                       0.2578125,
                       0.3828125,
                       0.5 };
// ********************************************************************************

// ********************************************************************************
//  Sinpix_pi:      a_k,b_k des Kettenbruchs K2 in A1 = [2^-27,0.05078125]:
// ********************************************************************************
  const real q_sin1_a[3] = {0.0,
                            -7408124450506663.0 / 4503599627370496.0,
                            +4595816915236340.0 / 36028797018963968.0 };

  const real q_sin1_b[3] = {+1.0,
                            +8889749340277122.0 / 18014398509481984.0,
                            -6103596185201565.0 / 36028797018963968.0};
// ********************************************************************************

// ********************************************************************************
//  Sinpix_pi:     a_k,b_k des Kettenbruchs K4 in A2 = [0.05078125,0.1015625]:
// ********************************************************************************
   const real q_sin2_a[5] = {0.0,
                             -4829370543434630.0 / 18014398509481984.0,
                             +5229154385037560.0 / 140737488355328.0,
                             -6000119407941526.0 / 576460752303423488.0,
                             +4592842702840413.0 / 1125899906842624.0 };

  const real q_sin2_b[5] = {-6359670668682560.0 / 576460752303423488.0,
                            -6771299865719208.0 / 1125899906842624.0,
                            +6861350950372367.0 / 1125899906842624.0,
                            -4669728680615660.0 / 2251799813685248.0,
                            +4607746187351458.0 / 2251799813685248.0 };
// ********************************************************************************

// ********************************************************************************
//   Sinpix_pi:    a_k,b_k des Kettenbruchs K4 in A3 = [0.1015625,0.15234375]:
// ********************************************************************************
  const real q_sin3_a[5] = {0.0,
	                         -7954137925457403.0 / 18014398509481984.0,
								    +7534721222201852.0 / 562949953421312.0,
				                -8481413201661564.0 / 288230376151711744.0,
						          +6484759126571114.0 / 4503599627370496.0};
  const real q_sin3_b[5] = {-8780869688499296.0 / 288230376151711744.0,
	                         -7929689932517363.0 / 2251799813685248.0,
								    +8223190551048856.0 / 2251799813685248.0,
				                -5789057670884168.0 / 4503599627370496.0,
						          +5586166285064551.0 / 4503599627370496.0 };
// ********************************************************************************

// ********************************************************************************
//   Sinpix_pi:    a_k,b_k des Kettenbruchs K5 in A4 = [0.15234375,0.2578125]:
// ********************************************************************************
  const real q_sin4_a[6] = {+0.0,
	                         -5469389903424399.0 / 9007199254740992.0,
								    +7704524171392577.0 / 1125899906842624.0,
				                -8513632445057187.0 / 144115188075855872.0,
						          +6447776473700967.0 / 9007199254740992.0,
					             -7362728796775514.0 / 288230376151711744.0 };

  const real q_sin4_b[6] = {-8529339040415872.0 / 144115188075855872.0,
	                         -5452397600582528.0 / 2251799813685248.0,
								    +5848879160432416.0 / 2251799813685248.0,
				                -8609222354677976.0 / 9007199254740992.0,
						          +8031029418839465.0 / 9007199254740992.0,
					             -4904826237544497.0 / 9007199254740992.0 };
// ********************************************************************************

// ********************************************************************************
//   Sinpix_pi:     a_k,b_k des Kettenbruchs K5 in A5 = [0.2578125,0.3828125]:
// ********************************************************************************
  const real q_sin5_a[6] = {0.0,
                            4*(5276398025110506.0 / 9007199254740992.0),
                            +7169490078506431.0 / 1125899906842624.0,
                            -7877917296929161.0 / 36028797018963968.0,
                            +5470058796541645.0 / 9007199254740992.0,
                            -7233501679212276.0 / 144115188075855872.0 };

  const real q_sin5_b[6] = {+4598161403289824.0 / 144115188075855872.0,
                            +4893639363223553.0 / 2251799813685248.0,
                            -5525704702280927.0 / 2251799813685248.0,
                            +4608444394217766.0 / 4503599627370496.0,
                            -8018266172128683.0 / 9007199254740992.0,
                            +5822429161405618.0 / 9007199254740992.0 };
// ********************************************************************************

// ********************************************************************************
//   Sinpix_pi:   a_k,b_k des Kettenbruchs K5 in A6 = [0.3828125,0.5]:
// ********************************************************************************
  const real q_sin6_a[6] = {0.0,
                            +7461973733147728.0 / 36028797018963968.0,
                            +7979710108639590.0 / 140737488355328.0,
                            -6289496525317218.0 / 288230376151711744.0,
                            +6964694082726429.0 / 1125899906842624.0,
                            -5522978779702360.0 / 1152921504606846976.0 };

  const real q_sin6_b[6] = {+5609829449006976.0 / 18014398509481984.0,
                            +8354019229473476.0 / 1125899906842624.0,
                            -8475200820952730.0 / 1125899906842624.0,
                            +5829377160898302.0 / 2251799813685248.0,
                            -5714036688191727.0 / 2251799813685248.0,
                            +7900583259891052.0 / 4503599627370496.0 };
// ********************************************************************************

// interne Funktion, wird fuer Gammafunktion benoetigt
int int_no(real *a,
           const int n,    // n: Anzahl der Elemente des Feldes a
           const real& x)  // x: Eine real-Zahl
// Ein Intervall [A,B] wird durch ein Feld a mit n Elementen in
// n-1 Teilintervalle unterteilt. Fuer ein x aus [A,B] wird eine
// Intervall-Nummer zurueckgegeben, wobei das erste Teilintervall
// die Nr. 0 und das letzte Teilintervall die Nr. n-2 erhaelt.
// Fuer  x<A  ist Nr. = -1, und fuer  x>=B  gilt  Nr. = n-1; 
{
    int i,j,k;
    i = 0;
    j = n-1;
    do 
    {
		k = (i+j)/2;
		if (x<a[k]) j = k-1;
		else i = k+1;
    } 
    while (i <= j);

    return j;
}
// *****************************************************************************
//    Sinpix_pi:         Approximation in A1 = [2^-27,0.05078125]:
// *****************************************************************************

real sinpi_A1(const real& x)
{
	real y,v;
	
	v = 1/(x*x);
	// Approximation:  sin(pi*x)/pi = x*K2(v);
	y = q_sin1_a[2] / (  v + q_sin1_b[2]);
	y = q_sin1_a[1] / ( (v + q_sin1_b[1]) + y);
	y *= x;
	y += x;
	
	return y;
}

// ********************************************************************************
//    Sinpix_pi:        Approximation in A2 = [0.05078125,0.1015625]:
// ********************************************************************************
real sinpi_A2(const real& x)
{
	real y,v;
	
	if (x == 0.08203125) 
		y = q_sin2_b[0];
	else
	{
		v = 1/(x-0.08203125);
		// Approximation:  sin(pi*x)/pi = x*K4(v);
		y = q_sin2_a[4] / (  v + q_sin2_b[4]);
		y = q_sin2_a[3] / ( (v + q_sin2_b[3]) + y);
		y = q_sin2_a[2] / ( (v + q_sin2_b[2]) + y);
		y = q_sin2_a[1] / ( (v + q_sin2_b[1]) + y) + q_sin2_b[0];
	}
	y *= x;
	y += x;
	
	return y;
}

// ********************************************************************************
//    Sinpix_pi:        Approximation in A3 = [0.1015625,0.15234375]:
// ********************************************************************************
real sinpi_A3(const real& x)
{
	real y,v;
	
	if (x == 0.13671875)
		y = q_sin3_b[0];
	else
	{
		v = 1/(x-0.13671875);
		// Approximation:  sin(pi*x)/pi = x*K4(v);
		y = q_sin3_a[4] / (  v + q_sin3_b[4]);
		y = q_sin3_a[3] / ( (v + q_sin3_b[3]) + y);
		y = q_sin3_a[2] / ( (v + q_sin3_b[2]) + y);
		y = q_sin3_a[1] / ( (v + q_sin3_b[1]) + y) + q_sin3_b[0];
	}
	y *= x;
	y += x;
		
	return y;
}

// ********************************************************************************
//   Sinpix_pi:        Approximation in A4 = [0.15234375,0.2578125]:
// ********************************************************************************
real sinpi_A4(const real& x)
{
	real y,v;
	
	if (x == 0.19140625)
		y = q_sin4_b[0];
	else
	{
		v = 1/(x-0.19140625);
		// Approximation:  sin(pi*x)/pi = x*K5(v);
		y = q_sin4_a[5] / (  v + q_sin4_b[5]);
		y = q_sin4_a[4] / ( (v + q_sin4_b[4]) + y);
		y = q_sin4_a[3] / ( (v + q_sin4_b[3]) + y);
		y = q_sin4_a[2] / ( (v + q_sin4_b[2]) + y);
		y = q_sin4_a[1] / ( (v + q_sin4_b[1]) + y) + q_sin4_b[0];
	}
	y *= x;
	y += x;
	
	return y;
}

// ********************************************************************************
//   Sinpix_pi:        Approximation in A5 = [0.2578125,0.3828125]:
// ********************************************************************************
real sinpi_A5(const real& x)
{
	real y,v;
	
	if (x == 0.30078125)
		y = q_sin5_b[0];
	else
	{
		v = 1/(x-0.30078125);
		// Approximation:  sin(pi*x)/pi = K5(v);
		y = q_sin5_a[5] / (  v + q_sin5_b[5]);
		y = q_sin5_a[4] / ( (v + q_sin5_b[4]) + y);
		y = q_sin5_a[3] / ( (v + q_sin5_b[3]) + y);
		y = q_sin5_a[2] / ( (v + q_sin5_b[2]) + y);
		y = q_sin5_a[1] / ( (v + q_sin5_b[1]) + y) + q_sin5_b[0];
	}
	y +=1.0;
	times2pown(y,-2); // Division durch 4
	
	return y;
}

// ********************************************************************************
//    Sinpix_pi:         Approximation in A6 = [0.3828125,0.5]:
// ********************************************************************************
real sinpi_A6(const real& x)
{
	real y,v;
	
	if (x == 0.43359375)
		y = q_sin6_b[0];
	else
	{
		v = 1/(x-0.43359375);
		// Approximation:  sin(pi*x)/pi = K5(v);
		y = q_sin6_a[5] / (  v + q_sin6_b[5]);
		y = q_sin6_a[4] / ( (v + q_sin6_b[4]) + y);
		y = q_sin6_a[3] / ( (v + q_sin6_b[3]) + y);
		y = q_sin6_a[2] / ( (v + q_sin6_b[2]) + y);
		y = q_sin6_a[1] / ( (v + q_sin6_b[1]) + y) + q_sin6_b[0];
	}
	
	return y;
}

real sinpix_pi(const real& x)
{  // Sin(Pi*x)/Pi;  rel. Fehlerschranke: 4.870167e-16;
	real res(0),xr;
	int m,nr;
	bool neg;

	m = Round(x);
	if (m == -2147483648)
		cxscthrow(STD_FKT_OUT_OF_DEF("real sinpix_pi(const real&)"));
	xr = x - m;  // xr: reduced argument, exactly calculated!
   // |xr| <= 0.5;
	neg = xr<0;
	if (neg) xr = -xr;
   // 0 <= xr <= 0.5;
	nr = int_no(a_sinpix_pi,8,xr);
	
	switch(nr)
	{
		case 0: res = xr;           break;
		case 1: res = sinpi_A1(xr); break;
		case 2: res = sinpi_A2(xr); break;
		case 3: res = sinpi_A3(xr); break;
		case 4: res = sinpi_A4(xr); break;
		case 5: res = sinpi_A5(xr); break;
		case 6: res = sinpi_A6(xr); break;		
		case 7: res = 5734161139222659.0 / 18014398509481984.0;
	}
	
	if (neg) res = -res;
	if (m%2 != 0) res = -res;
	
	return res;
}

// ********************************************************************************
// ***************  Konstanten bez. Gamma(x) und 1/Gamma(x) ********************
// ********************************************************************************

// Interval bounds of the neighbour intervals with respect to the
// intervals S_k for Gamma(x):
real gam_f85[19] =
{-0.5, 8.5, 16.5, 24.5, 35.5, 46.5, 57.5, 68.5, 79.5, 90.5,
 101.5, 112.5, 122.5, 132.5, 142.5, 150.5, 157.5, 164.5, 171.5 };

 // 1/Gamma(x): S0=[1.5,2.5], x0 = 2.0;
  const real q_gams0_a[8] = 
  {0.0, 
   -7616205496030159.0 / 18014398509481984.0,
   6808968414757234.0 / 9007199254740992.0,
  -7179656013820698.0 / 144115188075855872.0,
  7646522333236288.0 / 144115188075855872.0,
  -7301751514589296.0 / 144115188075855872.0,
  6062274112599236.0 / 288230376151711744.0,
  -7580146497393670.0 / 576460752303423488.0 };
	 
  const real q_gams0_b[8] = 
  {1.0,
  -4965940208012024.0 / 9007199254740992.0,
  8627036189024519.0 / 9007199254740992.0,
  -7819065539096245.0 / 72057594037927936.0,
  7433219784613049.0 / 36028797018963968.0,
  7389378817674059.0 / 72057594037927936.0,
  -6060163955665826.0 / 144115188075855872.0,
  5769165265609039.0 / 18014398509481984.0 };
	 
  // Gamma(x): S1=[10.5,11.5], x0 = 11.125;
  const real q_gams1_a[7] = 
  {0.0, 
  5262165366811426.0 / 72057594037927936.0,
  5574236285326272.0 / 9007199254740992.0,
  -8823729115230752.0 / 9223372036854775808.0,
  4639548867796070.0 / 72057594037927936.0,
  -4784760610310013.0 / 4611686018427387904.0,
  5900435828568799.0 / 288230376151711744.0 };
	 
  const real q_gams1_b[7] = 
  {7108556377252296.0 / 36028797018963968.0,
  -7219014979973550.0 / 9007199254740992.0,
  7224885284258947.0 / 9007199254740992.0,
  -8569030245421939.0 / 36028797018963968.0,
  4838653997792226.0 / 18014398509481984.0,
  -4567015707194151.0 / 36028797018963968.0,
  5705105513288442.0 / 36028797018963968.0 };
	 
  // Gamma(x): S2=[18.5,19.5], x0 = 18.96875;
  const real q_gams2_a[6] = 
  {0.0, 
  7322606177885925.0 / 72057594037927936.0,
  5856515731454725.0 / 144115188075855872.0,
  -5420981868254561.0 / 1152921504606846976.0,
  6209135328051357.0 / 1152921504606846976.0,
  -8261990134869242.0 / 2305843009213693952.0 };
	 
  const real q_gams2_b[6] = 
  {-5267325780498998.0 / 18014398509481984.0,
  -4688643295042637.0 / 18014398509481984.0,
  6865435581764388.0 / 36028797018963968.0,
  -5147410677045000.0 / 144115188075855872.0,
  5878944809110092.0 / 144115188075855872.0,
  5315270486582075.0 / 576460752303423488.0 };
  
  // Gamma(x): S3=[29.5,30.5], x0 = 29.9375;
  const real q_gams3_a[6] = 
  {0.0, 
  4608894695736103.0 / 9007199254740992.0,
  4622056591672246.0 / 144115188075855872.0,
  7936321632956658.0 / 1152921504606846976.0,
  7533691379187725.0 / 2305843009213693952.0,
  8317076627520632.0 / 4611686018427387904.0 };
	 
  const real q_gams3_b[6] = 
  {-5793090370845400.0 / 36028797018963968.0,
  -5993724737769479.0 / 18014398509481984.0,
  -7355055148590989.0 / 288230376151711744.0,
  -6811362569175008.0 / 288230376151711744.0,
  -5643865380181690.0 / 288230376151711744.0,
  -7257737463164532.0 / 288230376151711744.0 };
  
  // Gamma(x): S4=[40.5,41.5], x0 = 41.140625;
  const real q_gams4_a[6] = 
  {0.0, 
  -7504135817814391.0 / 18014398509481984.0,
  4990516833351222.0 / 576460752303423488.0,
  7409944094502278.0 / 4611686018427387904.0,
  8239437277348520.0 / 2305843009213693952.0,
  -5361398870955224.0 / 4611686018427387904.0 };
	 
  const real q_gams4_b[6] = 
  {7405548010914168.0 / 18014398509481984.0,
  6819421126036941.0 / 36028797018963968.0,
  8474998846370093.0 / 288230376151711744.0,
  7733944242652532.0 / 144115188075855872.0,
  -6547033063052812.0 / 144115188075855872.0,
  -6636555235967309.0 / 576460752303423488.0 };
  
  // Gamma(x): S5=[51.5,52.5], x0 = 52.015625;
  const real q_gams5_a[6] = 
  {0.0, 
  -7282435522169731.0 / 144115188075855872.0,
  7812862759300002.0 / 288230376151711744.0,
  -7591789799457117.0 / 9223372036854775808.0,
  7290376444377792.0 / 2305843009213693952.0,
  -7051827298983325.0 / 9223372036854775808.0 };
	 
  const real q_gams5_b[6] = 
  {-4692552597358020.0 / 36028797018963968.0,
  7065248316566075.0 / 36028797018963968.0,
  -5667667511562837.0 / 36028797018963968.0,
  7741409329322298.0 / 144115188075855872.0,
  -6371025742549657.0 / 144115188075855872.0,
  8045016330906311.0 / 288230376151711744.0 };
  
  // Gamma(x): S6=[62.5,63.5], x0 = 63.015625;
  const real q_gams6_a[6] = 
  {0.0, 
  6719861857699617.0 / 36028797018963968.0,
  6146108740657636.0 / 1152921504606846976.0,
  -8996491371873855.0 / 4611686018427387904.0,
  5082150623010284.0 / 2305843009213693952.0,
  -4650554434211955.0 / 18446744073709551616.0 };
	 
  const real q_gams6_b[6] = 
  {6795471792542416.0 / 18014398509481984.0,
  -4567367130077106.0 / 36028797018963968.0,
  8403023343856889.0 / 288230376151711744.0,
  5083390457138636.0 / 144115188075855872.0,
  -4614092542204999.0 / 144115188075855872.0,
  7104477795873427.0 / 72057594037927936.0 };
  
  // Gamma(x): S7=[73.5,74.5], x0 = 74.16015625;
  const real q_gams7_a[6] = 
  {0.0 / 1.0,
  5372276754422009.0 / 18014398509481984.0,
  4663470139182266.0 / 576460752303423488.0,
  8173857739398457.0 / 4611686018427387904.0,
  4932343261664244.0 / 4611686018427387904.0,
  -5822870305854791.0 / 295147905179352825856.0 };
	 
  const real q_gams7_b[6] = 
  {-4806355488125952.0 / 1152921504606846976.0,
  -6211400907258787.0 / 36028797018963968.0,
  -5429997192499743.0 / 288230376151711744.0,
  -5653057298064211.0 / 288230376151711744.0,
  -4746908442350821.0 / 1152921504606846976.0,
  8670733620601602.0 / 18014398509481984.0 };
  
  // Gamma(x): S8=[84.5,85.5], x0 = 85.1015625;
  const real q_gams8_a[5] = 
  { 0.0 / 1.0,
  -8754830070807807.0 / 72057594037927936.0,
  7931975738629914.0 / 2305843009213693952.0,
  6414800170661986.0 / 73786976294838206464.0,
  4558538137025832.0 / 18014398509481984.0};
	 
  const real q_gams8_b[5] = 
  {-4924952929015874.0 / 18014398509481984.0,
  8571263173618757.0 / 72057594037927936.0,
  7912939311540494.0 / 576460752303423488.0,
  8975599623306617.0 / 18014398509481984.0,
  -4569152133381960.0 / 9007199254740992.0 };
  
  // Gamma(x): S9=[95.5,96.5], x0 = 95.984375;
  const real q_gams9_a[5] = 
  {0.0 / 1.0,
  -6140046099729716.0 / 144115188075855872.0,
  7278997066304035.0 / 576460752303423488.0,
  -4759432412103962.0 / 9223372036854775808.0,
  6912711986126027.0 / 4611686018427387904.0 };
	 
  const real q_gams9_b[5] = 
  {-5611185402650416.0 / 72057594037927936.0,
  4915638659703209.0 / 36028797018963968.0,
  -7692469721238990.0 / 72057594037927936.0,
  5045754589592299.0 / 144115188075855872.0,
  -8290188163752846.0 / 288230376151711744.0 };
  
  // Gamma(x): S10=[106.5,107.5], x0 = 107.078125;
  const real q_gams10_a[5] = 
  {0.0 / 1.0,
  4717576521494070.0 / 72057594037927936.0,
  6906626643101502.0 / 1152921504606846976.0,
  -8147816895338065.0 / 9223372036854775808.0,
  7960901442281121.0 / 9223372036854775808.0 };
	 
  const real q_gams10_b[5] = 
  {7952058244493952.0 / 288230376151711744.0,
  -7601355771801470.0 / 72057594037927936.0,
  4923655207606594.0 / 72057594037927936.0,
  -7238499638430702.0 / 576460752303423488.0,
  5852764550899526.0 / 576460752303423488.0 };
  
  // Gamma(x): S11=[117.5,118.5], x0 = 117.8671875;
  const real q_gams11_a[5] = 
  {0.0 / 1.0,
  5000157748740957.0 / 36028797018963968.0,
  6733781271808263.0 / 2305843009213693952.0,
  4912448110743899.0 / 18446744073709551616.0,
  6917678469221762.0 / 576460752303423488.0 };
	 
  const real q_gams11_b[5] = 
  {-4805173503400964.0 / 36028797018963968.0,
  -7686561799456867.0 / 72057594037927936.0,
  -6649023763933472.0 / 576460752303423488.0,
  -7324214555543306.0 / 72057594037927936.0,
  8284786467509721.0 / 72057594037927936.0 };
  
  // Gamma(x): Stuetzinterval S12=[126.5,127.5], x0 = 126.7421875;
  const real q_gams12_a[5] = 
  {0.0 / 1.0,
  5226845396890227.0 / 36028797018963968.0,
  5602232597459763.0 / 1152921504606846976.0,
  4920783552645014.0 / 4611686018427387904.0,
  5570041179715558.0 / 9223372036854775808.0 };
	 
  const real q_gams12_b[5] = 
  {-6799658092196962.0 / 18014398509481984.0,
  -4810318281953418.0 / 36028797018963968.0,
  -8340064868499919.0 / 576460752303423488.0,
  -8375346434191481.0 / 576460752303423488.0,
  -6231299816839781.0 / 1152921504606846976.0 };
  
  // Gamma(x): S13=[137.5,138.5], x0 = 138.0390625;
  const real q_gams13_a[5] = 
  {0.0 / 1.0,
  5077602289066213.0 / 18014398509481984.0,
  4971389438544026.0 / 576460752303423488.0,
  8326326060413143.0 / 4611686018427387904.0,
  7582214946963028.0 / 9223372036854775808.0 };
	 
  const real q_gams13_b[5] = 
  {-8336661118182000.0 / 72057594037927936.0,
  -6152827173580884.0 / 36028797018963968.0,
  -6306415240260295.0 / 576460752303423488.0,
  -5949724688158565.0 / 576460752303423488.0,
  -5625482009313801.0 / 576460752303423488.0 };
  
  // Gamma(x): S14=[146.5,147.5], x0 = 146.94921875;
  const real q_gams14_a[6] = 
  {0.0 / 1.0,
  8624518348171320.0 / 18014398509481984.0,
  7049908374954842.0 / 576460752303423488.0,
  5769548496333642.0 / 2305843009213693952.0,
  5110826386599631.0 / 4611686018427387904.0,
  5902099580937297.0 / 9223372036854775808.0 };
	 
  const real q_gams14_b[6] = 
  {4591842870321392.0 / 18014398509481984.0,
  -7195103997153274.0 / 36028797018963968.0,
  -5062133254875404.0 / 576460752303423488.0,
  -4899776496053167.0 / 576460752303423488.0,
  -4703735605739871.0 / 576460752303423488.0,
  -8717622510435264.0 / 1152921504606846976.0 };
  
  // Gamma(x): Stuetzinterval S15=[153.5,154.5], x0 = 153.90234375;
  const real q_gams15_a[6] = 
  {0.0 / 1.0,
  5047093325633351.0 / 9007199254740992.0,
  8838568428828895.0 / 576460752303423488.0,
  7169875770591066.0 / 2305843009213693952.0,
  6277116142443993.0 / 4611686018427387904.0,
  7160380168510556.0 / 9223372036854775808.0 };
	 
  const real q_gams15_b[6] = 
  {5575895624898616.0 / 18014398509481984.0,
  -7982725134478240.0 / 36028797018963968.0,
  -8683131154525203.0 / 1152921504606846976.0,
  -8504376590080638.0 / 1152921504606846976.0,
  -8272396567408614.0 / 1152921504606846976.0,
  -7370684458128734.0 / 1152921504606846976.0 }; 
   
  // Gamma(x): S16=[160.5,161.5], x0 = 161.08984375;
  const real q_gams16_a[6] = 
  {0.0 / 1.0,
  8928052213955455.0 / 18014398509481984.0,
  5405751921217612.0 / 288230376151711744.0,
  8725861134087760.0 / 2305843009213693952.0,
  7583718370137025.0 / 4611686018427387904.0,
  8577302812898546.0 / 9223372036854775808.0 };
	 
  const real q_gams16_b[6] = 
  {6669445708872192.0 / 144115188075855872.0,
  -8769965967955460.0 / 36028797018963968.0,
  -7524042392635917.0 / 1152921504606846976.0,
  -7422547685779587.0 / 1152921504606846976.0,
  -7284571995868819.0 / 1152921504606846976.0,
  -7797018373646746.0 / 1152921504606846976.0 }; 
  
  // Gamma(x): S17=[167.5,168.5], x0 = 168.0;
  const real q_gams17_a[6] = 
  {0.0 / 1.0,
  4642907317603488.0 / 9007199254740992.0,
  6403634366308572.0 / 288230376151711744.0,
  5153510478021874.0 / 1152921504606846976.0,
  8919203526515161.0 / 4611686018427387904.0,
  5015558226948384.0 / 4611686018427387904.0 };
	 
  const real q_gams17_b[6] = 
  {-6229620043098112.0 / 9223372036854775808.0,
  -4750296278401558.0 / 18014398509481984.0,
  -6639785613322219.0 / 1152921504606846976.0,
  -6577910457226093.0 / 1152921504606846976.0,
  -6491061443859352.0 / 1152921504606846976.0,
  -6372399542597081.0 / 1152921504606846976.0 }; 									 

// ****************************************************************************
// ******************** Gamma(x), 1/Gamma(x) **********************************
// ****************************************************************************

 real gam_S0(const real& x)
// Calculating approximations for 1/Gamma(x) in S0 = [1.5,2.5];
// Rel. error bound:  3.930138e-16;
// Blomquist, 19.05.2009;
{
	real y(0),v;

	 // Continued fraction:  K_7(v), v = 1/(x-x0), x0=2;
	if (x==2)
		y = q_gams0_b[0];
	else
	{
		v = 1/(x-2);
	 	
		y = q_gams0_a[7] / (  v + q_gams0_b[7]);
		y = q_gams0_a[6] / ( (v + q_gams0_b[6]) + y );
		y = q_gams0_a[5] / ( (v + q_gams0_b[5]) + y );	
		y = q_gams0_a[4] / ( (v + q_gams0_b[4]) + y );
		y = q_gams0_a[3] / ( (v + q_gams0_b[3]) + y );
		y = q_gams0_a[2] / ( (v + q_gams0_b[2]) + y );
		y = q_gams0_a[1] / ( (v + q_gams0_b[1]) + y ) + q_gams0_b[0];
	}	
	return y;
}
 
 real gam_S0_n0(const real& x)
// Calculating approximations for K_7(v),  v=1/x,  x in [-0.5,+0.5];
// Rel. error bound:  3.930138e-16;
// Blomquist, 19.05.2009;
{
	real y(0),v;

	 // Continued fraction:  K_7(v), v = 1/x, x0=0;
	if (x==0)
		y = q_gams0_b[0];
	else
	{
		v = 1/x;
		 
		y = q_gams0_a[7] / (  v + q_gams0_b[7]);
		y = q_gams0_a[6] / ( (v + q_gams0_b[6]) + y );
		y = q_gams0_a[5] / ( (v + q_gams0_b[5]) + y );	
		y = q_gams0_a[4] / ( (v + q_gams0_b[4]) + y );
		y = q_gams0_a[3] / ( (v + q_gams0_b[3]) + y );
		y = q_gams0_a[2] / ( (v + q_gams0_b[2]) + y );
		y = q_gams0_a[1] / ( (v + q_gams0_b[1]) + y ) + q_gams0_b[0];
	}	
		
	return y;
}

real gam_S0_n1(const real& x)
// Calculating approximations for K_7(v),  v=1/(x-1),  x in [0.5,1.5];
// Rel. error bound:  3.930138e-16;
// Blomquist, 21.05.2009;
{
	real y(0),v;

	 // Continued fraction:  K_7(v), v = 1/(x-1), x0=1;
	if (x==1)
		y = q_gams0_b[0];
	else
	{
		v = 1/(x-1);
		 
		y = q_gams0_a[7] / (  v + q_gams0_b[7]);
		y = q_gams0_a[6] / ( (v + q_gams0_b[6]) + y );
		y = q_gams0_a[5] / ( (v + q_gams0_b[5]) + y );	
		y = q_gams0_a[4] / ( (v + q_gams0_b[4]) + y );
		y = q_gams0_a[3] / ( (v + q_gams0_b[3]) + y );
		y = q_gams0_a[2] / ( (v + q_gams0_b[2]) + y );
		y = q_gams0_a[1] / ( (v + q_gams0_b[1]) + y ) + q_gams0_b[0];
	}	
		
	return y;
}

real gammar_S0(const real& x)
// Calculating approximations for 1/Gamma(x), x in (-0.5,+8.5];
// eps = 1.725284e-15;
// Blomquist, 21.05.2009;
{
	real y,Ne;
	int n,p;
	
	n = Round(x);
	
	switch(n)
	{
		case 0: if (expo(x)<=-52) y = x;
		else y = x*(x+1)*gam_S0_n0(x); break; // x in (-0.5,+0.5);
		case 1: y = x*gam_S0_n1(x);            break; // x in [0.5,1.5);
		case 2: y = gam_S0(x);                 break; // x in [1.5,2.5);

		default: // n=3,4,...,8; x in [2.5,8.5];
			p = n-2; Ne = x-1;
			for (int k=2; k<=p; k++) Ne *= (x-k);
			y = gam_S0(x-p) / Ne;
	}
	
	return y;
}

real gamma_S0(const real& x)
// Calculating approximations for Gamma(x) in S0 = (-0.5,+8.5];
// eps = 1.947330E-15;
// Blomquist, 21.05.2009;
{
	real y = 1/gammar_S0(x);
	return y;
}

real gam_S1(const real& x)
// Calculating approximations for Gamma(x) in  S1 = [10.5,11.5];
// Rel. error bound:  7.696082e-16;
// Blomquist, 22.05.2009;		
{
	real y(0),v;

	// Continued fraction:  K_6(v), v = 1/(x-x0), x0=11.125;
	if (x==11.125)
		y = q_gams1_b[0];
	else
	{
		v = 1/(x-11.125);
	 	
		y = q_gams1_a[6] / (  v + q_gams1_b[6]);
		y = q_gams1_a[5] / ( (v + q_gams1_b[5]) + y );	
		y = q_gams1_a[4] / ( (v + q_gams1_b[4]) + y );
		y = q_gams1_a[3] / ( (v + q_gams1_b[3]) + y );
		y = q_gams1_a[2] / ( (v + q_gams1_b[2]) + y );
		y = q_gams1_a[1] / ( (v + q_gams1_b[1]) + y ) + q_gams1_b[0];
	}	
	y += 1;
	y *= q_ex10(x);
	times2pown(y,-15);
		
	return y;
}	
	
real gamma_S1(const real& x)
// Calculating approximations for Gamma(x), x in [8.5,16.5];
// eps = 1.879833e-15;
// Blomquist, 24.05.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-11;

	if (n>11) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S1(x-p);
	}
	else // neighbour intervals to the left and S1 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S1(x) : gam_S1(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S2(const real& x)
// Calculating approximations for Gamma(x) in S2 = [18.5,19.5];
// Rel. error bound:  7.411454e-16;
// Blomquist, 25.05.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=18.96875;
	if (x==18.96875)
		y = q_gams2_b[0];
	else
	{
		v = 1/(x-18.96875);

		y = q_gams2_a[5] / (  v + q_gams2_b[5]);
		y = q_gams2_a[4] / ( (v + q_gams2_b[4]) + y );
		y = q_gams2_a[3] / ( (v + q_gams2_b[3]) + y );
		y = q_gams2_a[2] / ( (v + q_gams2_b[2]) + y );
		y = q_gams2_a[1] / ( (v + q_gams2_b[1]) + y ) + q_gams2_b[0];
	}	
	y += 1;
	y *= q_exp2(4*x);
	times2pown(y,-23);
		
	return y;
}	

real gamma_S2(const real& x)
// Calculating approximations for Gamma(x), x in [16.5,24.5];
// eps = 1.851370e-15;
// Blomquist, 25.05.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-19;

	if (n>19) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S2(x-p);
	}
	else // neighbour intervals to the left and S2 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S2(x) : gam_S2(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S3(const real& x)
// Calculating approximations for Gamma(x) in  S3 = [29.5,30.5];
// Rel. error bound:  9.724029e-16;
// Blomquist, 26.05.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=29.9375;
	if (x==29.9375)
		y = q_gams3_b[0];
	else
	{
		v = 1/(x-29.9375);

		y = q_gams3_a[5] / (  v + q_gams3_b[5]);
		y = q_gams3_a[4] / ( (v + q_gams3_b[4]) + y );
		y = q_gams3_a[3] / ( (v + q_gams3_b[3]) + y );
		y = q_gams3_a[2] / ( (v + q_gams3_b[2]) + y );
		y = q_gams3_a[1] / ( (v + q_gams3_b[1]) + y ) + q_gams3_b[0];
	}	
	y += 1;
	y *= q_exp2(4*x);
	times2pown(y,-17);
		
	return y;
}	

real gamma_S3(const real& x)
// Calculating approximations for Gamma(x), x in [24.5,35.5];
// eps = 2.082627e-15;
// Blomquist, 26.05.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-30;

	if (n>30) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S3(x-p);
	}
	else // neighbour intervals to the left and S3 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S3(x) : gam_S3(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S4(const real& x)
// Calculating approximations for Gamma(x) in  S4 = [40.5,41.5];
// Rel. error bound:  8.170628e-16;
// Blomquist, 26.05.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=41.140625;
	if (x==41.140625)
		y = q_gams4_b[0];
	else
	{
		v = 1/(x-41.140625);

		y = q_gams4_a[5] / (  v + q_gams4_b[5]);
		y = q_gams4_a[4] / ( (v + q_gams4_b[4]) + y );
		y = q_gams4_a[3] / ( (v + q_gams4_b[3]) + y );
		y = q_gams4_a[2] / ( (v + q_gams4_b[2]) + y );
		y = q_gams4_a[1] / ( (v + q_gams4_b[1]) + y ) + q_gams4_b[0];
	}	
	y += 1;
	y *= exp(4*x);
	times2pown(y,-78);
		
	return y;
}	

real gamma_S4(const real& x)
// Calculating approximations for Gamma(x), x in [35.5,46.5];
// eps = 1.927286e-15;
// Blomquist, 28.05.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-41;

	if (n>41) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S4(x-p);
	}
	else // neighbour intervals to the left and S4 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S4(x) : gam_S4(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S5(const real& x)
// Calculating approximations for Gamma(x) in S5 = [51.5,52.5];
// Rel. error bound:  6.351369e-16;
// Blomquist, 28.05.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=52.015625;
	if (x==52.015625)
		y = q_gams5_b[0];
	else
	{
		v = 1/(x-52.015625);

		y = q_gams5_a[5] / (  v + q_gams5_b[5]);
		y = q_gams5_a[4] / ( (v + q_gams5_b[4]) + y );
		y = q_gams5_a[3] / ( (v + q_gams5_b[3]) + y );
		y = q_gams5_a[2] / ( (v + q_gams5_b[2]) + y );
		y = q_gams5_a[1] / ( (v + q_gams5_b[1]) + y ) + q_gams5_b[0];
	}	
	y += 1;
	y *= exp(4*x);
	times2pown(y,-80);
		
	return y;
}	

real gamma_S5(const real& x)
// Calculating approximations for Gamma(x), x in [46.5,57.5];
// eps = 1.745360e-15;
// Blomquist, 29.05.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-52;

	if (n>52) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S5(x-p);
	}
	else // neighbour intervals to the left and S5 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S5(x) : gam_S5(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S6(const real& x)
// Calculating approximations for Gamma(x) in S6 = [62.5,63.5];
// Rel. error bound:  8.310104e-16;
// Blomquist, 29.05.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=63.015625;
	if (x==63.015625)
		y = q_gams6_b[0];
	else
	{
		v = 1/(x-63.015625);

		y = q_gams6_a[5] / (  v + q_gams6_b[5]);
		y = q_gams6_a[4] / ( (v + q_gams6_b[4]) + y );
		y = q_gams6_a[3] / ( (v + q_gams6_b[3]) + y );
		y = q_gams6_a[2] / ( (v + q_gams6_b[2]) + y );
		y = q_gams6_a[1] / ( (v + q_gams6_b[1]) + y ) + q_gams6_b[0];
	}	
	y += 1;
	y *= exp(4*x);
	times2pown(y,-80);
		
	return y;
}	

real gamma_S6(const real& x)
// Calculating approximations for Gamma(x), x in [57.5,68.5];
// eps = 1.941234e-15;
// Blomquist, 06.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-63;

	if (n>63) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S6(x-p);
	}
	else // neighbour intervals to the left and S6 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S6(x) : gam_S6(x+p)/Pr;		
	}
	
	return y;
}	


real gam_S7(const real& x)
// Calculating approximations for Gamma(x) in S7 = [73.5,74.5];
// Rel. error bound:  7.999142e-16;
// Blomquist, 05.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=74.16015625;
	if (x==74.16015625)
		y = q_gams7_b[0];
	else
	{
		v = 1/(x-74.16015625);

		y = q_gams7_a[5] / (  v + q_gams7_b[5]);
		y = q_gams7_a[4] / ( (v + q_gams7_b[4]) + y );
		y = q_gams7_a[3] / ( (v + q_gams7_b[3]) + y );
		y = q_gams7_a[2] / ( (v + q_gams7_b[2]) + y );
		y = q_gams7_a[1] / ( (v + q_gams7_b[1]) + y ) + q_gams7_b[0];
	}	
	y += 1;
	y *= exp(4*x);
	times2pown(y,-76);
		
	return y;
}	

real gamma_S7(const real& x)
// Calculating approximations for Gamma(x), x in [68.5,79.5];
// eps = 1.910138e-15;
// Blomquist, 08.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-74;

	if (n>74) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S7(x-p);
	}
	else // neighbour intervals to the left and S7 itself 
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S7(x) : gam_S7(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S8(const real& x)
// Calculating approximations for Gamma(x) in  S8 = [84.5,85.5];
// Rel. error bound:  7.192334e-16;
// Blomquist, 08.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_4(v), v = 1/(x-x0), x0=85.1015625;
	if (x==85.1015625)
		y = q_gams8_b[0];
	else
	{
		v = 1/(x-85.1015625);

		y = q_gams8_a[4] / (  v + q_gams8_b[4]);
		y = q_gams8_a[3] / ( (v + q_gams8_b[3]) + y );
		y = q_gams8_a[2] / ( (v + q_gams8_b[2]) + y );
		y = q_gams8_a[1] / ( (v + q_gams8_b[1]) + y ) + q_gams8_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-144);
		
	return y;
}	

real gamma_S8(const real& x)
// Calculating approximations for Gamma(x), x in [79.5,90.5];
// eps = 1.829457e-15;
// Blomquist, 09.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-85;

	if (n>85) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S8(x-p);
	}
	else // neighbour intervals to the left and S8 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S8(x) : gam_S8(x+p)/Pr;		
	}
	
	return y;
}	


real gam_S9(const real& x)
// Calculating approximations for Gamma(x) in S9 = [95.5,96.5];
// Rel. error bound:  6.622145e-16;
// Blomquist, 09.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_4(v), v = 1/(x-x0), x0=95.984375;
	if (x==95.984375)
		y = q_gams9_b[0];
	else
	{
		v = 1/(x-95.984375);

		y = q_gams9_a[4] / (  v + q_gams9_b[4]);
		y = q_gams9_a[3] / ( (v + q_gams9_b[3]) + y );
		y = q_gams9_a[2] / ( (v + q_gams9_b[2]) + y );
		y = q_gams9_a[1] / ( (v + q_gams9_b[1]) + y ) + q_gams9_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-146);
		
	return y;
}	

real gamma_S9(const real& x)
// Calculating approximations for Gamma(x), x in [90.5,101.5];
// eps = 1.772438e-15;
// Blomquist, 10.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-96;

	if (n>96) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S9(x-p);
	}
	else // neighbour intervals to the left and S9 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S9(x) : gam_S9(x+p)/Pr;		
	}
	
	return y;
}	


real gam_S10(const real& x)
// Calculating approximations for Gamma(x) in S10 = [106.5,107.5];
// Rel. error bound:  8.254965e-16;
// Blomquist, 10.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_4(v), v = 1/(x-x0), x0=107.078125;
	if (x==107.078125)
		y = q_gams10_b[0];
	else
	{
		v = 1/(x-107.078125);

		y = q_gams10_a[4] / (  v + q_gams10_b[4]);
		y = q_gams10_a[3] / ( (v + q_gams10_b[3]) + y );
		y = q_gams10_a[2] / ( (v + q_gams10_b[2]) + y );
		y = q_gams10_a[1] / ( (v + q_gams10_b[1]) + y ) + q_gams10_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-146);
		
	return y;
}	

real gamma_S10(const real& x)
// Calculating approximations for Gamma(x), x in [101.5,112.5];
// eps = 1.935720e-15;
// Blomquist, 12.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-107;

	if (n>107) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S10(x-p);
	}
	else // neighbour intervals to the left and S10 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S10(x) : gam_S10(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S11(const real& x)
// Calculating approximations for Gamma(x) in S11 = [117.5,118.5];
// Rel. error bound:  6.766734e-16;
// Blomquist, 11.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_4(v), v = 1/(x-x0), x0=117.8671875;
	if (x==117.8671875)
		y = q_gams11_b[0];
	else
	{
		v = 1/(x-117.8671875);

		y = q_gams11_a[4] / (  v + q_gams11_b[4]);
		y = q_gams11_a[3] / ( (v + q_gams11_b[3]) + y );
		y = q_gams11_a[2] / ( (v + q_gams11_b[2]) + y );
		y = q_gams11_a[1] / ( (v + q_gams11_b[1]) + y ) + q_gams11_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-144);
		
	return y;
}	

real gamma_S11(const real& x)
// Calculating approximations for Gamma(x) in [112.5,122.5];
// eps = 1.786897e-15;
// Blomquist, 14.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-118;
	
	if (n>118) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S11(x-p);
	}
	else // left neighbour intervals and S11 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S11(x) : gam_S11(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S12(const real& x)
// Calculating approximations for Gamma(x) in S12 = [126.5,127.5];
// Rel. error bound:  8.175777e-16;
// Blomquist, 13.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_4(v), v = 1/(x-x0), x0=126.7421875;
	if (x==126.7421875)
		y = q_gams12_b[0];
	else
	{
		v = 1/(x-126.7421875);

		y = q_gams12_a[4] / (  v + q_gams12_b[4]);
		y = q_gams12_a[3] / ( (v + q_gams12_b[3]) + y );
		y = q_gams12_a[2] / ( (v + q_gams12_b[2]) + y );
		y = q_gams12_a[1] / ( (v + q_gams12_b[1]) + y ) + q_gams12_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-141);
		
	return y;
}	

real gamma_S12(const real& x)
// Calculating approximations for Gamma(x) in [122.5,132.5];
// eps = 1.927801e-15;
// Blomquist, 15.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-127;
	if (n>127) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S12(x-p);
	}
	else // left neighbour intervals and S12 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S12(x) : gam_S12(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S13(const real& x)
// Calculating approximations for Gamma(x) in S13 = [137.5,138.5];
// Rel. error bound:  8.563188e-16;
// Blomquist, 14.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_4(v), v = 1/(x-x0), x0=138.0390625;
	if (x==138.0390625)
		y = q_gams13_b[0];
	else
	{
		v = 1/(x-138.0390625);

		y = q_gams13_a[4] / (  v + q_gams13_b[4]);
		y = q_gams13_a[3] / ( (v + q_gams13_b[3]) + y );
		y = q_gams13_a[2] / ( (v + q_gams13_b[2]) + y );
		y = q_gams13_a[1] / ( (v + q_gams13_b[1]) + y ) + q_gams13_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-137);
		
	return y;
}	

real gamma_S13(const real& x)
// Calculating approximations for Gamma(x) in [132.5,142.5];
// eps = 1.966542e-15;
// Blomquist, 15.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-138;
	if (n>138) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S13(x-p);
	}
	else // left neighbour intervals and S13 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S13(x) : gam_S13(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S14(const real& x)
// Calculating approximations for Gamma(x) in S14 = [146.5,147.5];
// Rel. error bound:  8.570174e-16;
// Blomquist, 15.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=146.94921875;
	if (x==146.94921875)
		y = q_gams14_b[0];
	else
	{
		v = 1/(x-146.94921875);
		
		y = q_gams14_a[5] / (  v + q_gams14_b[5]);
		y = q_gams14_a[4] / ( (v + q_gams14_b[4]) + y );
		y = q_gams14_a[3] / ( (v + q_gams14_b[3]) + y );
		y = q_gams14_a[2] / ( (v + q_gams14_b[2]) + y );
		y = q_gams14_a[1] / ( (v + q_gams14_b[1]) + y ) + q_gams14_b[0];
	}	
	y += 1;
	y *= q_ex10(2*x);
	times2pown(y,-133);
		
	return y;
}	

real gamma_S14(const real& x)
// Calculating approximations for Gamma(x) in [142.5,150.5];
// eps = 1.745196e-15;
// Blomquist, 16.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-147;
	if (n>147) // neighbour intervals to the right:
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S14(x-p);
	}
	else // left neighbour intervals and S14 itself:
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S14(x) : gam_S14(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S15(const real& x)
// Calculating approximations for Gamma(x) in S15 = [153.5,154.5];
// Rel. error bound:  1.279047e-15;
// Blomquist, 16.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=153.90234375;
	if (x==153.90234375)
		y = q_gams15_b[0];
	else
	{
		v = 1/(x-153.90234375);
		
		y = q_gams15_a[5] / (  v + q_gams15_b[5]);
		y = q_gams15_a[4] / ( (v + q_gams15_b[4]) + y );
		y = q_gams15_a[3] / ( (v + q_gams15_b[3]) + y );
		y = q_gams15_a[2] / ( (v + q_gams15_b[2]) + y );
		y = q_gams15_a[1] / ( (v + q_gams15_b[1]) + y ) + q_gams15_b[0];
	}	
	y += 1;
	v = q_ex10(x);
	times2pown(y,-129);
	y *= v;
	y *= v;

	return y;
}

real gamma_S15(const real& x)
// Calculating approximations for Gamma(x) in [150.5,157.5];
// eps = 1.945181e-15;
// Blomquist, 19.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-154;
	if (n>154) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S15(x-p);
	}
	else // left neighbour intervals and S15 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S15(x) : gam_S15(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S16(const real& x)
// Calculating approximations for Gamma(x) in S16 = [160.5,161.5];
// Rel. error bound:  1.381209e-15;
// Blomquist, 19.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=161.08984375;
	if (x==161.08984375)
		y = q_gams16_b[0];
	else
	{
		v = 1/(x-161.08984375);
		
		y = q_gams16_a[5] / (  v + q_gams16_b[5]);
		y = q_gams16_a[4] / ( (v + q_gams16_b[4]) + y );
		y = q_gams16_a[3] / ( (v + q_gams16_b[3]) + y );
		y = q_gams16_a[2] / ( (v + q_gams16_b[2]) + y );
		y = q_gams16_a[1] / ( (v + q_gams16_b[1]) + y ) + q_gams16_b[0];
	}	
	y += 1;
	v = q_ex10(x);
	times2pown(v,-62);
	y *= sqr(v);
	return y;
}

real gamma_S16(const real& x)
// Calculating approximations for Gamma(x) in [157.5,164.5];
// eps = 2.047343e-15;
// Blomquist, 20.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-161;
	if (n>161) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S16(x-p);
	}
	else  // left neighbour intervals and S16 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S16(x) : gam_S16(x+p)/Pr;		
	}
	
	return y;
}	

real gam_S17(const real& x)
// Calculating approximations for Gamma(x) in:  S17 = [167.5,168.5];
// Rel. error bound:  1.390239e-15;
// Blomquist, 20.06.2009;		
{
	real y(0),v;

	// Continued fraction:  K_5(v), v = 1/(x-x0), x0=168.0;
	if (x==168.0)
		y = q_gams17_b[0];
	else
	{
		v = 1/(x-168.0);
		
		y = q_gams17_a[5] / (  v + q_gams17_b[5]);
		y = q_gams17_a[4] / ( (v + q_gams17_b[4]) + y );
		y = q_gams17_a[3] / ( (v + q_gams17_b[3]) + y );
		y = q_gams17_a[2] / ( (v + q_gams17_b[2]) + y );
		y = q_gams17_a[1] / ( (v + q_gams17_b[1]) + y ) + q_gams17_b[0];
	}	
	y += 1;
	v = q_ex10(x);
	times2pown(y,-119);
	y *= v;
	y *= v;	
	return y;
}

real gamma_S17(const real& x)
// Calculating approximations for Gamma(x), x in [164.5,171.5];
// eps = 2.056373e-15;
// Blomquist, 20.06.2009;
{
	real y,Pr;
	int n,p;
	
	n = Round(x);
	p = n-168;
	if (n>168) // neighbour intervals to the right
	{
		Pr = x-1;
		for (int k=2; k<=p; k++)
			Pr *= x-k;
		y = Pr*gam_S17(x-p);
	}
	else // left neighbour intervals and S17 itself
	{
		p = -p;  Pr = x;
		for (int k=1; k<=p-1; k++) 
			Pr *= x+k;
		y = (p==0)? gam_S17(x) : gam_S17(x+p)/Pr;		
	}
	
	return y;
}	

real gamma_05(const real& x)
// Calculating Gamma(x), x in (-0.5,+171.5]
// eps = 2.082627e-15;
// Blomquist, 21.06.2009;
{
	real y(0);
	int Nr;
	
	Nr = int_no(gam_f85,19,x);
	
	switch(Nr)
	{
		case 0:  y = gamma_S0(x);   break; // x in (-0.5,8.5);
		case 1:  y = gamma_S1(x);   break; // x in [8.5,16.5);
		case 2:  y = gamma_S2(x);   break; // x in [16.5,24.5);
		case 3:  y = gamma_S3(x);   break; // x in [24.5,35.5);
		case 4:  y = gamma_S4(x);   break; // x in [35.5,46.5);
		case 5:  y = gamma_S5(x);   break; // x in [46.5,57.5);
		case 6:  y = gamma_S6(x);   break; // x in [57.5,68.5);
		case 7:  y = gamma_S7(x);   break; // x in [68.5,79.5);
		case 8:  y = gamma_S8(x);   break; // x in [79.5,90.5);
		case 9:  y = gamma_S9(x);   break; // x in [90.5,101.5);
		case 10: y = gamma_S10(x);  break; // x in [101.5,112.5);
		case 11: y = gamma_S11(x);  break; // x in [112.5,122.5);
		case 12: y = gamma_S12(x);  break; // x in [122.5,132.5);
		case 13: y = gamma_S13(x);  break; // x in [132.5,142.5);
		case 14: y = gamma_S14(x);  break; // x in [142.5,150.5);
		case 15: y = gamma_S15(x);  break; // x in [150.5,157.5);
		case 16: y = gamma_S16(x);  break; // x in [157.5,164.5);
		case 17: y = gamma_S17(x);  break; // x in [164.5,171.5);

		default: // n=18; x in [171.5, ...];
			y = gamma_S17(x);
	}
	
	return y;
}

real gammar(const real& x)
// Calculating 1/Gamma(x) for x in [-170.0,+171.0];
// eps = 2.866906e-15;
// Blomquist, 21.06.2009;
{
	real y;
	
	if (x<-170.0 || x>171.0) 
		cxscthrow(STD_FKT_OUT_OF_DEF("real gammar(const real& x)"));
	
	if (x<=-0.5)
		y = -sinpix_pi(x)*x*gamma_05(-x);
	else
		if (x<=8.5)
			y = gammar_S0(x);
	else y = 1/gamma_05(x);
	
	return y;
}

real gamma(const real& x)
// Calculating Gamma(x) for x in [-170.0,+171.5]
// eps = 3.088951e-15;
// Blomquist, 21.06.2009;
{
	real y;
	
	if (x>171.5 || x<-170.0) 
		cxscthrow(STD_FKT_OUT_OF_DEF("real gamma(const real& x)"));
	
	if (x<=-0.5)
		y = -1/( sinpix_pi(x)*x*gamma_05(-x) );
	else
		y = gamma_05(x);
	
	return y;
}


extern "C" {
   void r_lfsr(void) {;} // Siehe real.hpp in real_ari...?!?!
}

} // namespace cxsc

