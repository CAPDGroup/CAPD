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
/* CVS $Id: l_complex.cpp,v 1.21 2014/01/30 17:23:46 cxsc Exp $ */

#include "l_complex.hpp"  // corresponding header file
#include "l_interval.hpp"

namespace cxsc {

// ---------------- Unary Operators -----------------------------------------
l_complex operator-(const l_complex& x)
{
    return l_complex(-x.re, -x.im);
}

l_complex operator+(const l_complex& x)
{
    return x;
}

// ----------------- cdotprecision +/- l_complex ----------------------------

cdotprecision operator+(const l_complex& lc, const cdotprecision& cd) throw() 
{ 
   return _cdotprecision(lc) + cd; 
}

cdotprecision operator+(const cdotprecision& cd, const l_complex& lc) throw() 
{ 
   return _cdotprecision(lc) + cd; 
}   

cdotprecision operator-(const l_complex& lc, const cdotprecision& cd) throw() 
{ 
   return _cdotprecision(lc) - cd; 
}

cdotprecision operator-(const cdotprecision& cd, const l_complex& lc) throw() 
{ 
   return cd - _cdotprecision(lc); 
}

// ------------------ dotprecision +/- l_complex ----------------------------

cdotprecision operator+(const l_complex& lc, const dotprecision& cd) throw() 
{ 
   return _cdotprecision(lc) + cd; 
}

cdotprecision operator+(const dotprecision& cd, const l_complex& lc) throw() 
{ 
   return _cdotprecision(lc) + cd; 
}   

cdotprecision operator-(const l_complex& lc, const dotprecision& cd) throw() 
{ 
   return _cdotprecision(lc) - cd; 
}

cdotprecision operator-(const dotprecision& cd, const l_complex& lc) throw() 
{ 
   return cd - _cdotprecision(lc); 
}

// ---------------- l_complex * l_complex -----------------------------------
l_complex operator * (const l_complex& a, const l_complex& b)
    throw()
{
    l_real r, i; 
    dotprecision dot(0.0);

    accumulate(dot,a.re,b.re);
    accumulate(dot,-a.im,b.im);
    r = dot;

    dot = 0.0;
    accumulate(dot,a.im,b.re);
    accumulate(dot,a.re,b.im);
    i = dot;

    return( l_complex(r,i) );            
}

// ------------------ l_complex * complex -----------------------------------
l_complex operator * (const l_complex& a, const complex& b)
    throw()
{
    l_real r, i; 
    dotprecision dot(0.0);
    
    accumulate(dot,a.re,Re(b));
    accumulate(dot,-a.im,Im(b));
    r = dot;  // r is real-part

    dot = 0.0;
    accumulate(dot,a.im,Re(b));
    accumulate(dot,a.re,Im(b));
    i = dot;  // i is imaginary part

    return( l_complex(r,i) );
}

l_complex operator * ( const complex& b, const l_complex& a )
    throw()
{
    l_real r, i; 
    dotprecision dot(0.0);

    accumulate(dot,a.re,Re(b));
    accumulate(dot,-a.im,Im(b));
    r = dot;  // r is real part

    dot = 0.0;
    accumulate(dot,a.im,Re(b));
    accumulate(dot,a.re,Im(b));
    i = dot;  // i is imaginary part

    return( l_complex(r,i) );
}

// ---------------- l_complex / l_complex -----------------------------------

static const int maxexpo  = 1020;

void skale_down_exp(int ex1, int ex2, int D, int& d1, int& d2)
// l_interval x1,x2 sind Punktintervalle mit den Zweierexponenten
// ex1,ex2 > -maxint, d.h. x1,x2 <> 0.0;
// x1*x2 ist mit 2^D, D<0, so zu skalieren, dass fuer das Produkt
// x1*x2 moeglichst wenig Stellen verloren gehen! x1 ist dazu mit
// 2^d1 und x2 ist mit 2^d2 zu multiplizieren. Es muss also gelten:
// d1+d2 = D<0;
// Diese Funktion berechnet zu den gegebenen Zweierexponenten
// ex1 und ex2 die Exponenten d1 und d2.
// ex1 berechnet sich z.B. durch:  ex1 = expo(x1).
// Blomquist, 25.10.2006;
  
{
    bool change(false);
    int Diff, D1, c;
    d1 = 0;  d2 = 0;

    if (D<0)
    {
	if (ex2>ex1)
	{ c = ex1;  ex1 = ex2;  ex2 = c;  change = true; }
        // ab jetzt gilt:  ex1 >= ex2;
	c = ex1 + D; 
	if (c<ex2)
	{
	    Diff = ex2 - c;   D1 = Diff/2;
	    d1 = D + D1;      d2 = D - d1; // d1+d2 == D<0;
	}
	else d1 = D;  // d2 = 0, s.o.
	if (change)
	{ c = d1;  d1 = d2;  d2 = c; }
    }

} // skale_down_exp(...)

void skale_up_exp1(int ex1, int ex2, int& fillin, int& d1, int& d2)
// l_interval x1,x2 sind Punktintervalle mit den Zweierexponenten
// ex1,ex2 > -maxint, d.h. x1,x2 <> 0.0;
// Das Maximum (ex1+ex2) der beiden gegebenen Exponentensummen
// ist <= 1022; Die Punktintervalle x1,x2 sind so mit 2^d1 bzw. mit
// 2^d2 zu multiplizieren, dass gilt.
//     I.   2*|(x1*x2)| < MaxReal;
//     II.  Bei dem Produkt x1*x2 duerfen moeglichst wenig Stellen 
//          verloren gehen.
// Blomquist, 25.10.2006;
{
    bool change(false);
    int c, pot2;
    d1 = 0;  d2 = 0;
    fillin = maxexpo - (ex1+ex2);
    if (fillin>0)
    {
	if (ex2>ex1)
	{ c = ex1;  ex1 = ex2;  ex2 = c;  change = true; }
        // ab jetzt gilt:  ex1 >= ex2;
	pot2 = maxexpo - ex2;  
        // Um pot2 kann x2 ohne Overflow hochskaliert werden.
	if (fillin <= pot2) d2 = fillin;
	else { d2 = pot2;  d1 = fillin - pot2; }
	if (change)
	{ c = d1;  d1 = d2;  d2 = c; }
    }
} // skale_up_exp1(...)

void skale_up_exp2(int ex1, int ex2, int fillin, int& d1, int& d2)
// l_interval x1,x2 sind Punktintervalle mit den Zweierexponenten
// ex1,ex2 > -maxint, d.h. x1,x2 <> 0.0;
// Das Minimum (ex1+ex2) der beiden gegebenen Exponentensummen
// ist <= 1022; Die Punktintervalle x1,x2 sind so mit 2^d1 bzw. mit
// 2^d2 zu multiplizieren, dass gilt.
//     I.   2*|(x1*x2)| < MaxReal;
//     II.  Bei dem Produkt x1*x2 duerfen moeglichst wenig Stellen 
//          verloren gehen.
// Blomquist, 25.10.2006;
{
    bool change(false);
    int c, pot2;
    d1 = 0;  d2 = 0;
    if (fillin>0)
    {
	if (ex2>ex1)
	{ c = ex1;  ex1 = ex2;  ex2 = c;  change = true; }
        // ab jetzt gilt:  ex1 >= ex2;
	pot2 = 1022 - ex2; // 1022 ??
        // Um pot2 kann x2 ohne Overflow hochskaliert werden.
	if (fillin <= pot2) d2 = fillin;
	else { d2 = pot2;  d1 = fillin - pot2; }
	if (change)
	{ c = d1;  d1 = d2;  d2 = c; }
    }
} // skale_up_exp2(...)

void product(const l_real& a, const l_real& b, const l_real& c, 
              const l_real& d, int& ex, l_interval& res)
// Calulation of an inclusion of: a*b + c*d
{
    l_real a1(a), b1(b), c1(c), d1(d);

    a1 += 0.0;   b1 += 0.0;   c1 += 0.0;   d1 += 0.0;
    int ex_a1( expo(a1[1]) ), ex_b1( expo(b1[1]) ),
	ex_c1( expo(c1[1]) ), ex_d1( expo(d1[1]) ), m,p,D1,D2;
    l_interval a_i(a1),b_i(b1),c_i(c1),d_i(d1),li;  // point intervals!
    idotprecision Akku(0.0);
    ex = expo(0.0);  res = 0.0;  // Initialization for a*b + c*d == 0;

    if ( ex_a1 == ex || ex_b1 == ex ) // a*b == 0;
	if ( ex_c1 == ex || ex_d1 == ex ) ;  // a*b == c*d == 0.0;
	else 
	{// a*b == 0;  c*d != 0;
	    m = ex_c1 + ex_d1;
	    if (m > maxexpo) 
	    {   
		p = maxexpo - m;
		skale_down_exp(ex_c1, ex_d1, p, D1, D2); 
		Times2pown(c_i,D1);
		Times2pown(d_i,D2);
		Akku = 0.0;
		accumulate(Akku,c_i,d_i);
		res = Akku;
		ex = -p;
	    }
	    else
	    {
		skale_up_exp1(ex_c1, ex_d1, p, D1, D2);
		Times2pown(c_i,D1);
		Times2pown(d_i,D2);
		Akku = 0.0;
		accumulate(Akku,c_i,d_i);
		res = Akku;
		ex = -p;
	    }   
	}
    else // a*b != 0.0;
	if ( ex_c1 == ex || ex_d1 == ex )
	{// a*b != 0.0;  c*d == 0.0;
	    m = ex_a1 + ex_b1;
	    if (m > maxexpo)
	    {
		p = maxexpo - m;
		skale_down_exp(ex_a1, ex_b1, p, D1, D2); 
		Times2pown(a_i,D1);
		Times2pown(b_i,D2);
		Akku = 0.0;
		accumulate(Akku,a_i,b_i);
		res = Akku;
		ex = -p;
	    }
	    else
	    {
		skale_up_exp1(ex_a1, ex_b1, p, D1, D2);
		Times2pown(a_i,D1);
		Times2pown(b_i,D2);
		Akku = 0.0;
		accumulate(Akku,a_i,b_i);
		res = Akku;
		ex = -p;
	    }
	}
	else // a*b != 0.0  and  c*d != 0.0;  
	{
	    if ( (ex_c1+ex_d1) > (ex_a1+ex_b1) )
	    {
		li = a_i;  a_i = c_i;  c_i = li;  // a_i <--> c_i
		li = b_i;  b_i = d_i;  d_i = li;  // b_i <--> d_i
		m = ex_a1;  ex_a1 = ex_c1;  ex_c1 = m;
		m = ex_b1;  ex_b1 = ex_d1;  ex_d1 = m;
	    }
	    m = ex_a1 + ex_b1;
            // ab jetzt gilt:  m = (ex_a1+ex_b1) >= (ex_c1+ex_d1);
	    if (m > maxexpo)
	    {   
		p = maxexpo - m;
		skale_down_exp(ex_a1,ex_b1,p,D1,D2);
		Times2pown(a_i,D1);
		Times2pown(b_i,D2);
		skale_down_exp(ex_c1,ex_d1,p,D1,D2);
		Times2pown(c_i,D1);
		Times2pown(d_i,D2);
		Akku = 0.0;
		accumulate(Akku,a_i,b_i);
		accumulate(Akku,c_i,d_i);
		res = Akku;
		ex = -p; 
	    }
	    else
	    {
		skale_up_exp1(ex_a1, ex_b1, p, D1, D2);
		Times2pown(a_i,D1);
		Times2pown(b_i,D2);
		skale_up_exp2(ex_c1, ex_d1, p, D1, D2);
		Times2pown(c_i,D1);
		Times2pown(d_i,D2);
		Akku = 0.0;
		accumulate(Akku,a_i,b_i);
		accumulate(Akku,c_i,d_i);
		res = Akku;
		ex = -p; 
	    }
	}
} // product(...)

void product(const l_real& c, const l_real& d, int& ex, l_interval& res)
// Calulation of an inclusion of: c*c + d*d
{
    l_real c1(c), d1(d);

    c1 += 0.0;   d1 += 0.0;
    int ex_c1( expo(c1[1]) ), ex_d1( expo(d1[1]) ), m;

    l_interval c_i(c1),d_i(d1),li;  // point intervals!
    idotprecision Akku(0.0);
    ex = expo(0.0);  res = 0.0;  // Initialization for c*c + d*d == 0;

    if (ex_c1 == ex) // c*c == 0;
	if (ex_d1 == ex) ; // c*c == d*d == 0; 
	else // c*c == 0;  d*d != 0;
	{
	    times2pown(d_i,-ex_d1); // d_i about 1.0;
	    Akku = 0.0;
	    accumulate(Akku,d_i,d_i);
	    res = Akku;
	    ex = 2*ex_d1;
	}
    else // c*c != 0;
	if (ex_d1 == ex) 
	{ // c*c != 0;   d*d == 0;
	    times2pown(c_i,-ex_c1); // c_i about 1.0;
	    Akku = 0.0;
	    accumulate(Akku,c_i,c_i);
	    res = Akku;
	    ex = 2*ex_c1;
	}
	else // c*c != 0;   d*d != 0; 
	{
	    if (ex_d1 > ex_c1)
	    { 
		li = c_i;  c_i = d_i;  d_i = li; // c_i <--> d_i
		m = ex_c1; ex_c1 = ex_d1;  ex_d1 = m; 
	    }  // c*c >= d*d:
	    times2pown(c_i,-ex_c1); // c_i about 1.0;
	    times2pown(d_i,-ex_c1);
	    Akku = 0.0;
	    accumulate(Akku,c_i,c_i);
	    accumulate(Akku,d_i,d_i);
	    res = Akku;
	    ex = 2*ex_c1;
	}

} // product(...)

l_real quotient(const l_interval& z, const l_interval& n, int round, 
		int ex_z, int ex_n)
// z is an inclusion of a numerator.
// n is an inclusion of a denominator.
// quotient(...) calculates with q1 an approximation of z/n
// using staggered arithmetic.
// Rounding with round (-1,0,+1) is considered.

{ 
    l_real q1; // return value
    int ex_diff; 
    l_interval res;

    if (0.0<=n) // 0 in denominator n leads to error message:
    {
	std::cerr << "quotient1(const l_interval& z, const l_interval& n, int round, int ex_z, int ex_n):  Division by zero" << std::endl;
	exit(1);
    }

    if ( zero_(z) ) { q1=0;  return q1; };

    ex_diff = ex_z - ex_n;
    res = z/n;
    Times2pown(res,ex_diff);
    switch(round)
    {
	case RND_DOWN:
	    q1 = Inf(res);
	    break;
	case RND_NEXT:
	    q1 = mid(res);
	    break;
	case RND_UP:
	    q1 = Sup(res);
	    break;
    } // switch
    return q1;
} // quotient

l_complex _c_division(l_complex a, l_complex b, int round)
{
    int ex1, ex2;
    l_interval z, n;
    l_complex tmp;

    product(Re(b),Im(b),ex2,n);
    product(Re(a),Re(b),Im(a),Im(b),ex1,z);
    SetRe(tmp, quotient(z,n,round,ex1,ex2));
    product(Im(a),Re(b),-Re(a),Im(b),ex1,z);
    SetIm(tmp, quotient(z,n,round,ex1,ex2));
    return tmp;
} // _c_division

l_complex divn (const l_complex & a, const l_complex & b)  
{  
   return _c_division(a,b,RND_NEXT);
}

l_complex divd (const l_complex & a, const l_complex & b) 
{  
   return _c_division(a,b,RND_DOWN);
}

l_complex divu (const l_complex & a, const l_complex & b)  
{  
   return _c_division(a,b,RND_UP);
}

l_complex operator / (const l_complex &a, const l_complex &b) throw()
{
   return divn(a,b);
}

int StagPrec(const l_complex& lc) throw()
{
    return StagPrec(lc.re);
}      

void accumulate(cdotprecision& cd, const l_complex& lc1, 
                                   const l_complex& lc2) throw()
{
    accumulate(Re(cd),lc1.re,lc2.re); 
    accumulate(Re(cd),-lc1.im,lc2.im);
    accumulate(Im(cd),lc1.im,lc2.re);
    accumulate(Im(cd),lc1.re,lc2.im);
}

void accumulate(cdotprecision& cd, const l_complex& lc, 
                                   const complex& c) throw()
{
    accumulate(Re(cd),lc.re,Re(c)); 
    accumulate(Re(cd),-lc.im,Im(c));
    accumulate(Im(cd),lc.im,Re(c));
    accumulate(Im(cd),lc.re,Im(c));
}

void accumulate(cdotprecision& cd, const l_complex& lc, 
                                   const real& r) throw()
{
    accumulate(Re(cd),lc.re,r); 
    accumulate(Im(cd),lc.im,r);
}

void accumulate(cdotprecision& cd, const l_complex& lc, 
                                   const l_real& lr) throw()
{
    accumulate(Re(cd),lc.re,lr); 
    accumulate(Im(cd),lc.im,lr);
}

l_real abs2(const l_complex &a) throw()
{
   dotprecision dot(0.0);
   accumulate(dot,a.re,a.re);
   accumulate(dot,a.im,a.im);
   return l_real(dot);
}

l_real abs(const l_complex& z) throw()
// Calculation of an approximation of |z|.
// In general the maximum precision is stagprec=19, predifined by the used
// sqrt-function declared in l_rmath.hpp.
// If the difference |exa-exb| of the exponents (base 2) is sufficient high,
// precision and accuracy can be choosen greater 19.
{
    return sqrtx2y2(Re(z),Im(z));
} // abs(...)

l_complex & SetIm(l_complex & a,const l_real & b) 
        { a.im=b; return a; } // See SetRe(...);

l_complex & SetRe(l_complex & a,const l_real & b)
	{ a.re=b; return a; } // The real part of a is substituted by b.

l_real & Re(l_complex& a) { return a.re; }
l_real   Re(const l_complex& a) { return a.re; }
l_real & Im(l_complex& a) { return a.im; }
l_real   Im(const l_complex& a) { return a.im; }

complex & complex::operator = (const l_complex& a) throw()
{
	real x(Re(a)), y(Im(a));
	return *this = complex(x,y);
}

} // end namespace cxsc






