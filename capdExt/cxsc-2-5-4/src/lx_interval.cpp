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

/* CVS $Id: lx_interval.cpp,v 1.13 2014/01/30 17:23:47 cxsc Exp $ */


/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

#include "lx_interval.hpp"

//class lx_real;

namespace cxsc {

// ----------------------------------------------------------------------
// ------------ Functions related to type lx_interval: ------------------
// ----------------------------------------------------------------------

lx_interval::lx_interval(const real& n, const string &s) throw()
// Constructor:
// (10^n, string) --- > (ex,li) vom Typ lx_interval
// string s must be of the form:  
//        s = "[1.234567,1.234567]";   or:  s = "[0,0]";
// or:    s = "[1.234567e+20,1.234568e+20]";
// or:    s = "[1.234567e-20,1.234568e-20]";
// or:    s = "[0,1]";   or: s = "[-1,0]";  or: s = "[-1,+1]";
// (ex,li) is an optimal inclusion of  10^n*s;

// Not allowed strings:
//        s = "1.234567";     or:   s = "1.234567e+20";
// or:    s = "[1.2,1.2";     or:   s = "1.2,1.2]";
// or:    s = "+[1.2,1.2]";   or:   s = "-[1.2,1.2]";

{
	const real c1 = 3.3219280948873623; // c1 = ln(10)/ln(2)
	const real c2 = 2.711437152598000E+015;
	const real c3 = 10000.0;
	l_interval x(0);
	lx_interval u,v;
	real y,ay,r;
	int p,q,stagsave;
	bool neg;
	
	stagsave = stagprec;

	if (!Is_Integer(n))
		cxscthrow(REAL_NOT_ALLOWED(
					 "lx_interval(const real&, const string& )"));
   // Check for brackets [...] in string s:
	string str1 = "[", str2 = "]";
	q = s.find_first_of(str1);
	p = s.find_first_of(str2);
	if (q == -1 || p == -1)
		cxscthrow(NO_BRACKETS_IN_STRING(
					 "lx_interval(const real&, const string& )"));
	stagprec = stagprec +1;	
	x = x + 0; // To get StagPrec(x) = actual value of stagprec;
	s >> x;    // [a,b] <= x; x of type l_interval;

	if (x==0)
	{
		li = 0;  ex = 0;
	}
	else
	{
		neg = (n<0);
		y = n*c1;
		y = std::floor(_double(y));
		ay = abs(y);  r = abs(n);
		if (ay>=Max_Int_R)
			cxscthrow(REAL_NOT_ALLOWED(
					 "lx_interval(const real&, const string& )"));
		u = power( lx_interval(0,l_interval(10)),r );
		v = lx_interval(ay,l_interval(1));
		if (r > c2)
		{
			v = lx_interval(ay-c3,l_interval(1));
			times2pown(u,-c3);
		}
      if (neg) u = v/u;
		else u = u/v;
		u = u*lx_interval(0,x);
		li = li_part(u);
		r = expo(u);
		ex = add_real(y,r);
	}
	stagprec = stagsave;
	li = adjust(li);
}

void scale_down(lx_interval &a)
// scale_down skaliert a moeglichst weit nach unten, ohne
// dass dabei binaere Stellen verloren gehen.
// Der mathematische Wert von a bleibt dabei erhalten, d.h.
// 2^(a.ex)*a.li  hat vor und nach der Skalierung stets den
// gleichen mathematischen Wert.
{
    int d;
    d = -1021 - expo_sm(a.li);
    // In case of a.li=0 we get d>0, and so nothing is done! 
    if (d<0)
    {
	Times2pown(a.li,d); // a.li = a.li * 2^d
	a.ex = sub_real(a.ex,d);
    }
}

void scale_up(lx_interval &a)
// scale_up skaliert a.li moeglichst weit nach oben.  
// Nach der Skalierung wird jedoch  a.li + a.li  keinen Overflow 
// erzeugen.
// 2^(a.ex)*a.li  hat fuer d>0 vor und nach der Skalierung stets
// den gleichen Wert. Nur im Fall d<0 kann es durch die 
// Skalierung nach unten nur zu einer minimalen Aufblaehung kommen,
// denn im Fall  d<0  nimmt d nur die Werte d = -1,-2 an.   
{
    int d;
    d = 1022 - expo_gr(a.li);
    if (d>-4) // if (d>-3) is also possible!
    {
	    Times2pown(a.li,d);      // a.li = a.li * 2^d
	    a.ex = sub_real(a.ex,d); // to keep 2^(a.ex)*a.li constant.
    }
}

l_interval & l_interval::operator = (const lx_interval &a) throw()
{
    int exa;
    real p(expo(a)),exr;
    l_interval al(li_part(a));
	 l_real lrI(Inf(al)), lrS(Sup(al));  // Neu !!!!!
    exa = expo_gr(al);
    if (exa<-100000) return (*this) = al; // al = 0;
    exr = p + exa;
    if (exr > 1024)
	cxscthrow(OVERFLOW_ERROR("LX_INTERVAL_UNREALIZABLE_AS_L_INTERVAL: l_interval & operator = (const lx_interval &a)"));
    Times2pown(al,p);
    if (Inf(al)<0 && lrI>=0)
       al = SetInf(al,l_real(0));
	 if (Sup(al)>0 && lrS<=0)
       al = SetSup(al,l_real(0));
    return (*this) = al;
}

interval & interval::operator = (const lx_interval &a) throw()
{
    l_interval al;
    interval z;
    al = a;
    z = al;
    return (*this) = z;
}

l_interval times2powr(const l_real &a, const real &r)
// b = times2powr(a,r); returns an inclusion of a*2^r <= b;  
// This inclusion holds only, if the following condition is
// fulfilled:
// ************************************************************** 
// *     For r in [-2100,+2100] r must be an integer value;     *
// **************************************************************
// The above condition is NOT tested in this function !!
// If r lies outside [-2100,+2100] then it holds a*2^r <= b for
// any real number, unless an overflow occurs.
// The function is for the internal use only!
// Blomquist, 25.02,2008;
{
    l_interval res(0),ai;
    double dbl;
    int ex;

    if ( a != 0 )
    {
	ai = l_interval(a);

	if (abs(r)<2147483647)
	{
	    dbl = _double(r);
	    ex = (int) dbl;
	    Times2pown(ai,ex);
	    res = ai;
	}
	else
	    if (r < -2147483645)
	    {
		    Times2pown(ai,-3000);
		    res = ai;
	    }
	    else
		   Times2pown(ai,3000); // produces overflow!
    }
    return res;
}

lx_interval expo2zero(const lx_interval &a) throw(OVERFLOW_ERROR)
// The exponent to base 2 of the object a is set to zero.
// The return value res is an inclusion of the object a
// unless an overflow occurs. 
{
    lx_interval res(0);
    l_interval la(li_part(a));
    int exa(expo_gr(la));
    real na(expo(a)), r;

    if (exa > -100000) // a != 0;
	{
	    r = add_real(exa,na);
	    if (r > 1024)
		cxscthrow(OVERFLOW_ERROR(
			      "lx_interval expo2zero(lx_interval &)"));
	    Times2pown(la,na);
	    res = lx_interval(0.0,la);
	}

    return res;
}

lx_interval::lx_interval(const lx_real& a, const lx_real& b) throw()
{
    lx_real a_(a),b_(b);
    l_real lr;
    l_interval h;
    bool a_zero(eq_zero(a)), b_zero(eq_zero(b));
    real na,nb,d;
    if (a_zero || b_zero)
	if (a_zero) 
	{
	    ex = expo(b);
	    li = l_interval(0,lr_part(b));
	} else
	{
	    ex = expo(a);
	    li = l_interval(lr_part(a),0);
	}
    else
    {  // Now we have:  a,b != 0:
	scale_up(a_);
	scale_up(b_);
	na = expo(a_);  nb = expo(b_);
	if (na >= nb)
	{
	    if (na==nb) // the most common case!
	    {
		ex = na;
		li = l_interval(lr_part(a_),lr_part(b_));
	    }
	    else // na > nb:
	    {
		d = na - nb;
		if (d>Max_Int_R) 
		    d = Max_Int_R;
		ex = na;
		lr = lr_part(b_);
		h = times2powr(lr,-d);
		li = l_interval(lr_part(a_),Sup(h)); 
	    }
	}
	else // na < nb:
	{
	    d = nb - na;
	    if (d>Max_Int_R) 
		d = Max_Int_R;
	    ex = nb;
	    lr = lr_part(a_);
	    h = times2powr(lr,-d);
	    li = l_interval(Inf(h),lr_part(b_));
	}
    }
}

lx_interval operator + (const lx_interval &a, const lx_interval &b) throw()
{
    const real c1 = 10000.0;
    int stagsave = stagprec,
	stagmax = 39,
	exa,exb;
    real sa,sb,p,n;

    if (stagprec>stagmax) stagprec = stagmax;

    l_interval a1,b1;
    lx_interval res,A(a),B(b);
    interval z;

    A = Blow(A);  // If |a| is too small, a is inflated. 
    B = Blow(B);  // A <= Blow(A) is guaranteed!
    // With these inflations a+b <= A+B is guaranteed and
    // A+B is calculated without an error message unless
    // |A| or |B| are too great. 
    a1 = li_part(A);  b1 = li_part(B);
    exa = expo_gr(a1);
    exb = expo_gr(b1);

    if (exa<-100000) return b;  // summand a = 0
    if (exb<-100000) return a;  // summand b = 0

    // a,b <> 0:
    sa = add_real(exa,expo(A)); 
    sb = add_real(exb,expo(B));
    if (sa > sb) // |a| >= |b|
    {
	p = 1022 - exa; // p >= 0
	Times2pown(a1,p);
	n = sub_real(expo(A),p);
	if (n>0 && expo(B)<sub_real(n,c1))
	{
	    z = b1;
	    if (Sup(z)<=0.0)
		b1 = l_interval(-MinReal,0);
	    else
		if (Inf(z)>=0.0)
		b1 = l_interval(0,MinReal);
		else
		    b1 = l_interval(-MinReal,MinReal);
	}
	else
	{
	    p = sub_real(expo(B),n); 
	    Times2pown(b1,p);
	}
    } else
    {
	p = 1022 - exb; // p >= 0
	Times2pown(b1,p);
	n = sub_real(expo(B),p);
	if (n>0 && expo(A)<sub_real(n,c1))
	{
	    z = a1;
	    if (Sup(z)<=0.0)
		a1 = l_interval(-MinReal,0);
	    else
		if (Inf(z)>=0.0)
		a1 = l_interval(0,MinReal);
		else
		    a1 = l_interval(-MinReal,MinReal);
	}
	else
	{
	    p = sub_real(expo(A),n);
	    Times2pown(a1,p);
	}
    }
    a1 = a1 + b1; // normal staggered addition
    res = lx_interval(n,a1);

    stagprec = stagsave;
    res = adjust(res);
    return res;
}

lx_interval operator * (const lx_interval &a, const lx_interval &b) throw()
{
    int stagsave = stagprec,
	stagmax = 39,
	exa,exb,d,D,k;
    real na,nb,diff;

    if (stagprec > stagmax) stagprec = stagmax;

    l_interval al,bl;
    lx_interval a_(a),b_(b),res(0);

    a_ = Blow(a_); // a_ is inflated, if |a_| is too small
    b_ = Blow(b_); // a_ = [0,0] remains unchanged 

    scale_down(a_);  scale_down(b_);

    al = li_part(a_);  bl = li_part(b_);
    exa = expo_gr(al);
    exb = expo_gr(bl);

    if (exa>-100000 && exb>-100000) // Beide Faktoren sind <> 0:
    {
	if (exa+exb <= 1022)
	{
	    if (exa<0)
	    {
		Times2pown(al,-exa);  
		na = add_real(expo(a_),exa);  // a = 2^na * al;
		d = 1022 - exb;
	    }
	    else // exa >= 0
	    {
		na = expo(a_);
		d = 1022 - (exa+exb);
	    }
	    Times2pown(bl,d);
	    nb = sub_real(expo(b_),d);  // b = 2^nb * bl;
	}
	else // exa+exb > 1022
	{
	    d = exa + exb - 1022;  // d > 0;
	    if (exa>exb)
	    {
		D = exa - exb; // D > 0;
		if (d <= D)
		{
		    Times2pown(al,-d);
		    na = add_real(expo(a_),d); // a = 2^na * al;
		    nb = expo(b_);
		}
		else // d > D;
		{
		    k = d - D;
		    if (k%2 != 0) k++;
		    k = k/2;
		    Times2pown(bl,-k);
		    nb = add_real(expo(b_),k); // b = 2^nb * bl;
		    k = k + D;
		    Times2pown(al,-k);
		    na = add_real(expo(a_),k); // a = 2^na * al;
		}
	    }
	    else // exb >= exa
	    {
		D = exb - exa; // D > 0;
		if (d <= D)
		{
		    Times2pown(bl,-d);
		    nb = add_real(expo(b_),d); // b = 2^nb * bl;
		    na = expo(a_);
		}
		else // d > D
		{
		    k = d - D;
		    if (k%2 != 0) k++;
		    k = k/2;
		    Times2pown(al,-k);
		    na = add_real(expo(a_),k);  // a = 2^na * al;
		    k = k + D;
		    Times2pown(bl,-k);
		    nb = add_real(expo(b_),k);  // b = 2^nb * bl;
		}
	    }
	}
	al = al*bl;      // normal staggered multiplication

	if (na+nb<-Max_Int_R) // Underflow:
	{
	    diff = (Max_Int_R + na) + nb; // diff < 0;
	    Times2pown(al,diff);

	    if (Inf(al)>0)
		SetInf(al,0.0);
	    else
		if (Sup(al)<0)
		    SetSup(al,0.0);

	    res = lx_interval(Max_Int_N,al);
	}
	else // Without underflow
	{
	    na = add_real(na,nb); // Hier normale Rechnung!
	    res = lx_interval(na,al);
            // Caution: This else-block can deliver an interval res,
            // with a so small |res| value , that for example  
            // res * res produces an error message. To avoide this 
            // message, res will be sufficiently inflated:

	    if (na<Max_Int_N) // Only in extreme cases!
	    {
		if (Inf(al)>0)
		    SetInf(al,0.0);
		else
		    if (Sup(al)<0)
			SetSup(al,0.0);

		res = lx_interval(Max_Int_N,al);
	    }
	}

    }
    stagprec = stagsave;
    res = adjust(res);
    return res;
} // operator *

lx_interval operator / (const lx_interval &a, const lx_interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{
    int stagsave = stagprec,
	stagmax = 39,
	exa,exb,d;
    real na,nb,diff;
    if (stagprec>stagmax) stagprec = stagmax;

    interval z;
    l_interval al,bl;
    lx_interval a_(a),b_(b),res;

    a_ = Blow(a_); // a_ is inflated, if |a_| is too small
    b_ = Blow(b_); // a_ = [0,0] remains unchanged 

    scale_down(a_);  
    scale_down(b_);
    al = li_part(a_);  bl = li_part(b_);

    exa = expo_gr(al);
    exb = expo_gr(bl);

    if (exb < -100000) 
	cxscthrow(ERROR_LINTERVAL_DIV_BY_ZERO("lx_interval operator/(const lx_interval &a, const lx_interval &b)"));

    if (exa < -100000) return a;
    // Now it holds:  a,b <> 0:

    bool bl_point( point_intv(bl) );

    // First: scaling the nominator al to 1022:
    d = 1022 - exa;
    Times2pown(al,d);   // expo_gr(al) = 1022;
    na = sub_real(expo(a_),d);  // a = 2^na * al;

    if (bl_point) 
    {
	if (0<=exb && exb<=511) 
	    nb = expo(b_);
	else 
	    if (exb<0)
	    {  // scaling the denominator bl to 0:
		Times2pown(bl,-exb);
		nb = add_real(expo(b_),exb);  // b = 2^nb * bl;
	    }
	    else // exb > 511;
                 // Now scaling the denominator bl to 511:
	    {
		d = 511 - exb;
		Times2pown(bl,d);
		nb = sub_real(expo(b_),d);
	    }
    }
    else // scaling the denominator bl to 511:
    {
	d = 511 - exb;
	Times2pown(bl,d);
	nb = sub_real(expo(b_),d);
    }

    z = al;
    al = al / bl; // Normal staggered division; 
    // However, with too wide interval operands an overflow is possible:
    if ( (al==0 && z!=0) || Inf(al)>Sup(al) )
	cxscthrow(OVERFLOW_ERROR("TOO_WIDE_INTERVAL_OPERANDS_IN: lx_interval operator / (const lx_interval &a, const lx_interval &b)"));

    if (na-nb<-Max_Int_R) // Underflow:
    {
	diff = (Max_Int_R - nb) + na; // diff < 0;
	Times2pown(al,diff);

	if (Inf(al)>0)
	    SetInf(al,0.0);
	else
	    if (Sup(al)<0)
		SetSup(al,0.0);

	res = lx_interval(Max_Int_N,al);
        // +2100 delivers an additional inflation
        // to avoide too small results
    }
    else
    {   // Without underflow:
	na = sub_real(na,nb); // Hier normale Rechnung!
	res = lx_interval(na,al);
	// Caution: This else-block can deliver an interval res,
	// with a so small |res| value , that for example  
	// res * res produces an error message. To avoide this 
	// message, res will be sufficiently inflated:

	if (na<Max_Int_N) // Only in extreme cases!
	{
	    if (Inf(al)>0)
		SetInf(al,0.0);
	    else
		if (Sup(al)<0)
		    SetSup(al,0.0);

	    res = lx_interval(Max_Int_N,al);
	}
    }

    stagprec = stagsave;
    res = adjust(res);
    return res;
}

// ----------------------- Input --------------------------------------------

std::string & operator >> (std::string &s, lx_interval &a) throw()
// Writes string s to variable a of type lx_interval 
// and returns an empty string s;
// Example:  s = "{-4000,[2,2]}" delivers an interval a
// with:    10^(-4000)*[2,2] <= a;
// Whitespaces are NOT allowed in s ! 
{
    l_interval la;
    real exr;
    int i;

    s = skipwhitespacessinglechar (s, '{'); 
    s >> exr;
    s = skipwhitespacessinglechar (s, ','); 
    i = s.find("]");
    s.erase(i+1); 
    a = lx_interval(exr,s);
    s = "";

    return s;
}

void operator >> (const std::string &s, lx_interval &a) throw()
{
// Writes strings s to variable a of type lx_interval;
    std::string r(s);
    r >> a;
}

void operator >> (const char *s, lx_interval &a) throw()
{
    std::string r(s);
    r >> a;
}

bool StrContains(const string &s, const char &a, const char &b)
// returns 1, if s contains char a or char b, otherwise 0;
{
    int ia,ib;
    ia = s.find(a,0);
    ib = s.find(b,0);
    return (ia>=0 || ib>=0);
}

std::istream & operator >> (std::istream &s, lx_interval &a) throw()
// s must have the form: { p,[r1,r2] }, where r1<=r2 are real values
// and the exponent p of type real must be an integer value to base 10;
// The integer value condition and r1<=r2 are tested in the constructor
// call:   a = lx_interval(exr,str);
// Blomquist, 27.01.2008; 
{
    char c;
    real exr;
    string str,sh;

    skipeolnflag = inpdotflag = true;
    c = skipwhitespacessinglechar (s, '{');
    if (inpdotflag) s.putback(c);
    s >> SaveOpt >> exr; 
    c = skipwhitespacessinglechar (s, ',');
    if (inpdotflag) s.putback(c);
    skipeolnflag = inpdotflag = true;
    s >> str >> RestoreOpt;

    while ( !StrContains(str,']','}') )
    {
	c = skipwhitespaces (s);
	if (inpdotflag && c != '}') s.putback(c);
	if (c == '}') break;
	if (c == ' ') break;
	s >> sh;
	str = str + sh;
    }

    a = lx_interval(exr,str);

    if (!waseolnflag)
    {
	skipeolnflag = false; inpdotflag = true;
	c = skipwhitespaces (s);
	if (inpdotflag && c != '}')
	    s.putback(c);
    }

    return s;
}

// ----------------------- Output --------------------------------------------

void Bin2Dec(const lx_interval& a, real& p, l_interval& m)
// Calculates for a given interval a of type lx_interval the
// exponent p and the classic staggered interval m of type
// l_interval, so that the inclusion 
//        a <= 10^p * m        is guaranteed.
// Blomquist, 26.09.2008;
{
	const real c1 = 0.301029995663981195;
	const real c2 = 9.007199254738000E+015;
	const real c3 = 10000.0;
	l_interval x;
	lx_interval u,v;
	real y,k,l;
	int stagsave;
	bool neg;
	
	stagsave = stagprec;
	stagprec = stagprec + 1;
		
	x = li_part(a);
	y = expo(a);
	neg = y<0;
	if (x==0)
	{
		p = 0;  m = 0;
	}
	else
	{
		p = y*c1;
		p = std::floor(_double(p));
		y = abs(y);
		u = lx_interval(y,l_interval(1));
		v = power( lx_interval(0,l_interval(10)),abs(p) );
		if (y > c2)
		{
			u = lx_interval(y-c3,l_interval(1));
			times2pown(v,-c3);
		}
		if (neg) u = v/u;
		else u = u/v;
		u = u*lx_interval(0,x);
		m = li_part(u);
		k = expo(u);
		l = k*c1;
		l = std::floor(_double(l)) + 1;
		p = p+l;
		u = lx_interval(k,l_interval(1)) /
				power( lx_interval(0,l_interval(10)),l );
		u = u * lx_interval(0,m);
		k = expo(u);
		m = li_part(u);
		Times2pown(m,k);
	}
	stagprec = stagsave;
	m = adjust(m);
}  // Bin2Dec(...)

std::ostream& operator << (std::ostream& s,const lx_interval& a) throw()
// An interval a of type lx_interval is written to the 
// output channel in the decimal form: 
// { 10**(p), l_interval m } = 10^p * m;
// The following inclusion is realized:   a <= 10^p * m;
// Blomquist, 18.01.2008;
{
	real p;
	l_interval m;

	Bin2Dec(a,p,m);
	
   s << "{ " 
     << "10**("
     << SaveOpt << SetPrecision(0,0) << Fixed << p << RestoreOpt
     << ")" 
     << "*" 
     << m       
     << " }";
   return s;
}

std::string & operator << (std::string &s, const lx_interval &a) throw()
// The value of a variable a of type lx_interval is copied to a string s.
// s has the form:  {2**(ex), li} = 2^ex * li;
{  
    std::stringstream ss;
    string str;
    s+="{2**(";
    ss << SaveOpt << SetPrecision(0,0) << Fixed << a.ex << RestoreOpt;
    ss >> str;
    s += str;
    s+=")*";
    s << a.li; 
    s+='}';
    return s;
}


// ----- Help Functions, declared in lx_interval.hpp outside ---------------
// --------------------- the class lx_interval------------------------------

l_interval point_max(void)
// returns a staggered point interval with maximum exponent 1020,
// whereby nearly all mantissa bits are set. 
{
    l_real lr;
    l_interval li = sqrt(l_interval(3.995)), res(li);
    times2pown(li,1019);
    lr = Inf(li);
    lr = ((lr + Inf(res)) + MinReal) + minreal;
    res = l_interval(lr);
    return res;
}

l_interval point_any(int n)
// returns a staggered point interval with exponent n,
// whereby nearly all mantissa bits are set.
// -1074 <= n <= +1020; 
{
    int n_(n);
    l_interval res;

    if (n_>1020)  n_ = 1020;
    if (n_<-1074) n_ = -1074;
    res = point_max();
    times2pown(res,n_ - 1020);
    return l_interval(Inf(res+MinReal));
}

l_interval wide_max(void)
// returns a staggered interval a with maximum exponent 1020
// and diam(a)>0, whereby nearly all mantissa bits are set. 
{
    l_interval a = point_max();
    return a + interval(MinReal,2*MinReal);
}

l_interval wide_any(int n)
// returns a staggered point interval a with exponent n,
// and diam(a)>0, whereby nearly all mantissa bits are set.
// -1074 <= n <= +1020;
{
    l_interval a = point_any(n);
    return a + interval(MinReal,2*MinReal);
}

// ----------------------------------------------------------------------------
// ---------------- Functions related to class lx_real ------------------------
// ----------------------------------------------------------------------------

	lx_real::lx_real(const real& n, const string &s) throw()
// Constructor:
// (10^n, string) --- > (ex,lr) vom Typ lx_real
// string s must be of the form:  
//        s = "1.234567";   or:  s = "0";
// or:    s = "1.234567e+20";
// or:    s = "1.234567e-20";
// (ex,li) is an approximation of  10^n*s;
{
	const real c1 = 3.3219280948873623; // c1 = ln(10)/ln(2)
	const real c2 = 2.711437152598000E+015;
	const real c3 = 10000.0;
	l_real x(0);
	lx_interval u,v;
	real y,ay,r;
	int stagsave;
	bool neg;
	
	stagsave = stagprec;

	if (!Is_Integer(n))
		cxscthrow(REAL_NOT_ALLOWED(
					 "lx_real(const real&, const string& )"));
	stagprec = stagprec +1;	
	x = x + 0; // To get StagPrec(x) = actual value of stagprec;
	s >> x;    // x ~ s; x of type l_real;
	if (x==0)
	{
		lr = 0;  ex = 0;
	}
	else
	{
		neg = (n<0);
		y = n*c1;
		y = std::floor(_double(y));
		ay = abs(y);  r = abs(n);
		if (ay>=Max_Int_R)
			cxscthrow(REAL_NOT_ALLOWED(
						 "lx_real(const real&, const string& )"));
		u = power( lx_interval(0,l_interval(10)),r );
		v = lx_interval(ay,l_interval(1));
		if (r > c2)
		{
			v = lx_interval(ay-c3,l_interval(1));
			times2pown(u,-c3);
		}
		if (neg) u = v/u;
		else u = u/v;
		u = u*lx_interval(0,l_interval(x));
		lr = mid(li_part(u));
		r = expo(u);
		ex = add_real(y,r);
	}
	stagprec = stagsave;
	lr = adjust(lr);
}

	std::string & operator >> (std::string &s, lx_real &a) throw()
// Writes string s to variable a of type lx_real 
// and returns an empty string s;
// Example:  s = "{-4000,2}" delivers a value a
// with:    10^(-4000)*2 ~ a;
{
	real exr;
	int i;

	s = skipwhitespacessinglechar (s, '{'); 
	s >> exr;
	s = skipwhitespacessinglechar (s, ','); 
	i = s.find("}");
	s.erase(i+1); 
	a = lx_real(exr,s);
	s = "";

	return s;
}

	void operator >> (const std::string &s, lx_real &a) throw()
{
   // Writes string s to variable a of type lx_real;
	std::string r(s);
	r >> a;
}

	void operator >> (const char *s, lx_real &a) throw()
{
	std::string r(s);
	r >> a;
}	
	
	bool Str_Contains(const string &s, const char &a, const char &b)
// returns 1, if s contains char a or char b, otherwise 0;
{
	int ia,ib;
	ia = s.find(a,0);
	ib = s.find(b,0);
	return (ia>=0 || ib>=0);
}
	
	std::istream & operator >> (std::istream &s, lx_real &a) throw()
// s must have the form: { p,r }, where r is a real value
// and the exponent p of type real must be an integer value to base 10;
// The integer value condition is tested in the constructor call:
//    a = lx_interval(exr,str);
// Blomquist, 06.101.2008; 
{
	char c;
	real exr;
	string str,sh;
	skipeolnflag = inpdotflag = true;
	c = skipwhitespacessinglechar (s, '{');
	if (inpdotflag) s.putback(c);
	s >> SaveOpt >> exr; 
	c = skipwhitespacessinglechar (s, ',');
	if (inpdotflag) s.putback(c);
	skipeolnflag = inpdotflag = true;
	s >> str >> RestoreOpt;
	while ( !Str_Contains(str,']','}') )
	{
		c = skipwhitespaces (s);
		if (inpdotflag && c != '}') s.putback(c);
		if (c == '}') break;
		if (c == ' ') break;
		s >> sh;
		str = str + sh;
	}
		
	a = lx_real(exr,str);

	if (!waseolnflag)
	{
		skipeolnflag = false; inpdotflag = true;
		c = skipwhitespaces (s);
		if (inpdotflag && c != '}')
			s.putback(c);
	}
	return s;
}


bool operator == (const lx_real &a, const lx_real &b) throw()
{
	const real c1 = 1e20;
	l_real ar(lr_part(a)), br(lr_part(b));
	real na(expo(a)), nb(expo(b));
	int exa(expo_gr(ar)), exb(expo_gr(br)),d;
	bool a_0(exa<-100000), b_0(exb<-100000), res;
	interval z(0);

	if (a_0 || b_0) res = (a_0 == b_0);
    else // a,b <> 0:
	 {
		 d = exa-exb;
		 if (d>0)
		 {  // br nach oben skalieren: (Overflow dabei nicht moeglich!)
			 Times2pown(br,z,d);
			 nb = nb - d;
			 if (abs(nb)>Max_Int_R)
				 nb = c1;        // equality not possible
		 } else // d<=0:
		 {  // ar nach oben skalieren: (Overflow dabei nicht moeglich!)
			 Times2pown(ar,z,-d);
			 na = na + d;
			 if (abs(na)>Max_Int_R)
				 na = c1;        // equality not possible
		 }
		 res = (na==nb && ar==br);
	 }
	 return res;
}

bool operator > (const lx_real &a, const lx_real &b) throw()
{
	l_real lra(lr_part(a)), lrb(lr_part(b));
	bool zero_a(eq_zero(a)), zero_b(eq_zero(b)), bl(false);
	real na,nb,d;
	double dbl;
	interval z(0);
	int D;

	if (zero_a) return sm_zero(b);
	if (zero_b) return gr_zero(a);
    // Now we have:  a,b != 0:
	na = expo(a);  nb = expo(b);
	d = na - nb;

	if (na==nb)
		bl = (lra >lrb);
	else 
		if (na>nb)
	{
		if (d > 1024-expo_gr(lra))
            // Scaling not possible
			bl = (sign(lra)>0) ? true : false;
		else
		{
			dbl = _double(d);
			D = (int) dbl;
			Times2pown(lra,z,D);
			bl = lra > lrb;
		}
	}
	else // nb > na
	{
		if (d < expo_gr(lrb)-1024)
	    // Scaling not possible
			bl = (sign(lrb)<0) ? true : false;
		else
		{
			dbl = _double(-d);
			D = (int) dbl;
			Times2pown(lrb,z,D);
			bl = lra > lrb;
		}
	}

	return bl;
}

void scale_up(lx_real &a) throw()
// scale_up(a) scales a.lr upwardly as far as possible. 
// Notice:    a.lr must absolutely be sorted, i.e. for
//            example: a.lr[1]=1e50, a.lr[2]=1e20, a.lr[3]=1;   
// In case of   (d>0 && a.ex-d >= -Max_Int_R)   the values 
// of 2^(a.ex)*a.lr are equal before and after the scaling.
// In case of (d<=0 || a.ex-d < -Max_Int_R) no scaling is done.
{
	interval z(0);
	int d;
	d = 1023 - expo_gr(a.lr);  // d >= 0, if a.lr <> 0;
	if (d>0 && a.ex >= -Max_Int_R+d) // If (a==0) then d<0 and nothing is done!
	{
		Times2pown(a.lr,z,d);    // a.lr wird mit 2^d multipliziert.
		a.ex = a.ex - d; // damit 2^(a.ex)*a.li konstant bleibt.
	}
}

void scale_down(lx_real &a) throw()
// scale_down skaliert a moeglichst weit nach unten, ohne
// dass dabei binaere Stellen verloren gehen.
// Der mathematische Wert von a bleibt dabei erhalten, d.h.
// 2^(a.ex)*a.li  hat vor und nach der Skalierung stets den
// gleichen mathematischen Wert.
{
	interval z(0);
	int d;
	d = -1021 - expo_sm(a.lr);
    // In case of a.lr=0 we get d>0, and so nothing is done! 
	if (d<0 && a.ex <= Max_Int_R+d)
	{
		Times2pown(a.lr,z,d); // a.lr = a.lr * 2^d
		a.ex = a.ex - d;
	}
}


lx_real upper_bnd(const lx_real& x) throw()
// y = upper_bnd(x) calculates an upper bound y of  x < y;
// lr = lr_part(x); (See the following code.)
// The components lr[1],lr[2], ... ,lr[StagPrec(x)] must be sorted!
{
	int stagsave = stagprec,
 p(StagPrec(x)),ex_gr;
 stagprec = p;

 lx_real res;
 l_real lr(lr_part(x));
 real ex(expo(x));

 lr = lr + 0; // sorting the lr[i]
 res = lx_real(ex,lr);
 if (p>1)
	 scale_up(res);
 lr = lr_part(res);
 ex = expo(res);
 ex_gr = expo_gr(lr);
 if (ex_gr>-10000000)
 {   // x<>0:
	 if (lr[1]==MaxReal)
	 {   
		 times2pown(lr,-1);
		 ex = add_real(ex,1);
	 }
	 lr[p] = succ(lr[p]);
	 lr = lr + 0; // sorting lr !
	 res = lx_real(ex,lr);
 }
    else // x = 0:
	 {
		 lr = minreal;
		 lr = adjust(lr);
		 lr = lr + 0; // sorting lr !
		 res = lx_real(-Max_Int_R,lr);
	 }
	 stagprec = stagsave;

	 return res;
} // upper_bnd(...)

lx_real lower_bnd(const lx_real& x) throw()
// y = lower_bnd(x) calculates a rather great lower bound y of  x > y;
// lr = lr_part(x); (See the following code.)
// The components lr[1],lr[2], ... ,lr[StagPrec(x)] must be sorted!
{
	int stagsave = stagprec,
 p(StagPrec(x)),ex_gr;

 stagprec = p;
 lx_real res;
 l_real lr(lr_part(x));
 real ex(expo(x));

 lr = lr + 0; // sorting the lr[i]
 res = lx_real(ex,lr);
 if (p>1)
	 scale_up(res);
 lr = lr_part(res);
 ex = expo(res);
 ex_gr = expo_gr(lr);
 if (ex_gr>-10000000)
 {
	 if (lr[1]==-MaxReal)
	 {
		 times2pown(lr,-1);
		 ex = add_real(ex,1);
	 }
	 lr[p] = pred(lr[p]);
	 lr = lr + 0; // sorting lr !
	 res = lx_real(ex,lr);
 }
 else // x = 0;
 {
	 lr = -minreal;
	 lr = adjust(lr);
	 lr = lr + 0; // sorting lr !
	 res = lx_real(-Max_Int_R,lr);
 }
 stagprec = stagsave;
 return res;
} // lower_bnd(...)

lx_real operator + (const lx_real &a, const lx_real &b) throw()
// An approximation of (a+b) is calculated. In most of the cases
// a maximum of accuracy is achieved.
{
   //const real c1 = 10000.0;
   const int  c2 = 1022,
              c3 = -100000;
 int stagsave = stagprec,
 stagmax = 39,
 exa,exb;
 real sa,sb,p,n;
	
 if (stagprec>stagmax) stagprec = stagmax;
 l_real a1,b1;
 lx_real res(0),A(a),B(b);
	
 a1 = lr_part(A);  b1 = lr_part(B);
 exa = expo_gr(a1);
 exb = expo_gr(b1);
	
 if (exa<c3) return b;  // summans a = 0;
 if (exb<c3) return a;  // summans b = 0;
	
	// From now on: a,b <> 0
 sa = add_real(exa,expo(A));
 sb = add_real(exb,expo(B));
 if (sa > sb)  // |a| >= |b|
 {
	 p = c2 - exa;
	 Times2pown(a1,p);
	 n = sub_real(expo(A),p);
	 p = sub_real(expo(B),n);
	 Times2pown(b1,p);
 }
 else
 {
	 p = c2 - exb;
	 Times2pown(b1,p);
	 n = sub_real(expo(B),p);
	 p = sub_real(expo(A),n);
	 Times2pown(a1,p);
 }
 a1 = a1 + b1;
 res = lx_real(n,a1);
	
 stagprec = stagsave;
 res = adjust(res);
 return res;
}

lx_real operator * (const lx_real& a, const lx_real& b) throw()
// An approximation of (a*b) is calculated. In most of the cases
// a maximum of accuracy is achieved.
{
	int stagsave = stagprec,
 stagmax = 39,
 exa,exb,d,D,k;
 real na,nb,diff;

 if (stagprec > stagmax) stagprec = stagmax;

 l_real al,bl;
 lx_real a_(a),b_(b),res(0);

 scale_down(a_);  scale_down(b_);

 al = lr_part(a_);  bl = lr_part(b_);
 exa = expo_gr(al);
 exb = expo_gr(bl);

   if (exa>-100000 && exb>-100000) // Beide Faktoren sind <> 0:
	{
		if (exa+exb <= 1022)
		{
			if (exa<0)
			{
				Times2pown(al,-exa);  
				na = add_real(expo(a_),exa);  // a = 2^na * al;
				d = 1022 - exb;
			}
			else // exa >= 0
			{
				na = expo(a_);
				d = 1022 - (exa+exb);
			}
			Times2pown(bl,d);
			nb = sub_real(expo(b_),d);  // b = 2^nb * bl;
		}
		else // exa+exb > 1022
		{
			d = exa + exb - 1022;  // d > 0;
			if (exa>exb)
			{
				D = exa - exb; // D > 0;
				if (d <= D)
				{
					Times2pown(al,-d);
					na = add_real(expo(a_),d); // a = 2^na * al;
					nb = expo(b_);
				}
				else // d > D;
				{
					k = d - D;
					if (k%2 != 0) k++;
					k = k/2;
					Times2pown(bl,-k);
					nb = add_real(expo(b_),k); // b = 2^nb * bl;
					k = k + D;
					Times2pown(al,-k);
					na = add_real(expo(a_),k); // a = 2^na * al;
				}
			}
			else // exb >= exa
			{
				D = exb - exa; // D > 0;
				if (d <= D)
				{
					Times2pown(bl,-d);
					nb = add_real(expo(b_),d); // b = 2^nb * bl;
					na = expo(a_);
				}
				else // d > D
				{
					k = d - D;
					if (k%2 != 0) k++;
					k = k/2;
					Times2pown(al,-k);
					na = add_real(expo(a_),k);  // a = 2^na * al;
					k = k + D;
					Times2pown(bl,-k);
					nb = add_real(expo(b_),k);  // b = 2^nb * bl;
				}
			}
		}
		al = al*bl;      // normal staggered multiplication

	if (na+nb<-Max_Int_R) // Underflow:
	{
		diff = (Max_Int_R + na) + nb; // diff < 0;
		Times2pown(al,diff);
		res = lx_real(-Max_Int_R,al);
	}
	else // Without underflow
	{
		na = add_real(na,nb); // Hier normale Rechnung!
		res = lx_real(na,al);
	}

	}
	stagprec = stagsave;
	res = adjust(res);
	return res;
} // operator *

lx_real operator / (const lx_real &a, const lx_real &b) throw(DIV_BY_ZERO)
{
	int stagsave = stagprec,
 stagmax = 39,
 exa,exb,d;
 real na,nb,diff;
 if (stagprec>stagmax) stagprec = stagmax;

 l_real al,bl;
 lx_real a_(a),b_(b),res;

 scale_down(a_);  
 scale_down(b_);
 al = lr_part(a_);  bl = lr_part(b_);

 exa = expo_gr(al);
 exb = expo_gr(bl);

 if (exb < -100000) 
	 cxscthrow(DIV_BY_ZERO("lx_real operator/(const lx_real &a, const lx_real &b)"));

 if (exa < -100000) return a;
   // Now it holds:  a,b <> 0:
   // First: scaling the nominator al to 1022:
 d = 1022 - exa;
 Times2pown(al,d);   // expo_gr(al) = 1022;
 na = sub_real(expo(a_),d);  // a = 2^na * al;
 if (0<=exb && exb<=511) 
	 nb = expo(b_);
 else 
	 if (exb<0)
 {  // scaling the denominator bl to 0:
	 Times2pown(bl,-exb);
	 nb = add_real(expo(b_),exb);  // b = 2^nb * bl;
 }
 else // exb > 511;
           // Now scaling the denominator bl to 511:
 {
	 d = 511 - exb;
	 Times2pown(bl,d);
	 nb = sub_real(expo(b_),d);
 }

 al = al / bl; // Normal staggered division;
    if (na-nb<-Max_Int_R) // Underflow:
	 {
		 diff = (Max_Int_R - nb) + na; // diff < 0;
		 Times2pown(al,diff);
		 res = lx_real(-Max_Int_R,al);
	 }
	 else
	 {   // Without underflow:
		 na = sub_real(na,nb); // Hier normale Rechnung!
		 res = lx_real(na,al);
	 }

	 stagprec = stagsave;
	 res = adjust(res);

	 return res;
} // operator /

l_real & l_real::operator = (const lx_real& a) throw()
{
	int exar;
	real p(expo(a)), exr;
	l_real ar(lr_part(a));
	exar = expo_gr(ar);
	if (exar<-100000) return (*this) = ar; // ar = 0;
	exr = p + exar;
	if (exr > 1024)
		cxscthrow(OVERFLOW_ERROR(
					 "LX_REAL_UNREALIZABLE_AS_L_REAL: l_real & operator = (const lx_real& a)"));
	Times2pown(ar,p);
	
	return *this = ar;
}

real & real::operator = (const lx_real& a) throw()
{
	l_real lr;
	real x;
	
	lr = a;
	x = lr;
	
	return *this = x;
}

// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------

real expo_RelDiam(const l_interval &x)
// rd := RelDiam(x) is defined as usual;
// rd '=' 2^n; n is returned by expo_RelDiam();
// |n| is a good approximation of the first
// common binary digits of the interval x.
// x is a point interval --> n <= -2147483647 
// x = [4,5] --> n = -1;
// x = [20000,20001] --> n = -14; 
// This function is only for the internal use!
// Blomquist, 11.02.2008; 
{
	lx_real R;
	real n;
	lx_interval t;
	l_real r1(Inf(x)),r2(Sup(x));
	t = lx_interval(r2) - lx_interval(r1);

	if ( in(0.0,x) )
		R = Sup(t);
	else
	{
		R = Sup( t / lx_interval( AbsMin(x) ) );
	}
	n = expo(R) + expo_gr(lr_part(R));

	return n;
} // expo_RelDiam()

// --------------------------------------------------------------------------
// -------------- Elementary functions of type lx_interval ------------------
// --------------------------------------------------------------------------

lx_interval sqrt(const lx_interval &a) throw()
{
	int stagsave = stagprec,
 stagmax = 30; // l_imath.cpp: sqrt() uses stagmax = 30;
 if (stagprec>stagmax) stagprec = stagmax;

 lx_interval res;
 l_interval la(li_part(a));
 real na(expo(a)),r(0);
 int exa(expo_gr(la)),d;

 if (exa<-100000) res = 0.0; // la = 0;
 else // la <> 0;
 {
	 if ( !Is_Integer(na/2) )
	 {
		 r = -1.0;
		 if (na>0)
			 na = (na-1)/2 + 1;
		 else na = (na+1)/2;
	 }
	 else // na ist gerade
		 na = na/2;

	 d = 1024 - exa; // d >= 0;
	 r = r - d;
	 if ( !Is_Integer(r/2) )
	 {
		 r = r + 1;
		 d--;
	 }
	 Times2pown(la,d);
	 na = na + r/2;
	 la = sqrt(la);
	 res = lx_interval(na,la);
 }

 stagprec = stagsave;
 res = adjust(res);
 return res;
} // sqrt()

lx_interval sqr(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 39;
 if (stagprec>stagmax) stagprec = stagmax;

 lx_interval y;
 l_interval xl(li_part(x));

 if (Inf(xl) >= 0.0)
	 y = x*x;
 else if (Sup(xl) <= 0.0)
 {
	 y = lx_interval(expo(x),l_interval(-Sup(xl),-Inf(xl)));
	 y = (y) * (y);
 } else
 {
	 if (abs(Inf(xl)) >= abs(Sup(xl)))
		 y = lx_interval(expo(x),l_interval(0.0,abs(Inf(xl))));
	 else y = lx_interval(expo(x),l_interval(0.0,abs(Sup(xl))));
	 y = (y) * (y);
 }

 stagprec = stagsave;
 y = adjust(y);
 return y;
}

lx_interval Lnp1(const lx_interval &x) throw()
// Calculating an inclusion of ln(1+x) for
// not too wide intervals,    |x| <= 1e-7;
// This function is for the internal use only,
// see the functions ln(...), lnp1(...);
// Blomquist, 27.02.2008;
{
	lx_interval res(0),z2,zeta,Ri,Two;
	l_interval xli;
	int N,ex,p;
	real m,expox;

	p = stagprec;
	xli = li_part(x);
	ex = expo_gr(xli);
	if (ex>-100000) // x <> 0
	{
		expox = expo(x);
		if (expox < 3-27*p-ex) N = 0;
		else
		{
			m = ex + expox;  // m = ex + expo(x);
			N = (int) _double( (53*p-4)/(2*(1-m)) ); 
		}
        // N: Polynomial degree, 0<=N<=42;
		zeta = x / (2+x);
		Two = lx_interval(2);
		Ri = lx_interval( Sup(abs(zeta)) );
		Ri = sqr(Ri); 
		if (N==0)
			res = Two;
	else // N >= 1:
	{
		z2 = sqr(zeta); // z2 = zeta^2
	    // Evaluating the polynomial:
		res = Two / (2*N+1);  // res = a_(N)
		for (int i=N-1; i>=0; --i)
			res = res*z2 + Two/(2*i+1);
	    // Calculating the absolute approximation error: 
		Ri = power(Ri,N+1)/(N+1);
	}
	// Implementing the approximation error:
	res += lx_interval(lx_real(0),Sup(Ri));
	res *= zeta;
	} // x <> 0;

	return res;
} // Lnp1(...)

lx_interval Ln_(const lx_interval &x) throw() 
{
	lx_interval T(0.0);
	l_interval lx(li_part(x));
	interval z;
	real eps;
	int exa(expo_gr(lx)),k(0);
	real nx(expo(x)),n;
    // we have: Inf(x)>0 --> exa <> -2147483648:
	if (abs(nx)<5000)
	{
		n = exa + nx; // nx is a guaranteed integer!
		if (n==0 || n==1)
		{
			T = x - 1.0;
			lx = T;
			if (Sup(abs(lx))<1e-7)
			{
				T = Lnp1(T);
				goto fertig;
			}
		}
	}
	T = x;
	times2pown(T,-nx);
	times2pown(T,-exa); // T: 1st reduced argument
	lx = T;
	z = lx;
	eps = Sup(1-z);
	if (eps>1e-12)
	{
		eps = (ln(eps) + 6*ln_N[8])/ln_N[0]; 
		k = (int) _double(eps);
		if (k<=0) k = 0;
	}
	for (int j=1; j<=k; j++) 
		T = sqrt(T); // T: 2nd reduced argument,
	             // T = [2^k]sqrt(T);
	T = Lnp1(T-1.0);
	if (k>0) times2pown(T,k);
	T = nx*Ln2_lx_interval() + exa*Ln2_lx_interval() + T;

 fertig:
		 return T;
} // Ln_(...)

lx_interval ln(const lx_interval &x) throw() 
{
	int stagsave = stagprec,
 stagmax = 39;
 if (stagprec>stagmax) stagprec = stagmax;

 lx_interval res,a;

 if (Inf(li_part(x))<=0) cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF(
	  "lx_interval ln(const lx_interval &)"));
    // Now we have: Inf(x)>0:
 real r;
 r = expo_RelDiam(li_part(x));
 if (r > -107) // If the rel. diam of x is too great
 {   // the ln function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is:          [Inf(res),Sup(a)]
	 res = lx_interval(Inf(x));
	 res = Ln_(res);
	 a = lx_interval(Sup(x));
	 a = Ln_(a);
	 res = lx_interval(Inf(res),Sup(a));
 }
 else 
	 res = Ln_(x);

 stagprec = stagsave;
 res = adjust(res);

 return res;
} // ln(...)

lx_interval log2(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 39;
 if (stagprec>stagmax) stagprec = stagmax;
 lx_interval res;
	
 res = ln(x)/Ln2_lx_interval();
 stagprec = stagsave;
 res = adjust(res);
	
 return res;
}

lx_interval log10(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 39;
 if (stagprec>stagmax) stagprec = stagmax;
 lx_interval res;
	
 res = ln(x)/Ln10_lx_interval();
 stagprec = stagsave;
 res = adjust(res);
	
 return res;	
}

lx_interval Lnp1_(const lx_interval &x) throw()
{
	const real r = 1e-7;
	lx_interval res;

	res = (Inf(x)>-r && Sup(x)<r)? Lnp1(x) : Ln_(1.0+x); 

	return res;
} // Lnp1_(...)

lx_interval lnp1(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 39;
 if (stagprec>stagmax) stagprec = stagmax;

 lx_interval res,a;

 if (Inf(x)<= -1) cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF(
	  "lx_interval lnp1(const lx_interval &)"));
    // Now we have: Inf(x)>-1:
 real r;

 r = expo_RelDiam(li_part(x));
 if (r > -107) // If the rel. diam of x is too great
 {   // the lnp1 function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is:          [Inf(res),Sup(a)]
	 res = lx_interval(Inf(x));
	 res = Lnp1_(res);
	 a = lx_interval(Sup(x));
	 a = Lnp1_(a);
	 res = lx_interval(Inf(res),Sup(a));
 }
 else 
	 res = Lnp1_(x);

 stagprec = stagsave;
 res = adjust(res);

 return res;
}

// -------------------------------------------------------------------------
// --------------------- power function ------------------------------------
// -------------------------------------------------------------------------

lx_interval Power_(const lx_interval &x, const real &n) throw()
// Calculates the inclusion of [x]^n for integer values n.
// The rel. diam of x should not be too great.
// This function is only for the internal use! ( power(...) ) 
// Blomquist, 11.02.2008;
{
	real one(1.0),zhi(2.0),N(n),r; 
	double dbl;

	lx_interval y,neu,X(x);

	if (x == one) y = x;
	else 
		if (N == 0.0) y = one; 
	else 
	{
		if (N == 1) y = x; 
		else 
			if (N == 2) y = sqr(x);
		else 
		{
			if (N < 0) 
			{
				X = 1.0/X;
				N = -N;
			}
		    // Initialisierung
			if ( !Is_Integer(N/2) )
				y = X;
			else y = one; 
			neu = sqr(X);   // neu = X*X;
			do {
				dbl = _double(N/zhi);
				dbl = std::floor(dbl);
				r = (real) dbl;
				if ( !Is_Integer( r/2 ) )
					y *= neu;
				zhi += zhi;
				if (zhi <= N) 
					neu = sqr(neu); // neu = neu * neu;
			} while (zhi <= N);
		}
	}

	return y;
} // Power_()

lx_interval power(const lx_interval &x, const real &n) throw()
{
	int stagsave = stagprec,
 stagmax = 38;
 if ( !(Is_Integer(n)) ) 
	 cxscthrow(REAL_NOT_ALLOWED("lx_interval power(const lx_interval&, const real&)"));

 if (stagprec > stagmax) 
	 stagprec = stagmax;

 lx_interval res,a;
 real r;

 r = expo_RelDiam(li_part(x));
 if (r > -107) // If the rel. diam of x is too great
 {   // the power function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	 res = lx_interval(Inf(x));
	 res = Power_(res,n);
	 a = lx_interval(Sup(x));
	 a = Power_(a,n);
	 res = res | a; // Konvex hull
 }
 else 
	 res = Power_(x,n);

 stagprec = stagsave;
 res = adjust(res);

 return res;
}


// -------------------- Exponential Function -----------------------------

int Poly_N_Exp(const lx_interval &x)
// Calculation of the polynomial degee N:
// Only for internal use in exp function !
{
	lx_interval absx(abs(x));
	int N(0),
			m;
			real S(-53*StagPrec(x)*ln_N[0]),
					  c,
		 D; // S: -53*(precision of x)*ln(2);

		 m = expo_gr(li_part(absx)); 
		 if (m>-1000000)
		 {   // x != 0:
			 c = add_real(m,expo(absx))*ln_N[0];
			 D = c;  // D(0) = m*ln(2);
			 while (D > S) {
				 N++;
				 D = D + c - ln_N[N-1];
			 }
		 }

		 return N;
} // Poly_N_Exp()


lx_interval Exp_(const lx_interval &x) throw()
{
	int ex,n,N;

	real m;
	lx_interval res(1.0),absx,T;
	lx_real S;

	absx = abs(x);
	ex = expo_gr(li_part(absx));
	if (ex>-1000000)
	{ // x != 0
		if (Sup(x) <= -6243314768166103.0)
			res = lx_interval(lx_real(0),lx_real(-Max_Int_R,l_real(minreal)));
		else
		{
			m = add_real(expo(absx),ex); // absx '=' 2^m
			n = (int) _double( (9*ln_N[8] + m*ln_N[0])/ln_N[0]);
			T = x;
			T = adjust(T);
			if (n>0)
				times2pown(T,-n); // T '=' 10^-9
	    // T: reduced argument
			N =  Poly_N_Exp(T); // N: Polynomial degree
	    // Calculation of the N+1 polynomial coefficients
	    // Koff[0], ... , Koff[N]:
			lx_interval *Koff = new lx_interval[N+1];
			Koff[0] = 1.0;
			for (int k=1; k<=N; k++)
				Koff[k] = Koff[k-1]/k;
	    // Horner-evaluation of the polynomial P_N(T):
			res = Koff[N];
			for (int k=N-1; k>=0; k--)
				res = res*T + Koff[k];
	    // Calculation of the absolute approximation error:
			Koff[0] = 1.0;
			for (int k=1; k<=N+1; k++)
				Koff[0] = Koff[0]*k;  // Koff[0] = (N+1)!
			T = lx_interval( Sup(abs(T)) );
			T = (power(T,N+1)/Koff[0])/(1-T);
	    // T: Inclusion of the absolute approximation error
			S = Sup(T);
			res = res + lx_interval(-S,S); // We now considered 
	    // the absolute approximation error of P_N(T).
	    // the following loop is only used if n = 1,2,3,...
			for (int k=1; k<=n; k++)
				res = sqr(res);
			delete[] Koff;
		}
	}
	return res;
} // Exp_(...)   

lx_interval exp(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 40;

 if (stagprec > stagmax) 
	 stagprec = stagmax;
 if (stagprec<3) stagprec = 3;

 lx_interval res,a;
 real r;

 if (Sup(x)>6243314768166065.0)
	 cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF(
				  "lx_interval exp(const lx_interval &x)"));
 r = expo_RelDiam(li_part(x));
 if (r > -107) // If the rel. diam of x is too great
 {   // the exp function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	 res = lx_interval(Inf(x));
	 res = Exp_(res);
	 a = lx_interval(Sup(x));
	 a = Exp_(a);
	 res = lx_interval(Inf(res),Sup(a));
 }
 else 
	 res = Exp_(x);

 stagprec = stagsave;
 res = adjust(res);

 return res;
} // exp()

lx_interval exp2(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 40;
 if (stagprec>stagmax) stagprec = stagmax;
 if (stagprec<3) stagprec = 3;

 lx_interval res;
 res = exp( x*Ln2_lx_interval() );

 stagprec = stagsave;
 res = adjust(res);

 return res;
} // exp2(...)  

lx_interval exp10(const lx_interval &x) throw()
{
	int stagsave = stagprec,
 stagmax = 40;
 if (stagprec>stagmax) stagprec = stagmax;
 if (stagprec<3) stagprec = 3;

 lx_interval res;
 res = exp( x*Ln10_lx_interval() );

 stagprec = stagsave;
 res = adjust(res);

 return res;
} // exp10(...)   

lx_interval pow(const lx_interval &x, const lx_interval &e) throw()
// Calculating an inclusion of x^y;
// If y is an integer value n with 
//      -9007199254740991.0 <= n <= +9007199254740991.0,
// then power(x,n) is used, otherwise y = exp( e*ln(x) ), Inf(x)>0;
// Blomquist, 03.03.2008;
{
	int stagsave = stagprec,
 stagmax = 40;
 lx_interval y;
 lx_real supabse = Sup(abs(e));
 interval z;
 real r;
 bool fertig(false);

 if (point_intv(e))
	 if (supabse < Max_Int_R)
 {
	 z = e;
	 r = Inf(z);
	 if ( Is_Integer(r) )
	 {
		 if (r==1) y = x;
		 else 
			 if (r==0) y = 1;
		 else y = power(x,r);
		 fertig = true;
	 }
 }
 if (!fertig)
 {
	 if (stagprec < stagmax) stagprec++;
	 else stagprec = stagmax;
	 y = exp( e*ln(x) );
	 stagprec = stagsave;
	 y = adjust(y);
 }

 return y;
} // pow(...)

lx_interval xp1_pow_y(const lx_interval &x, const lx_interval &y) throw()
{
	int stagsave = stagprec,
 stagmax = 40;
 if (stagprec>stagmax) stagprec = stagmax;
 lx_interval res, su(1.0+x);

 if (point_intv(su) && Sup(su) == 1.0)
	 res = 1.0;
 else 
	 if (point_intv(y) && Sup(y) == 1.0)
		 res = su;
 else
	 if (point_intv(y) && Sup(y) == 0.0)
		 res = 1.0;
 else
 {
	 if (stagprec < stagmax) stagprec++;
	 else stagprec = stagmax;
	 res = exp(y*lnp1(x));
	 stagprec = stagsave;
	 res = adjust(res);
 }

 return res;
} // xp1_pow_y(...)

int Poly_N_Expm1(const lx_interval &x) 
// Calculation of the polynomial degee N:
// Only for internal use in Expm1(...) 
{
	lx_interval absx(abs(x));
	int N(0), m;
	real S(-53*stagprec*ln_N[0]),
			  c,
	  D; // S: -53*(Praezision von x)*ln(2);
	  m = expo_gr(li_part(absx)); // Neu
	  if (m>-1000000)
	  {
		  D = m + expo(absx); // Calculation of D is uncritical!
		  c = D*ln_N[0];
		  D = (D-1)*ln_N[0];  // D = (m-1)*ln(2);
		  while (D > S) 
		  {
			  N++;
			  D = D + c - ln_N[N];
		  }
	  }
	  return N;
}

lx_interval Expm1(const lx_interval &x)
//  e^x - 1;  |x| <= 1e-7;
//  Only for internal use and for
//  not too wide intervals! 
{
	lx_interval res(0),Fak,D;
	lx_real S;
	int N;
	N = Poly_N_Expm1(x);
    // Calculation of the (N+1) Taylor coefficients
	lx_interval *Koff = new lx_interval[N+1];
	Koff[0] = 1.0;
	for (int k=1; k<=N; k++)
		Koff[k] = Koff[k-1]/(k+1);
    // Evaluation of the Taylor polynimial P_N(x) (Horner)
	res = Koff[N];
	for (int k=N-1; k>=0; k--)
		res = res*x + Koff[k];
    // Calculation of the absolute approximation error:
	Fak = 1.0;
	for (int k=1; k<=N+2; k++)
		Fak = Fak * k;  // Fak = (N+2)!
	D = lx_interval( Sup(abs(x)) );
	D = (power(D,N+1)/Fak)/(1-D); // D: inclusion of the absolute
                                  // approximation error;
	S = Sup(D);
	res = res + lx_interval(-S,S); // Considering the approx. error;
	res = res * x;
	delete[] Koff;
	return res;
}

lx_interval EXPm1(const lx_interval &x) throw()
//  e^x - 1;
//  Only for internal use and for
//  not too wide intervals! 
{
	int stagsave = stagprec,
 stagmax = 40;
 lx_interval res;
 if (Sup(abs(x))<1e-7)
	 res = Expm1(x);
 else 
 {
	 if (stagprec < stagmax) stagprec++;
	 res = Exp_(x)-1;
 }
 stagprec = stagsave;
 res = adjust(res);

 return res;
} // EXPm1(...)

lx_interval expm1(const lx_interval &x) throw()
//  e^x - 1;
{
	int stagsave = stagprec,
 stagmax = 40;
 if (stagprec > stagmax) 
	 stagprec = stagmax;

 lx_interval res,a;
 real r;

 r = expo_RelDiam(li_part(x));
 if (r > -107) // If the rel. diam of x is too great
 {   // the expm1 function is calculated at the
        // boundary points of the interval x. Thus, due to
        // the monotony, the inclusion is given by: 
        // res = [Inf(res),Sup(a)];
	 res = lx_interval(Inf(x));
	 res = EXPm1(res);
	 a = lx_interval(Sup(x));
	 a = EXPm1(a);
	 res = lx_interval(Inf(res),Sup(a));
 }
 else 
	 res = EXPm1(x);

 stagprec = stagsave;
 res = adjust(res);

 return res;
}

// -------------------------- sin(...) --------------------------------------

real pot3_n[20] = {3.0, 9.0, 27.0, 81.0, 243.0, 729.0, 2187.0, 6561.0,
	19683.0, 59049.0, 177147.0, 531441.0, 1594323.0, 4782969.0,
 14348907.0, 43046721.0, 129140163.0, 387420489.0,
 1162261467.0, 3486784401.0 }; 
// pot3_n[k] = 3^(k+1); k = 0,1,2, ... ,19;
// pot3_n[0] = 3.0; ...  pot3_n[19] = 3^20 = 3486784401.0; 

 lx_interval Sin_Rek(const lx_interval &x, int n)
// Only for internal use in sin-function.
 {
	 lx_interval u,v;
	 u = x;
	 for (int k=1; k<=n; k++)
	 { // Rekursion:  u = u*(3-4*u^2);
		 v = sqr(u);
		 times2pown(v,2); // v = 4*sqr(u);
		 u = u*(3-v);
	 }
	 return u;
 }

 int Poly_N_Sin(const lx_interval &x)
// Calculation of the polynomial degee N:
// Only for internal use in sin-function.
 {
	 lx_interval absx(abs(x));
	 int N(0),
			 k;
			 real S(-53*stagprec*ln_N[0]), // S: -53*(actual pecision)*ln(2);
						c,kr,D; 

						k = expo_gr(li_part(absx));
						if (k>-1000000) 
						{   // x != 0
							kr = k + expo(absx); 
							c = (2*kr)*ln_N[0];
							D = c - ln_N[4];  // D = D(0) = 2*kr*ln(2) - ln(3!);
							while (D > S) 
							{
								N++;
								D = D + c - ln_N[2*N] - ln_N[2*N+1]; // Recursion formula
							}
						}  // N=0 if x==0;
						return N;
 }

 l_real floor(const l_real &x)
// Rounding to the next integer number z, with z <= x;
 {
	 l_real y(x),z(0);
	 int p(StagPrec(y)),k;

	 y = y + 0; // Sorting the components.

	 if (expo(y[1])>-1000000)
    // x != 0
	 {
		 k=1;
		 while(expo(y[k])>=53 && k<=p)
		 {
			 z += y[k]; // Addition of the integer parts
			 k++;
		 }
		 if (k<=p)
			 z += std::floor(_double(y[k]));  // Addition of the last
	 }                                        // floor part of y[k].

	 return z;
 } // floor

 lx_interval sin(const lx_interval &x) throw()
// Inclusion of sin(x)
// Blomquist, 11.03.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39,
  n,N,exl;
  const real ln2 = 0.6931471806,
  ln10 = 2.302585093,
  ln3r = 0.9102392264, // 1/ln(3)
  MaxR = 1.71e308,     // < MaxReal
  pi2  = 6.283185309;  // pi2 > 2*pi   
  l_interval xl;
  l_real m;
  real r,ex;
  interval z;
  bool bl(false);

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res(0.0),T,Pi2;
  lx_real S,Si;

  if (Sup(abs(x))>MaxR)
	  res = lx_interval(0,l_interval(-1,1));
  else 
	  if (diam(x)>pi2) // Pi2 > 2*pi;
		  res = lx_interval(0,l_interval(-1,1));
  else // diam(x) <= 6.283185309 && Sup(abs(x)) <= MaxReal/1.05
  {
	  ex = expo(x);
	  exl = expo_gr(li_part(x));
	  if (exl>-1000000)
	  {  // x != 0
		  if (ex<-1080-exl) // sin(x) approx x
			  if (0<=x)
				  res = x;
		  else 
			  if (Inf(x)>0)
		  {
			  Pi2 = lx_interval(Inf(x)); 
			  T = Pi2 * One_m_lx_interval();
			  res = lx_interval(Inf(T),Sup(x));
		  }
		  else 
		  {
			  Pi2 = lx_interval(Sup(x));
			  T = Pi2 * One_m_lx_interval();
			  res = lx_interval(Inf(x),Sup(T));
		  }
		  else
		  {
			  xl = x;
			  r = expo_RelDiam(li_part(x));
			  if (r > -107) // If the rel. diam of x is too great
			  {   // the sin function is calculated using 
                        // the class l_interval
				  xl = sin(xl);
				  res = xl;
			  }
			  else
			  {
				  Pi2 = Pi_lx_interval();
				  times2pown(Pi2,1);  // Pi2 = 2*Pi
				  m = floor( Inf(xl)/Inf(Pi2_l_interval()) ); 
				  T = x - m*Pi2;  // T: reduced argument
				  xl = T;  z = xl;
				  if (diam(z)>Pi2_real)
					  return lx_interval(0,l_interval(-1,1));

				  r = expo( Sup(abs(z)) );
				  n = (int) _double((r*ln2 + 8*ln10)*ln3r);
                        // if (n<=0) then |z|<=10^(-8), 
		        // i.e. a second argument reduction is not necessary!
		        // A second argument reduction exist  <==>  n>0.
				  if (n>0)
					  T = T / pot3_n[n-1]; // T = T / 3^n;
			    // T is the reduced argument;
				  N = Poly_N_Sin(T);  // N: Polynomial degree;

			// Calculation of the N+1 polynomial coefficients
			// Koff[0], ... , Koff[N]: 
				  lx_interval *Koff = new lx_interval[N+1];
				  Koff[0] = 1.0;
				  for (int k=1; k<=N; k++)
				  {
					  m = 2*k;
					  Koff[k] = -Koff[k-1]/(m*(m+1));
				  }
			// Horner evaluation of the polynomial P_N(T^2):
				  res = Koff[N];  Pi2 = sqr(T);
				  for (int k=N-1; k>=0; k--)
					  res = res*Pi2 + Koff[k];
			// Calculation of the absolute approximation error: 
				  Koff[0] = 1.0;
				  for (int k=1; k<=N+1; k++)
				  {
					  m = 2*k;
					  Koff[0] = Koff[0]*m*(m+1);
				  }  // Koff[0] = (2N+3)!
				  Pi2 = lx_interval( Sup(abs(T)) );
				  Pi2 = sqr(Pi2);
				  Pi2 = power(Pi2,N+1)/Koff[0];
			// Pi2: Inclusion of the absolute approx. error
				  S = Sup(Pi2);
				  res = res + lx_interval(-S,S); // We now considered 
			// the absolute approximation error of P_N(T^2).
				  res = res*T; // res: Inclusion of sin(T/3^n)
			// Now calculating an inclusion of sin(T) using
			// the recurrence:
				  res = Sin_Rek(res,n);
				  delete[] Koff;
			  }  
		  }
	  }
  }
  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res); 
  S = Sup(res);  Si = Inf(res);
  if (S>1.0)   { S =  1.0;  bl = true; }
  if (Si<-1.0) { Si = -1.0; bl = true; }
  if (bl)
	  res = lx_interval(Si,S);
    // Now   -1 <= res <= +1   is guaranteed

  return res;
 } // sin()

 lx_interval sin_n(const lx_interval &x, const real& n) throw()
// Inclusion of sin(n*Pi + x)
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;
  lx_interval res;

  if ( !(Is_Integer(n)) ) 
	  cxscthrow(REAL_NOT_ALLOWED("lx_interval sin_n(const lx_interval&, const real&)"));

  res = sin(x);
//    if (n%2 != 0) res = -res;
  if ( !(Is_Integer(n/2)) ) 
	  res = -res;

  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res);   

  return res;
 } // sin_n(...)

 lx_interval cos(const lx_interval &x) throw()
// Inclusion of cos(x)
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;
  lx_interval res;
  lx_real S,Si;
  bool bl(false);

  if (li_part(x)==0) res = 1.0;
  else
  {
	  res = Pi_lx_interval();
	  times2pown(res,-1); // res = Pi/2;
	  res = x + res;
	  res = sin(res);
  }

  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res);   

  S = Sup(res);  Si = Inf(res);
  if (S>1.0)   { S =  1.0;  bl = true; }
  if (Si<-1.0) { Si = -1.0; bl = true; }
  if (bl)
	  res = lx_interval(Si,S);
    // Now   -1 <= res <= +1   is realized

  return res;
 } // cos(...)

 lx_interval cos_n(const lx_interval &x, const real& n) throw()
// Inclusion of cos((n+0.5)*Pi + x)
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;
  lx_interval res;

  if ( !(Is_Integer(n)) ) 
	  cxscthrow(REAL_NOT_ALLOWED("lx_interval cos_n(const lx_interval&, const real&)"));
  res = sin(x);
//    if (n%2 == 0) res = -res;
  if ( (Is_Integer(n/2)) ) 
	  res = -res;

  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res);   

  return res;
 } // cos_n(...)

 lx_interval tan(const lx_interval& x) throw() 
 {
	 lx_interval c,y;

	 if (li_part(x) == 0)
		 y = lx_interval(0,l_interval(0.0));
	 else 
	 {
		 c = cos(x);
		 if (0.0 <= c) 
		 {
			 cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("lx_interval tan(const lx_interval &x)"));
		 }
		 y = sin(x)/c;
	 }
	 return y;
 } // tan()

 lx_interval cot(const lx_interval & x) throw() 
 {
	 lx_interval s,y;

	 s = sin(x);
	 if (0.0 <= s) 
	 {
		 cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF("lx_interval cot(const lx_interval &x)"));
	 }
	 y = cos(x)/s;

	 return y;
 } // cot()

 lx_interval sqrt1px2(const lx_interval &x) throw() 
 {   // Inclusion of sqrt(1+x^2)
	 int ext,c(3210);
	 real ex;
	 lx_interval y(abs(x));
	 ext = expo_gr(li_part(y));
	 ex = expo(y);

	 if (ext>-1000000 && ex>c)
		 y = lx_interval(Inf(x),upper_bnd(Sup(x)));
	 else 
		 if (ext<-1000000) 
			 y = lx_interval(0,l_interval(1));
	 else // x != 0;
		 y = (ex<-c)? 1.0 + lx_interval(lx_real(0),Sup(y)) : sqrt(1+sqr(y));

	 return y;
 } // sqrt1px2()

 int Poly_N_Atan(const lx_interval &x)
// Calculation of the polynomial degee N:
// Only for the internal use in the Atan_-function.
 {
	 lx_interval absx(abs(x));
	 int N(0),
			 ex;
			 real S(-53*stagprec*ln_N[0]), // S: -53*(actual precision)*ln(2);
						c,D;
						ex = expo_gr(li_part(absx));
						if (ex>-1000000) 
						{   // x != 0
							D = ex + expo(absx); // Calculation of D is uncritical!
							c = D*ln_N[0];
							D = 2*D*ln_N[0] - ln_N[1];  // D = D(0) = 2ex*ln(2)-ln(3);
							while (D > S) {
								N++;
								D = (2*N+2)*c - ln_N[2*N+1]; 
							}
						}  // N=0 if x==0;
						return N;
 } // Poly_N_Atan()

 lx_interval Atan_(const lx_interval &x) throw() 
// Inclusion of atan(x) for rather tight intervals x only!
// Only for the internal use in atan(x)!
 { 
	 const int c(4); 
	 int exl, m(0), N;
	 real ex(expo(x));

	 lx_interval res(0.0),T,T2;
	 l_interval xl(li_part(x));
	 bool neg;
	 lx_real Xl;

	 exl = expo_gr(xl);
	 if (exl>-100000)
	 { // x != 0
		 res = x; 
		 if ( abs(Inf(res))>lx_real(2092,8.567562) )
		 { 
			 neg = Sup(xl)<0;
			 res = Pi_lx_interval();
			 times2pown(res,-1); // Pi/2 enclosed by res;
			 if (neg) res = -res;
		 }
		 else 
			 if (ex<-1080-exl) // atan(x) approx x
				 if (Inf(x)>0)
		 {
			 T2 = lx_interval(Inf(x));
			 T = T2 * One_m_lx_interval();
			 res = lx_interval(Inf(T),Sup(x));
		 }
		 else // Sup(x)<0, because x is a rather tight interval
		 {
			 T2 = lx_interval(Sup(x));
			 T = T2 * One_m_lx_interval();
			 res = lx_interval(Inf(x),Sup(T));
		 }
		 else
		 {   // Now with argument reduction (Rothmaier)
			 ex = ex + exl;
			 if (ex>1) 
				 m = (int) _double( (ln_N[0]+c*ln_N[8])/ln_N[0] );
			 else
				 if (ex>-10000)
					 m = (int) _double((c*ln_N[8] + ex*ln_N[0])/ln_N[0]);
			 T = x;
		// Now with m >= 0 transformations:
			 for (int k=1; k<=m; k++) 
				 T = T /(1 + sqrt1px2(T)); // T: reduced argument
		// Now calculating the polynomial degree:
			 N = Poly_N_Atan(T);
		// Calculation of the N+1 polynomial coefficients
		// Koff[0], ... , Koff[N]:
			 lx_interval *Koff = new lx_interval[N+1];
			 Koff[0] = 1.0;
			 for (int k=1; k<=N; k++)
			 {
				 Koff[k] = lx_interval(0,l_interval(1))/ (2*k+1);
				 if (k%2 != 0) Koff[k] = -Koff[k];
			 }
		// Horner-evaluation of the polynomial P_N(T2):
			 res = Koff[N]; 
			 T2 = sqr(T);
			 for (int k=N-1; k>=0; k--)
				 res = res*T2 + Koff[k]; // res = P_N(T2)
		// Calculation of the absolute approximation error:
			 T2 = lx_interval( Sup(abs(T)) );
			 T2 = sqr(T2);
			 T2 = power(T2,N+1)/(2*N+3);
		// T2: Inclusion of the absolute approximation error
		// Implementing the approximation error:
			 Xl = Sup(T2);
			 res += lx_interval(-Xl,Xl);
			 res = T*res;
			 if (m>0) times2pown(res,m);
			 delete[] Koff;
		 }
	 }

	 return res;
 } // Atan_()

 lx_interval atan(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the atan function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Atan_(res);
	  a = lx_interval(Sup(x));
	  a = Atan_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Atan_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // atan()

// ---------------------- sqrt(1-x^2) -----------------------------------

 lx_interval Sqrt1mx2_(const lx_interval &x) throw() 
// Inclusion of sqrt(1-x^2) for rather tight intervals x;
// For the internal use only!!
// Blomquist, 08.04.2008;
 { 
	 int exl;

	 lx_interval res(1);
	 l_interval li(li_part(x));
	 real ex(expo(x));

	 exl = expo_gr(li);
	 if (exl>-100000)
	 { // x != 0
		 if (ex<=-2098-exl)
		 {   // Construction of the inclusion 
            // interval [1 - 2^-2097,1], only possible if p>=2; 
			 res = lx_interval(-1023,l_interval(comp(0.5,1024)));
			 li = li_part(res);
			 li += 0; // To get the actual precision and sorting the li[k];
			 li[StagPrec(li)] = -minreal; // minreal = 2^-1074;
			 res = lx_interval(-1023,li); 
            // res = 2^-1023*[2^1023 - 2^-1074,2^1023] = [1-2^-2097,1];
		 }
		 else 
			 res = (ex>=-exl)? sqrt( (1-x)*(1+x) ) : sqrt( 1-sqr(x) );
	 }
	 return res;
 } // Sqrt(1-x^2) 

 lx_interval sqrt1mx2(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) stagprec = stagmax;
  if (stagprec==1) stagprec++; // stagprec >= 2;

  lx_interval res(abs(x)),a;
  real r;

  r = expo_RelDiam(li_part(res));
  if (r > -107) // If the rel. diam of x is too great
  {   // the sqrt1mx2 function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  a = lx_interval(Inf(res));
	  a = Sqrt1mx2_(a);

	  res = lx_interval(Sup(res));
	  res = Sqrt1mx2_(res);

	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Sqrt1mx2_(res);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // sqrt1mx2

// ------------------- sqrt(x^2-1) ------------------------------------

 lx_interval sqrtx2m1(const lx_interval &x) throw() 
// Inclusion of sqrt(x^2-1)
// Blomquist, 22.04.08;
 { 
	 int stagsave = stagprec,
  stagmax = 39;
  real ex(expo(x));
  if (stagprec>stagmax) stagprec = stagmax;

  lx_interval u(abs(x)),res;
  l_interval li(li_part(u));
  lx_real Infu(Inf(u));

  if (Infu < 1)
  {
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF(
					"lx_interval sqrtx2m1(const lx_interval &x)"));
  }
    // Now: Inf(u) >= 1;
  if (ex > 1604-expo_gr(li))
	  res = u + lx_interval(lx_real(-1),lx_real(0));
  else 
  {
	  res = sqrt( (u-1)*(u+1) );
	  if (Inf(res)<=0 && Infu>1)
	  {
		  u = lx_interval(Infu);
		  u = u + 1;
		  times2pown(u,-2097);
		  u = sqrt(u);
		  res = lx_interval( Inf(u), Sup(res) );
	  }
  }

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // sqrt(x^2-1)

// ------------------- asin(x) ----------------------------------------

 lx_interval Asin_(const lx_interval& x) throw()
// Inclusion of asin(x) for rather tight intervals x only!
// Only for the internal use in asin(x)!
 {
	 lx_interval y,pihalbe,h;
	 bool neg;

	 l_interval xl;
	 interval z;
	 xl = x;
	 z = xl;

	 real supabsz = Sup(abs(z)),
		infz = Inf(z),
		supz = Sup(z);

	pihalbe = Pi_lx_interval();
	times2pown(pihalbe,-1); // pihalbe = pi/2;

	if (infz == supz && supz == 0)
		y = 0;
	else 
		if (infz == supz && supabsz == 1)
	{
		if (supz == 1.0)
			y = pihalbe;
		else y = -pihalbe;
	}
	else
	{
		if (0<=z)
			y = atan(x/sqrt1mx2(x));
		else 
		{
			y = x;
			neg = Sup(z)<0;
			if (neg) y = -y;
			if (supabsz < 0.75)
				y = atan(y/sqrt1mx2(y));
			else 
				y = pihalbe - atan(sqrt1mx2(y)/y);
			if (neg) y = -y;
		}
	}

	return y;
} // Asin_(x)

 lx_interval asin(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  real r;

  if (Inf(x)<-1 || Sup(x)>1) 
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval asin(const lx_interval& x)"));

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the asin function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Asin_(res);
	  a = lx_interval(Sup(x));
	  a = Asin_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Asin_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // asin()

// ------------------- acos(x) ----------------------------------------

 lx_interval Acos_(const lx_interval& x) throw()
// Inclusion of acos(x) for rather tight intervals x only!
// Only for the internal use in acos(x)!
// Blomquist, 24.04.2008;
 {
	 lx_interval y,pihalbe,pi;
	 bool neg;
	 l_interval xl;
	 interval z;

	 xl = x;
	 z = xl;

	 real supabsz = Sup(abs(z)),
		infz = Inf(z),
	supz = Sup(z);

	pi = Pi_lx_interval();
	pihalbe = pi;
	times2pown(pihalbe,-1); // pihalbe = pi/2;

	if (infz == supz && supz == 0)
		y = pihalbe;
	else 
		if (infz == supz && supabsz == 1)
	{
		if (supz == 1.0)
			y = 0;
		else y = pi;
	}
	else
	{
		if (0<=z)
			y = pihalbe - atan(x/sqrt1mx2(x));
		else 
		{
			y = x;
			neg = Sup(z)<0;
			if (neg) y = -y;
			if (supabsz < 0.25)
				y = pihalbe - atan(y/sqrt1mx2(y));
			else 
				y = atan(sqrt1mx2(y)/y);
			if (neg) y = pi - y;
		}
	}

	return y;
} // Acos_(x)

 lx_interval acos(const lx_interval &x) throw()
// Blomquist, 24.04.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  real r;

  if (Inf(x)<-1 || Sup(x)>1) 
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval acos(const lx_interval& x)"));

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the acos function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Acos_(res);
	  a = lx_interval(Sup(x));
	  a = Acos_(a);
	  res = lx_interval(Inf(a),Sup(res));
  }
  else 
	  res = Acos_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // acos()


// ------------------- acot(x) ----------------------------------------

 lx_interval Acot_(const lx_interval& x) throw()
// codomain of acot(x): (0,pi);
// acot(x) is continuous at x=0;
// Blomquist, 25.04.2008; 
 {
	 int exl;
	 real ex(expo(x));
	 lx_interval y,pihalbe,pi;
	 bool neg;

	 l_interval xl(li_part(x));
	 exl = expo_gr(xl);

	 pi = Pi_lx_interval();
	 pihalbe = pi;
	 times2pown(pihalbe,-1); // pihalbe = pi/2;

	 if (0<=xl)
		 y = pihalbe - atan(x);
	 else 
	 {
		 y = x;
		 neg = Sup(xl) < 0;
		 if (neg) y = -y;
		 if (ex<=-exl) 
			 if (ex<-exl-1580)
				 y = pihalbe;
		 else 
			 y = pihalbe - atan(y);
		 else y = atan(1/y);
		 if (neg) y = pi - y;
	 }

	 return y;
 } // Acot_(x)

 lx_interval acot(const lx_interval &x) throw()
// Blomquist, 24.04.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the acot function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Acot_(res);
	  a = lx_interval(Sup(x));
	  a = Acot_(a);
	  res = lx_interval(Inf(a),Sup(res));
  }
  else 
	  res = Acot_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // acot()

// -------------------------- sinh(...) --------------------------------------

 lx_interval Sinh_(const lx_interval &x) throw()
// Inclusion of sinh(x) for sufficiently tight intervals x
 {
	 int exl;
	 real ex;
	 lx_interval res(0.0),T;
	 l_interval xl(li_part(x));
	 bool neg;

	 ex = expo(x);
	 exl = expo_gr(xl);

	 if (0<=xl) // Null in x enthalten!
	 {
		 T = expm1(-x);
		 res = -T*(1+1/(T+1));
		 times2pown(res,-1);
	 }
	 else // Null ist nicht in x enthalten!
		 if (ex<=-exl)
	 {   // |abs(x)| < 1;
		 T = expm1(-x);
		 res = -T*(1+1/(T+1));
		 times2pown(res,-1);
	 }
	 else 
	 {   // |abs(x)| >= 1;
		 T = x;
		 neg = Inf(xl)<0;
		 if (neg) T = -T;

		 if (ex>12-exl)
		 {
			 res = exp(T);
			 times2pown(res,-1);
			 res += l_interval(-0.5,0);
		 }
		 else
		 {
			 res = exp(T);
			 res = res - 1/res;
			 times2pown(res,-1);
		 }
		 if (neg) res = -res;
	 }

	 return res;
 } // Sinh_(...)

 lx_interval sinh(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the sinh function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Sinh_(res);
	  a = lx_interval(Sup(x));
	  a = Sinh_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Sinh_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // sinh()

// -------------------------- cosh(...) ------------------------------------

 lx_interval Cosh_(const lx_interval& x) throw()
// codomain of cosh(x): [1,+infty);
// This function is only used for tight intervals x,
// with Inf(x)>=0;
// Blomquist, 26.04.2008; 
 {
	 lx_interval y;
	 lx_real Si;
	 l_interval xl;
	 interval z;

	 xl = x;
	 z = xl;

	 y = exp(x);
	 if (Inf(z)>4096)
	 {
		 times2pown(y,-1);
		 y += l_interval(0,0.5);
	 }
	 else 
	 {
		 y = y + 1/y;
		 times2pown(y,-1);
	 }
	 Si = Inf(y);
	 if (Si<1.0)  y = lx_interval(lx_real(1.0),Sup(y));

	 return y;
 } // Cosh_(x)

 lx_interval cosh(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a,y(abs(x));
  real r;

  r = expo_RelDiam(li_part(y));
  if (r > -107) // If the rel. diam of x is too great
  {   // the cosh function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(y));
	  res = Cosh_(res);
	  a = lx_interval(Sup(y));
	  a = Cosh_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Cosh_(y);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // cosh()

// -------------------------- tanh(...) ------------------------------------

 lx_interval Tanh_(const lx_interval& x) throw()
 {
// Blomquist, 27.04.2008;
	 int exl;
	 real ex(expo(x));
	 lx_interval y(x),t;
	 lx_real Iy;
	 bool neg;
	 l_interval xl(li_part(x));
	 exl = expo_gr(xl);

	 if (0<=xl)
	 {
		 times2pown(y,1); // y = 2x;
		 y = expm1(-y);   // y = exp(-2x) - 1;
		 y = -y/(2+y);
	 }
	 else 
	 {   // 0 not included by x.
		 if (ex<=-exl) // 0<|x|<1;
		 {
			 times2pown(y,1); // y = 2x;
			 y = expm1(-y);   // y = exp(-2x) - 1;
			 y = -y/(2+y);
		 }
		 else
		 {   // |x| >= 1;
			 neg = Sup(xl)<0;
			 if (neg) y = -y; // Inf(y)>0, bzw. y = abs(x);
			 Iy = Inf(y);
			 ex = expo(Iy);  exl = expo_gr(lr_part(Iy));
			 if (ex >= -exl+13)
				 y = lx_interval(Inf(One_m_lx_interval()),lx_real(1));
			 else
			 {
				 times2pown(y,1); // y = 2x;
				 y = exp(-y);     // y = exp(-2x);
				 y = (1-y) / (1+y);
			 }
			 if (neg) y = - y;
		 }
	 }
 
	 return y;
 } // Tanh_(x)

 lx_interval tanh(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the tanh function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Tanh_(res);
	  a = lx_interval(Sup(x));
	  a = Tanh_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Tanh_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // tanh()

// -------------------------- coth(...) ------------------------------------

 lx_interval Coth_(const lx_interval& x) throw()
// This function is only used for tight intervals x,
// with Inf(x) > 0;
// Blomquist, 27.04.2008;
 {
	 int exl;
	 lx_interval y(x);
	 real ex(expo(x)); 
	 lx_real Iy;
	 l_interval xl(li_part(x));

	 exl = expo_gr(xl);

	 if (ex<=-exl) // 0<|x|<1;
	 {
		 times2pown(y,1); // y = 2x;
		 y = expm1(-y);   // y = exp(-2x) - 1;
		 y = -(2+y)/y;
	 }
	 else
	 {   // |x| >= 1;
		 Iy = Inf(y);
		 ex = expo(Iy);  exl = expo_gr(lr_part(Iy));

		 if (ex>=-exl+13)
			 y = lx_interval(lx_real(1),Sup(One_p_lx_interval()));
		 else
		 {
			 times2pown(y,1); // y = 2x;
			 y = exp(-y);     // y = exp(-2x);
			 y = (1+y) / (1-y);
		 }
	 }
 
	 return y;
 } // coth(x)

 lx_interval coth(const lx_interval &x) throw()
// Blomquist, 27.04.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  l_interval xl(li_part(x));
  if (0<=xl) 
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval coth(const lx_interval& x)"));

  lx_interval res,a,y(x);
  real r;
  bool neg(Sup(xl)<0);

  if (neg) y = -y; // Inf(y)>0;

  r = expo_RelDiam(xl);
  if (r > -107) // If the rel. diam of x is too great
  {   // the coth function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(y));
	  res = Coth_(res);
	  a = lx_interval(Sup(y));
	  a = Coth_(a);
	  res = lx_interval(Inf(a),Sup(res));
  }
  else 
	  res = Coth_(y);
  if (neg) res = -res;

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // coth()

// -------------------------- sqrt(1+x)-1 ------------------------------------

 lx_interval Sqrtp1m1_(const lx_interval& x) throw()
// sqrtp1m1(x) calculates an inclusion of sqrt(x+1)-1;
// Blomquist, 25.08.07;
 {
	 int exl;
	 lx_interval y(0),tmp;
	 l_interval xl(li_part(x));
	 interval z;
	 real ex(expo(x));
	 const real c = 0.1;

	 exl = expo_gr(xl);
	 if (exl>-1000000) // x != 0;
	 {
		 if (ex>-exl+500) 
			 y = sqrt(1+x) - 1;
		 else
		 {
			 xl = x;  // xl is an inclusion of x,
			 z = xl;  // and z is an inclusion of x;
			 tmp = x+1;
			 y = (z<=interval(-c,c))? x/(sqrt(tmp)+1) : sqrt(tmp)-1;
		 }
	 }

	 return y; 
 } // Sqrtp1m1_

 lx_interval sqrtp1m1(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 30;

  stagprec++;
  if (stagprec>stagmax) stagprec = stagmax;

  lx_interval res,a;
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the sqrtp1m1 function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Sqrtp1m1_(res);
	  a = lx_interval(Sup(x));
	  a = Sqrtp1m1_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Sqrtp1m1_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // sqrtp1m1()


// -------------------------- asinh(x) ------------------------------------

 lx_interval Asinh_(const lx_interval &x) throw()
// Inclusion of asinh(x)
 {
	 int exl;
	 lx_interval res,u;
	 lx_real S1,S2;
	 l_interval xl(li_part(x));
	 interval z;
	 real ex;
	 bool neg;

	 ex = expo(x);
	 exl = expo_gr(xl);

	 if (0<=xl) // Null in x enthalten!
	 {
		 xl = x; // Falls hier Fehlermeldung, ist das OK,
		 z = xl; // denn sonst ist x viel zu breit!!
		 res = asinh(z); // Bis hier OK!!!
	 }
	 else // Null ist nicht in x enthalten!
	 {
		 neg = (Inf(xl)<0);
		 u = (neg)? -x : x; // Inf(u)>0
		 if (ex<=-exl-790)
		 {  
			 std::cout << "Hallo: (1-u)*(1+u)*res" << std::endl; 
			 res = u;
			 S2 = Sup(res);
			 S1 = Inf( (1-u)*(1+u)*res );
			 res = lx_interval(S1,S2);
		 }
		 else 
			 res = (ex<=-exl)? lnp1(u+sqrtp1m1(sqr(u))) : ln(u+sqrt1px2(u));

		 if (neg) res = -res;
	 }

	 return res;
 } // Asinh_(...)

 lx_interval asinh(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec>stagmax) 
	  stagprec = stagmax;

  lx_interval res,a; 
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the asinh function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Asinh_(res);
	  a = lx_interval(Sup(x));
	  a = Asinh_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Asinh_(x);
  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // asinh()

// -------------------------- acosh(...) -------------------------------

 lx_interval acosh(const lx_interval &x) throw()
// Inclusion of acosh(x)
 {
	 int stagsave = stagprec,
  stagmax = 39,
  exl;
  if (stagprec > stagmax) 
	  stagprec = stagmax;
  lx_interval res(0.0),x_1;
  lx_real S1,S2;
  real ex;

  if (Inf(x)<1) 
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval acosh(const lx_interval& x)"));
  x_1 = x-1;
  ex = expo(x_1);
  l_interval x_1l(li_part(x_1));
  exl = expo_gr(x_1l);

  if (exl>-100000) // x != [0,0]
  {
     if (ex <= -exl-1600) 
     {
	  res = x;
	  times2pown(res,1);   // res = 2x;
	  res = sqrt(res - 2); // res = sqrt(2x - 2);
	  S2 = Sup(res);                       
	  res = res * (2-x); // Dieser Trick spart Division durch 12 !!
	  S1 = Inf(res);
	  res = lx_interval(S1,S2);
     }
     else 
	  res = (ex<=-exl)? lnp1( x_1 + sqrt(x_1*(2+x_1)) ) :
			  ln( x+sqrtx2m1(x) );
  }
  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res);   

  return res;
 } // acosh(...)

// ------------------------- acoshp1(...) -----------------------------

 lx_interval Acoshp1_(const lx_interval &x) throw()
// This function is only used for tight intervals x,
// with Inf(x) >= 0;
// Blomquist, 03.05.2008;
 {
	 int exl;
	 lx_interval res(0.0);
	 lx_real S1,S2;
	 real ex;

	 ex = expo(x);
	 l_interval xl(li_part(x));
	 exl = expo_gr(xl);

	 if (exl>-100000) {
		 if (ex <= -exl-1600) 
	 {
		 res = x;
		 times2pown(res,1);   // res = 2x;
		 res = sqrt(res); // res = sqrt(2x);
		 S2 = Sup(res);                       
		 res = res * (1-x); 
		 S1 = Inf(res);
		 res = lx_interval(S1,S2);
	 }
	 else 
		 if (ex <= -exl)
			 res = lnp1( x + sqrt(x*(2+x)) );
	 else  // res = acosh(1+x);
	 {
		 res = 1+x;
		 res = ln(res + sqrtx2m1(res));
	 }
   }
	 return res;
 } // Acoshp1_(...)

 lx_interval acoshp1(const lx_interval &x) throw()
// Inclusion of acosh(1+x);
// Blomquist, 03.05.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec>stagmax) 
	  stagprec = stagmax;
  l_interval lix(li_part(x));
  if (Inf(lix)<0) 
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval acoshp1(const lx_interval& x)"));
  lx_interval res,a; 
  real r;

  r = expo_RelDiam(lix);
  if (r > -107) // If the rel. diam of x is too great
  {   // the acoshp1 function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Acoshp1_(res);
	  a = lx_interval(Sup(x));
	  a = Acoshp1_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Acoshp1_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // acoshp1()

// ------------------------ atanh(...) ----------------------------------

 lx_interval Atanh_(const lx_interval &x) throw()
// Inclusion of atanh(x) only for sufficiently small
// intervals x. Only for the internal use!
// Blomquist, 03.05.2008;
 {
	 int exl;
	 lx_interval res(0.0);
	 real ex(expo(x));
	 l_interval xl(li_part(x));

	 exl = expo_gr(xl);

	 if (0<=xl) // Zero included by x
	 {
		 res = x/(1-x);
		 times2pown(res,1);
		 res = lnp1(res);
	 }
	 else // Zero not included by x!
		 if (ex <= -2 - exl)
	 {  
		 res = x/(1-x);
		 times2pown(res,1);
		 res = lnp1(res);
	 }
	 else 
		 res = ln( (1.0+x)/(1.0-x) );
	 times2pown(res,-1);

	 res = adjust(res);   

	 return res;
 } // Atanh_(...)

 lx_interval atanh(const lx_interval &x) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec>stagmax) 
	  stagprec = stagmax;

  lx_interval res,a; 
  real r;

  if (Inf(x)<=-1 || Sup(x)>=1) 
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval atanh(const lx_interval& x)"));

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the asinh function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Atanh_(res);
	  a = lx_interval(Sup(x));
	  a = Atanh_(a);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Atanh_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // atanh()

// ------------------------ atanh1m(...) ---------------------------------

 lx_interval Atanh1m_(const lx_interval &x) throw()
// Inclusion of atanh(1-x); 0<x<2;
// Only for not too wide intervals x;
// For the internal use only!
// Blomquist, 07.05.2008;
 {
	 lx_interval res(0.0);
	 res = ln(2/x-1);
	 times2pown(res,-1); 

	 return res;
 }  // Atanh1m_(...)

 lx_interval atanh1m(const lx_interval &x) throw()
// Blomquist, 07.05.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax)
	  stagprec = stagmax;

  if (Inf(x)<=0 || Sup(x)>=2)
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval atanh1m(const lx_interval& x)"));

  lx_interval res,a;
  real r;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the acot function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Atanh1m_(res);
	  a = lx_interval(Sup(x));
	  a = Atanh1m_(a);
	  res = lx_interval(Inf(a),Sup(res));
  }
  else 
	  res = Atanh1m_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // atanh1m()

// ------------------------ atanhm1p(...) ------------------------------

 lx_interval atanhm1p(const lx_interval &x) throw()
// Inclusion of atanh(-1+x) , 0<x<2;
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;
  lx_interval res;

  res = -atanh1m(x);

  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res);   

  return res;
 } // atanhm1p

// ------------------------ acoth(...) ---------------------------------

 lx_interval Acoth_(const lx_interval& x) throw()
//  acoth(x), x=[x1,x2], x1>1; Calculating 
//  inclusions only for not too wide intervals;
//  Only for the internal use in acoth(x) !
//  Blomquist, 13.05.2008;
 {
	 lx_interval res;

	 res = lnp1(2/(x-1));
	 times2pown(res,-1); 

	 return res;
 } // Acoth_

 lx_interval acoth(const lx_interval &x) throw()
// Blomquist, 13.05.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) stagprec = stagmax;

  lx_interval res,a,u;
  l_interval xl(li_part(x));
  bool neg;
  real r;

  res = lx_interval(0,l_interval(-1,1));
  if ( (Inf(x) <= res) || ((Sup(x)) <= res) || res <= x)
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval acoth(const lx_interval& x)"));
  neg = Inf(xl)<0;
  if (neg) u = -x;
  else u = x; // u = [u1,u2], u1 > 1;

  r = expo_RelDiam(xl);
  if (r > -107) // If the rel. diam of x is too great
  {   // the acoth function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(u));
	  res = Acoth_(res);
	  a = lx_interval(Sup(u));
	  a = Acoth_(a);
	  res = lx_interval(Inf(a),Sup(res));
  }
  else 
	  res = Acoth_(u);
  if (neg) res = -res;

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // acoth()

// ------------------------ acoth(1+x) ---------------------------------

 lx_interval Acothp1_(const lx_interval& x) throw()
//  arcoth(1+x), x>0; 
//  only for not too wide intervals x and
//  only for the internal use in acothp1(...)
//  Blomquist, 14.05.2008;
 {
	 lx_interval res;

	 res = lnp1(2/x);
	 times2pown(res,-1); 
 
	 return res;
 } // Acothp1_

 lx_interval acothp1(const lx_interval &x) throw()
// Blomquist, 14.05.2008;
 {
	 int stagsave = stagprec,
  stagmax = 39;

  if (stagprec > stagmax) 
	  stagprec = stagmax;

  lx_interval res,a;
  l_interval xl(li_part(x));
  real r;

  if (Inf(xl) <= 0.0)
	  cxscthrow(ERROR_LINTERVAL_STD_FKT_OUT_OF_DEF
			  ("lx_interval acothp1(const lx_interval& x)"));

  r = expo_RelDiam(xl);
  if (r > -107) // If the rel. diam of x is too great
  {   // the acot function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Acothp1_(res);
	  a = lx_interval(Sup(x));
	  a = Acothp1_(a);
	  res = lx_interval(Inf(a),Sup(res));
  }
  else 
	  res = Acothp1_(x);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // acothp1()

// ------------------------ acothm1m(...) ------------------------------

 lx_interval acothm1m(const lx_interval &x) throw()
// Inclusion of acoth(-1-x) , x>0;
 {
	 int stagsave = stagprec,
  stagmax = 39;
  if (stagprec > stagmax) 
	  stagprec = stagmax;
  lx_interval res;

  res = -acothp1(x);

  stagprec = stagsave; // Restore the old stagprec value
  res = adjust(res);   

  return res;
 } // acothm1m

// ------------------------ sqrt(x^2+y^2) ------------------------------

 lx_interval Sqrtx2y2_(const lx_interval &x, const lx_interval &y) throw()
// Inclusion of sqrt(x^2+y^2);
// x = [x1,x2], x1>=0;  y = [y1,y2], y1>=0;  x2>=y2;
// For not too wide intervals x,y;
// Blomquist, 15.05.2008;
 {
	 const int c1=4000,
  c2=-100000; 
  int exxl,exyl;
  real exx,exy;
  lx_interval res,ax,ay;
  lx_real R;
  R = Sup(y);
  exy = expo(R);
  exyl = expo_gr( lr_part(R) );

  R = Inf(x);
  exx = expo(R);
  exxl = expo_gr( lr_part(R) );

  if (exyl > c2) // y2 > 0;
  {
	  exxl = exxl-exyl-1051;    // Without integer-Overflow 
	  if (exy < exx + (exxl-1)) // (exxl-1) for security
		  res = x * lx_interval(lx_real(1.0),Sup(One_p_lx_interval()));
	  else 
	  {
		  res = x;
		  ax = x;  ay = y;
		  if (exx<-9007199254735000.0) 
		  {
			  times2pown(ax,c1);
			  times2pown(ay,c1);
		  };
		  res = res * sqrt1px2(ay/ax);
	  }	
  }
  else // y2 = 0;
	  res = x;

  return res;
 }  // Sqrtx2y2_(...)

 lx_interval sqrtx2y2(const lx_interval &x, const lx_interval &y) throw()
 {
	 int stagsave = stagprec,
  stagmax = 30,
  exyl,exxl;
  if (stagprec>stagmax) 
	  stagprec = stagmax;

  const int c1=4000,
  c2=-100000;
  lx_interval res,ax,ay,u,v; 
  lx_real R;
  real rx,ry;

  ax = abs(x);    ay = abs(y);
  if ( Sup(ay) > Sup(ax) ) 
  {
	  res = ax; ax = ay; ay = res;
  }   // Sup(ax) >= Sup(ay);

  if (0<=ax) // 0 is included by ax
    // The inclusion res is calculated as the convex hull
    // of the inclusions at the boundary bounds of ax,ay.
  {
    // Inclusion of the upper bound:
	  R = Sup(ay);
	  v = lx_interval(R);
	  ry = expo(R);
	  exyl = expo_gr( lr_part(R) );

	  R = Sup(ax);
	  u = lx_interval(R);
	  rx = expo(R);
	  exxl = expo_gr( lr_part(R) );

	  if (exyl>c2)
	  {
		  exxl = exxl-exyl-1051;  // Without integer overflow
		  if (ry < rx + (exxl-1)) // (exxl-1) for security
			  res = u * lx_interval(lx_real(1.0),Sup(One_p_lx_interval()));
		  else
		  {
			  res = u;
			  if (rx < -9007199254735000.0)
			  {
				  times2pown(u,c1);
				  times2pown(v,c1);
			  }
			  res = res * sqrt1px2(v/u);
		  }
	  }
	  else // v=0;
		  res = u;
       // res is the inclusion of the function value
       // at the upper bounds of ax,ay;
       // Remember: Inf(ax)=0, so Inf(ay) is the lower bound
       // of the function value at the lower bounds of ax,ay;
       // Thus, the inclusion of sqrt(ax^2 + ay^2) is given by:
	  res = lx_interval(Inf(ay),Sup(res));
  }
    else // Inf(ax)>0:
	 {
		 rx = expo_RelDiam(li_part(ax));
		 ry = expo_RelDiam(li_part(ay));
		 if (rx > -107 || ry > -107) // If the rel. diam of ax or ay 
		 {   // is too great the sqrtx2y2 function is calculated at 
          // the boundary points of the intervals ax, ay. Thus, the
          // inclusion res is the convex hull of these intervals:
			 u = lx_interval(Inf(ax));
			 v = lx_interval(Inf(ay));
			 res = Sqrtx2y2_(u,v);
			 u = lx_interval(Sup(ax));
			 v = lx_interval(Sup(ay));
			 u = Sqrtx2y2_(u,v);
			 res = lx_interval(Inf(res),Sup(u));
		 }
		 else 
			 res = Sqrtx2y2_(ax,ay);
	 }

	 stagprec = stagsave;
	 res = adjust(res);
	 return res;
 } // sqrtx2y2

// ---------------------- ln(sqrt(x^2+y^2)) ----------------------------

 lx_interval Ln_sqrtx2y2_(const lx_interval &x, const lx_interval &y) throw()
// Inclusion of ln(sqrt(x^2+y^2));
// Blomquist, 19.05.2008;
 {
	 lx_interval res,ax,ay;
	 lx_real Rx,Ry;

	 ax = abs(x);    ay = abs(y);
	 Rx = Sup(ax);   Ry = Sup(ay);
	 if (Ry>Rx) 
	 {
		 Rx = Ry; // Rx ist maximaler Wert
		 res = ay;   ay = ax;   ax = res; // ax ist maximaler Wert
	 }

	 if (Rx>100) res = ln( sqrtx2y2(ax,ay) );
	 else 
		 if (Rx<0.25)
			 res = ln( sqrtx2y2(ax,ay) );
	 else
	 {
	    // Calculating:  0.5*ln(x*x + y*y)
            //             = 0.5*ln(1 + (x*x + y*y - 1))
            //             = 0.5*ln(1 + ((x-1)*(x+1) + y*y))
		 if (ay==0.0)
			 res = ln(ax);  // Without Overflow!
		 else
		 {
			 res = lnp1((ax-1)*(ax+1) + ay*ay);
			 times2pown(res,-1);
		 }
	 }

	 return res;
 }  // Ln_sqrtx2y2_(...)

 lx_interval ln_sqrtx2y2(const lx_interval &x, const lx_interval &y) throw()
 {
	 int stagsave = stagprec,
  stagmax = 30;
  if (stagprec>stagmax) stagprec = stagmax;

  lx_interval res,ax,ay,u,v; 
  real rx,ry;

  ax = abs(x);    ay = abs(y);
  if ( Sup(ay) > Sup(ax) ) 
  {
	  res = ax; ax = ay; ay = res;
  }   // Sup(ax) >= Sup(ay);

  rx = expo_RelDiam(li_part(ax));
  ry = expo_RelDiam(li_part(ay));
  if (rx > -107 || ry > -107) // If the rel. diam of ax or ay 
  {   // is too great the ln_sqrtx2y2 function is calculated at 
          // the boundary points of the intervals ax, ay. Thus, the
          // inclusion res is the convex hull of these intervals:
	  u = lx_interval(Inf(ax));
	  v = lx_interval(Inf(ay));
	  res = Ln_sqrtx2y2_(u,v);

	  u = lx_interval(Sup(ax));
	  v = lx_interval(Sup(ay));
	  u = Ln_sqrtx2y2_(u,v);

	  res = lx_interval(Inf(res),Sup(u));
  }
  else 
	  res = Ln_sqrtx2y2_(ax,ay);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // ln_sqrtx2y2

// -------------------------- n-th root --------------------------------

 lx_interval Sqrt_(const lx_interval& x, int n) throw()
 {   // Inclusion of the n-th root, based on the pow(...) function
    // n>=2; For not too wide intervals x;
    // and for the internal use in sqrt(x,n) only!
    // Blomquist, 20.05.2008;
	 lx_interval res;
	 lx_real Sx;

	 if (eq_zero(Inf(x)))
	 { 
		 Sx = Sup(x);
		 if (gr_zero(Sx))
		 {
			 res = pow(lx_interval(Sx),1/lx_interval(n));
			 res = lx_interval(lx_real(0),Sup(res));
		 }
		 else res = 0.0;
	 }
	 else
		 res = pow(x,1/lx_interval(n));

	 return res;
 } // Sqrt_(...)

 lx_interval sqrt(const lx_interval &x, int n) throw()
 {
	 int stagsave = stagprec,
  stagmax = 39,
  k(n);
  if (stagprec>stagmax) stagprec = stagmax;

  lx_interval res,a; 
  real r;

  if (k<2) k=2;

  r = expo_RelDiam(li_part(x));
  if (r > -107) // If the rel. diam of x is too great
  {   // the sqrt(x,n) function is calculated at the
        // boundary points of the interval x. Thus, the
        // inclusion res is the convex hull of these intervals:
	  res = lx_interval(Inf(x));
	  res = Sqrt_(res,k);
	  a = lx_interval(Sup(x));
	  a = Sqrt_(a,k);
	  res = lx_interval(Inf(res),Sup(a));
  }
  else 
	  res = Sqrt_(x,k);

  stagprec = stagsave;
  res = adjust(res);

  return res;
 } // sqrt()


// ---------------- Interval constants in high accuracy --------------------

 static real CXSC_Pi[40]; // CXSC_Pi[0], ... ,CXSC_Pi[39]
 static bool CXSC_Pi_initialized = false;

 lx_interval Pi_lx_interval() throw()
// Inclusion of Pi, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;  // alt
  if (!CXSC_Pi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1921FB54442D18e7FC";
	  str >> CXSC_Pi[0];
	  str = "+11A62633145C07e7C6";
	  str >> CXSC_Pi[1];
	  str = "-1F1976B7ED8FBCe78E";
	  str >> CXSC_Pi[2];
	  str = "+14CF98E804177De758";
	  str >> CXSC_Pi[3];
	  str = "+131D89CD9128A5e722";
	  str >> CXSC_Pi[4];
	  str = "+10F31C6809BBDFe6E8";
	  str >> CXSC_Pi[5];
	  str = "+1519B3CD3A431Be6B1";
	  str >> CXSC_Pi[6];
	  str = "+18158536F92F8Ae67A";
	  str >> CXSC_Pi[7];
	  str = "+1BA7F09AB6B6A9e642";
	  str >> CXSC_Pi[8];
	  str = "-1EDD0DBD2544CFe60A";
	  str >> CXSC_Pi[9];
	  str = "+179FB1BD1310BAe5D3";
	  str >> CXSC_Pi[10];
	  str = "+1A637ED6B0BFF6e59D";
	  str >> CXSC_Pi[11];
	  str = "-1A485FCA40908Ee566";
	  str >> CXSC_Pi[12];
	  str = "-1E501295D98169e52F";
	  str >> CXSC_Pi[13];
	  str = "-1160DBEE83B4E0e4F9";
	  str >> CXSC_Pi[14];
	  str = "-19B6D799AE131Ce4C1";
	  str >> CXSC_Pi[15];
	  str = "+16CF70801F2E28e48B";
	  str >> CXSC_Pi[16];
	  str = "+163BF0598DA483e455";
	  str >> CXSC_Pi[17];
	  str = "+1871574E69A459e41F";
	  str >> CXSC_Pi[18];
	  str = "-15C0B6CC281F27e3E3";
	  str >> CXSC_Pi[19];
	  str = "+15D23DCA3AD962e3AD";
	  str >> CXSC_Pi[20];
	  str = "-1CE8654EFBD56Ae376";
	  str >> CXSC_Pi[21];
	  str = "-1184AB5BE23DA6e33F";
	  str >> CXSC_Pi[22];
	  str = "+166D670C354E4Be309";
	  str >> CXSC_Pi[23];
	  str = "-10D9FEC3A2E4FEe2D3";
	  str >> CXSC_Pi[24];
	  str = "+1943042F86520Ce29C";
	  str >> CXSC_Pi[25];
	  str = "-1B9D1C931C41C6e265";
	  str >> CXSC_Pi[26];
	  str = "-188D3E7F179FC6e22D";
	  str >> CXSC_Pi[27];
	  str = "-1361F1744FE176e1F7";
	  str >> CXSC_Pi[28];
	  str = "+1F6B8ABBE0DE99e1C0";
	  str >> CXSC_Pi[29];
	  str = "-169B10EA1A04B5e18A";
	  str >> CXSC_Pi[30];
	  str = "-14FD1CF8CD56D0e154";
	  str >> CXSC_Pi[31];
	  str = "-18AB54A8D7516Fe11E";
	  str >> CXSC_Pi[32];
	  str = "+186263E8144056e0E7";
	  str >> CXSC_Pi[33];
	  str = "-1AE34AEAAA77A5e0B0";
	  str >> CXSC_Pi[34];
	  str = "+16998B8682283De07A";
	  str >> CXSC_Pi[35];
	  str = "+19D42A90D5EF8Ee042";
	  str >> CXSC_Pi[36];
	  str = "+174C9D9F70A08Be00C";
	  str >> CXSC_Pi[37];
	  str = "+100000000000DBe000";
	  str >> CXSC_Pi[38];
	  str = "+100000000000DCe000";
	  str >> CXSC_Pi[39];

	  CXSC_Pi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1020,y);
 } // Pi


 static real CXSC_Ln2[40]; // CXSC_Ln2[0], ... ,CXSC_Ln2[39]
 static bool CXSC_Ln2_initialized = false;

 lx_interval Ln2_lx_interval() throw()
// Inclusion of ln(2), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ln2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+162E42FEFA39EFe7FC";
	  str >> CXSC_Ln2[0];
	  str = "+1ABC9E3B39803Fe7C5";
	  str >> CXSC_Ln2[1];
	  str = "+17B57A079A1934e78E";
	  str >> CXSC_Ln2[2];
	  str = "-1ACE93A4EBE5D1e758";
	  str >> CXSC_Ln2[3];
	  str = "-123A2A82EA0C24e722";
	  str >> CXSC_Ln2[4];
	  str = "+1D881B7AEB2615e6EB";
	  str >> CXSC_Ln2[5];
	  str = "+19552FB4AFA1B1e6B5";
	  str >> CXSC_Ln2[6];
	  str = "+1DA5D5C6B82704e67C";
	  str >> CXSC_Ln2[7];
	  str = "+14427573B29117e645";
	  str >> CXSC_Ln2[8];
	  str = "-191F6B05A4D7A7e60F";
	  str >> CXSC_Ln2[9];
	  str = "-1DB5173AE53426e5D9";
	  str >> CXSC_Ln2[10];
	  str = "+11317C387EB9EBe5A1";
	  str >> CXSC_Ln2[11];
	  str = "-190F13B267F137e56B";
	  str >> CXSC_Ln2[12];
	  str = "+16FA0EC7657F75e535";
	  str >> CXSC_Ln2[13];
	  str = "-1234C5E1398A6Be4FF";
	  str >> CXSC_Ln2[14];
	  str = "+1195EBBF4D7A70e4C8";
	  str >> CXSC_Ln2[15];
	  str = "+18192432AFD0C4e492";
	  str >> CXSC_Ln2[16];
	  str = "-1A1BE38BA4BA4De45C";
	  str >> CXSC_Ln2[17];
	  str = "-1D7860151CFC06e422";
	  str >> CXSC_Ln2[18];
	  str = "+19423F6B7F720Ce3EC";
	  str >> CXSC_Ln2[19];
	  str = "+10D30F88FE551Ae3B5";
	  str >> CXSC_Ln2[20];
	  str = "+1772B4EB6FE0F8e37E";
	  str >> CXSC_Ln2[21];
	  str = "-17AA0B477087B0e347";
	  str >> CXSC_Ln2[22];
	  str = "+1672C2E8C0EEBBe30C";
	  str >> CXSC_Ln2[23];
	  str = "+1C4C872E4A1F3Ae2D6";
	  str >> CXSC_Ln2[24];
	  str = "-1A970C65986667e2A0";
	  str >> CXSC_Ln2[25];
	  str = "-18CD36365759DAe26A";
	  str >> CXSC_Ln2[26];
	  str = "+1A1E0BD1D6095De231";
	  str >> CXSC_Ln2[27];
	  str = "+12B34D999AB252e1FA";
	  str >> CXSC_Ln2[28];
	  str = "-1912AC700EB43De1C4";
	  str >> CXSC_Ln2[29];
	  str = "-1B8BEFC5924FF5e18E";
	  str >> CXSC_Ln2[30];
	  str = "-180C2AE79DBFADe156";
	  str >> CXSC_Ln2[31];
	  str = "-17D195E5A6D545e120";
	  str >> CXSC_Ln2[32];
	  str = "-1743270F423129e0EA";
	  str >> CXSC_Ln2[33];
	  str = "+189E6DB6303659e0B2";
	  str >> CXSC_Ln2[34];
	  str = "-1F0E11945C9A4Ae07C";
	  str >> CXSC_Ln2[35];
	  str = "+18DAFA85A8C283e046";
	  str >> CXSC_Ln2[36];
	  str = "+13062D3458B6CFe00F";
	  str >> CXSC_Ln2[37];
	  str = "-10000000000C9Be000";
	  str >> CXSC_Ln2[38];
	  str = "-10000000000C9Ae000";
	  str >> CXSC_Ln2[39];

	  CXSC_Ln2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ln2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // ln(2)

 static real CXSC_Ln10[40]; // CXSC_Ln10[0], ... ,CXSC_Ln10[39]
 static bool CXSC_Ln10_initialized = false;

 lx_interval Ln10_lx_interval() throw()
// Inclusion of ln(10), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ln10_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+126BB1BBB55516e7FD";
	  str >> CXSC_Ln10[0];
	  str = "-1F48AD494EA3E9e7C7";
	  str >> CXSC_Ln10[1];
	  str = "-19EBAE3AE0260Ce791";
	  str >> CXSC_Ln10[2];
	  str = "-12D10378BE1CF1e75B";
	  str >> CXSC_Ln10[3];
	  str = "+10403E05AE52C6e725";
	  str >> CXSC_Ln10[4];
	  str = "-1FA509CAFDF466e6ED";
	  str >> CXSC_Ln10[5];
	  str = "-1C79A1FE9D0795e6B7";
	  str >> CXSC_Ln10[6];
	  str = "+1058C448308218e681";
	  str >> CXSC_Ln10[7];
	  str = "-1D250470877BFDe64A";
	  str >> CXSC_Ln10[8];
	  str = "-1AE92987D3075De612";
	  str >> CXSC_Ln10[9];
	  str = "-1D5CDBB8626956e5DC";
	  str >> CXSC_Ln10[10];
	  str = "-13C4F27CE0410Ae5A6";
	  str >> CXSC_Ln10[11];
	  str = "+1B3AC12ACF1BE9e570";
	  str >> CXSC_Ln10[12];
	  str = "+1161BB49D219C8e53A";
	  str >> CXSC_Ln10[13];
	  str = "-110D6613293728e504";
	  str >> CXSC_Ln10[14];
	  str = "+142163A4CDA351e4CC";
	  str >> CXSC_Ln10[15];
	  str = "+1E2713D6C22C16e494";
	  str >> CXSC_Ln10[16];
	  str = "-15090EF85CB0ADe45B";
	  str >> CXSC_Ln10[17];
	  str = "-1C5B3E859F876Ee424";
	  str >> CXSC_Ln10[18];
	  str = "-1C0D54B14459D9e3EC";
	  str >> CXSC_Ln10[19];
	  str = "+1AB685CD44E2C3e3B5";
	  str >> CXSC_Ln10[20];
	  str = "+1A47ECB26978C6e37E";
	  str >> CXSC_Ln10[21];
	  str = "-15812716B8AD41e347";
	  str >> CXSC_Ln10[22];
	  str = "-16047E37E81868e311";
	  str >> CXSC_Ln10[23];
	  str = "+1E14126A45765De2DA";
	  str >> CXSC_Ln10[24];
	  str = "-10ECBE631205C0e2A3";
	  str >> CXSC_Ln10[25];
	  str = "-15A485363BE7D4e26C";
	  str >> CXSC_Ln10[26];
	  str = "-1DEDE455922FF8e234";
	  str >> CXSC_Ln10[27];
	  str = "-1C202C3AE8B719e1FE";
	  str >> CXSC_Ln10[28];
	  str = "-148E3DB9B96D03e1C7";
	  str >> CXSC_Ln10[29];
	  str = "+1E3795D1008FE3e191";
	  str >> CXSC_Ln10[30];
	  str = "-13C992BD5AD855e158";
	  str >> CXSC_Ln10[31];
	  str = "-152096175A0882e122";
	  str >> CXSC_Ln10[32];
	  str = "+1BB0274A1CB072e0EB";
	  str >> CXSC_Ln10[33];
	  str = "-1D6A3FC0087494e0B4";
	  str >> CXSC_Ln10[34];
	  str = "+1AD6BFBFFD821Ce07E";
	  str >> CXSC_Ln10[35];
	  str = "-17D6CD3EE64A79e048";
	  str >> CXSC_Ln10[36];
	  str = "-166DC44198DC68e010";
	  str >> CXSC_Ln10[37];
	  str = "-100000000012D2e000";
	  str >> CXSC_Ln10[38];
	  str = "-100000000012D1e000";
	  str >> CXSC_Ln10[39];

	  CXSC_Ln10_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ln10[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // ln(10)

 static real CXSC_Pir[40]; // CXSC_Pir[0], ... ,CXSC_Pir[39]
 static bool CXSC_Pir_initialized = false;

 lx_interval Pir_lx_interval() throw()
// Inclusion of 1/Pi, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Pir_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+145F306DC9C883e7FC";
	  str >> CXSC_Pir[0];
	  str = "-16B01EC5417056e7C6";
	  str >> CXSC_Pir[1];
	  str = "-16447E493AD4CEe790";
	  str >> CXSC_Pir[2];
	  str = "+1E21C820FF28B2e75A";
	  str >> CXSC_Pir[3];
	  str = "-1508510EA79237e723";
	  str >> CXSC_Pir[4];
	  str = "+1B8E909374B802e6EB";
	  str >> CXSC_Pir[5];
	  str = "-1B6D115F62E6DEe6B5";
	  str >> CXSC_Pir[6];
	  str = "-180F10A71A76B3e67E";
	  str >> CXSC_Pir[7];
	  str = "+1CFBA208D7D4BBe647";
	  str >> CXSC_Pir[8];
	  str = "-12EDEC598E3F65e60F";
	  str >> CXSC_Pir[9];
	  str = "-1741037D8CDC54e5D8";
	  str >> CXSC_Pir[10];
	  str = "+1CC1A99CFA4E42e5A2";
	  str >> CXSC_Pir[11];
	  str = "+17E2EF7E4A0EC8e56B";
	  str >> CXSC_Pir[12];
	  str = "-1DA00087E99FC0e52F";
	  str >> CXSC_Pir[13];
	  str = "-10D0EE74A5F593e4F9";
	  str >> CXSC_Pir[14];
	  str = "+1F6D367ECF27CBe4C1";
	  str >> CXSC_Pir[15];
	  str = "+136E9E8C7ECD3De488";
	  str >> CXSC_Pir[16];
	  str = "-100AE9456C229Ce452";
	  str >> CXSC_Pir[17];
	  str = "-141A0E84C2F8C6e419";
	  str >> CXSC_Pir[18];
	  str = "-10EB5ADA2B2809e3E0";
	  str >> CXSC_Pir[19];
	  str = "-10277039517BD5e3AA";
	  str >> CXSC_Pir[20];
	  str = "+198237E3DB5D60e36E";
	  str >> CXSC_Pir[21];
	  str = "-1E6087BECA1794e338";
	  str >> CXSC_Pir[22];
	  str = "+1DA9E391615EE6e301";
	  str >> CXSC_Pir[23];
	  str = "+1B086599855F15e2C9";
	  str >> CXSC_Pir[24];
	  str = "-17E5EFDC8009E0e293";
	  str >> CXSC_Pir[25];
	  str = "+135CC9CC418185e25B";
	  str >> CXSC_Pir[26];
	  str = "+156CA73A8C960Ee225";
	  str >> CXSC_Pir[27];
	  str = "+13DE04635A3E21e1EE";
	  str >> CXSC_Pir[28];
	  str = "-18F260C88C5FDBe1B7";
	  str >> CXSC_Pir[29];
	  str = "-157CA63B89746Ae181";
	  str >> CXSC_Pir[30];
	  str = "+1CA6DDAF44D157e149";
	  str >> CXSC_Pir[31];
	  str = "+19053EA5FF0705e111";
	  str >> CXSC_Pir[32];
	  str = "+1FBF19F419616Fe0DA";
	  str >> CXSC_Pir[33];
	  str = "+13E60C9F6EF0CFe0A3";
	  str >> CXSC_Pir[34];
	  str = "+126EF6B1E5EF8Ae06D";
	  str >> CXSC_Pir[35];
	  str = "-18BC1946A1B01Ce034";
	  str >> CXSC_Pir[36];
	  str = "-12780EDE6F8384e000";
	  str >> CXSC_Pir[37];
	  str = "+10000000000000e000";
	  str >> CXSC_Pir[38];
	  str = "+10000000000001e000";
	  str >> CXSC_Pir[39];

	  CXSC_Pir_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pir[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // 1/Pi

 static real CXSC_SqrtPi[40]; // CXSC_SqrtPi[0], ... ,CXSC_SqrtPi[39]
 static bool CXSC_SqrtPi_initialized = false;

 lx_interval SqrtPi_lx_interval() throw()
// Inclusion of sqrt(Pi), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_SqrtPi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1C5BF891B4EF6Be7FC";
	  str >> CXSC_SqrtPi[0];
	  str = "-1618F13EB7CA89e7C6";
	  str >> CXSC_SqrtPi[1];
	  str = "-1B1F0071B7AAE4e78E";
	  str >> CXSC_SqrtPi[2];
	  str = "-1389B5A46BDFE8e757";
	  str >> CXSC_SqrtPi[3];
	  str = "-160AF5C5C89448e721";
	  str >> CXSC_SqrtPi[4];
	  str = "-14835F07122994e6E5";
	  str >> CXSC_SqrtPi[5];
	  str = "+1CEC283C18EE8Fe6AF";
	  str >> CXSC_SqrtPi[6];
	  str = "-13ADEBB9223CA8e678";
	  str >> CXSC_SqrtPi[7];
	  str = "+1454912430D291e642";
	  str >> CXSC_SqrtPi[8];
	  str = "-1E8B2345020EF6e60C";
	  str >> CXSC_SqrtPi[9];
	  str = "-17262982556291e5D5";
	  str >> CXSC_SqrtPi[10];
	  str = "+1196FA9B140CABe59E";
	  str >> CXSC_SqrtPi[11];
	  str = "-175EEE59D91D39e568";
	  str >> CXSC_SqrtPi[12];
	  str = "+1789268B7D9D48e52D";
	  str >> CXSC_SqrtPi[13];
	  str = "+17162E2F06B89Ce4F7";
	  str >> CXSC_SqrtPi[14];
	  str = "+1EC9C08F40A3DBe4C0";
	  str >> CXSC_SqrtPi[15];
	  str = "+1B6048DD0729E2e48A";
	  str >> CXSC_SqrtPi[16];
	  str = "+1471CF4C33FF6Be453";
	  str >> CXSC_SqrtPi[17];
	  str = "+1D75FBD8B36F94e41D";
	  str >> CXSC_SqrtPi[18];
	  str = "+16BA59D137CC6Ee3E7";
	  str >> CXSC_SqrtPi[19];
	  str = "-1FDFA25FB4BFD8e3B1";
	  str >> CXSC_SqrtPi[20];
	  str = "-1699363F68227Be379";
	  str >> CXSC_SqrtPi[21];
	  str = "-1BDD2FD4684487e342";
	  str >> CXSC_SqrtPi[22];
	  str = "+1122B2D8015ED6e30B";
	  str >> CXSC_SqrtPi[23];
	  str = "-17EB1A81424DE5e2D2";
	  str >> CXSC_SqrtPi[24];
	  str = "+1C08B42B2EB0E0e29C";
	  str >> CXSC_SqrtPi[25];
	  str = "-1316DE24F93E9Fe266";
	  str >> CXSC_SqrtPi[26];
	  str = "+129354F2D42931e230";
	  str >> CXSC_SqrtPi[27];
	  str = "-1CB7B480D41490e1FA";
	  str >> CXSC_SqrtPi[28];
	  str = "+1608DE7786C4ABe1C3";
	  str >> CXSC_SqrtPi[29];
	  str = "-117732A85F48BCe18D";
	  str >> CXSC_SqrtPi[30];
	  str = "-18BFB034DC2D75e156";
	  str >> CXSC_SqrtPi[31];
	  str = "+155DAB8C4A398Ee120";
	  str >> CXSC_SqrtPi[32];
	  str = "+161C9A5BA77FF3e0E8";
	  str >> CXSC_SqrtPi[33];
	  str = "-1ECF0081DB503Ce0B2";
	  str >> CXSC_SqrtPi[34];
	  str = "-192FF4749E0FD8e07B";
	  str >> CXSC_SqrtPi[35];
	  str = "-1B84C9BCD51654e044";
	  str >> CXSC_SqrtPi[36];
	  str = "-11CF482677D72Fe00A";
	  str >> CXSC_SqrtPi[37];
	  str = "-10000000000025e000";
	  str >> CXSC_SqrtPi[38];
	  str = "-10000000000024e000";
	  str >> CXSC_SqrtPi[39];

	  CXSC_SqrtPi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_SqrtPi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // sqrt(Pi)

 static real CXSC_Sqrt2Pi[40]; // CXSC_Sqrt2Pi[0], ... ,CXSC_Sqrt2Pi[39]
 static bool CXSC_Sqrt2Pi_initialized = false;

 lx_interval Sqrt2Pi_lx_interval() throw()
// Inclusion of sqrt(2*Pi), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt2Pi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+140D931FF62706e7FD";
	  str >> CXSC_Sqrt2Pi[0];
	  str = "-1A6A0D6F814637e7C7";
	  str >> CXSC_Sqrt2Pi[1];
	  str = "-1311D073060ACEe791";
	  str >> CXSC_Sqrt2Pi[2];
	  str = "+16000B50DC2F41e758";
	  str >> CXSC_Sqrt2Pi[3];
	  str = "+16EF75CA45A834e721";
	  str >> CXSC_Sqrt2Pi[4];
	  str = "+19BDB2B4C39342e6E9";
	  str >> CXSC_Sqrt2Pi[5];
	  str = "+1F5582E2063EE6e6B2";
	  str >> CXSC_Sqrt2Pi[6];
	  str = "+183F879BEA150Ce679";
	  str >> CXSC_Sqrt2Pi[7];
	  str = "-1F1EA3CA289B00e641";
	  str >> CXSC_Sqrt2Pi[8];
	  str = "-1699CDA77736F9e60A";
	  str >> CXSC_Sqrt2Pi[9];
	  str = "-11A379D298B55Ee5D1";
	  str >> CXSC_Sqrt2Pi[10];
	  str = "-1A6DDB0152BA94e59B";
	  str >> CXSC_Sqrt2Pi[11];
	  str = "-1957E2E58A02FEe564";
	  str >> CXSC_Sqrt2Pi[12];
	  str = "-1D6160F18E604De52E";
	  str >> CXSC_Sqrt2Pi[13];
	  str = "+1311860CDF7215e4F5";
	  str >> CXSC_Sqrt2Pi[14];
	  str = "+12271F44C50274e4BE";
	  str >> CXSC_Sqrt2Pi[15];
	  str = "-100BF5C5497A21e487";
	  str >> CXSC_Sqrt2Pi[16];
	  str = "+1E94B6E6AD51E2e44F";
	  str >> CXSC_Sqrt2Pi[17];
	  str = "-1C910B5F3D27CEe416";
	  str >> CXSC_Sqrt2Pi[18];
	  str = "+1F266C0EA9E7FBe3E0";
	  str >> CXSC_Sqrt2Pi[19];
	  str = "+1D84A8782A175De3AA";
	  str >> CXSC_Sqrt2Pi[20];
	  str = "-1C75F70BEA4BE2e36F";
	  str >> CXSC_Sqrt2Pi[21];
	  str = "+1ABE52E82F3797e338";
	  str >> CXSC_Sqrt2Pi[22];
	  str = "-10A542CE87822Be302";
	  str >> CXSC_Sqrt2Pi[23];
	  str = "+1457EC878576D9e2CC";
	  str >> CXSC_Sqrt2Pi[24];
	  str = "-1E99158B23E861e296";
	  str >> CXSC_Sqrt2Pi[25];
	  str = "-178F7C0F0F2130e25F";
	  str >> CXSC_Sqrt2Pi[26];
	  str = "-12619EF8E11367e223";
	  str >> CXSC_Sqrt2Pi[27];
	  str = "-1F5C2604AB7BA5e1EC";
	  str >> CXSC_Sqrt2Pi[28];
	  str = "+13ED3039D91C88e1B2";
	  str >> CXSC_Sqrt2Pi[29];
	  str = "-19055B38C434DDe17A";
	  str >> CXSC_Sqrt2Pi[30];
	  str = "+1A7547A2C9CA85e144";
	  str >> CXSC_Sqrt2Pi[31];
	  str = "-1543DE12F9EBD5e10E";
	  str >> CXSC_Sqrt2Pi[32];
	  str = "+180AD772E1F1A0e0D8";
	  str >> CXSC_Sqrt2Pi[33];
	  str = "+1FEF40C66F715Ee0A2";
	  str >> CXSC_Sqrt2Pi[34];
	  str = "+14C244F317CD0De069";
	  str >> CXSC_Sqrt2Pi[35];
	  str = "-1A29C32DBBFCAFe030";
	  str >> CXSC_Sqrt2Pi[36];
	  str = "-1038D4217EB5FCe000";
	  str >> CXSC_Sqrt2Pi[37];
	  str = "-10000000000001e000";
	  str >> CXSC_Sqrt2Pi[38];
	  str = "-10000000000000e000";
	  str >> CXSC_Sqrt2Pi[39];


	  CXSC_Sqrt2Pi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt2Pi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // sqrt(2*Pi)

 static real CXSC_Sqrt2[40]; // CXSC_Sqrt2[0], ... ,CXSC_Sqrt2[39]
 static bool CXSC_Sqrt2_initialized = false;

 lx_interval Sqrt2_lx_interval() throw()
// Inclusion of sqrt(2), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+16A09E667F3BCDe7FC";
	  str >> CXSC_Sqrt2[0];
	  str = "-1BDD3413B26456e7C6";
	  str >> CXSC_Sqrt2[1];
	  str = "+157D3E3ADEC175e790";
	  str >> CXSC_Sqrt2[2];
	  str = "+12775099DA2F59e758";
	  str >> CXSC_Sqrt2[3];
	  str = "+160CCE64552BF2e71F";
	  str >> CXSC_Sqrt2[4];
	  str = "+1821D5C5161D46e6E6";
	  str >> CXSC_Sqrt2[5];
	  str = "-1C032046F8498Ee6B0";
	  str >> CXSC_Sqrt2[6];
	  str = "+1EE950BC8738F7e678";
	  str >> CXSC_Sqrt2[7];
	  str = "-1AC3FDBC64E103e642";
	  str >> CXSC_Sqrt2[8];
	  str = "+13B469101743A1e60A";
	  str >> CXSC_Sqrt2[9];
	  str = "+15E3E9CA60B38Ce5D4";
	  str >> CXSC_Sqrt2[10];
	  str = "+11BC337BCAB1BDe599";
	  str >> CXSC_Sqrt2[11];
	  str = "-1BBA5DEE9D6E7De563";
	  str >> CXSC_Sqrt2[12];
	  str = "-1438DD083B1CC4e52D";
	  str >> CXSC_Sqrt2[13];
	  str = "+1B56A28E2EDFA7e4F7";
	  str >> CXSC_Sqrt2[14];
	  str = "+1CCB2A634331F4e4C1";
	  str >> CXSC_Sqrt2[15];
	  str = "-1BD9056876F83Ee48A";
	  str >> CXSC_Sqrt2[16];
	  str = "-1234FA22AB6BEFe454";
	  str >> CXSC_Sqrt2[17];
	  str = "+19040CA4A81395e41D";
	  str >> CXSC_Sqrt2[18];
	  str = "-15249C0BC4082De3E7";
	  str >> CXSC_Sqrt2[19];
	  str = "+13A02CEBC93E0Ce3B1";
	  str >> CXSC_Sqrt2[20];
	  str = "+109936AF354A2Ee37B";
	  str >> CXSC_Sqrt2[21];
	  str = "-1AE4730CBE4908e345";
	  str >> CXSC_Sqrt2[22];
	  str = "+11B6380826E010e30E";
	  str >> CXSC_Sqrt2[23];
	  str = "-1CDCAD0CCD5A16e2D5";
	  str >> CXSC_Sqrt2[24];
	  str = "-1084BC28012BC8e29C";
	  str >> CXSC_Sqrt2[25];
	  str = "-1C035DDECF8216e265";
	  str >> CXSC_Sqrt2[26];
	  str = "+18907DEAA070B0e22B";
	  str >> CXSC_Sqrt2[27];
	  str = "+1FCBDDEA2F7DC3e1F5";
	  str >> CXSC_Sqrt2[28];
	  str = "+18C41C51757FB0e1BE";
	  str >> CXSC_Sqrt2[29];
	  str = "-189A5B616B1381e188";
	  str >> CXSC_Sqrt2[30];
	  str = "+165C417EFF0B88e152";
	  str >> CXSC_Sqrt2[31];
	  str = "-1627043F832999e11A";
	  str >> CXSC_Sqrt2[32];
	  str = "+105E5FCA017092e0E3";
	  str >> CXSC_Sqrt2[33];
	  str = "-187A16D6A8FDCAe0AD";
	  str >> CXSC_Sqrt2[34];
	  str = "-1838421AE0AE62e072";
	  str >> CXSC_Sqrt2[35];
	  str = "-165073EB433984e03C";
	  str >> CXSC_Sqrt2[36];
	  str = "+1F0A42F9DA4A6Ce006";
	  str >> CXSC_Sqrt2[37];
	  str = "+10000000000002e000";
	  str >> CXSC_Sqrt2[38];
	  str = "+10000000000003e000";
	  str >> CXSC_Sqrt2[39];

	  CXSC_Sqrt2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // sqrt(2)

 static real CXSC_Sqrt3[40]; // CXSC_Sqrt3[0], ... ,CXSC_Sqrt3[39]
 static bool CXSC_Sqrt3_initialized = false;

 lx_interval Sqrt3_lx_interval() throw()
// Inclusion of sqrt(3), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt3_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1BB67AE8584CAAe7FC";
	  str >> CXSC_Sqrt3[0];
	  str = "+1CEC95D0B5C1E3e7C6";
	  str >> CXSC_Sqrt3[1];
	  str = "-1F11DB689F2CCFe78E";
	  str >> CXSC_Sqrt3[2];
	  str = "+13DA4798C720A6e758";
	  str >> CXSC_Sqrt3[3];
	  str = "+121B9169B89243e722";
	  str >> CXSC_Sqrt3[4];
	  str = "-1813508751212Be6E9";
	  str >> CXSC_Sqrt3[5];
	  str = "-1B3D547B775C1Ee6B2";
	  str >> CXSC_Sqrt3[6];
	  str = "-19D986D92E2F0Ae679";
	  str >> CXSC_Sqrt3[7];
	  str = "+1A34334CE806B6e642";
	  str >> CXSC_Sqrt3[8];
	  str = "+1A383B9E122E61e60C";
	  str >> CXSC_Sqrt3[9];
	  str = "+1C61D736F2F6F2e5D5";
	  str >> CXSC_Sqrt3[10];
	  str = "-10AF49233F9250e59E";
	  str >> CXSC_Sqrt3[11];
	  str = "-1558A109EC0523e567";
	  str >> CXSC_Sqrt3[12];
	  str = "+1F799D4D4FF2BCe531";
	  str >> CXSC_Sqrt3[13];
	  str = "-1AD7B219E34EDBe4FB";
	  str >> CXSC_Sqrt3[14];
	  str = "+15AB940B6677E3e4C5";
	  str >> CXSC_Sqrt3[15];
	  str = "-1D9B2A8203B8F0e48E";
	  str >> CXSC_Sqrt3[16];
	  str = "-1DB0C8975A3834e458";
	  str >> CXSC_Sqrt3[17];
	  str = "-1BCAAB3F6BE884e422";
	  str >> CXSC_Sqrt3[18];
	  str = "+14C70ADB1EC1BBe3E8";
	  str >> CXSC_Sqrt3[19];
	  str = "-14E1EF77987E55e3AF";
	  str >> CXSC_Sqrt3[20];
	  str = "-19695FC6269D28e378";
	  str >> CXSC_Sqrt3[21];
	  str = "+10D0652AAC5936e342";
	  str >> CXSC_Sqrt3[22];
	  str = "-1BD0891D370824e30C";
	  str >> CXSC_Sqrt3[23];
	  str = "-129B4C6252D061e2D4";
	  str >> CXSC_Sqrt3[24];
	  str = "+1DC9B1A4C31275e29E";
	  str >> CXSC_Sqrt3[25];
	  str = "+11FF9B8422294Ee267";
	  str >> CXSC_Sqrt3[26];
	  str = "-1E4A6AA47F3A85e231";
	  str >> CXSC_Sqrt3[27];
	  str = "+17043E01AA3F3De1FA";
	  str >> CXSC_Sqrt3[28];
	  str = "+188EF377D2D5B6e1C0";
	  str >> CXSC_Sqrt3[29];
	  str = "-1735E8C815F031e185";
	  str >> CXSC_Sqrt3[30];
	  str = "-1B89330FD8417Ce14F";
	  str >> CXSC_Sqrt3[31];
	  str = "+16D1A627670F5Ce117";
	  str >> CXSC_Sqrt3[32];
	  str = "+1AF43BBA8154D3e0DB";
	  str >> CXSC_Sqrt3[33];
	  str = "+1DA9A969A91295e0A5";
	  str >> CXSC_Sqrt3[34];
	  str = "-1636594394C675e06E";
	  str >> CXSC_Sqrt3[35];
	  str = "+1064B9DA1A3185e037";
	  str >> CXSC_Sqrt3[36];
	  str = "-1CE514CF1825CCe001";
	  str >> CXSC_Sqrt3[37];
	  str = "+10000000000000e000";
	  str >> CXSC_Sqrt3[38];
	  str = "+10000000000001e000";
	  str >> CXSC_Sqrt3[39];

	  CXSC_Sqrt3_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt3[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // sqrt(3)

 static real CXSC_Ln2r[40]; // CXSC_Ln2r[0], ... ,CXSC_Ln2r[39]
 static bool CXSC_Ln2r_initialized = false;

 lx_interval Ln2r_lx_interval() throw()
// Inclusion of 1/ln(2), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ln2r_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+171547652B82FEe7FC";
	  str >> CXSC_Ln2r[0];
	  str = "+1777D0FFDA0D24e7C4";
	  str >> CXSC_Ln2r[1];
	  str = "-160BB8A5442AB9e78E";
	  str >> CXSC_Ln2r[2];
	  str = "-14B52D3BA6D74De756";
	  str >> CXSC_Ln2r[3];
	  str = "+19A342648FBC39e720";
	  str >> CXSC_Ln2r[4];
	  str = "-1E0455744994EEe6EA";
	  str >> CXSC_Ln2r[5];
	  str = "+1B25EEB82D7C16e6B4";
	  str >> CXSC_Ln2r[6];
	  str = "+1F5485CF306255e67E";
	  str >> CXSC_Ln2r[7];
	  str = "-1EC07680A1F958e648";
	  str >> CXSC_Ln2r[8];
	  str = "-106326680EB5B6e612";
	  str >> CXSC_Ln2r[9];
	  str = "-1B3D04C549BC98e5DC";
	  str >> CXSC_Ln2r[10];
	  str = "+1EABCEAD10305Be5A6";
	  str >> CXSC_Ln2r[11];
	  str = "-14440C57D7AB97e56D";
	  str >> CXSC_Ln2r[12];
	  str = "-17185D42A4E6D6e536";
	  str >> CXSC_Ln2r[13];
	  str = "-1F332B5BE48526e4FE";
	  str >> CXSC_Ln2r[14];
	  str = "+12CE4F199E108De4C8";
	  str >> CXSC_Ln2r[15];
	  str = "-18DAFCC6077F2Ae48F";
	  str >> CXSC_Ln2r[16];
	  str = "+19ABB71EC25E12e458";
	  str >> CXSC_Ln2r[17];
	  str = "-11473D7A3366BDe41F";
	  str >> CXSC_Ln2r[18];
	  str = "-125DF4E28B5ED4e3E8";
	  str >> CXSC_Ln2r[19];
	  str = "+1C64262D010330e3B2";
	  str >> CXSC_Ln2r[20];
	  str = "-17DCAE42742BDEe37C";
	  str >> CXSC_Ln2r[21];
	  str = "+109C8C7E7B896Fe346";
	  str >> CXSC_Ln2r[22];
	  str = "+10C470FE2464B9e310";
	  str >> CXSC_Ln2r[23];
	  str = "-1B5F6CFB7C34BEe2DA";
	  str >> CXSC_Ln2r[24];
	  str = "-125E5DBA4A1165e2A1";
	  str >> CXSC_Ln2r[25];
	  str = "-1FA683975309E6e26A";
	  str >> CXSC_Ln2r[26];
	  str = "-140C23C4E5CC64e233";
	  str >> CXSC_Ln2r[27];
	  str = "+117670EC70E797e1FD";
	  str >> CXSC_Ln2r[28];
	  str = "+1B2A04B8E7416Ce1C7";
	  str >> CXSC_Ln2r[29];
	  str = "+11D96159397087e18F";
	  str >> CXSC_Ln2r[30];
	  str = "+10E29D810B4C60e159";
	  str >> CXSC_Ln2r[31];
	  str = "+1D7442ECEFEFA1e123";
	  str >> CXSC_Ln2r[32];
	  str = "+1CE25B70026529e0ED";
	  str >> CXSC_Ln2r[33];
	  str = "-12CA24549E0811e0B7";
	  str >> CXSC_Ln2r[34];
	  str = "+1220755E0827AEe080";
	  str >> CXSC_Ln2r[35];
	  str = "+1086BCE30D4370e04A";
	  str >> CXSC_Ln2r[36];
	  str = "-16FF855E4293BCe011";
	  str >> CXSC_Ln2r[37];
	  str = "+10000000002A50e000";
	  str >> CXSC_Ln2r[38];
	  str = "+10000000002A51e000";
	  str >> CXSC_Ln2r[39];

	  CXSC_Ln2r_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ln2r[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // 1/ln(2)

 static real CXSC_Pid3[40]; // CXSC_Pid3[0], ... ,CXSC_Pid3[39]
 static bool CXSC_Pid3_initialized = false;

 lx_interval Pid3_lx_interval() throw()
// Inclusion of Pi/3, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Pid3_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+10C152382D7366e7FD";
	  str >> CXSC_Pid3[0];
	  str = "-1EE6913347C2A6e7C7";
	  str >> CXSC_Pid3[1];
	  str = "-14BBA47A9E5FD2e78F";
	  str >> CXSC_Pid3[2];
	  str = "-1CCAEF65529B02e759";
	  str >> CXSC_Pid3[3];
	  str = "+197CB7BCC18B87e722";
	  str >> CXSC_Pid3[4];
	  str = "-13EBBDA1FF3058e6EC";
	  str >> CXSC_Pid3[5];
	  str = "-11D10CB320F4D1e6B4";
	  str >> CXSC_Pid3[6];
	  str = "+1958EB892987ECe67D";
	  str >> CXSC_Pid3[7];
	  str = "+167C54B11CF247e647";
	  str >> CXSC_Pid3[8];
	  str = "+12C2E985923A44e60E";
	  str >> CXSC_Pid3[9];
	  str = "+1945484A2DD81Fe5D6";
	  str >> CXSC_Pid3[10];
	  str = "+1197A9E475D54Fe59E";
	  str >> CXSC_Pid3[11];
	  str = "-1E181FEE158585e568";
	  str >> CXSC_Pid3[12];
	  str = "+1047FCE7066A6Ee532";
	  str >> CXSC_Pid3[13];
	  str = "+1D1A8602EA0C85e4FC";
	  str >> CXSC_Pid3[14];
	  str = "+14430C5998BF34e4C6";
	  str >> CXSC_Pid3[15];
	  str = "+173BF40AAD43D9e48F";
	  str >> CXSC_Pid3[16];
	  str = "-137B014DDEDCF5e459";
	  str >> CXSC_Pid3[17];
	  str = "-1A5F1B210EE7C5e420";
	  str >> CXSC_Pid3[18];
	  str = "+151B536DDF9502e3EA";
	  str >> CXSC_Pid3[19];
	  str = "+10E4DB4F709CEEe3B4";
	  str >> CXSC_Pid3[20];
	  str = "+16841F78EC058Ee37E";
	  str >> CXSC_Pid3[21];
	  str = "+1D269E370AFA06e346";
	  str >> CXSC_Pid3[22];
	  str = "+119123BD75E37Be310";
	  str >> CXSC_Pid3[23];
	  str = "+1C7DBAADF64D9De2DA";
	  str >> CXSC_Pid3[24];
	  str = "+16CC595AEA086De2A4";
	  str >> CXSC_Pid3[25];
	  str = "+1942EC979DED29e26E";
	  str >> CXSC_Pid3[26];
	  str = "+1EFBE875957C10e238";
	  str >> CXSC_Pid3[27];
	  str = "+133B7D68BA4029e1FF";
	  str >> CXSC_Pid3[28];
	  str = "-11EB0DA382BF6Ce1C9";
	  str >> CXSC_Pid3[29];
	  str = "+1970EDF4B943FDe193";
	  str >> CXSC_Pid3[30];
	  str = "-11C6A6D14BBC74e15C";
	  str >> CXSC_Pid3[31];
	  str = "+1FBE371E3DC1D2e125";
	  str >> CXSC_Pid3[32];
	  str = "-1F34D225753A55e0EF";
	  str >> CXSC_Pid3[33];
	  str = "+1D0DA1E2E38EC1e0B7";
	  str >> CXSC_Pid3[34];
	  str = "-18C889B4CA7CA6e07E";
	  str >> CXSC_Pid3[35];
	  str = "+1B346B8DAF1FA8e048";
	  str >> CXSC_Pid3[36];
	  str = "+1326EDF35258AEe012";
	  str >> CXSC_Pid3[37];
	  str = "-1000000000DEDBe000";
	  str >> CXSC_Pid3[38];
	  str = "-1000000000DEDAe000";
	  str >> CXSC_Pid3[39];

	  CXSC_Pid3_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pid3[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // Pi/3

 static real CXSC_SqrtPir[40]; // CXSC_SqrtPir[0], ... ,CXSC_SqrtPir[39]
 static bool CXSC_SqrtPir_initialized = false;

 lx_interval SqrtPir_lx_interval() throw()
// Inclusion of 1/sqrt(Pi), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_SqrtPir_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+120DD750429B6De7FD";
	  str >> CXSC_SqrtPir[0];
	  str = "+11AE3A914FED80e7C5";
	  str >> CXSC_SqrtPir[1];
	  str = "-13CBBEBF65F145e78E";
	  str >> CXSC_SqrtPir[2];
	  str = "-1E0C574632F53Ee757";
	  str >> CXSC_SqrtPir[3];
	  str = "-1E6633BE9E7F15e721";
	  str >> CXSC_SqrtPir[4];
	  str = "+1CF859270F1141e6EA";
	  str >> CXSC_SqrtPir[5];
	  str = "-1FE4FB499C328Ae6B3";
	  str >> CXSC_SqrtPir[6];
	  str = "-10B82C446DC78De67C";
	  str >> CXSC_SqrtPir[7];
	  str = "-1878B089078800e646";
	  str >> CXSC_SqrtPir[8];
	  str = "-13DAEADA9E233Ee60E";
	  str >> CXSC_SqrtPir[9];
	  str = "+1137197A708BD2e5D8";
	  str >> CXSC_SqrtPir[10];
	  str = "-109009506D5BA2e59D";
	  str >> CXSC_SqrtPir[11];
	  str = "+17C9F0B5951E94e567";
	  str >> CXSC_SqrtPir[12];
	  str = "-1735F4949633A4e530";
	  str >> CXSC_SqrtPir[13];
	  str = "-146014DBC90D0Ee4FA";
	  str >> CXSC_SqrtPir[14];
	  str = "+1CAB0B222EEEA0e4C4";
	  str >> CXSC_SqrtPir[15];
	  str = "+1B1C750754B40Ae48E";
	  str >> CXSC_SqrtPir[16];
	  str = "-16B2CD2E72C16Ee456";
	  str >> CXSC_SqrtPir[17];
	  str = "-148C024FF194B2e420";
	  str >> CXSC_SqrtPir[18];
	  str = "+1CF866DD09628De3EA";
	  str >> CXSC_SqrtPir[19];
	  str = "-16CBF3DC0C536Ee3B3";
	  str >> CXSC_SqrtPir[20];
	  str = "+1EEDA3436CC85Fe37C";
	  str >> CXSC_SqrtPir[21];
	  str = "+16928D2C52986Ae346";
	  str >> CXSC_SqrtPir[22];
	  str = "+1EF82FD6BDE19De30E";
	  str >> CXSC_SqrtPir[23];
	  str = "+1218F360779AEEe2D8";
	  str >> CXSC_SqrtPir[24];
	  str = "+14F13DB14D3AC2e2A2";
	  str >> CXSC_SqrtPir[25];
	  str = "+1461BD3C8DC495e26C";
	  str >> CXSC_SqrtPir[26];
	  str = "-1CAE7855A13FF8e234";
	  str >> CXSC_SqrtPir[27];
	  str = "-18CEA8571A02F7e1FE";
	  str >> CXSC_SqrtPir[28];
	  str = "-1AC80A19057BBDe1C8";
	  str >> CXSC_SqrtPir[29];
	  str = "+1A1910D82DC198e191";
	  str >> CXSC_SqrtPir[30];
	  str = "+16B8F9198BC17Fe15B";
	  str >> CXSC_SqrtPir[31];
	  str = "+1D7B579C6CAF5De123";
	  str >> CXSC_SqrtPir[32];
	  str = "+1E97B95E80FE25e0ED";
	  str >> CXSC_SqrtPir[33];
	  str = "-188EF7630D4F86e0AF";
	  str >> CXSC_SqrtPir[34];
	  str = "+1988FBA498490Ae076";
	  str >> CXSC_SqrtPir[35];
	  str = "+178235EAEC9403e03F";
	  str >> CXSC_SqrtPir[36];
	  str = "+1067F03DFDFF93e005";
	  str >> CXSC_SqrtPir[37];
	  str = "+10000000000006e000";
	  str >> CXSC_SqrtPir[38];
	  str = "+10000000000007e000";
	  str >> CXSC_SqrtPir[39];

	  CXSC_SqrtPir_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_SqrtPir[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // 1/sqrt(Pi)

 static real CXSC_Sqrt2Pir[40]; // CXSC_Sqrt2Pir[0], ... ,CXSC_Sqrt2Pir[39]
 static bool CXSC_Sqrt2Pir_initialized = false;

 lx_interval Sqrt2Pir_lx_interval() throw()
// Inclusion of 1/sqrt(2*Pi), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt2Pir_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+19884533D43651e7FC";
	  str >> CXSC_Sqrt2Pir[0];
	  str = "-1CBC0D30EBFD15e7C6";
	  str >> CXSC_Sqrt2Pir[1];
	  str = "-1C7402C7D60CFBe78E";
	  str >> CXSC_Sqrt2Pir[2];
	  str = "+12706D8C0471B5e756";
	  str >> CXSC_Sqrt2Pir[3];
	  str = "-1FF6718B45881De720";
	  str >> CXSC_Sqrt2Pir[4];
	  str = "-13AABB82C248DCe6EA";
	  str >> CXSC_Sqrt2Pir[5];
	  str = "-1458A899162EE4e6B1";
	  str >> CXSC_Sqrt2Pir[6];
	  str = "-14EBD8868F41EBe67A";
	  str >> CXSC_Sqrt2Pir[7];
	  str = "+13278E993445F1e642";
	  str >> CXSC_Sqrt2Pir[8];
	  str = "-1CC019F5F4780Ae60C";
	  str >> CXSC_Sqrt2Pir[9];
	  str = "+147CE4B4ECDBD7e5D6";
	  str >> CXSC_Sqrt2Pir[10];
	  str = "-19A3DCC6A3534Be59E";
	  str >> CXSC_Sqrt2Pir[11];
	  str = "+11379A7BA8CB0Ae568";
	  str >> CXSC_Sqrt2Pir[12];
	  str = "-12D909C875312Ee531";
	  str >> CXSC_Sqrt2Pir[13];
	  str = "+1C1CEC4882C77Be4FA";
	  str >> CXSC_Sqrt2Pir[14];
	  str = "-14C4078263DF36e4C4";
	  str >> CXSC_Sqrt2Pir[15];
	  str = "+1AB3FC8D2AB243e48E";
	  str >> CXSC_Sqrt2Pir[16];
	  str = "+17B9172454310Ae458";
	  str >> CXSC_Sqrt2Pir[17];
	  str = "-1444B6B781B7F2e422";
	  str >> CXSC_Sqrt2Pir[18];
	  str = "-1DB5C6773B74B7e3EC";
	  str >> CXSC_Sqrt2Pir[19];
	  str = "-12D4CD4FAB8CF9e3B5";
	  str >> CXSC_Sqrt2Pir[20];
	  str = "-12C319ACF346DCe37D";
	  str >> CXSC_Sqrt2Pir[21];
	  str = "+193ED298857CB8e346";
	  str >> CXSC_Sqrt2Pir[22];
	  str = "+1AB87659565E92e30E";
	  str >> CXSC_Sqrt2Pir[23];
	  str = "-1AEB785019F78Ee2D6";
	  str >> CXSC_Sqrt2Pir[24];
	  str = "-17DAF38DE68CA0e29F";
	  str >> CXSC_Sqrt2Pir[25];
	  str = "-14D672D025580Ce265";
	  str >> CXSC_Sqrt2Pir[26];
	  str = "-17AA87F2ABB794e22F";
	  str >> CXSC_Sqrt2Pir[27];
	  str = "+16E85953CBD917e1F9";
	  str >> CXSC_Sqrt2Pir[28];
	  str = "-10AB555C9A9735e1C3";
	  str >> CXSC_Sqrt2Pir[29];
	  str = "+1020FEB8ED1EA2e18A";
	  str >> CXSC_Sqrt2Pir[30];
	  str = "+18282C79079F71e152";
	  str >> CXSC_Sqrt2Pir[31];
	  str = "+13D282FF699FC4e11C";
	  str >> CXSC_Sqrt2Pir[32];
	  str = "-1A547E139AE10Ce0E6";
	  str >> CXSC_Sqrt2Pir[33];
	  str = "+167FCF0E311B0De0AE";
	  str >> CXSC_Sqrt2Pir[34];
	  str = "+1AC5D5E32ED719e078";
	  str >> CXSC_Sqrt2Pir[35];
	  str = "+1B311EB9071956e042";
	  str >> CXSC_Sqrt2Pir[36];
	  str = "+1D237BE89A6494e00C";
	  str >> CXSC_Sqrt2Pir[37];
	  str = "+10000000000075e000";
	  str >> CXSC_Sqrt2Pir[38];
	  str = "+10000000000076e000";
	  str >> CXSC_Sqrt2Pir[39];

	  CXSC_Sqrt2Pir_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt2Pir[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // 1/sqrt(2*Pi)

 static real CXSC_LnPi[40]; // CXSC_LnPi[0], ... ,CXSC_LnPi[39]
 static bool CXSC_LnPi_initialized = false;

 lx_interval LnPi_lx_interval() throw()
// Inclusion of ln(Pi), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_LnPi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1250D048E7A1BDe7FD";
	  str >> CXSC_LnPi[0];
	  str = "+17ABF2AD8D5088e7C4";
	  str >> CXSC_LnPi[1];
	  str = "-16CCF43244818Ae78C";
	  str >> CXSC_LnPi[2];
	  str = "+1F9303719C0176e756";
	  str >> CXSC_LnPi[3];
	  str = "+15DF52611CB54Ee720";
	  str >> CXSC_LnPi[4];
	  str = "-1D9056E74F8C97e6EA";
	  str >> CXSC_LnPi[5];
	  str = "+100B095B6C2E1Ae6B3";
	  str >> CXSC_LnPi[6];
	  str = "-18C7557878A9E7e67D";
	  str >> CXSC_LnPi[7];
	  str = "+1B9BBBB4F4CEE7e646";
	  str >> CXSC_LnPi[8];
	  str = "+1B477FCC702F86e610";
	  str >> CXSC_LnPi[9];
	  str = "+141F1344A31799e5DA";
	  str >> CXSC_LnPi[10];
	  str = "+1B6740BE95CD58e5A4";
	  str >> CXSC_LnPi[11];
	  str = "-1F2C63904D27DBe56C";
	  str >> CXSC_LnPi[12];
	  str = "+1426F00B933976e534";
	  str >> CXSC_LnPi[13];
	  str = "+125703BE5FAA20e4FE";
	  str >> CXSC_LnPi[14];
	  str = "-1DADAE5397F95Be4C7";
	  str >> CXSC_LnPi[15];
	  str = "+17C9D110381543e48F";
	  str >> CXSC_LnPi[16];
	  str = "-1259230E627FCAe459";
	  str >> CXSC_LnPi[17];
	  str = "+191CEAB6B13A33e422";
	  str >> CXSC_LnPi[18];
	  str = "+109D49A13CB595e3EB";
	  str >> CXSC_LnPi[19];
	  str = "-12C574CDCD41C2e3B4";
	  str >> CXSC_LnPi[20];
	  str = "+1D4141476C3E9De37E";
	  str >> CXSC_LnPi[21];
	  str = "+1D26D892B64467e344";
	  str >> CXSC_LnPi[22];
	  str = "-16BAF3B607E1ADe30B";
	  str >> CXSC_LnPi[23];
	  str = "+165A768D8BC5ADe2D5";
	  str >> CXSC_LnPi[24];
	  str = "+12185364B32BD1e29D";
	  str >> CXSC_LnPi[25];
	  str = "-114D72550F0B90e266";
	  str >> CXSC_LnPi[26];
	  str = "+1E586BAEEB8BF4e230";
	  str >> CXSC_LnPi[27];
	  str = "-1F4B9322D4506Fe1F8";
	  str >> CXSC_LnPi[28];
	  str = "+16D32BAA9A4FCCe1C2";
	  str >> CXSC_LnPi[29];
	  str = "+1A12D8CF8CC6DAe18B";
	  str >> CXSC_LnPi[30];
	  str = "-1215CF3BD682CAe155";
	  str >> CXSC_LnPi[31];
	  str = "-184CA7D8873E45e11D";
	  str >> CXSC_LnPi[32];
	  str = "-1F02ECC3E58C6Ee0E7";
	  str >> CXSC_LnPi[33];
	  str = "-11F6EC8ED0D92Be0B1";
	  str >> CXSC_LnPi[34];
	  str = "-199F29ACE1FC18e077";
	  str >> CXSC_LnPi[35];
	  str = "+119F3673AA919Ae041";
	  str >> CXSC_LnPi[36];
	  str = "-1A8359A2831626e00A";
	  str >> CXSC_LnPi[37];
	  str = "-10000000000072e000";
	  str >> CXSC_LnPi[38];
	  str = "-10000000000071e000";
	  str >> CXSC_LnPi[39];

	  CXSC_LnPi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_LnPi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // ln(Pi)

 static real CXSC_Ln2Pi[40]; // CXSC_Ln2Pi[0], ... ,CXSC_Ln2Pi[39]
 static bool CXSC_Ln2Pi_initialized = false;

 lx_interval Ln2Pi_lx_interval() throw()
// Inclusion of ln(2*Pi), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ln2Pi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1D67F1C864BEB5e7FC";
	  str >> CXSC_Ln2Pi[0];
	  str = "-165B5A1B7FF5DFe7C6";
	  str >> CXSC_Ln2Pi[1];
	  str = "-1B7F70C13DC1CCe78F";
	  str >> CXSC_Ln2Pi[2];
	  str = "+13458B4DDEC6A3e759";
	  str >> CXSC_Ln2Pi[3];
	  str = "+133DAA155D2130e721";
	  str >> CXSC_Ln2Pi[4];
	  str = "-18A007FC5E501Be6EB";
	  str >> CXSC_Ln2Pi[5];
	  str = "-15406FA3AA9644e6B1";
	  str >> CXSC_Ln2Pi[6];
	  str = "-13E8D52A392CC9e67B";
	  str >> CXSC_Ln2Pi[7];
	  str = "-1A43099131E88De645";
	  str >> CXSC_Ln2Pi[8];
	  str = "-114835B6623C4De60F";
	  str >> CXSC_Ln2Pi[9];
	  str = "-1ABB7858CF827Ae5D9";
	  str >> CXSC_Ln2Pi[10];
	  str = "+1D8D7045A5A495e5A3";
	  str >> CXSC_Ln2Pi[11];
	  str = "+1A26094B3F6FC5e56C";
	  str >> CXSC_Ln2Pi[12];
	  str = "-1EF27932D0E3D0e534";
	  str >> CXSC_Ln2Pi[13];
	  str = "-12128804136AB6e4FD";
	  str >> CXSC_Ln2Pi[14];
	  str = "+15F8A4AC0BEE17e4C4";
	  str >> CXSC_Ln2Pi[15];
	  str = "+1892F2A5B69B5Fe48E";
	  str >> CXSC_Ln2Pi[16];
	  str = "+1CC7C09477ADCEe458";
	  str >> CXSC_Ln2Pi[17];
	  str = "-116DD579AF074Ae41F";
	  str >> CXSC_Ln2Pi[18];
	  str = "+190E43C1DCCD69e3E7";
	  str >> CXSC_Ln2Pi[19];
	  str = "-11F55BBD0978D3e3AF";
	  str >> CXSC_Ln2Pi[20];
	  str = "+167EC65B83F29Be378";
	  str >> CXSC_Ln2Pi[21];
	  str = "-14C0D466FC8C7Ae33C";
	  str >> CXSC_Ln2Pi[22];
	  str = "-1D56DE4860435Ce306";
	  str >> CXSC_Ln2Pi[23];
	  str = "-10C7B15DFFBDFCe2D0";
	  str >> CXSC_Ln2Pi[24];
	  str = "-15007E40803B52e299";
	  str >> CXSC_Ln2Pi[25];
	  str = "+1DF2A457B56D15e261";
	  str >> CXSC_Ln2Pi[26];
	  str = "-16B7CAD686151De22B";
	  str >> CXSC_Ln2Pi[27];
	  str = "-11F972F1A61CA1e1F5";
	  str >> CXSC_Ln2Pi[28];
	  str = "+1443CF52FBF6B4e1BE";
	  str >> CXSC_Ln2Pi[29];
	  str = "-12652AE82DC678e187";
	  str >> CXSC_Ln2Pi[30];
	  str = "-11712858901127e151";
	  str >> CXSC_Ln2Pi[31];
	  str = "+124D51F4842F1Fe11B";
	  str >> CXSC_Ln2Pi[32];
	  str = "+1B8F6B0823A92Ae0E4";
	  str >> CXSC_Ln2Pi[33];
	  str = "-1742244E0D8F40e0AB";
	  str >> CXSC_Ln2Pi[34];
	  str = "+127A8F1E2AEAD3e074";
	  str >> CXSC_Ln2Pi[35];
	  str = "-180BC6B9E8F00Ce03D";
	  str >> CXSC_Ln2Pi[36];
	  str = "+190933A24F0ECEe007";
	  str >> CXSC_Ln2Pi[37];
	  str = "-10000000000007e000";
	  str >> CXSC_Ln2Pi[38];
	  str = "-10000000000006e000";
	  str >> CXSC_Ln2Pi[39];

	  CXSC_Ln2Pi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ln2Pi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // ln(2*Pi)

 static real CXSC_E[40]; // CXSC_E[0], ... ,CXSC_E[39]
 static bool CXSC_E_initialized = false;

 lx_interval E_lx_interval() throw()
// Inclusion of e, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_E_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+15BF0A8B145769e7FD";
	  str >> CXSC_E[0];
	  str = "+14D57EE2B1013Ae7C7";
	  str >> CXSC_E[1];
	  str = "-1618713A31D3E2e78F";
	  str >> CXSC_E[2];
	  str = "+1C5A6D2B53C26De759";
	  str >> CXSC_E[3];
	  str = "-1F75CDE60219B6e723";
	  str >> CXSC_E[4];
	  str = "-188C76D93041A1e6EC";
	  str >> CXSC_E[5];
	  str = "+12FE363630C75Ee6B6";
	  str >> CXSC_E[6];
	  str = "-1C25F937F544EEe680";
	  str >> CXSC_E[7];
	  str = "-1E852C20E12A2Ae64A";
	  str >> CXSC_E[8];
	  str = "-14D4F6DE605705e60F";
	  str >> CXSC_E[9];
	  str = "-1F3225EF539355e5D5";
	  str >> CXSC_E[10];
	  str = "-16109728625547e59F";
	  str >> CXSC_E[11];
	  str = "-194301506D94CFe569";
	  str >> CXSC_E[12];
	  str = "-1879C78F8CBA44e533";
	  str >> CXSC_E[13];
	  str = "-1D5976250C1018e4FA";
	  str >> CXSC_E[14];
	  str = "+1C877C56284DABe4C4";
	  str >> CXSC_E[15];
	  str = "+1E73530ACCA4F5e48E";
	  str >> CXSC_E[16];
	  str = "-1F161A150FD53Ae458";
	  str >> CXSC_E[17];
	  str = "+159927DB0E8845e41F";
	  str >> CXSC_E[18];
	  str = "+12976591C7F773e3E9";
	  str >> CXSC_E[19];
	  str = "-1525489F280B98e3B2";
	  str >> CXSC_E[20];
	  str = "+1D4F42A3DE394Ee37A";
	  str >> CXSC_E[21];
	  str = "-16A3522431391Be341";
	  str >> CXSC_E[22];
	  str = "+1D8C8583D3E477e30B";
	  str >> CXSC_E[23];
	  str = "+14DAE13C05F9C4e2D1";
	  str >> CXSC_E[24];
	  str = "-19040E899FE5FEe29B";
	  str >> CXSC_E[25];
	  str = "+19A50685EC322Ee265";
	  str >> CXSC_E[26];
	  str = "+17F4E74C2C1FFCe22F";
	  str >> CXSC_E[27];
	  str = "+1C9E2465DDE503e1F9";
	  str >> CXSC_E[28];
	  str = "+1E1FF1D8DA637De1BF";
	  str >> CXSC_E[29];
	  str = "+1AE6776BF9785Ee189";
	  str >> CXSC_E[30];
	  str = "-1EEFFD1D38873Ee153";
	  str >> CXSC_E[31];
	  str = "-105D2F89A72197e11D";
	  str >> CXSC_E[32];
	  str = "+11360D977FD443e0E7";
	  str >> CXSC_E[33];
	  str = "+168470C23F9FBAe0B1";
	  str >> CXSC_E[34];
	  str = "-10E552624D737Ee07B";
	  str >> CXSC_E[35];
	  str = "-148879616420CCe045";
	  str >> CXSC_E[36];
	  str = "-1FEE3CF25C81B1e00F";
	  str >> CXSC_E[37];
	  str = "-10000000001233e000";
	  str >> CXSC_E[38];
	  str = "-10000000001232e000";
	  str >> CXSC_E[39];

	  CXSC_E_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_E[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // e

 static real CXSC_EpPi[40]; // CXSC_EpPi[0], ... ,CXSC_EpPi[39]
 static bool CXSC_EpPi_initialized = false;

 lx_interval EpPi_lx_interval() throw()
// Inclusion of e^Pi, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_EpPi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1724046EB0933Ae7FD";
	  str >> CXSC_EpPi[0];
	  str = "-184C962DD81952e7C7";
	  str >> CXSC_EpPi[1];
	  str = "-12D659C0BCD22Ee790";
	  str >> CXSC_EpPi[2];
	  str = "+117496B8A92F91e75A";
	  str >> CXSC_EpPi[3];
	  str = "+16A8C4203E5FCDe724";
	  str >> CXSC_EpPi[4];
	  str = "-166B11F99A663Be6EE";
	  str >> CXSC_EpPi[5];
	  str = "-118EC2076DABB1e6B8";
	  str >> CXSC_EpPi[6];
	  str = "+19776E5BEB18A5e682";
	  str >> CXSC_EpPi[7];
	  str = "+1AD4091E84B051e64C";
	  str >> CXSC_EpPi[8];
	  str = "+1E89AA12909B40e616";
	  str >> CXSC_EpPi[9];
	  str = "+1ACE3C0DDBB994e5DD";
	  str >> CXSC_EpPi[10];
	  str = "+141EC9379CBBFEe5A6";
	  str >> CXSC_EpPi[11];
	  str = "+1FC4E78D00A016e56D";
	  str >> CXSC_EpPi[12];
	  str = "+1608BE35B9A409e537";
	  str >> CXSC_EpPi[13];
	  str = "-1A0D8AA90EB6B9e4FD";
	  str >> CXSC_EpPi[14];
	  str = "+106FE8AFD21ACFe4C7";
	  str >> CXSC_EpPi[15];
	  str = "+1C072FEA1BFCAFe48F";
	  str >> CXSC_EpPi[16];
	  str = "+1915B9F352EC68e455";
	  str >> CXSC_EpPi[17];
	  str = "-13FA07C37897E9e41E";
	  str >> CXSC_EpPi[18];
	  str = "-1EC01C89B8DDFFe3E8";
	  str >> CXSC_EpPi[19];
	  str = "+1EDC3A402336AFe3B0";
	  str >> CXSC_EpPi[20];
	  str = "-12677080620EA5e378";
	  str >> CXSC_EpPi[21];
	  str = "-1C63FD21D891DEe340";
	  str >> CXSC_EpPi[22];
	  str = "-1FB6165FFF8730e309";
	  str >> CXSC_EpPi[23];
	  str = "-177AB93E2523EFe2D3";
	  str >> CXSC_EpPi[24];
	  str = "+16D78E0B522E2Ce29C";
	  str >> CXSC_EpPi[25];
	  str = "-17473D7DD61EBEe266";
	  str >> CXSC_EpPi[26];
	  str = "-1F082665C53E27e22A";
	  str >> CXSC_EpPi[27];
	  str = "-17CDF823ACB5D6e1F3";
	  str >> CXSC_EpPi[28];
	  str = "-1D95D856C4BF74e1BC";
	  str >> CXSC_EpPi[29];
	  str = "-1327665D26E23Ae186";
	  str >> CXSC_EpPi[30];
	  str = "-1EC5804BDCCA81e150";
	  str >> CXSC_EpPi[31];
	  str = "-1C73760E976CC3e117";
	  str >> CXSC_EpPi[32];
	  str = "-10B1DCE92BE86Ce0E0";
	  str >> CXSC_EpPi[33];
	  str = "-17372866D0A1CCe0AA";
	  str >> CXSC_EpPi[34];
	  str = "-15510B0AF58D1Ee074";
	  str >> CXSC_EpPi[35];
	  str = "+1B9820D80B02D9e03C";
	  str >> CXSC_EpPi[36];
	  str = "-17765315D853BAe002";
	  str >> CXSC_EpPi[37];
	  str = "-10000000000001e000";
	  str >> CXSC_EpPi[38];
	  str = "-10000000000000e000";
	  str >> CXSC_EpPi[39];

	  CXSC_EpPi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_EpPi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1018,y);
 } // e^Pi

 static real CXSC_EulerGamma[40]; // CXSC_EulerGamma[0],...,CXSC_EulerGamma[39]
 static bool CXSC_EulerGamma_initialized = false;

 lx_interval EulerGamma_lx_interval() throw()
// Inclusion of EulerGamma, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_EulerGamma_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+12788CFC6FB619e7FD";
	  str >> CXSC_EulerGamma[0];
	  str = "-16CB90701FBFABe7C4";
	  str >> CXSC_EulerGamma[1];
	  str = "-134A95E3133C51e78E";
	  str >> CXSC_EulerGamma[2];
	  str = "+19730064300F7De758";
	  str >> CXSC_EulerGamma[3];
	  str = "-171ECA0084E369e721";
	  str >> CXSC_EulerGamma[4];
	  str = "-1302FE2B078898e6EB";
	  str >> CXSC_EulerGamma[5];
	  str = "+192732D88415F4e6B4";
	  str >> CXSC_EulerGamma[6];
	  str = "+11056AE9132136e67E";
	  str >> CXSC_EulerGamma[7];
	  str = "-17DC6F12E630A3e648";
	  str >> CXSC_EulerGamma[8];
	  str = "+175FD4B1BD70F2e611";
	  str >> CXSC_EulerGamma[9];
	  str = "-19BC9466120C20e5DB";
	  str >> CXSC_EulerGamma[10];
	  str = "-18FD5699260EADe5A5";
	  str >> CXSC_EulerGamma[11];
	  str = "-12EA987665551Fe56E";
	  str >> CXSC_EulerGamma[12];
	  str = "-1FB159BA4A423De537";
	  str >> CXSC_EulerGamma[13];
	  str = "+1FA543D43BCC60e501";
	  str >> CXSC_EulerGamma[14];
	  str = "-1E6F04E0F639F6e4C8";
	  str >> CXSC_EulerGamma[15];
	  str = "-1A23768654F43De490";
	  str >> CXSC_EulerGamma[16];
	  str = "-14F1C5CB4F55EBe457";
	  str >> CXSC_EulerGamma[17];
	  str = "+1E71DF52EDAA7Fe41F";
	  str >> CXSC_EulerGamma[18];
	  str = "+1C398F9B427E3Fe3E8";
	  str >> CXSC_EulerGamma[19];
	  str = "+1432C7402B3D24e3AF";
	  str >> CXSC_EulerGamma[20];
	  str = "-1810CF88C5F0D1e377";
	  str >> CXSC_EulerGamma[21];
	  str = "-1E9610AE5B38C5e340";
	  str >> CXSC_EulerGamma[22];
	  str = "+18220365594965e30A";
	  str >> CXSC_EulerGamma[23];
	  str = "-11F19DA40D550De2D4";
	  str >> CXSC_EulerGamma[24];
	  str = "+1936632B38107Ee29D";
	  str >> CXSC_EulerGamma[25];
	  str = "+11569F72C8E198e267";
	  str >> CXSC_EulerGamma[26];
	  str = "-1534E2EDC357EAe230";
	  str >> CXSC_EulerGamma[27];
	  str = "-126925B08D68FFe1F7";
	  str >> CXSC_EulerGamma[28];
	  str = "-1F62F86B1E6840e1C0";
	  str >> CXSC_EulerGamma[29];
	  str = "+133239C12DF171e188";
	  str >> CXSC_EulerGamma[30];
	  str = "+17E60989B189E2e14D";
	  str >> CXSC_EulerGamma[31];
	  str = "+1F6FB023E2AA98e117";
	  str >> CXSC_EulerGamma[32];
	  str = "+16911C6C6708F0e0DF";
	  str >> CXSC_EulerGamma[33];
	  str = "-1B7B61C9327A74e0A9";
	  str >> CXSC_EulerGamma[34];
	  str = "-16E397C5D924C3e073";
	  str >> CXSC_EulerGamma[35];
	  str = "-1435DBEDE3F382e03C";
	  str >> CXSC_EulerGamma[36];
	  str = "+1C0AEA1C4BB6FFe006";
	  str >> CXSC_EulerGamma[37];
	  str = "-1000000000000Ee000";
	  str >> CXSC_EulerGamma[38];
	  str = "-1000000000000De000";
	  str >> CXSC_EulerGamma[39];

	  CXSC_EulerGamma_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_EulerGamma[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // EulerGamma

 static real CXSC_Catalan[40]; // CXSC_Catalan[0],...,CXSC_Catalan[39]
 static bool CXSC_Catalan_initialized = false;

 lx_interval Catalan_lx_interval() throw()
// Inclusion of Catalan, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Catalan_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1D4F9713E8135De7FC";
	  str >> CXSC_Catalan[0];
	  str = "+11485608B8DF4De7C3";
	  str >> CXSC_Catalan[1];
	  str = "-12F39C13BC1EC8e78D";
	  str >> CXSC_Catalan[2];
	  str = "+1C2FF8094A263Ee755";
	  str >> CXSC_Catalan[3];
	  str = "+168F335DBE5370e71F";
	  str >> CXSC_Catalan[4];
	  str = "+16291BBB16163Ee6E7";
	  str >> CXSC_Catalan[5];
	  str = "+124D663F739C43e6B1";
	  str >> CXSC_Catalan[6];
	  str = "-136A0725ED0E94e679";
	  str >> CXSC_Catalan[7];
	  str = "-1D3A26F9C06FCEe63E";
	  str >> CXSC_Catalan[8];
	  str = "-164E42486BFCD2e607";
	  str >> CXSC_Catalan[9];
	  str = "+14F358CFDEC843e5D1";
	  str >> CXSC_Catalan[10];
	  str = "-11EB82210976ABe59B";
	  str >> CXSC_Catalan[11];
	  str = "-17D31F6DF5E801e565";
	  str >> CXSC_Catalan[12];
	  str = "+13FD19CE3E396Ae52F";
	  str >> CXSC_Catalan[13];
	  str = "-1C8CBB3852FF3Fe4F6";
	  str >> CXSC_Catalan[14];
	  str = "+1A86EB34EAD01Ae4C0";
	  str >> CXSC_Catalan[15];
	  str = "+1C68C37800513Be485";
	  str >> CXSC_Catalan[16];
	  str = "+1D46EBB334D7C9e44E";
	  str >> CXSC_Catalan[17];
	  str = "-1944C5E2711625e417";
	  str >> CXSC_Catalan[18];
	  str = "-17885C649BB92Fe3E1";
	  str >> CXSC_Catalan[19];
	  str = "+1A2A0CEE24DD91e3A9";
	  str >> CXSC_Catalan[20];
	  str = "+159AEC52EB2869e372";
	  str >> CXSC_Catalan[21];
	  str = "-1D26976389F1E1e339";
	  str >> CXSC_Catalan[22];
	  str = "+1E9AF9FF2E2FB1e302";
	  str >> CXSC_Catalan[23];
	  str = "+1E8B66677323FEe2CC";
	  str >> CXSC_Catalan[24];
	  str = "+164BBD8A306F6Ae296";
	  str >> CXSC_Catalan[25];
	  str = "-1EE36A15C4872Be25F";
	  str >> CXSC_Catalan[26];
	  str = "+1C3A35B39DC2FFe228";
	  str >> CXSC_Catalan[27];
	  str = "-1CCAF1572CDFC2e1F1";
	  str >> CXSC_Catalan[28];
	  str = "+15FB135902BFEEe1BA";
	  str >> CXSC_Catalan[29];
	  str = "+19FDE4873721BAe183";
	  str >> CXSC_Catalan[30];
	  str = "-17ABB7B5115456e14D";
	  str >> CXSC_Catalan[31];
	  str = "+1458F7F79FA825e117";
	  str >> CXSC_Catalan[32];
	  str = "-1416ED1E24CEFDe0E1";
	  str >> CXSC_Catalan[33];
	  str = "+15A6293C127A02e0A9";
	  str >> CXSC_Catalan[34];
	  str = "-1F1AABC6E5593Ce073";
	  str >> CXSC_Catalan[35];
	  str = "+1A7A8FB50B3479e03D";
	  str >> CXSC_Catalan[36];
	  str = "+1B853813268EF2e005";
	  str >> CXSC_Catalan[37];
	  str = "+10000000000003e000";
	  str >> CXSC_Catalan[38];
	  str = "+10000000000004e000";
	  str >> CXSC_Catalan[39];

	  CXSC_Catalan_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Catalan[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // Catalan

 static real CXSC_sqrt5[40]; // CXSC_sqrt5[0],...,CXSC_sqrt5[39]
 static bool CXSC_sqrt5_initialized = false;

 lx_interval sqrt5_lx_interval() throw()
// Inclusion of sqrt(5), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_sqrt5_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+11E3779B97F4A8e7FD";
	  str >> CXSC_sqrt5[0];
	  str = "-1F506319FCFD19e7C6";
	  str >> CXSC_sqrt5[1];
	  str = "+1B906821044ED8e790";
	  str >> CXSC_sqrt5[2];
	  str = "-18BB1B5C0F272Ce758";
	  str >> CXSC_sqrt5[3];
	  str = "+11D0C18E952768e721";
	  str >> CXSC_sqrt5[4];
	  str = "-1E9D585B0901F9e6E8";
	  str >> CXSC_sqrt5[5];
	  str = "-1C7DD252073EC0e6B2";
	  str >> CXSC_sqrt5[6];
	  str = "-1FCEF21EDAF7FAe67C";
	  str >> CXSC_sqrt5[7];
	  str = "+160EB25D20799Be63E";
	  str >> CXSC_sqrt5[8];
	  str = "-1C90F95285168Fe605";
	  str >> CXSC_sqrt5[9];
	  str = "+1E1DFA160E75BCe5CF";
	  str >> CXSC_sqrt5[10];
	  str = "-10A08E66CB368Ce593";
	  str >> CXSC_sqrt5[11];
	  str = "+1C5371682CADD1e55D";
	  str >> CXSC_sqrt5[12];
	  str = "-1998100220F4EDe526";
	  str >> CXSC_sqrt5[13];
	  str = "+1C6771A0968663e4F0";
	  str >> CXSC_sqrt5[14];
	  str = "+1DFB9E3C86CA7Ce4BA";
	  str >> CXSC_sqrt5[15];
	  str = "-18AE38ED5304B1e483";
	  str >> CXSC_sqrt5[16];
	  str = "+182A5FEC507706e44D";
	  str >> CXSC_sqrt5[17];
	  str = "-1B5191A18C5647e415";
	  str >> CXSC_sqrt5[18];
	  str = "+1F3AA4DB287AE4e3DD";
	  str >> CXSC_sqrt5[19];
	  str = "+13CCCBA48E9CF3e3A7";
	  str >> CXSC_sqrt5[20];
	  str = "-1BA6DC6F5C5A29e36F";
	  str >> CXSC_sqrt5[21];
	  str = "-1C87E14FE4B628e339";
	  str >> CXSC_sqrt5[22];
	  str = "+19D3E854210678e303";
	  str >> CXSC_sqrt5[23];
	  str = "+1D54807166A5B3e2CC";
	  str >> CXSC_sqrt5[24];
	  str = "+1D6987A242DB8De296";
	  str >> CXSC_sqrt5[25];
	  str = "+11DEC88A541BA8e260";
	  str >> CXSC_sqrt5[26];
	  str = "+1B91971F345A99e228";
	  str >> CXSC_sqrt5[27];
	  str = "-1758726EB11AAFe1F2";
	  str >> CXSC_sqrt5[28];
	  str = "+115613474B3372e1BA";
	  str >> CXSC_sqrt5[29];
	  str = "+150CCED228BEBFe184";
	  str >> CXSC_sqrt5[30];
	  str = "-13BFAB17A99623e14E";
	  str >> CXSC_sqrt5[31];
	  str = "-17397D0B07EFCEe115";
	  str >> CXSC_sqrt5[32];
	  str = "+157D728069CE1Ee0DC";
	  str >> CXSC_sqrt5[33];
	  str = "+10B046030D58E1e0A6";
	  str >> CXSC_sqrt5[34];
	  str = "-1A8870A118821Ce06E";
	  str >> CXSC_sqrt5[35];
	  str = "-1BD96261F224FBe038";
	  str >> CXSC_sqrt5[36];
	  str = "+1FDBEC758D20C8e002";
	  str >> CXSC_sqrt5[37];
	  str = "+10000000000000e000";
	  str >> CXSC_sqrt5[38];
	  str = "+10000000000001e000";
	  str >> CXSC_sqrt5[39];

	  CXSC_sqrt5_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_sqrt5[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // sqrt(5)

 static real CXSC_sqrt7[40]; // CXSC_sqrt7[0],...,CXSC_sqrt7[39]
 static bool CXSC_sqrt7_initialized = false;

 lx_interval sqrt7_lx_interval() throw()
// Inclusion of sqrt(7), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_sqrt7_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+152A7FA9D2F8EAe7FD";
	  str >> CXSC_sqrt7[0];
	  str = "-121C62B033C079e7C7";
	  str >> CXSC_sqrt7[1];
	  str = "-177CAAD6200612e78E";
	  str >> CXSC_sqrt7[2];
	  str = "-1EFA880DC72D64e756";
	  str >> CXSC_sqrt7[3];
	  str = "-171D206D5B1A4Ce71C";
	  str >> CXSC_sqrt7[4];
	  str = "+119392FA9B0494e6E3";
	  str >> CXSC_sqrt7[5];
	  str = " +17BB8A64890057e6AA";
	  str >> CXSC_sqrt7[6];
	  str = "-17E89300383DDEe674";
	  str >> CXSC_sqrt7[7];
	  str = "+130FB7AF68A6FBe63E";
	  str >> CXSC_sqrt7[8];
	  str = "+1322281D303D36e606";
	  str >> CXSC_sqrt7[9];
	  str = "+1996109A16D3B1e5D0";
	  str >> CXSC_sqrt7[10];
	  str = "+1F239C301DFBB4e599";
	  str >> CXSC_sqrt7[11];
	  str = "-1B5CA40AB771A2e560";
	  str >> CXSC_sqrt7[12];
	  str = "-1675711487FEAAe527";
	  str >> CXSC_sqrt7[13];
	  str = "+122CB7FA26ABA5e4F1";
	  str >> CXSC_sqrt7[14];
	  str = "+1059211B7D5398e4BA";
	  str >> CXSC_sqrt7[15];
	  str = "-10F15BFA46EB7Fe484";
	  str >> CXSC_sqrt7[16];
	  str = "+15AB71566CE72Be44E";
	  str >> CXSC_sqrt7[17];
	  str = "-1386BDCA3845C7e417";
	  str >> CXSC_sqrt7[18];
	  str = "+158978E1A883F0e3E1";
	  str >> CXSC_sqrt7[19];
	  str = "+1F8A772604CAF1e3AB";
	  str >> CXSC_sqrt7[20];
	  str = "-169730DF195425e374";
	  str >> CXSC_sqrt7[21];
	  str = "+121341C1C58312e33D";
	  str >> CXSC_sqrt7[22];
	  str = "+1117A6A8D5A423e307";
	  str >> CXSC_sqrt7[23];
	  str = "-12D6F266C4AE34e2D0";
	  str >> CXSC_sqrt7[24];
	  str = "-1B8D3340BC4497e299";
	  str >> CXSC_sqrt7[25];
	  str = "-1E0712BA5C43CFe263";
	  str >> CXSC_sqrt7[26];
	  str = "+14BF59F22FDCEDe22D";
	  str >> CXSC_sqrt7[27];
	  str = "-189028E9147846e1F3";
	  str >> CXSC_sqrt7[28];
	  str = "-1BBE8044991D23e1BB";
	  str >> CXSC_sqrt7[29];
	  str = "+1B07A12AB80C0Ae185";
	  str >> CXSC_sqrt7[30];
	  str = "+1970A29C5148A7e14E";
	  str >> CXSC_sqrt7[31];
	  str = "-18DBF939495F77e117";
	  str >> CXSC_sqrt7[32];
	  str = "-1865ECB3F101B4e0E1";
	  str >> CXSC_sqrt7[33];
	  str = "-1DCDFCFE6050DFe0A4";
	  str >> CXSC_sqrt7[34];
	  str = "+1E635021BBF848e06B";
	  str >> CXSC_sqrt7[35];
	  str = "+15C53CF344AF94e035";
	  str >> CXSC_sqrt7[36];
	  str = "+124CB182B20E69e000";
	  str >> CXSC_sqrt7[37];
	  str = "-10000000000001e000";
	  str >> CXSC_sqrt7[38];
	  str = "-10000000000000e000";
	  str >> CXSC_sqrt7[39];

	  CXSC_sqrt7_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_sqrt7[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // sqrt(7)

 static real CXSC_One_m[3]; // CXSC_One_m[0], ... ,CXSC_One_m[2]
 static bool CXSC_One_m_initialized = false;

 lx_interval One_m_lx_interval() throw()
// Inclusion of 2^(-1023)*(2^(+1023) - 2^(-1074)) = 1-2^(-2097), 
// Blomquist, 26.07.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 2;

  if (!CXSC_One_m_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+10000000000000e7FE";
	  str >> CXSC_One_m[0];
	  str = "-10000000000001e000";
	  str >> CXSC_One_m[1];
	  str = "-10000000000001e000";
	  str >> CXSC_One_m[2];
	  CXSC_One_m_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_One_m[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // One_m_lx_interval()

 static real CXSC_One_p[3]; // CXSC_One_p[0], ... ,CXSC_One_p[2]
 static bool CXSC_One_p_initialized = false;

 lx_interval One_p_lx_interval() throw()
// Inclusion of 2^(-1023)*(2^(+1023) + 2^(-1074)) = 1+2^(-2097), 
// Blomquist, 27.07.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 2;

  if (!CXSC_One_p_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+10000000000000e7FE";
	  str >> CXSC_One_p[0];
	  str = "+10000000000001e000";
	  str >> CXSC_One_p[1];
	  str = "+10000000000001e000";
	  str >> CXSC_One_p[2];
	  CXSC_One_p_initialized = true;
	  std::cout << RestoreOpt;
  }
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_One_p[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // One_p_lx_interval()

 static real CXSC_Ln10r[40]; // CXSC_Ln10r[0], ... ,CXSC_Ln10r[39]
 static bool CXSC_Ln10r_initialized = false;

 lx_interval Ln10r_lx_interval() throw()
// Inclusion of 1/ln(10), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ln10r_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1BCB7B1526E50Ee7FC";
	  str >> CXSC_Ln10r[0];
	  str = "+195355BAAAFAD3e7C5";
	  str >> CXSC_Ln10r[1];
	  str = "+1EE191F71A3012e78E";
	  str >> CXSC_Ln10r[2];
	  str = "+17268808E8FCB5e757";
	  str >> CXSC_Ln10r[3];
	  str = "+13DE3A94F1D509e71F";
	  str >> CXSC_Ln10r[4];
	  str = "+1DF42805E7E524e6E8";
	  str >> CXSC_Ln10r[5];
	  str = "+11AAC96323250Be6B2";
	  str >> CXSC_Ln10r[6];
	  str = "-1CE63884C058E4e67C";
	  str >> CXSC_Ln10r[7];
	  str = "-1A1C82EA3969BAe646";
	  str >> CXSC_Ln10r[8];
	  str = "+1B4F6686AD7A33e610";
	  str >> CXSC_Ln10r[9];
	  str = "-1B97C8035FFC70e5DA";
	  str >> CXSC_Ln10r[10];
	  str = "+1630771369962Ee59F";
	  str >> CXSC_Ln10r[11];
	  str = "-1E15BD37B295AFe569";
	  str >> CXSC_Ln10r[12];
	  str = "-132484B432318Be533";
	  str >> CXSC_Ln10r[13];
	  str = "+15430212AE68C0e4FD";
	  str >> CXSC_Ln10r[14];
	  str = "+1351923B322731e4C7";
	  str >> CXSC_Ln10r[15];
	  str = "+11F934D794D64Fe491";
	  str >> CXSC_Ln10r[16];
	  str = "+13E4B475D9FF20e45A";
	  str >> CXSC_Ln10r[17];
	  str = "+185D9B63ED9A24e424";
	  str >> CXSC_Ln10r[18];
	  str = "+1ADC650C65E948e3ED";
	  str >> CXSC_Ln10r[19];
	  str = "-1149FBC70C04EAe3B6";
	  str >> CXSC_Ln10r[20];
	  str = "+1056270A8CDF9Ce380";
	  str >> CXSC_Ln10r[21];
	  str = "-1D339476D1076Fe34A";
	  str >> CXSC_Ln10r[22];
	  str = "-1635343A8B5E85e314";
	  str >> CXSC_Ln10r[23];
	  str = "+1DB78151377249e2DC";
	  str >> CXSC_Ln10r[24];
	  str = "+14DD9BD4601639e2A6";
	  str >> CXSC_Ln10r[25];
	  str = "+1D545BF0E9F470e26F";
	  str >> CXSC_Ln10r[26];
	  str = "+17CC4CE204C9F6e239";
	  str >> CXSC_Ln10r[27];
	  str = "+10CED2851AF1ABe200";
	  str >> CXSC_Ln10r[28];
	  str = "-16C9D9EB4EF234e1C9";
	  str >> CXSC_Ln10r[29];
	  str = "-1D67966CCC4205e192";
	  str >> CXSC_Ln10r[30];
	  str = "-1A83D3193C0A22e15C";
	  str >> CXSC_Ln10r[31];
	  str = "-17066494F5F3BEe126";
	  str >> CXSC_Ln10r[32];
	  str = "+12A7753E2FBCACe0F0";
	  str >> CXSC_Ln10r[33];
	  str = "+1FE8D8E367317Be0BA";
	  str >> CXSC_Ln10r[34];
	  str = "+12BB8A7F6B3745e084";
	  str >> CXSC_Ln10r[35];
	  str = "+1D906BB4F052BDe04B";
	  str >> CXSC_Ln10r[36];
	  str = "+1FEED057798219e013";
	  str >> CXSC_Ln10r[37];
	  str = "+100000000010E1e000";
	  str >> CXSC_Ln10r[38];
	  str = "+100000000010E2e000";
	  str >> CXSC_Ln10r[39];

	  CXSC_Ln10r_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ln10r[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // 1/ln(10)

// **************************************************************************

 static real CXSC_Pid4[40]; // CXSC_Pid4[0], ... ,CXSC_Pid4[39]
 static bool CXSC_Pid4_initialized = false;

 lx_interval Pid4_lx_interval() throw()
// Inclusion of 1/ln(10), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Pid4_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1921FB54442D18e7FC";
	  str >> CXSC_Pid4[0];
	  str = "+11A62633145C07e7C6";
	  str >> CXSC_Pid4[1];
	  str = "-1F1976B7ED8FBCe78E";
	  str >> CXSC_Pid4[2];
	  str = "+14CF98E804177De758";
	  str >> CXSC_Pid4[3];
	  str = "+131D89CD9128A5e722";
	  str >> CXSC_Pid4[4];
	  str = "+10F31C6809BBDFe6E8";
	  str >> CXSC_Pid4[5];
	  str = "+1519B3CD3A431Be6B1";
	  str >> CXSC_Pid4[6];
	  str = "+18158536F92F8Ae67A";
	  str >> CXSC_Pid4[7];
	  str = "+1BA7F09AB6B6A9e642";
	  str >> CXSC_Pid4[8];
	  str = "-1EDD0DBD2544CFe60A";
	  str >> CXSC_Pid4[9];
	  str = "+179FB1BD1310BAe5D3";
	  str >> CXSC_Pid4[10];
	  str = "+1A637ED6B0BFF6e59D";
	  str >> CXSC_Pid4[11];
	  str = "-1A485FCA40908Ee566";
	  str >> CXSC_Pid4[12];
	  str = "-1E501295D98169e52F";
	  str >> CXSC_Pid4[13];
	  str = "-1160DBEE83B4E0e4F9";
	  str >> CXSC_Pid4[14];
	  str = "-19B6D799AE131Ce4C1";
	  str >> CXSC_Pid4[15];
	  str = "+16CF70801F2E28e48B";
	  str >> CXSC_Pid4[16];
	  str = "+163BF0598DA483e455";
	  str >> CXSC_Pid4[17];
	  str = "+1871574E69A459e41F";
	  str >> CXSC_Pid4[18];
	  str = "-15C0B6CC281F27e3E3";
	  str >> CXSC_Pid4[19];
	  str = "+15D23DCA3AD962e3AD";
	  str >> CXSC_Pid4[20];
	  str = "-1CE8654EFBD56Ae376";
	  str >> CXSC_Pid4[21];
	  str = "-1184AB5BE23DA6e33F";
	  str >> CXSC_Pid4[22];
	  str = "+166D670C354E4Be309";
	  str >> CXSC_Pid4[23];
	  str = "-10D9FEC3A2E4FEe2D3";
	  str >> CXSC_Pid4[24];
	  str = "+1943042F86520Ce29C";
	  str >> CXSC_Pid4[25];
	  str = "-1B9D1C931C41C6e265";
	  str >> CXSC_Pid4[26];
	  str = "-188D3E7F179FC6e22D";
	  str >> CXSC_Pid4[27];
	  str = "-1361F1744FE176e1F7";
	  str >> CXSC_Pid4[28];
	  str = "+1F6B8ABBE0DE99e1C0";
	  str >> CXSC_Pid4[29];
	  str = "-169B10EA1A04B5e18A";
	  str >> CXSC_Pid4[30];
	  str = "-14FD1CF8CD56D0e154";
	  str >> CXSC_Pid4[31];
	  str = "-18AB54A8D7516Fe11E";
	  str >> CXSC_Pid4[32];
	  str = "+186263E8144056e0E7";
	  str >> CXSC_Pid4[33];
	  str = "-1AE34AEAAA77A5e0B0";
	  str >> CXSC_Pid4[34];
	  str = "+16998B8682283De07A";
	  str >> CXSC_Pid4[35];
	  str = "+19D42A90D5EF8Ee042";
	  str >> CXSC_Pid4[36];
	  str = "+174C9D9F70A08Be00C";
	  str >> CXSC_Pid4[37];
	  str = "+100000000000DBe000";
	  str >> CXSC_Pid4[38];
	  str = "+100000000000DCe000";
	  str >> CXSC_Pid4[39];

	  CXSC_Pid4_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pid4[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // pi/4

// *************************************************************************

 static real CXSC_Pid2[40]; // CXSC_Pid2[0], ... ,CXSC_Pid2[39]
 static bool CXSC_Pid2_initialized = false;

 lx_interval Pid2_lx_interval() throw()
// Inclusion of 1/ln(10), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Pid2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1921FB54442D18e7FC";
	  str >> CXSC_Pid2[0];
	  str = "+11A62633145C07e7C6";
	  str >> CXSC_Pid2[1];
	  str = "-1F1976B7ED8FBCe78E";
	  str >> CXSC_Pid2[2];
	  str = "+14CF98E804177De758";
	  str >> CXSC_Pid2[3];
	  str = "+131D89CD9128A5e722";
	  str >> CXSC_Pid2[4];
	  str = "+10F31C6809BBDFe6E8";
	  str >> CXSC_Pid2[5];
	  str = "+1519B3CD3A431Be6B1";
	  str >> CXSC_Pid2[6];
	  str = "+18158536F92F8Ae67A";
	  str >> CXSC_Pid2[7];
	  str = "+1BA7F09AB6B6A9e642";
	  str >> CXSC_Pid2[8];
	  str = "-1EDD0DBD2544CFe60A";
	  str >> CXSC_Pid2[9];
	  str = "+179FB1BD1310BAe5D3";
	  str >> CXSC_Pid2[10];
	  str = "+1A637ED6B0BFF6e59D";
	  str >> CXSC_Pid2[11];
	  str = "-1A485FCA40908Ee566";
	  str >> CXSC_Pid2[12];
	  str = "-1E501295D98169e52F";
	  str >> CXSC_Pid2[13];
	  str = "-1160DBEE83B4E0e4F9";
	  str >> CXSC_Pid2[14];
	  str = "-19B6D799AE131Ce4C1";
	  str >> CXSC_Pid2[15];
	  str = "+16CF70801F2E28e48B";
	  str >> CXSC_Pid2[16];
	  str = "+163BF0598DA483e455";
	  str >> CXSC_Pid2[17];
	  str = "+1871574E69A459e41F";
	  str >> CXSC_Pid2[18];
	  str = "-15C0B6CC281F27e3E3";
	  str >> CXSC_Pid2[19];
	  str = "+15D23DCA3AD962e3AD";
	  str >> CXSC_Pid2[20];
	  str = "-1CE8654EFBD56Ae376";
	  str >> CXSC_Pid2[21];
	  str = "-1184AB5BE23DA6e33F";
	  str >> CXSC_Pid2[22];
	  str = "+166D670C354E4Be309";
	  str >> CXSC_Pid2[23];
	  str = "-10D9FEC3A2E4FEe2D3";
	  str >> CXSC_Pid2[24];
	  str = "+1943042F86520Ce29C";
	  str >> CXSC_Pid2[25];
	  str = "-1B9D1C931C41C6e265";
	  str >> CXSC_Pid2[26];
	  str = "-188D3E7F179FC6e22D";
	  str >> CXSC_Pid2[27];
	  str = "-1361F1744FE176e1F7";
	  str >> CXSC_Pid2[28];
	  str = "+1F6B8ABBE0DE99e1C0";
	  str >> CXSC_Pid2[29];
	  str = "-169B10EA1A04B5e18A";
	  str >> CXSC_Pid2[30];
	  str = "-14FD1CF8CD56D0e154";
	  str >> CXSC_Pid2[31];
	  str = "-18AB54A8D7516Fe11E";
	  str >> CXSC_Pid2[32];
	  str = "+186263E8144056e0E7";
	  str >> CXSC_Pid2[33];
	  str = "-1AE34AEAAA77A5e0B0";
	  str >> CXSC_Pid2[34];
	  str = "+16998B8682283De07A";
	  str >> CXSC_Pid2[35];
	  str = "+19D42A90D5EF8Ee042";
	  str >> CXSC_Pid2[36];
	  str = "+174C9D9F70A08Be00C";
	  str >> CXSC_Pid2[37];
	  str = "+100000000000DBe000";
	  str >> CXSC_Pid2[38];
	  str = "+100000000000DCe000";
	  str >> CXSC_Pid2[39];

	  CXSC_Pid2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pid2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // pi/2

// **************************************************************************

 static real CXSC_Pi2[40]; // CXSC_Pi2[0], ... ,CXSC_Pi2[39]
 static bool CXSC_Pi2_initialized = false;

 lx_interval Pi2_lx_interval() throw()
// Inclusion of 2*pi, Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;  // alt
  if (!CXSC_Pi2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1921FB54442D18e7FC";
	  str >> CXSC_Pi2[0];
	  str = "+11A62633145C07e7C6";
	  str >> CXSC_Pi2[1];
	  str = "-1F1976B7ED8FBCe78E";
	  str >> CXSC_Pi2[2];
	  str = "+14CF98E804177De758";
	  str >> CXSC_Pi2[3];
	  str = "+131D89CD9128A5e722";
	  str >> CXSC_Pi2[4];
	  str = "+10F31C6809BBDFe6E8";
	  str >> CXSC_Pi2[5];
	  str = "+1519B3CD3A431Be6B1";
	  str >> CXSC_Pi2[6];
	  str = "+18158536F92F8Ae67A";
	  str >> CXSC_Pi2[7];
	  str = "+1BA7F09AB6B6A9e642";
	  str >> CXSC_Pi2[8];
	  str = "-1EDD0DBD2544CFe60A";
	  str >> CXSC_Pi2[9];
	  str = "+179FB1BD1310BAe5D3";
	  str >> CXSC_Pi2[10];
	  str = "+1A637ED6B0BFF6e59D";
	  str >> CXSC_Pi2[11];
	  str = "-1A485FCA40908Ee566";
	  str >> CXSC_Pi2[12];
	  str = "-1E501295D98169e52F";
	  str >> CXSC_Pi2[13];
	  str = "-1160DBEE83B4E0e4F9";
	  str >> CXSC_Pi2[14];
	  str = "-19B6D799AE131Ce4C1";
	  str >> CXSC_Pi2[15];
	  str = "+16CF70801F2E28e48B";
	  str >> CXSC_Pi2[16];
	  str = "+163BF0598DA483e455";
	  str >> CXSC_Pi2[17];
	  str = "+1871574E69A459e41F";
	  str >> CXSC_Pi2[18];
	  str = "-15C0B6CC281F27e3E3";
	  str >> CXSC_Pi2[19];
	  str = "+15D23DCA3AD962e3AD";
	  str >> CXSC_Pi2[20];
	  str = "-1CE8654EFBD56Ae376";
	  str >> CXSC_Pi2[21];
	  str = "-1184AB5BE23DA6e33F";
	  str >> CXSC_Pi2[22];
	  str = "+166D670C354E4Be309";
	  str >> CXSC_Pi2[23];
	  str = "-10D9FEC3A2E4FEe2D3";
	  str >> CXSC_Pi2[24];
	  str = "+1943042F86520Ce29C";
	  str >> CXSC_Pi2[25];
	  str = "-1B9D1C931C41C6e265";
	  str >> CXSC_Pi2[26];
	  str = "-188D3E7F179FC6e22D";
	  str >> CXSC_Pi2[27];
	  str = "-1361F1744FE176e1F7";
	  str >> CXSC_Pi2[28];
	  str = "+1F6B8ABBE0DE99e1C0";
	  str >> CXSC_Pi2[29];
	  str = "-169B10EA1A04B5e18A";
	  str >> CXSC_Pi2[30];
	  str = "-14FD1CF8CD56D0e154";
	  str >> CXSC_Pi2[31];
	  str = "-18AB54A8D7516Fe11E";
	  str >> CXSC_Pi2[32];
	  str = "+186263E8144056e0E7";
	  str >> CXSC_Pi2[33];
	  str = "-1AE34AEAAA77A5e0B0";
	  str >> CXSC_Pi2[34];
	  str = "+16998B8682283De07A";
	  str >> CXSC_Pi2[35];
	  str = "+19D42A90D5EF8Ee042";
	  str >> CXSC_Pi2[36];
	  str = "+174C9D9F70A08Be00C";
	  str >> CXSC_Pi2[37];
	  str = "+100000000000DBe000";
	  str >> CXSC_Pi2[38];
	  str = "+100000000000DCe000";
	  str >> CXSC_Pi2[39];

	  CXSC_Pi2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pi2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1019,y);
 } // 2pi

// **************************************************************************

 static real CXSC_Pi2r[40]; // CXSC_Pi2r[0], ... ,CXSC_Pi2r[39]
 static bool CXSC_Pi2r_initialized = false;

 lx_interval Pi2r_lx_interval() throw()
// Inclusion of 1/Pi, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Pi2r_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+145F306DC9C883e7FC";
	  str >> CXSC_Pi2r[0];
	  str = "-16B01EC5417056e7C6";
	  str >> CXSC_Pi2r[1];
	  str = "-16447E493AD4CEe790";
	  str >> CXSC_Pi2r[2];
	  str = "+1E21C820FF28B2e75A";
	  str >> CXSC_Pi2r[3];
	  str = "-1508510EA79237e723";
	  str >> CXSC_Pi2r[4];
	  str = "+1B8E909374B802e6EB";
	  str >> CXSC_Pi2r[5];
	  str = "-1B6D115F62E6DEe6B5";
	  str >> CXSC_Pi2r[6];
	  str = "-180F10A71A76B3e67E";
	  str >> CXSC_Pi2r[7];
	  str = "+1CFBA208D7D4BBe647";
	  str >> CXSC_Pi2r[8];
	  str = "-12EDEC598E3F65e60F";
	  str >> CXSC_Pi2r[9];
	  str = "-1741037D8CDC54e5D8";
	  str >> CXSC_Pi2r[10];
	  str = "+1CC1A99CFA4E42e5A2";
	  str >> CXSC_Pi2r[11];
	  str = "+17E2EF7E4A0EC8e56B";
	  str >> CXSC_Pi2r[12];
	  str = "-1DA00087E99FC0e52F";
	  str >> CXSC_Pi2r[13];
	  str = "-10D0EE74A5F593e4F9";
	  str >> CXSC_Pi2r[14];
	  str = "+1F6D367ECF27CBe4C1";
	  str >> CXSC_Pi2r[15];
	  str = "+136E9E8C7ECD3De488";
	  str >> CXSC_Pi2r[16];
	  str = "-100AE9456C229Ce452";
	  str >> CXSC_Pi2r[17];
	  str = "-141A0E84C2F8C6e419";
	  str >> CXSC_Pi2r[18];
	  str = "-10EB5ADA2B2809e3E0";
	  str >> CXSC_Pi2r[19];
	  str = "-10277039517BD5e3AA";
	  str >> CXSC_Pi2r[20];
	  str = "+198237E3DB5D60e36E";
	  str >> CXSC_Pi2r[21];
	  str = "-1E6087BECA1794e338";
	  str >> CXSC_Pi2r[22];
	  str = "+1DA9E391615EE6e301";
	  str >> CXSC_Pi2r[23];
	  str = "+1B086599855F15e2C9";
	  str >> CXSC_Pi2r[24];
	  str = "-17E5EFDC8009E0e293";
	  str >> CXSC_Pi2r[25];
	  str = "+135CC9CC418185e25B";
	  str >> CXSC_Pi2r[26];
	  str = "+156CA73A8C960Ee225";
	  str >> CXSC_Pi2r[27];
	  str = "+13DE04635A3E21e1EE";
	  str >> CXSC_Pi2r[28];
	  str = "-18F260C88C5FDBe1B7";
	  str >> CXSC_Pi2r[29];
	  str = "-157CA63B89746Ae181";
	  str >> CXSC_Pi2r[30];
	  str = "+1CA6DDAF44D157e149";
	  str >> CXSC_Pi2r[31];
	  str = "+19053EA5FF0705e111";
	  str >> CXSC_Pi2r[32];
	  str = "+1FBF19F419616Fe0DA";
	  str >> CXSC_Pi2r[33];
	  str = "+13E60C9F6EF0CFe0A3";
	  str >> CXSC_Pi2r[34];
	  str = "+126EF6B1E5EF8Ae06D";
	  str >> CXSC_Pi2r[35];
	  str = "-18BC1946A1B01Ce034";
	  str >> CXSC_Pi2r[36];
	  str = "-12780EDE6F8384e000";
	  str >> CXSC_Pi2r[37];
	  str = "+10000000000000e000";
	  str >> CXSC_Pi2r[38];
	  str = "+10000000000001e000";
	  str >> CXSC_Pi2r[39];

	  CXSC_Pi2r_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pi2r[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1024,y);
 } // 1/(2*pi)

// **************************************************************************

 static real CXSC_Pip2[40]; // CXSC_Pip2[0], ... ,CXSC_Pip2[39]
 static bool CXSC_Pip2_initialized = false;

 lx_interval Pip2_lx_interval() throw()
// Inclusion of 1/Pi, Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Pip2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+13BD3CC9BE45DEe7FD";
	  str >> CXSC_Pip2[0];
	  str = "+1692B71366CC04e7C7";
	  str >> CXSC_Pip2[1];
	  str = "+18358E10ACD480e791";
	  str >> CXSC_Pip2[2];
	  str = "-1F2F5DD7997DDFe75A";
	  str >> CXSC_Pip2[3];
	  str = "+129E39B47B884Ee71F";
	  str >> CXSC_Pip2[4];
	  str = "-12CF7459DD5DAFe6E9";
	  str >> CXSC_Pip2[5];
	  str = "-11842F87B5FE0Fe6B3";
	  str >> CXSC_Pip2[6];
	  str = "+1FFD8A79616A21e67D";
	  str >> CXSC_Pip2[7];
	  str = "+12492A6663E899e647";
	  str >> CXSC_Pip2[8];
	  str = "-1A15F4352CC511e610";
	  str >> CXSC_Pip2[9];
	  str = "-1301AA1792FF3Ce5D9";
	  str >> CXSC_Pip2[10];
	  str = "+122B6F31626EFEe5A3";
	  str >> CXSC_Pip2[11];
	  str = "+1B317FA13BDD8Fe56D";
	  str >> CXSC_Pip2[12];
	  str = "+16F83B49040075e537";
	  str >> CXSC_Pip2[13];
	  str = "-1B1890A945FE17e501";
	  str >> CXSC_Pip2[14];
	  str = "+12DCD389B96CDBe4CB";
	  str >> CXSC_Pip2[15];
	  str = "-1743F5DDE2F157e492";
	  str >> CXSC_Pip2[16];
	  str = "-153F96FFD4AEB5e45B";
	  str >> CXSC_Pip2[17];
	  str = "+13CD6F5847D569e423";
	  str >> CXSC_Pip2[18];
	  str = "+1471E79A7B0882e3EC";
	  str >> CXSC_Pip2[19];
	  str = "-14C5022456E37Ae3B6";
	  str >> CXSC_Pip2[20];
	  str = "-1471463BD938A3e380";
	  str >> CXSC_Pip2[21];
	  str = "+13EABA147FEB41e349";
	  str >> CXSC_Pip2[22];
	  str = "-1D7FBCA9B23073e312";
	  str >> CXSC_Pip2[23];
	  str = "-17B06B8196DD15e2DC";
	  str >> CXSC_Pip2[24];
	  str = "-13A91786954EA1e2A6";
	  str >> CXSC_Pip2[25];
	  str = "+1C841C914201E8e26F";
	  str >> CXSC_Pip2[26];
	  str = "-1BDD5D6465807De239";
	  str >> CXSC_Pip2[27];
	  str = "-1BD8C694B35945e202";
	  str >> CXSC_Pip2[28];
	  str = "+181914426DA9A7e1CC";
	  str >> CXSC_Pip2[29];
	  str = "-1C83BE2430C1FEe192";
	  str >> CXSC_Pip2[30];
	  str = "+16530E2CE920C0e157";
	  str >> CXSC_Pip2[31];
	  str = "+103B8F2850B82Ee121";
	  str >> CXSC_Pip2[32];
	  str = "-1F116B01D43595e0EB";
	  str >> CXSC_Pip2[33];
	  str = "-13AF6BD210759Fe0B5";
	  str >> CXSC_Pip2[34];
	  str = "+1F82CDE6A1FFF3e07E";
	  str >> CXSC_Pip2[35];
	  str = "+1C8EF198F1ACD2e048";
	  str >> CXSC_Pip2[36];
	  str = "+18077590C18251e011";
	  str >> CXSC_Pip2[37];
	  str = "-10000000002C7Ee000";
	  str >> CXSC_Pip2[38];
	  str = "-10000000002C7De000";
	  str >> CXSC_Pip2[39];

	  CXSC_Pip2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Pip2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1019,y);
 } // pi^2

// **************************************************************************

 static real CXSC_Sqrt2r[40]; // CXSC_Sqrt2r[0], ... ,CXSC_Sqrt2r[39]
 static bool CXSC_Sqrt2r_initialized = false;

 lx_interval Sqrt2r_lx_interval() throw()
// Inclusion of sqrt(2), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt2r_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+16A09E667F3BCDe7FC";
	  str >> CXSC_Sqrt2r[0];
	  str = "-1BDD3413B26456e7C6";
	  str >> CXSC_Sqrt2r[1];
	  str = "+157D3E3ADEC175e790";
	  str >> CXSC_Sqrt2r[2];
	  str = "+12775099DA2F59e758";
	  str >> CXSC_Sqrt2r[3];
	  str = "+160CCE64552BF2e71F";
	  str >> CXSC_Sqrt2r[4];
	  str = "+1821D5C5161D46e6E6";
	  str >> CXSC_Sqrt2r[5];
	  str = "-1C032046F8498Ee6B0";
	  str >> CXSC_Sqrt2r[6];
	  str = "+1EE950BC8738F7e678";
	  str >> CXSC_Sqrt2r[7];
	  str = "-1AC3FDBC64E103e642";
	  str >> CXSC_Sqrt2r[8];
	  str = "+13B469101743A1e60A";
	  str >> CXSC_Sqrt2r[9];
	  str = "+15E3E9CA60B38Ce5D4";
	  str >> CXSC_Sqrt2r[10];
	  str = "+11BC337BCAB1BDe599";
	  str >> CXSC_Sqrt2r[11];
	  str = "-1BBA5DEE9D6E7De563";
	  str >> CXSC_Sqrt2r[12];
	  str = "-1438DD083B1CC4e52D";
	  str >> CXSC_Sqrt2r[13];
	  str = "+1B56A28E2EDFA7e4F7";
	  str >> CXSC_Sqrt2r[14];
	  str = "+1CCB2A634331F4e4C1";
	  str >> CXSC_Sqrt2r[15];
	  str = "-1BD9056876F83Ee48A";
	  str >> CXSC_Sqrt2r[16];
	  str = "-1234FA22AB6BEFe454";
	  str >> CXSC_Sqrt2r[17];
	  str = "+19040CA4A81395e41D";
	  str >> CXSC_Sqrt2r[18];
	  str = "-15249C0BC4082De3E7";
	  str >> CXSC_Sqrt2r[19];
	  str = "+13A02CEBC93E0Ce3B1";
	  str >> CXSC_Sqrt2r[20];
	  str = "+109936AF354A2Ee37B";
	  str >> CXSC_Sqrt2r[21];
	  str = "-1AE4730CBE4908e345";
	  str >> CXSC_Sqrt2r[22];
	  str = "+11B6380826E010e30E";
	  str >> CXSC_Sqrt2r[23];
	  str = "-1CDCAD0CCD5A16e2D5";
	  str >> CXSC_Sqrt2r[24];
	  str = "-1084BC28012BC8e29C";
	  str >> CXSC_Sqrt2r[25];
	  str = "-1C035DDECF8216e265";
	  str >> CXSC_Sqrt2r[26];
	  str = "+18907DEAA070B0e22B";
	  str >> CXSC_Sqrt2r[27];
	  str = "+1FCBDDEA2F7DC3e1F5";
	  str >> CXSC_Sqrt2r[28];
	  str = "+18C41C51757FB0e1BE";
	  str >> CXSC_Sqrt2r[29];
	  str = "-189A5B616B1381e188";
	  str >> CXSC_Sqrt2r[30];
	  str = "+165C417EFF0B88e152";
	  str >> CXSC_Sqrt2r[31];
	  str = "-1627043F832999e11A";
	  str >> CXSC_Sqrt2r[32];
	  str = "+105E5FCA017092e0E3";
	  str >> CXSC_Sqrt2r[33];
	  str = "-187A16D6A8FDCAe0AD";
	  str >> CXSC_Sqrt2r[34];
	  str = "-1838421AE0AE62e072";
	  str >> CXSC_Sqrt2r[35];
	  str = "-165073EB433984e03C";
	  str >> CXSC_Sqrt2r[36];
	  str = "+1F0A42F9DA4A6Ce006";
	  str >> CXSC_Sqrt2r[37];
	  str = "+10000000000002e000";
	  str >> CXSC_Sqrt2r[38];
	  str = "+10000000000003e000";
	  str >> CXSC_Sqrt2r[39];

	  CXSC_Sqrt2r_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt2r[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // 1/sqrt(2)

// ************************************************************************

 static real CXSC_Sqrt3d2[40]; // CXSC_Sqrt3d2[0], ... ,CXSC_Sqrt3d2[39]
 static bool CXSC_Sqrt3d2_initialized = false;

 lx_interval Sqrt3d2_lx_interval() throw()
// Inclusion of sqrt(3), Blomquist, 15.06.2007;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt3d2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1BB67AE8584CAAe7FC";
	  str >> CXSC_Sqrt3d2[0];
	  str = "+1CEC95D0B5C1E3e7C6";
	  str >> CXSC_Sqrt3d2[1];
	  str = "-1F11DB689F2CCFe78E";
	  str >> CXSC_Sqrt3d2[2];
	  str = "+13DA4798C720A6e758";
	  str >> CXSC_Sqrt3d2[3];
	  str = "+121B9169B89243e722";
	  str >> CXSC_Sqrt3d2[4];
	  str = "-1813508751212Be6E9";
	  str >> CXSC_Sqrt3d2[5];
	  str = "-1B3D547B775C1Ee6B2";
	  str >> CXSC_Sqrt3d2[6];
	  str = "-19D986D92E2F0Ae679";
	  str >> CXSC_Sqrt3d2[7];
	  str = "+1A34334CE806B6e642";
	  str >> CXSC_Sqrt3d2[8];
	  str = "+1A383B9E122E61e60C";
	  str >> CXSC_Sqrt3d2[9];
	  str = "+1C61D736F2F6F2e5D5";
	  str >> CXSC_Sqrt3d2[10];
	  str = "-10AF49233F9250e59E";
	  str >> CXSC_Sqrt3d2[11];
	  str = "-1558A109EC0523e567";
	  str >> CXSC_Sqrt3d2[12];
	  str = "+1F799D4D4FF2BCe531";
	  str >> CXSC_Sqrt3d2[13];
	  str = "-1AD7B219E34EDBe4FB";
	  str >> CXSC_Sqrt3d2[14];
	  str = "+15AB940B6677E3e4C5";
	  str >> CXSC_Sqrt3d2[15];
	  str = "-1D9B2A8203B8F0e48E";
	  str >> CXSC_Sqrt3d2[16];
	  str = "-1DB0C8975A3834e458";
	  str >> CXSC_Sqrt3d2[17];
	  str = "-1BCAAB3F6BE884e422";
	  str >> CXSC_Sqrt3d2[18];
	  str = "+14C70ADB1EC1BBe3E8";
	  str >> CXSC_Sqrt3d2[19];
	  str = "-14E1EF77987E55e3AF";
	  str >> CXSC_Sqrt3d2[20];
	  str = "-19695FC6269D28e378";
	  str >> CXSC_Sqrt3d2[21];
	  str = "+10D0652AAC5936e342";
	  str >> CXSC_Sqrt3d2[22];
	  str = "-1BD0891D370824e30C";
	  str >> CXSC_Sqrt3d2[23];
	  str = "-129B4C6252D061e2D4";
	  str >> CXSC_Sqrt3d2[24];
	  str = "+1DC9B1A4C31275e29E";
	  str >> CXSC_Sqrt3d2[25];
	  str = "+11FF9B8422294Ee267";
	  str >> CXSC_Sqrt3d2[26];
	  str = "-1E4A6AA47F3A85e231";
	  str >> CXSC_Sqrt3d2[27];
	  str = "+17043E01AA3F3De1FA";
	  str >> CXSC_Sqrt3d2[28];
	  str = "+188EF377D2D5B6e1C0";
	  str >> CXSC_Sqrt3d2[29];
	  str = "-1735E8C815F031e185";
	  str >> CXSC_Sqrt3d2[30];
	  str = "-1B89330FD8417Ce14F";
	  str >> CXSC_Sqrt3d2[31];
	  str = "+16D1A627670F5Ce117";
	  str >> CXSC_Sqrt3d2[32];
	  str = "+1AF43BBA8154D3e0DB";
	  str >> CXSC_Sqrt3d2[33];
	  str = "+1DA9A969A91295e0A5";
	  str >> CXSC_Sqrt3d2[34];
	  str = "-1636594394C675e06E";
	  str >> CXSC_Sqrt3d2[35];
	  str = "+1064B9DA1A3185e037";
	  str >> CXSC_Sqrt3d2[36];
	  str = "-1CE514CF1825CCe001";
	  str >> CXSC_Sqrt3d2[37];
	  str = "+10000000000000e000";
	  str >> CXSC_Sqrt3d2[38];
	  str = "+10000000000001e000";
	  str >> CXSC_Sqrt3d2[39];

	  CXSC_Sqrt3d2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt3d2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // sqrt(3)/2

// **********************************************************************

 static real CXSC_Sqrt3r[40]; // CXSC_Sqrt3r[0], ... ,CXSC_Sqrt3r[39]
 static bool CXSC_Sqrt3r_initialized = false;

 lx_interval Sqrt3r_lx_interval() throw()
// Inclusion of 1/sqrt(3), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Sqrt3r_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1279A74590331Ce7FC";
	  str >> CXSC_Sqrt3r[0];
	  str = "+134863E0792BEDe7C6";
	  str >> CXSC_Sqrt3r[1];
	  str = "-1A82F9E6C53222e790";
	  str >> CXSC_Sqrt3r[2];
	  str = "-1CB0F41134253Ae75A";
	  str >> CXSC_Sqrt3r[3];
	  str = "+1859ED919EC30Be724";
	  str >> CXSC_Sqrt3r[4];
	  str = "+1454874FB1F3F4e6ED";
	  str >> CXSC_Sqrt3r[5];
	  str = "-1DE69C6D3D2741e6B7";
	  str >> CXSC_Sqrt3r[6];
	  str = "+17EEC450C48BE1e681";
	  str >> CXSC_Sqrt3r[7];
	  str = "-16F743EEE65D53e64B";
	  str >> CXSC_Sqrt3r[8];
	  str = "-1887B505D7E7C2e613";
	  str >> CXSC_Sqrt3r[9];
	  str = "-1484D2E10C1161e5DC";
	  str >> CXSC_Sqrt3r[10];
	  str = "-1A0B1F86177FB7e5A6";
	  str >> CXSC_Sqrt3r[11];
	  str = "+1FE389D3F2C54Ee56E";
	  str >> CXSC_Sqrt3r[12];
	  str = "+1F29F77C671544e538";
	  str >> CXSC_Sqrt3r[13];
	  str = "-16CE74ED77D9BEe502";
	  str >> CXSC_Sqrt3r[14];
	  str = "-1E38708FF0CCB5e4CC";
	  str >> CXSC_Sqrt3r[15];
	  str = "-1F13BCC70157D1e496";
	  str >> CXSC_Sqrt3r[16];
	  str = "+17EC34CF9B1930e460";
	  str >> CXSC_Sqrt3r[17];
	  str = "-117A638EFF3A8Be429";
	  str >> CXSC_Sqrt3r[18];
	  str = "-16A8EF69C312C5e3F3";
	  str >> CXSC_Sqrt3r[19];
	  str = "-1835C4B4FD2883e3BC";
	  str >> CXSC_Sqrt3r[20];
	  str = "+178E66E7009A44e386";
	  str >> CXSC_Sqrt3r[21];
	  str = "-136FD3299CE38Ae350";
	  str >> CXSC_Sqrt3r[22];
	  str = "+10CE607E925CC3e31A";
	  str >> CXSC_Sqrt3r[23];
	  str = "-1AC2B71232EC37e2E4";
	  str >> CXSC_Sqrt3r[24];
	  str = "-1ACAD8486879A7e2AD";
	  str >> CXSC_Sqrt3r[25];
	  str = "+1B7C6155125817e277";
	  str >> CXSC_Sqrt3r[26];
	  str = "-1CE6286338DB54e240";
	  str >> CXSC_Sqrt3r[27];
	  str = "-13819B52815439e20A";
	  str >> CXSC_Sqrt3r[28];
	  str = "-159653496BB604e1D3";
	  str >> CXSC_Sqrt3r[29];
	  str = "+1E3CEA8BB81EF5e19C";
	  str >> CXSC_Sqrt3r[30];
	  str = "+1C0A89B6922280e165";
	  str >> CXSC_Sqrt3r[31];
	  str = "+1A7F02B9E11970e12F";
	  str >> CXSC_Sqrt3r[32];
	  str = "-197D70AA62CA0Be0F9";
	  str >> CXSC_Sqrt3r[33];
	  str = "-18FF1CC85B90E7e0C3";
	  str >> CXSC_Sqrt3r[34];
	  str = "+19C2DC378988CFe08B";
	  str >> CXSC_Sqrt3r[35];
	  str = "-1A1BB3467EF366e055";
	  str >> CXSC_Sqrt3r[36];
	  str = "+1166CBADB2F273e01F";
	  str >> CXSC_Sqrt3r[37];
	  str = "+1000000B453C22e000";
	  str >> CXSC_Sqrt3r[38];
	  str = "+1000000B453C23e000";
	  str >> CXSC_Sqrt3r[39];

	  CXSC_Sqrt3r_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Sqrt3r[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1022,y);
 } // 1/sqrt(3)

// **********************************************************************

 static real CXSC_Er[40]; // CXSC_Er[0], ... ,CXSC_Er[39]
 static bool CXSC_Er_initialized = false;

 lx_interval Er_lx_interval() throw()
// Inclusion of 1/sqrt(3), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Er_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+178B56362CEF38e7FC";
	  str >> CXSC_Er[0];
	  str = "-1CA8A4270FADF5e7C5";
	  str >> CXSC_Er[1];
	  str = "-1837912B3FD2AAe78F";
	  str >> CXSC_Er[2];
	  str = "-152711999FB68Ce759";
	  str >> CXSC_Er[3];
	  str = "-17AD7C1289274Ee723";
	  str >> CXSC_Er[4];
	  str = "+17E8E56842B705e6E5";
	  str >> CXSC_Er[5];
	  str = "-1D24CB13796C2De6AF";
	  str >> CXSC_Er[6];
	  str = "-1456AABDA5C8F2e678";
	  str >> CXSC_Er[7];
	  str = "+1229F03C6276DDe642";
	  str >> CXSC_Er[8];
	  str = "-1569CFC4F53109e60C";
	  str >> CXSC_Er[9];
	  str = "-155B63C9B68091e5D4";
	  str >> CXSC_Er[10];
	  str = "+1580CF14DC087Ce59E";
	  str >> CXSC_Er[11];
	  str = "+1F9FF222313669e567";
	  str >> CXSC_Er[12];
	  str = "+15BC9CB1A22487e531";
	  str >> CXSC_Er[13];
	  str = "-1857E415C89B13e4FA";
	  str >> CXSC_Er[14];
	  str = "+13DF75706E3643e4C4";
	  str >> CXSC_Er[15];
	  str = "+13BDF5B7646234e48C";
	  str >> CXSC_Er[16];
	  str = "+1C956A5A3BE55De456";
	  str >> CXSC_Er[17];
	  str = "-167243FE9CD95Ee41F";
	  str >> CXSC_Er[18];
	  str = "+1798666D9D76F9e3E9";
	  str >> CXSC_Er[19];
	  str = "-195BC96299ED95e3B3";
	  str >> CXSC_Er[20];
	  str = "-1962287D82F280e37D";
	  str >> CXSC_Er[21];
	  str = "+1C3CF6DDC027D8e347";
	  str >> CXSC_Er[22];
	  str = "-182A3C09F5C0B7e310";
	  str >> CXSC_Er[23];
	  str = "+181C26FE7F6AB1e2DA";
	  str >> CXSC_Er[24];
	  str = "+19F6D7E4825294e2A4";
	  str >> CXSC_Er[25];
	  str = "+1BBC423BEA892Fe26D";
	  str >> CXSC_Er[26];
	  str = "+1342C7A3A14AB4e237";
	  str >> CXSC_Er[27];
	  str = "+12A70DFB042173e201";
	  str >> CXSC_Er[28];
	  str = "-10325653502352e1CB";
	  str >> CXSC_Er[29];
	  str = "-10AD4492DE41FFe191";
	  str >> CXSC_Er[30];
	  str = "-1E4529AB93CDA1e156";
	  str >> CXSC_Er[31];
	  str = "+1128833F39DF0Ae11E";
	  str >> CXSC_Er[32];
	  str = "-1E7EDF8F9B8A50e0E8";
	  str >> CXSC_Er[33];
	  str = "+1A42CBDB5BB8D0e0B0";
	  str >> CXSC_Er[34];
	  str = "+1973F3BD8250A1e07A";
	  str >> CXSC_Er[35];
	  str = "+116AF9EF0E6C71e040";
	  str >> CXSC_Er[36];
	  str = "-1786993285AA7Ae00A";
	  str >> CXSC_Er[37];
	  str = "-1000000000007De000";
	  str >> CXSC_Er[38];
	  str = "-1000000000007Ce000";
	  str >> CXSC_Er[39];

	  CXSC_Er_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Er[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1023,y);
 } // 1/e

// **********************************************************************

 static real CXSC_Ep2[40]; // CXSC_Ep2[0], ... ,CXSC_Ep2[39]
 static bool CXSC_Ep2_initialized = false;

 lx_interval Ep2_lx_interval() throw()
// Inclusion of e^2, Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ep2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1D8E64B8D4DDAEe7FC";
	  str >> CXSC_Ep2[0];
	  str = "-19E62E22EFCA4Ce7C5";
	  str >> CXSC_Ep2[1];
	  str = "+1577508F5CF5EDe78F";
	  str >> CXSC_Ep2[2];
	  str = "-186EF0294C2511e759";
	  str >> CXSC_Ep2[3];
	  str = "+177D109F148782e722";
	  str >> CXSC_Ep2[4];
	  str = "+166BBC354AB700e6EB";
	  str >> CXSC_Ep2[5];
	  str = "-1273AEC0115969e6B5";
	  str >> CXSC_Ep2[6];
	  str = "-1C5AE00D3BEEF1e67F";
	  str >> CXSC_Ep2[7];
	  str = "+15ACA3FDC9595Fe647";
	  str >> CXSC_Ep2[8];
	  str = "-113FCDFE2B1F0Ce610";
	  str >> CXSC_Ep2[9];
	  str = "+10EEDFD1AE90C9e5DA";
	  str >> CXSC_Ep2[10];
	  str = "+1D2CB8EDC7078Be5A4";
	  str >> CXSC_Ep2[11];
	  str = "+11827A19F175F8e56E";
	  str >> CXSC_Ep2[12];
	  str = "-10267512A9BFB2e537";
	  str >> CXSC_Ep2[13];
	  str = "-19A1E2FC413AE3e500";
	  str >> CXSC_Ep2[14];
	  str = "+1170C7A5981ADBe4CA";
	  str >> CXSC_Ep2[15];
	  str = "-1FC991480067CFe494";
	  str >> CXSC_Ep2[16];
	  str = "-12E9A54CF5CFB5e45D";
	  str >> CXSC_Ep2[17];
	  str = "-166FA6C468910Ae425";
	  str >> CXSC_Ep2[18];
	  str = "+10FA9B7050AF8De3EE";
	  str >> CXSC_Ep2[19];
	  str = "+198127CED41761e3B7";
	  str >> CXSC_Ep2[20];
	  str = "+107FD1EB487B65e380";
	  str >> CXSC_Ep2[21];
	  str = "+1B63EE064187DBe348";
	  str >> CXSC_Ep2[22];
	  str = "+13C943324AF1B5e311";
	  str >> CXSC_Ep2[23];
	  str = "+16AAE6F376094Ee2DA";
	  str >> CXSC_Ep2[24];
	  str = "+15DBB3D45B5A29e2A4";
	  str >> CXSC_Ep2[25];
	  str = "-181BC5BF587296e26E";
	  str >> CXSC_Ep2[26];
	  str = "-1819FC0B42A502e235";
	  str >> CXSC_Ep2[27];
	  str = "-1E06AE15A2D879e1FF";
	  str >> CXSC_Ep2[28];
	  str = "+171395ABE3E6CEe1C8";
	  str >> CXSC_Ep2[29];
	  str = "-1B76514AE69513e192";
	  str >> CXSC_Ep2[30];
	  str = "-1707F6C56433B7e15C";
	  str >> CXSC_Ep2[31];
	  str = "+108C1FADE66FE9e126";
	  str >> CXSC_Ep2[32];
	  str = "+1FB253285CC9E2e0F0";
	  str >> CXSC_Ep2[33];
	  str = "-16B3E49A6C1691e0B9";
	  str >> CXSC_Ep2[34];
	  str = "-12B135E875C44Ae080";
	  str >> CXSC_Ep2[35];
	  str = "+1385B11510C48Ce04A";
	  str >> CXSC_Ep2[36];
	  str = "-184E63EB2E35F1e010";
	  str >> CXSC_Ep2[37];
	  str = "+1000000000063Ee000";
	  str >> CXSC_Ep2[38];
	  str = "+1000000000064Ee000";
	  str >> CXSC_Ep2[39];

	  CXSC_Ep2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ep2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1019,y);
 } // e^2

// **********************************************************************

 static real CXSC_Ep2r[40]; // CXSC_Ep2r[0], ... ,CXSC_Ep2r[39]
 static bool CXSC_Ep2r_initialized = false;

 lx_interval Ep2r_lx_interval() throw()
// Inclusion of 1/e^2, Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ep2r_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+1152AAA3BF81CCe7FD";
	  str >> CXSC_Ep2r[0];
	  str = "-1809224547B4BFe7C7";
	  str >> CXSC_Ep2r[1];
	  str = "-16A8E079134F13e791";
	  str >> CXSC_Ep2r[2];
	  str = "+14564CACF0994Ee759";
	  str >> CXSC_Ep2r[3];
	  str = "+1B796438129AF8e723";
	  str >> CXSC_Ep2r[4];
	  str = "-1ACFED57EF2AE5e6ED";
	  str >> CXSC_Ep2r[5];
	  str = "-1A968CBDBB5D9De6B6";
	  str >> CXSC_Ep2r[6];
	  str = "+1A7238CBD97B71e67D";
	  str >> CXSC_Ep2r[7];
	  str = "-146C53DB77BB01e646";
	  str >> CXSC_Ep2r[8];
	  str = "-1EEC161C3EBBD7e60D";
	  str >> CXSC_Ep2r[9];
	  str = "-12D084DC157ACEe5D6";
	  str >> CXSC_Ep2r[10];
	  str = "+12A61F46883347e5A0";
	  str >> CXSC_Ep2r[11];
	  str = "+1993BAF10CAE0Be565";
	  str >> CXSC_Ep2r[12];
	  str = "+1F9224351178FFe52F";
	  str >> CXSC_Ep2r[13];
	  str = "-1C366D1C7BA64Ae4F8";
	  str >> CXSC_Ep2r[14];
	  str = "-17D9938EFA4657e4C1";
	  str >> CXSC_Ep2r[15];
	  str = "+1B6668DF0C1286e48B";
	  str >> CXSC_Ep2r[16];
	  str = "+1F7A4FFC9B48C6e451";
	  str >> CXSC_Ep2r[17];
	  str = "+1F3E3AF6F17591e41B";
	  str >> CXSC_Ep2r[18];
	  str = "+1B1E0C69EE87BEe3E4";
	  str >> CXSC_Ep2r[19];
	  str = "+12DF95CABED5A7e3AD";
	  str >> CXSC_Ep2r[20];
	  str = "+1665EF8CBE7C05e377";
	  str >> CXSC_Ep2r[21];
	  str = "-1B3C417ABEAD6Be340";
	  str >> CXSC_Ep2r[22];
	  str = "-1A19C35B2B0C58e30A";
	  str >> CXSC_Ep2r[23];
	  str = "+18607193ADB301e2D2";
	  str >> CXSC_Ep2r[24];
	  str = "-1BB2016A08F428e298";
	  str >> CXSC_Ep2r[25];
	  str = "+1154B5B7FEDB3Ee262";
	  str >> CXSC_Ep2r[26];
	  str = "-13AE5C2DBEA451e22C";
	  str >> CXSC_Ep2r[27];
	  str = "-10F60BC60CDBFCe1F6";
	  str >> CXSC_Ep2r[28];
	  str = "+1F540F667B3746e1BE";
	  str >> CXSC_Ep2r[29];
	  str = "-1B205D40167EC7e187";
	  str >> CXSC_Ep2r[30];
	  str = "+1C8A0A08DE85F9e151";
	  str >> CXSC_Ep2r[31];
	  str = "+1856ED169F0183e11B";
	  str >> CXSC_Ep2r[32];
	  str = "-147D787794462Ce0E4";
	  str >> CXSC_Ep2r[33];
	  str = "+1516AB11B003F3e0AB";
	  str >> CXSC_Ep2r[34];
	  str = "-1D652A20732EA1e075";
	  str >> CXSC_Ep2r[35];
	  str = "-12EB9F673FF3EDe03F";
	  str >> CXSC_Ep2r[36];
	  str = "+1A5FE9239C5237e008";
	  str >> CXSC_Ep2r[37];
	  str = "-1000000000002Ce000";
	  str >> CXSC_Ep2r[38];
	  str = "-1000000000002Be000";
	  str >> CXSC_Ep2r[39];

	  CXSC_Ep2r_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ep2r[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1025,y);
 } // 1/e^2

// **********************************************************************

 static real CXSC_Ep2Pi[40]; // CXSC_Ep2Pi[0], ... ,CXSC_Ep2Pi[39]
 static bool CXSC_Ep2Pi_initialized = false;

 lx_interval Ep2Pi_lx_interval() throw()
// Inclusion of e^(2pi), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_Ep2Pi_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+10BBEEE9177E19e7FD";
	  str >> CXSC_Ep2Pi[0];
	  str = "+1C7DD9272526B1e7C5";
	  str >> CXSC_Ep2Pi[1];
	  str = "+15200F57AB89EDe78F";
	  str >> CXSC_Ep2Pi[2];
	  str = "-1FCCB6EDBE9C36e758";
	  str >> CXSC_Ep2Pi[3];
	  str = "+1BEA0BF179A589e722";
	  str >> CXSC_Ep2Pi[4];
	  str = "-1F3AD5A6B77F9Ee6EC";
	  str >> CXSC_Ep2Pi[5];
	  str = "-1622F702B57637e6B5";
	  str >> CXSC_Ep2Pi[6];
	  str = "-100C09AE818734e67C";
	  str >> CXSC_Ep2Pi[7];
	  str = "+10DA7ADA79EFE6e642";
	  str >> CXSC_Ep2Pi[8];
	  str = "+1FF9BF48B72959e60B";
	  str >> CXSC_Ep2Pi[9];
	  str = "-17AD7A3F6D2A14e5D5";
	  str >> CXSC_Ep2Pi[10];
	  str = "+1FCD4B0FA971E4e59E";
	  str >> CXSC_Ep2Pi[11];
	  str = "+193A2CDC04526Be567";
	  str >> CXSC_Ep2Pi[12];
	  str = "-18CBE5FDFAF25Fe531";
	  str >> CXSC_Ep2Pi[13];
	  str = "+1D47EEE171DA93e4FA";
	  str >> CXSC_Ep2Pi[14];
	  str = "-15B0F8DA29DB32e4C4";
	  str >> CXSC_Ep2Pi[15];
	  str = "-19207AD7E637D8e48C";
	  str >> CXSC_Ep2Pi[16];
	  str = "+191CA743F265A6e456";
	  str >> CXSC_Ep2Pi[17];
	  str = "+1A15069182EF28e41F";
	  str >> CXSC_Ep2Pi[18];
	  str = "-1D58BF80B501B6e3E9";
	  str >> CXSC_Ep2Pi[19];
	  str = "+1435920A849065e3B3";
	  str >> CXSC_Ep2Pi[20];
	  str = "-11931903C826FBe37C";
	  str >> CXSC_Ep2Pi[21];
	  str = "+169B0688CF564Ee346";
	  str >> CXSC_Ep2Pi[22];
	  str = "-12539A43ACDD10e309";
	  str >> CXSC_Ep2Pi[23];
	  str = "+172B8963B0CE58e2D3";
	  str >> CXSC_Ep2Pi[24];
	  str = "-13E6A7B1E3A306e29D";
	  str >> CXSC_Ep2Pi[25];
	  str = "-17F20768EDB9E7e267";
	  str >> CXSC_Ep2Pi[26];
	  str = "+130F006E28050Fe22F";
	  str >> CXSC_Ep2Pi[27];
	  str = "+149C245E1C5FEFe1F9";
	  str >> CXSC_Ep2Pi[28];
	  str = "-102CDEE5CA2F95e1C2";
	  str >> CXSC_Ep2Pi[29];
	  str = "+1040AABBBB0BFBe18B";
	  str >> CXSC_Ep2Pi[30];
	  str = "+18D7DB731409F2e154";
	  str >> CXSC_Ep2Pi[31];
	  str = "-1868ADF8479A20e11A";
	  str >> CXSC_Ep2Pi[32];
	  str = "+1BCB4CE8F6AF6Ae0E2";
	  str >> CXSC_Ep2Pi[33];
	  str = "-1A6BA8B081A793e0AC";
	  str >> CXSC_Ep2Pi[34];
	  str = "-1DE2841143A816e075";
	  str >> CXSC_Ep2Pi[35];
	  str = "+1CB5B248339C0Ee03F";
	  str >> CXSC_Ep2Pi[36];
	  str = "-1B1B84E7980944e007";
	  str >> CXSC_Ep2Pi[37];
	  str = "-10000000000003e000";
	  str >> CXSC_Ep2Pi[38];
	  str = "-10000000000002e000";
	  str >> CXSC_Ep2Pi[39];

	  CXSC_Ep2Pi_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_Ep2Pi[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1013,y);
 } // e^(2pi)

// **********************************************************************

 static real CXSC_EpPid2[40]; // CXSC_EpPid2[0], ... ,CXSC_EpPid2[39]
 static bool CXSC_EpPid2_initialized = false;

 lx_interval EpPid2_lx_interval() throw()
// Inclusion of e^(pi/2), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_EpPid2_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+133DEDC855935Fe7FD";
	  str >> CXSC_EpPid2[0];
	  str = "+13E45A768FB73Ce7C7";
	  str >> CXSC_EpPid2[1];
	  str = "-1FB31CF300FF3Ce791";
	  str >> CXSC_EpPid2[2];
	  str = "-1E80D8BEB83F79e75B";
	  str >> CXSC_EpPid2[3];
	  str = "-14A3DE039142DDe722";
	  str >> CXSC_EpPid2[4];
	  str = "-18792D7A37282Be6E7";
	  str >> CXSC_EpPid2[5];
	  str = "-19DF43A5980C28e6B1";
	  str >> CXSC_EpPid2[6];
	  str = "-1C6F0F641C0D67e67B";
	  str >> CXSC_EpPid2[7];
	  str = "-1779C86C2DB5ACe645";
	  str >> CXSC_EpPid2[8];
	  str = "+168521EE91B16Fe60D";
	  str >> CXSC_EpPid2[9];
	  str = "+12530F905D97BDe5D7";
	  str >> CXSC_EpPid2[10];
	  str = "+13498112CB7585e5A1";
	  str >> CXSC_EpPid2[11];
	  str = "+1BA4546B13A434e56B";
	  str >> CXSC_EpPid2[12];
	  str = "+14FF791C56421Ce534";
	  str >> CXSC_EpPid2[13];
	  str = "-1F375C223A2152e4FE";
	  str >> CXSC_EpPid2[14];
	  str = "-126AB0C8C77412e4C8";
	  str >> CXSC_EpPid2[15];
	  str = "-1B39C9C0B8C54Ae490";
	  str >> CXSC_EpPid2[16];
	  str = "-167741414E31E3e459";
	  str >> CXSC_EpPid2[17];
	  str = "+1DEFB4462546C1e421";
	  str >> CXSC_EpPid2[18];
	  str = "-10F7B89CC30514e3E9";
	  str >> CXSC_EpPid2[19];
	  str = "+1E87D3145A3CEEe3B3";
	  str >> CXSC_EpPid2[20];
	  str = "+18AA09D5CD3B7Be37D";
	  str >> CXSC_EpPid2[21];
	  str = "+1E738C390E548Be347";
	  str >> CXSC_EpPid2[22];
	  str = "+147542CC36F28Be30E";
	  str >> CXSC_EpPid2[23];
	  str = "+1B217FFE679632e2D8";
	  str >> CXSC_EpPid2[24];
	  str = "+1A8F3962771086e2A0";
	  str >> CXSC_EpPid2[25];
	  str = "-187231F1E3EFC2e26A";
	  str >> CXSC_EpPid2[26];
	  str = "-15010B009CF001e233";
	  str >> CXSC_EpPid2[27];
	  str = "-1F22E68271119Fe1FB";
	  str >> CXSC_EpPid2[28];
	  str = "+11CA8D2164A3BAe1C5";
	  str >> CXSC_EpPid2[29];
	  str = "+1C20B237A324D7e18F";
	  str >> CXSC_EpPid2[30];
	  str = "-18C70E40461930e157";
	  str >> CXSC_EpPid2[31];
	  str = "+1025F32E109A37e120";
	  str >> CXSC_EpPid2[32];
	  str = "-12087D5EA8F469e0EA";
	  str >> CXSC_EpPid2[33];
	  str = "-14E1EE796B734Ae0B4";
	  str >> CXSC_EpPid2[34];
	  str = "-176EBB3BB1E41Ce07E";
	  str >> CXSC_EpPid2[35];
	  str = "+1374F617B0FF49e048";
	  str >> CXSC_EpPid2[36];
	  str = "+1D28C408575ECEe011";
	  str >> CXSC_EpPid2[37];
	  str = "-10000000006878e000";
	  str >> CXSC_EpPid2[38];
	  str = "-10000000006877e000";
	  str >> CXSC_EpPid2[39];

	  CXSC_EpPid2_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_EpPid2[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1020,y);
 } // e^(pi/2)

// **********************************************************************

 static real CXSC_EpPid4[40]; // CXSC_EpPid4[0], ... ,CXSC_EpPid4[39]
 static bool CXSC_EpPid4_initialized = false;

 lx_interval EpPid4_lx_interval() throw()
// Inclusion of e^(pi/4), Blomquist, 27.11.2008;
 {
	 l_interval y;
	 int stagsave = stagprec,
  stagmax = 39;
  if (!CXSC_EpPid4_initialized)
  {
	  std::string str;
	  std::cout << SaveOpt;
	  std::cout << Hex;
	  str = "+118BD669471CAAe7FD";
	  str >> CXSC_EpPid4[0];
	  str = "+1F0ED609715756e7C7";
	  str >> CXSC_EpPid4[1];
	  str = "-1B9C7B871FE1DBe791";
	  str >> CXSC_EpPid4[2];
	  str = "+15C0FECE98F209e75A";
	  str >> CXSC_EpPid4[3];
	  str = "+18C9FACC5DF3CEe724";
	  str >> CXSC_EpPid4[4];
	  str = "+15EDE838B4A399e6EC";
	  str >> CXSC_EpPid4[5];
	  str = "-1C7EFACA363051e6B6";
	  str >> CXSC_EpPid4[6];
	  str = "-1A1EBEA1646411e680";
	  str >> CXSC_EpPid4[7];
	  str = "+1AEF54E68CE03Be649";
	  str >> CXSC_EpPid4[8];
	  str = "-11250CB97FDDBFe60F";
	  str >> CXSC_EpPid4[9];
	  str = "-169ADC0E65B8A7e5D8";
	  str >> CXSC_EpPid4[10];
	  str = "+198A501DB90EDDe5A2";
	  str >> CXSC_EpPid4[11];
	  str = "-1586909A3F6365e56B";
	  str >> CXSC_EpPid4[12];
	  str = "+1BE542410F8CE7e535";
	  str >> CXSC_EpPid4[13];
	  str = "+1E7EEC51889EECe4FF";
	  str >> CXSC_EpPid4[14];
	  str = "-1913C9FC19333Ce4C9";
	  str >> CXSC_EpPid4[15];
	  str = "+1112C71EA1E6F0e492";
	  str >> CXSC_EpPid4[16];
	  str = "-1C4CCF0F5D1E14e45B";
	  str >> CXSC_EpPid4[17];
	  str = "+1AC4A72310FA27e425";
	  str >> CXSC_EpPid4[18];
	  str = "-13EC6A07AD7C15e3EE";
	  str >> CXSC_EpPid4[19];
	  str = "+1114CC16D255A3e3B6";
	  str >> CXSC_EpPid4[20];
	  str = "+17FA54DD584C6Ee380";
	  str >> CXSC_EpPid4[21];
	  str = "-1BBFE0D94FA881e34A";
	  str >> CXSC_EpPid4[22];
	  str = "-1AACE950D75AB5e314";
	  str >> CXSC_EpPid4[23];
	  str = "-1C197FC5941652e2DD";
	  str >> CXSC_EpPid4[24];
	  str = "-1FC33183282E6Ee2A7";
	  str >> CXSC_EpPid4[25];
	  str = "+15ADA6B92D282Ee26A";
	  str >> CXSC_EpPid4[26];
	  str = "-10CC6C2A9B1995e233";
	  str >> CXSC_EpPid4[27];
	  str = "+1767402CD3F07Be1FC";
	  str >> CXSC_EpPid4[28];
	  str = "-134D3B2AED1AACe1C5";
	  str >> CXSC_EpPid4[29];
	  str = "+1C87E322B76BC8e18D";
	  str >> CXSC_EpPid4[30];
	  str = "+1CBF921AE01812e157";
	  str >> CXSC_EpPid4[31];
	  str = "+16B04C1CCCDEAEe11D";
	  str >> CXSC_EpPid4[32];
	  str = "+11D2A5DAE175A6e0E5";
	  str >> CXSC_EpPid4[33];
	  str = "-1D6F0482D56D0Ee0AF";
	  str >> CXSC_EpPid4[34];
	  str = "-1E1B89C6DE0660e079";
	  str >> CXSC_EpPid4[35];
	  str = "-1202BA3792C129e042";
	  str >> CXSC_EpPid4[36];
	  str = "+1B23FE2BAFDF56e00A";
	  str >> CXSC_EpPid4[37];
	  str = "-100000000000DEe000";
	  str >> CXSC_EpPid4[38];
	  str = "-100000000000DDe000";
	  str >> CXSC_EpPid4[39];

	  CXSC_EpPid4_initialized = true;
	  std::cout << RestoreOpt;
  } 
  stagprec = stagmax;
  y = adjust(l_interval(0));

  for (int i=0; i<=stagmax; i++)
	  y[i+1] = CXSC_EpPid4[i];

  stagprec = stagsave;
  y = adjust(y);

  return lx_interval(-1021,y);
 } // e^(pi/4)


// ------------------------------------------------------------------------------
// ---------------- lx_real constants in high accuracy --------------------
// ------------------------------------------------------------------------------

 lx_real Pi_lx_real() throw()
 { return mid(Pi_lx_interval()); }
 lx_real Pip2_lx_real() throw()
 { return mid(Pip2_lx_interval()); }
 lx_real Pi2r_lx_real() throw()
 { return mid(Pi2r_lx_interval()); }
 lx_real Pi2_lx_real() throw()
 { return mid(Pi2_lx_interval()); }
 lx_real Pid4_lx_real() throw()
 { return mid(Pid4_lx_interval()); }
 lx_real Pid2_lx_real() throw()
 { return mid(Pid2_lx_interval()); }
 lx_real Ln2_lx_real() throw()
 { return mid(Ln2_lx_interval()); }
 lx_real Ln10_lx_real() throw()
 { return mid(Ln10_lx_interval()); }
 lx_real Ln10r_lx_real() throw()
 { return mid(Ln10r_lx_interval()); }
 lx_real Pir_lx_real() throw()
 { return mid(Pir_lx_interval()); }
 lx_real SqrtPi_lx_real() throw()
 { return mid(SqrtPi_lx_interval()); }
 lx_real Sqrt2Pi_lx_real() throw()
 { return mid(Sqrt2Pi_lx_interval()); }
 lx_real Sqrt2_lx_real() throw()
 { return mid(Sqrt2_lx_interval()); }
 lx_real Sqrt2r_lx_real() throw()
 { return mid(Sqrt2r_lx_interval()); }
 lx_real Sqrt3_lx_real() throw()
 { return mid(Sqrt3_lx_interval()); }
 lx_real Sqrt3r_lx_real() throw()
 { return mid(Sqrt3r_lx_interval()); }
 lx_real Sqrt3d2_lx_real() throw()
 { return mid(Sqrt3d2_lx_interval()); }
 lx_real Ln2r_lx_real() throw()
 { return mid(Ln2r_lx_interval()); }
 lx_real Pid3_lx_real() throw()
 { return mid(Pid3_lx_interval()); }
 lx_real SqrtPir_lx_real() throw()
 { return mid(SqrtPir_lx_interval()); }
 lx_real Sqrt2Pir_lx_real() throw()
 { return mid(Sqrt2Pir_lx_interval()); }
 lx_real LnPi_lx_real() throw()
 { return mid(LnPi_lx_interval()); }
 lx_real Ln2Pi_lx_real() throw()
 { return mid(Ln2Pi_lx_interval()); }
 lx_real E_lx_real() throw()
 { return mid(E_lx_interval()); }
 lx_real Ep2r_lx_real() throw()
 { return mid(Ep2r_lx_interval()); }
 lx_real Ep2_lx_real() throw()
 { return mid(Ep2_lx_interval()); }
 lx_real Er_lx_real() throw()
 { return mid(Er_lx_interval()); }
 lx_real EpPi_lx_real() throw()
 { return mid(EpPi_lx_interval()); }
 lx_real EpPid2_lx_real() throw()
 { return mid(EpPid2_lx_interval()); }
 lx_real EpPid4_lx_real() throw()
 { return mid(EpPid4_lx_interval()); }
 lx_real Ep2Pi_lx_real() throw()
 { return mid(Ep2Pi_lx_interval()); }
 lx_real EulerGamma_lx_real() throw()
 { return mid(EulerGamma_lx_interval()); }
 lx_real Catalan_lx_real() throw()
 { return mid(Catalan_lx_interval()); }
 lx_real sqrt5_lx_real() throw()
 { return mid(sqrt5_lx_interval()); }
 lx_real sqrt7_lx_real() throw()
 { return mid(sqrt7_lx_interval()); }
 lx_real One_m_lx_real() throw()
 { return mid(One_m_lx_interval()); }
 lx_real One_p_lx_real() throw()
 { return mid(One_p_lx_interval()); }

// --------------------------------------------------------------------------
// -------Elementary functions related to type lx_real ----------------------
// --------------------------------------------------------------------------

 lx_real sqrt(const lx_real& x) throw()
 { return mid(sqrt(lx_interval(x))); }
 lx_real sqr(const lx_real& x) throw()
 { return mid(sqr(lx_interval(x))); }
 lx_real ln(const lx_real& x) throw()
 { return mid(ln(lx_interval(x))); }
 lx_real log2(const lx_real& x) throw()
 { return mid(log2(lx_interval(x))); }
 lx_real log10(const lx_real& x) throw()
 { return mid(log10(lx_interval(x))); }
 lx_real lnp1(const lx_real& x) throw()
 { return mid(lnp1(lx_interval(x))); }
 lx_real exp(const lx_real& x) throw()
 { return mid(exp(lx_interval(x))); }
 lx_real exp2(const lx_real& x) throw()
 { return mid(exp2(lx_interval(x))); }
 lx_real exp10(const lx_real& x) throw()
 { return mid(exp10(lx_interval(x))); }
 lx_real expm1(const lx_real& x) throw()
 { return mid(expm1(lx_interval(x))); }
 lx_real power(const lx_real& x, const real& r) throw()
 { return mid(power(lx_interval(x),r)); }
 lx_real pow(const lx_real& x, const lx_real& y) throw()
 { return mid(pow(lx_interval(x),lx_interval(y))); }
 lx_real xp1_pow_y(const lx_real& x, const lx_real& y) throw()
 { return mid(xp1_pow_y(lx_interval(x),lx_interval(y))); }
 lx_real sin(const lx_real& x) throw()
 { return mid(sin(lx_interval(x))); }
 lx_real sin_n(const lx_real& x, const real& n) throw()
 { return mid(sin_n(lx_interval(x),n)); }
 lx_real cos(const lx_real& x) throw()
 { return mid(cos(lx_interval(x))); }
 lx_real cos_n(const lx_real& x, const real& n) throw()
 { return mid(cos_n(lx_interval(x),n)); }
 lx_real tan(const lx_real& x) throw()
 { return mid(tan(lx_interval(x))); }
 lx_real cot(const lx_real& x) throw()
 { return mid(cot(lx_interval(x))); }
 lx_real sqrt1px2(const lx_real& x) throw()
 { return mid(sqrt1px2(lx_interval(x))); }
 lx_real atan(const lx_real& x) throw()
 { return mid(atan(lx_interval(x))); }
 lx_real sqrt1mx2(const lx_real& x) throw()
 { return mid(sqrt1mx2(lx_interval(x))); }
 lx_real sqrtx2m1(const lx_real& x) throw()
 { return mid(sqrtx2m1(lx_interval(x))); }
 lx_real asin(const lx_real& x) throw()
 { return mid(asin(lx_interval(x))); }
 lx_real acos(const lx_real& x) throw()
 { return mid(acos(lx_interval(x))); }
 lx_real acot(const lx_real& x) throw()
 { return mid(acot(lx_interval(x))); }
 lx_real sinh(const lx_real& x) throw()
 { return mid(sinh(lx_interval(x))); }
 lx_real cosh(const lx_real& x) throw()
 { return mid(cosh(lx_interval(x))); }
 lx_real tanh(const lx_real& x) throw()
 { return mid(tanh(lx_interval(x))); }
 lx_real coth(const lx_real& x) throw()
 { return mid(coth(lx_interval(x))); }
 lx_real sqrtp1m1(const lx_real& x) throw()
 { return mid(sqrtp1m1(lx_interval(x))); }
 lx_real asinh(const lx_real& x) throw()
 { return mid(asinh(lx_interval(x))); }
 lx_real acosh(const lx_real& x) throw()
 { return mid(acosh(lx_interval(x))); }
 lx_real acoshp1(const lx_real& x) throw()
 { return mid(acoshp1(lx_interval(x))); }
 lx_real atanh(const lx_real& x) throw()
 { return mid(atanh(lx_interval(x))); }
 lx_real atanh1m(const lx_real& x) throw()
 { return mid(atanh1m(lx_interval(x))); }
 lx_real atanhm1p(const lx_real& x) throw()
 { return mid(atanhm1p(lx_interval(x))); }
 lx_real acoth(const lx_real& x) throw()
 { return mid(acoth(lx_interval(x))); }
 lx_real acothp1(const lx_real& x) throw()
 { return mid(acothp1(lx_interval(x))); }
 lx_real acothm1m(const lx_real& x) throw()
 { return mid(acothm1m(lx_interval(x))); }
 lx_real sqrtx2y2(const lx_real& x, const lx_real& y) throw()
 { return mid(sqrtx2y2(lx_interval(x),lx_interval(y))); }
 lx_real ln_sqrtx2y2(const lx_real& x, const lx_real& y) throw()
 { return mid(ln_sqrtx2y2(lx_interval(x),lx_interval(y))); }
 lx_real sqrt(const lx_real& x, int n) throw()
 { return mid(sqrt(lx_interval(x),n)); }

} // namespace cxsc
