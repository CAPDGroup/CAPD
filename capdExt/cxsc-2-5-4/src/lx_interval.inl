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

/* CVS $Id: lx_interval.inl,v 1.8 2014/01/30 17:23:47 cxsc Exp $ */

/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

namespace cxsc {
	

// --------------------------------------------------------------------------
// ------ Inline functions and operators related to type lx_interval --------
// --------------------------------------------------------------------------

inline int StagPrec(const lx_interval &a) throw() 
{ return StagPrec(a.li); }

inline real expo(const lx_interval &a) throw()
{ return (a.ex); }

inline l_interval li_part(const lx_interval &a) throw()
{ return (a.li); }

inline int Disjoint(const lx_interval& a, const lx_interval& b)
{ return (Inf(a)>Sup(b) || Inf(b)>Sup(a)); }

inline int in(const lx_interval& x, const lx_interval& y) 
{ return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) ); }

inline int in(const l_interval& x, const lx_interval& y) 
{ return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) ); }

inline int in(const interval& x, const lx_interval& y) 
{ return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) ); }

inline int in(const lx_real& x, const lx_interval& y ) 
{ return ( (Inf(y) <= x) && (x <= Sup(y)) ); }

inline int in(const l_real& x, const lx_interval& y ) 
{ return ( (Inf(y) <= x) && (x <= Sup(y)) ); }

inline int in(const real& x, const lx_interval& y ) 
{ return ( (Inf(y) <= x) && (x <= Sup(y)) ); }

inline lx_interval Blow(const lx_interval& a) throw()
// This function returns an inflated interval res,
// if the following condition is complied:
//      expo(a) < -Max_Int_R + 2100 = Max_Int_N;
// After the inflation it holds: 
//      expo(a) = -Max_Int_R + 2100;
// a = [0,0] is not modified.
// Blomquist, 31.01.2008;
{
    lx_interval res(a);
    l_interval al(a.li);

    if (a.ex<Max_Int_N)
    {
	if (Inf(al)>0)
	    SetInf(al,0.0);
	else
	    if (Sup(al)<0)
		SetSup(al,0.0);
	res = lx_interval(Max_Int_N,al);
    }

    return res;
}

inline lx_interval Blow( const lx_interval& x, const real& eps )
{
    l_interval y(li_part(x));
    return lx_interval( expo(x), Blow(y,eps) );
}

inline lx_real AbsMin (const lx_interval &x)
{
   if ( in(0.0,x) ) return lx_real(0.0);
   else 
   { 
       lx_real y(Inf(x));
       if (y > 0.0) return y;
       else return -Sup(x);
   }
}
  
inline lx_real AbsMax (const lx_interval &x)
{ 
    lx_real a, b;

    a = abs(Inf(x));
    b = abs(Sup(x));

    if (a > b)
	return a;
    else
	return b;
}

inline lx_interval & lx_interval::operator = (const lx_interval &a) throw()
{
    ex = a.ex;
    li = a.li;
    return *this;
} 

inline lx_interval & lx_interval::operator = (const lx_real &a) throw()
{
    ex = expo(a);
    li = lr_part(a);
    return *this;
} 

inline lx_interval & lx_interval::operator = (const l_interval &a) throw()
{
    ex = 0;
    li = a;
    return *this;
}

inline lx_interval & lx_interval::operator = (const l_real & a) throw()
{
    ex = 0;
    li = a;
    return *this;
}

inline lx_interval & lx_interval::operator = (const real & a) throw()
{
    ex = 0;
    li = a;
    return *this;
}

inline lx_interval & lx_interval::operator = (const interval & a) throw()
{
    ex = 0;
    li = a;
    return *this;
} 

inline lx_interval operator+(const lx_interval & a) throw() { return a; }
inline lx_interval operator-(const lx_interval & a) throw()
{  return lx_interval(a.ex,-a.li);  }

inline lx_interval operator + (const lx_interval &a, const l_interval &b) throw()
{  return a + lx_interval(b); }

inline lx_interval operator + (const l_interval &a, const lx_interval &b) throw()
{  return lx_interval(a) + b; }

inline lx_interval operator + (const lx_interval &a, const l_real &b) throw()
{  return a + lx_interval(b); }

inline lx_interval operator + (const l_real &a, const lx_interval &b) throw()
{  return lx_interval(a) + b; }

inline lx_interval operator + (const lx_interval &a, const lx_real &b) throw()
{  return a + lx_interval(b); }

inline lx_interval operator + (const lx_real &a, const lx_interval &b) throw()
{  return lx_interval(a) + b; }

inline lx_interval operator + (const lx_interval &a, const real &b) throw()
{  return a + lx_interval(b); }

inline lx_interval operator + (const real &a, const lx_interval &b) throw()
{  return lx_interval(a) + b; }

inline lx_interval operator + (const lx_interval &a, const interval &b) throw()
{  return a + lx_interval(b); }

inline lx_interval operator + (const interval &a, const lx_interval &b) throw()
{  return lx_interval(a) + b; }

inline lx_interval & operator +=(lx_interval &a, const lx_interval &b) throw()
{  return a = a+b; }

inline lx_interval & operator +=(lx_interval &a, const l_interval &b) throw()
{  return a = a+b; }

inline lx_interval & operator +=(lx_interval &a, const l_real &b) throw()
{  return a = a+b; }

inline lx_interval & operator +=(lx_interval &a, const lx_real &b) throw()
{  return a = a+b; }

inline lx_interval & operator +=(lx_interval &a, const real &b) throw()
{  return a = a+b; }

inline lx_interval & operator +=(lx_interval &a, const interval &b) throw()
{  return a = a+b; }

inline lx_interval operator - (const lx_interval &a, const lx_interval &b) throw()
{  return a + lx_interval(-b); }

inline lx_interval operator - (const lx_interval &a, const l_interval &b) throw()
{  return a + lx_interval(-b); }

inline lx_interval operator - (const l_interval &a, const lx_interval &b) throw()
{  return lx_interval(a) + (-b); }
 
inline lx_interval operator - (const lx_interval &a, const l_real &b) throw()
{  return a + lx_interval(-b); }

inline lx_interval operator - (const l_real &a, const lx_interval &b) throw()
{  return lx_interval(a) + (-b); }

inline lx_interval operator - (const lx_interval &a, const lx_real &b) throw()
{  return a + lx_interval(-b); }

inline lx_interval operator - (const lx_real &a, const lx_interval &b) throw()
{  return lx_interval(a) + (-b); }
 
inline lx_interval operator - (const lx_interval &a, const real &b) throw()
{  return a + lx_interval(-b); }

inline lx_interval operator - (const real &a, const lx_interval &b) throw()
{  return lx_interval(a) + (-b); }

inline lx_interval operator - (const lx_interval &a, const interval &b) throw()
{  return a + lx_interval(-b); }

inline lx_interval operator - (const interval &a, const lx_interval &b) throw()
{  return lx_interval(a) + (-b); }

inline lx_interval & operator -=(lx_interval &a, const lx_interval &b) throw()
{  return a = a-b; }
inline lx_interval & operator -=(lx_interval &a, const l_interval &b) throw()
{  return a = a-b; }
inline lx_interval & operator -=(lx_interval &a, const l_real &b) throw()
{  return a = a-b; }
inline lx_interval & operator -=(lx_interval &a, const lx_real &b) throw()
{  return a = a-b; }
inline lx_interval & operator -=(lx_interval &a, const real &b) throw()
{  return a = a-b; }
inline lx_interval & operator -=(lx_interval &a, const interval &b) throw()
{  return a = a-b; }

inline lx_interval operator * (const lx_interval &a, const l_interval &b) 
    throw()
{ return a * lx_interval(b); }

inline lx_interval operator * (const l_interval &a, const lx_interval &b) 
    throw()
{ return lx_interval(a) * b; }

inline lx_interval operator * (const lx_interval &a, const l_real &b) throw()
{ return a * lx_interval(b); }

inline lx_interval operator * (const l_real &a, const lx_interval &b) throw()
{ return lx_interval(a) * b; }

inline lx_interval operator * (const lx_interval &a, const lx_real &b) throw()
{ return a * lx_interval(b); }

inline lx_interval operator * (const lx_real &a, const lx_interval &b) throw()
{ return lx_interval(a) * b; }

inline lx_interval operator * (const lx_interval &a, const real &b) throw()
{ return a * lx_interval(b); }

inline lx_interval operator * (const real &a, const lx_interval &b) throw()
{ return lx_interval(a) * b; }

inline lx_interval operator * (const lx_interval &a, const interval &b) throw()
{ return a * lx_interval(b); }

inline lx_interval operator * (const interval &a, const lx_interval &b) throw()
{ return lx_interval(a) * b; }

inline lx_interval & operator *=(lx_interval &a, const lx_interval &b) throw()
{  return a = a*b; }
inline lx_interval & operator *=(lx_interval &a, const l_interval &b) throw()
{  return a = a*b; }
inline lx_interval & operator *=(lx_interval &a, const l_real &b) throw()
{  return a = a*b; }
inline lx_interval & operator *=(lx_interval &a, const lx_real &b) throw()
{  return a = a*b; }
inline lx_interval & operator *=(lx_interval &a, const real &b) throw()
{  return a = a*b; }
inline lx_interval & operator *=(lx_interval &a, const interval &b) throw()
{  return a = a*b; } 

inline lx_interval operator / (const lx_interval &a, const l_interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return a / lx_interval(b); }

inline lx_interval operator / (const l_interval &a, const lx_interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return lx_interval(a) / b; }

inline lx_interval operator / (const lx_interval &a, const l_real &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return a / lx_interval(b); }

inline lx_interval operator / (const l_real &a, const lx_interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return lx_interval(a) / b; }

inline lx_interval operator / (const lx_interval &a, const real &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return a / lx_interval(b); }

inline lx_interval operator / (const real &a, const lx_interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return lx_interval(a) / b; }

inline lx_interval operator / (const lx_interval &a, const interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return a / lx_interval(b); }

inline lx_interval operator / (const interval &a, const lx_interval &b) 
                              throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ return lx_interval(a) / b; }

inline lx_interval operator / (const lx_interval &a, const lx_real &b) 
                                      throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{  return a / lx_interval(b); }

inline lx_interval operator / (const lx_real &a, const lx_interval &b) 
                                      throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{  return lx_interval(a) / b; }


inline lx_interval & operator /=(lx_interval &a, const lx_interval &b) throw()
{  return a = a/b; }
inline lx_interval & operator /=(lx_interval &a, const l_interval &b) throw()
{  return a = a/b; }
inline lx_interval & operator /=(lx_interval &a, const l_real &b) throw()
{  return a = a/b; }
inline lx_interval & operator /=(lx_interval &a, const real &b) throw()
{  return a = a/b; }
inline lx_interval & operator /=(lx_interval &a, const interval &b) throw()
{  return a = a/b; }
inline lx_interval & operator /=(lx_interval &a, const lx_real &b) throw()
{  return a = a/b; } 

// ------------------ Vergleichsoperatoren ----------------------------------

inline bool operator ! (const lx_interval& a) throw()
{ return !a.li; }

inline bool operator == (const lx_interval &a, const lx_interval &b) throw()
{
    l_interval al(li_part(a)), bl(li_part(b));
    int exa(expo_gr(al)), exb(expo_gr(bl));
    real na(expo(a)), nb(expo(b)), d;
    bool a_0(exa<-100000), b_0(exb<-100000), res;

    if (a_0 || b_0) res = (a_0 == b_0);
    else // a,b <> 0:
    {
	d = exa-exb;
	if (d>0)
	{  // bl nach oben skalieren:
	    Times2pown(bl,d);
	    d = nb - d;
	    nb = (abs(d) > Max_Int_R)? MaxReal : d; // b = 2^nb * bl;
	} else // d<=0:
	{  // al nach oben skalieren:
	    Times2pown(al,-d);
	    d = na + d;
	    na = (abs(d) > Max_Int_R)? MaxReal : d; // a = 2^na * al;
	}
	res = (na==nb && al==bl);
    }
    return res;
}

inline bool operator == (const lx_interval &a, const l_interval &b) throw()
{  return (a == lx_interval(b)); }

inline bool operator == (const l_interval &a, const lx_interval &b) throw()
{  return ( lx_interval(a) == b); }

inline bool operator == (const lx_interval &a, const interval &b) throw()
{  return (a == lx_interval(b)); }

inline bool operator == (const interval &a, const lx_interval &b) throw()
{  return ( lx_interval(a) == b); }

inline bool operator == (const lx_interval &a, const real &b) throw()
{  return ( a == lx_interval(b)); }

inline bool operator == (const lx_interval &a, const l_real &b) throw()
{  return ( a == lx_interval(b)); }

inline bool operator == (const real &a, const lx_interval &b) throw()
{  return ( lx_interval(a) == b); }

inline bool operator == (const l_real &a, const lx_interval &b) throw()
{  return ( lx_interval(a) == b); }

inline bool operator == (const lx_interval &a, const lx_real &b) throw()
{  return ( a == lx_interval(b)); }

inline bool operator == (const lx_real &a, const lx_interval &b) throw()
{  return ( lx_interval(a) == b); }

inline bool operator != (const lx_interval &a, const lx_interval &b) throw()
{ return !(a==b); }

inline bool operator != (const lx_interval &a, const l_interval &b) throw()
{ return !(a==lx_interval(b)); }

inline bool operator != (const l_interval &a, const lx_interval &b) throw()
{ return !(lx_interval(a) == b); }

inline bool operator != (const lx_interval &a, const interval &b) throw()
{ return !(a==lx_interval(b)); }

inline bool operator != (const interval &a, const lx_interval &b) throw()
{ return !(lx_interval(a) == b); }

inline bool operator != (const lx_interval &a, const real &b) throw()
{ return !(a == lx_interval(b)); }

inline bool operator != (const real &a, const lx_interval &b) throw()
{ return !(lx_interval(a) == b); }

inline bool operator != (const lx_interval &a, const l_real &b) throw()
{ return !(a == lx_interval(b)); }

inline bool operator != (const l_real &a, const lx_interval &b) throw()
{ return !(lx_interval(a) == b); }

inline bool operator != (const lx_interval &a, const lx_real &b) throw()
{ return !(a == lx_interval(b)); }

inline bool operator != (const lx_real &a, const lx_interval &b) throw()
{ return !(lx_interval(a) == b); }

// ------------------------------------------------------------------------
// --------------------------- set comparisons ----------------------------
// ------------------------------------------------------------------------

// ---- lx_interval--lx_interval
inline bool operator <  (const lx_interval &a, const lx_interval &b) throw()
{ return( Inf(a)>Inf(b) && Sup(a)<Sup(b) ); }
inline bool operator <= (const lx_interval &a, const lx_interval &b) throw()
{ return( Inf(a)>=Inf(b) && Sup(a)<=Sup(b) ); }
inline bool operator > (const lx_interval& a, const lx_interval& b) throw()
{ return(Inf(a)<Inf(b) && Sup(a)>Sup(b)); }
inline bool operator >= (const lx_interval& a, const lx_interval& b) throw()
{ return(Inf(a)<=Inf(b) && Sup(a)>=Sup(b)); }

// ---- lx_interval--l_interval
inline bool operator <  (const lx_interval &a, const l_interval &b) throw()
{ return( Inf(a)>Inf(b) && Sup(a)<Sup(b) ); }
inline bool operator <= (const lx_interval &a, const l_interval &b) throw()
{ return( Inf(a)>=Inf(b) && Sup(a)<=Sup(b) ); }
inline bool operator <  (const l_interval &a, const lx_interval &b) throw()
{ return( Inf(a)>Inf(b) && Sup(a)<Sup(b) ); }
inline bool operator <= (const l_interval &a, const lx_interval &b) throw()
{ return( Inf(a)>=Inf(b) && Sup(a)<=Sup(b) ); }
inline bool operator > (const lx_interval& a, const l_interval& b) throw()
{ return(Inf(a)<Inf(b) && Sup(a)>Sup(b)); }
inline bool operator >= (const lx_interval& a, const l_interval& b) throw()
{ return(Inf(a)<=Inf(b) && Sup(a)>=Sup(b)); }
inline bool operator > (const l_interval& a, const lx_interval& b) throw()
{ return(Inf(a)<Inf(b) && Sup(a)>Sup(b)); }
inline bool operator >= (const l_interval& a, const lx_interval& b) throw()
{ return(Inf(a)<=Inf(b) && Sup(a)>=Sup(b)); }

// ---- lx_interval--interval
inline bool operator <  (const lx_interval &a, const interval &b) throw()
{ return( Inf(a)>Inf(b) && Sup(a)<Sup(b) ); }
inline bool operator <= (const lx_interval &a, const interval &b) throw()
{ return( Inf(a)>=Inf(b) && Sup(a)<=Sup(b) ); }
inline bool operator <  (const interval &a, const lx_interval &b) throw()
{ return( Inf(a)>Inf(b) && Sup(a)<Sup(b) ); }
inline bool operator <= (const interval &a, const lx_interval &b) throw()
{ return( Inf(a)>=Inf(b) && Sup(a)<=Sup(b) ); }
inline bool operator > (const lx_interval& a, const interval& b) throw()
{ return(Inf(a)<Inf(b) && Sup(a)>Sup(b)); }
inline bool operator >= (const lx_interval& a, const interval& b) throw()
{ return(Inf(a)<=Inf(b) && Sup(a)>=Sup(b)); }
inline bool operator > (const interval& a, const lx_interval& b) throw()
{ return(Inf(a)<Inf(b) && Sup(a)>Sup(b)); }
inline bool operator >= (const interval& a, const lx_interval& b) throw()
{ return(Inf(a)<=Inf(b) && Sup(a)>=Sup(b)); }

// ---- lx_interval--real
inline bool operator <  (const real &a, const lx_interval &b) throw()
{ return( a>Inf(b) && a<Sup(b) ); }
inline bool operator <= (const real &a, const lx_interval &b) throw()
{ return( a>=Inf(b) && a<=Sup(b) ); }
inline bool operator >  (const lx_interval &a, const real &b) throw()
{ return( b>Inf(a) && b<Sup(a) ); }
inline bool operator >=  (const lx_interval &a, const real &b) throw()
{ return( b>=Inf(a) && b<=Sup(a) ); }

// ---- lx_interval--l_real
inline bool operator <  (const l_real &a, const lx_interval &b) throw()
{ return( a>Inf(b) && a<Sup(b) ); }
inline bool operator <= (const l_real &a, const lx_interval &b) throw()
{ return( a>=Inf(b) && a<=Sup(b) ); }
inline bool operator >  (const lx_interval &a, const l_real &b) throw()
{ return( b>Inf(a) && b<Sup(a) ); }
inline bool operator >=  (const lx_interval &a, const l_real &b) throw()
{ return( b>=Inf(a) && b<=Sup(a) ); }

// ---- lx_interval--lx_real
inline bool operator <  (const lx_real &a, const lx_interval &b) throw()
{ return( a>Inf(b) && a<Sup(b) ); }
inline bool operator <= (const lx_real &a, const lx_interval &b) throw()
{ return( a>=Inf(b) && a<=Sup(b) ); }
inline bool operator >  (const lx_interval &a, const lx_real &b) throw()
{ return( b>Inf(a) && b<Sup(a) ); }
inline bool operator >=  (const lx_interval &a, const lx_real &b) throw()
{ return( b>=Inf(a) && b<=Sup(a) ); }


inline lx_interval adjust(const lx_interval &a) throw()
{  return lx_interval(a.ex,adjust(a.li));  }

inline lx_interval abs(const lx_interval &a) throw()
{  return lx_interval(a.ex,abs(a.li));  }

inline bool point_intv(const lx_interval &a)
{  return point_intv(a.li);  }

inline void times2pown(lx_interval &a, const real& n) throw()
{  a = lx_interval(add_real(n,a.ex),a.li); }

inline void times2pown_neg(lx_interval& a, const real& n) throw()
// Calculating an inclusion of  a*2^n, n = 0,-1,-2,...,-9007199254740991.0;
// n MUST be an integer and n MUST not be positive! 
// These conditions are not tested in this function!
// Blomquist, 09.06.2008;
{
    int exal(expo_gr(a.li));
    real exa,d,n_d;
    l_interval lia(a.li);
    int k;

    if (exal>-100000) // a != [0,0]
    {
	exa = a.ex;
	if (exa < -Max_Int_R - n) // exa+n < -9007199254740991;
	{   // It holds: -Max_Int_R - n in {-Max_Int_R, ...,-1,0},
            // Thus, -Max_Int_R - n is always in integer value.
            // Furthermore it holds: exa in {-Max_Int_R,...,-2,-1}.
            d = -Max_Int_R - exa; // d in {-Max_Int_R+1,..., -1,0} 
            n_d = n-d;
            // With exa+n < -Max_Int_R  and with  exa+d = -Max_Int_R
            // it follows:  n-d < 0, and:
            // n-d in {-Max_Int_R,-Max_Int_R+1, ... , -1};
            // Thus, n-d is a negative and integer value.
            if (n_d < -2147483647)
            {
               if (Inf(lia)>=0)
                  lia = l_interval(0,minreal);
               else 
                  if (Sup(lia)<=0) 
                      lia = l_interval(-minreal,0);
                  else lia = l_interval(-minreal,minreal);
            }
            else  // n_d >=-2147483647:
            {
                k = (int) _double(n_d);
                Times2pown(lia,k);
            }
            a = lx_interval(-Max_Int_R,lia);
	}
	else // n+a.ex >= -9007199254740991, so an integer overflow
             // is not possible here!
	    a = lx_interval(n+a.ex,lia);
    }
} // times2pown_neg(...)

inline lx_real Inf(const lx_interval &a) throw()
{  return lx_real(a.ex,Inf(a.li));  }

inline lx_real Sup(const lx_interval &a) throw()
{  return lx_real(a.ex,Sup(a.li));  }

inline lx_real RelDiam( const lx_interval &a )
{ 
	lx_real x;
	if (0<=a.li)
		x = lx_real(a.ex,RelDiam(a.li));
	else x = RelDiam(a.li);
	return x;
}

inline lx_real diam(const lx_interval &a) throw()
{  return lx_real(a.ex,diam(a.li)); }

inline lx_real mid(const lx_interval& a) throw()
{  return lx_real(a.ex,mid(a.li)); }

inline bool IsEmpty(const lx_interval& a) throw()
{ return  Inf(a.li) > Sup(a.li); }

// ----------------------------- Convex hull -------------------------------

inline lx_interval operator |(const lx_interval &a, const lx_interval &b) throw() 
{
   return lx_interval( (Inf(a)<Inf(b)) ? Inf(a) : Inf(b),
                      (Sup(a)>Sup(b)) ? Sup(a) : Sup(b) );
}

inline lx_interval operator |(const lx_interval &a, const l_interval &b) throw() 
{
    lx_interval Lb(0.0,b);
    return lx_interval( (Inf(a)<Inf(Lb)) ? Inf(a) : Inf(Lb),
		       (Sup(a)>Sup(Lb)) ? Sup(a) : Sup(Lb) );
}

inline lx_interval operator |(const l_interval &a, const lx_interval &b) throw() 
{
    lx_interval La(0.0,a);
    return lx_interval( (Inf(La)<Inf(b)) ? Inf(La) : Inf(b),
		       (Sup(La)>Sup(b)) ? Sup(La) : Sup(b) );
}

inline lx_interval operator |(const lx_interval &a, const interval &b) throw() 
{
    lx_interval Lb(0.0,l_interval(b));
    return lx_interval( (Inf(a)<Inf(Lb)) ? Inf(a) : Inf(Lb),
		       (Sup(a)>Sup(Lb)) ? Sup(a) : Sup(Lb) );
}

inline lx_interval operator |(const interval &a, const lx_interval &b) throw() 
{
    lx_interval La(0.0,l_interval(a));
    return lx_interval( (Inf(La)<Inf(b)) ? Inf(La) : Inf(b),
		       (Sup(La)>Sup(b)) ? Sup(La) : Sup(b) );
}

inline lx_interval & operator |= (lx_interval &a,const lx_interval &b) throw() 
{
   Inf(a)=(Inf(a)<Inf(b))?Inf(a):Inf(b),Sup(a)=(Sup(a)>Sup(b))?Sup(a):Sup(b);
   return a;
}

inline lx_interval & operator |= (lx_interval &a,const l_interval &b) throw() 
{
    lx_interval Lb(0,b);
    Inf(a) = (Inf(a)<Inf(Lb)) ? Inf(a) : Inf(Lb),
	Sup(a) = (Sup(a)>Sup(Lb)) ? Sup(a) : Sup(Lb);
    return a;
}

inline lx_interval & operator |= (lx_interval &a,const interval &b) throw() 
{
    lx_interval Lb(0,l_interval(b));
    Inf(a) = (Inf(a)<Inf(Lb)) ? Inf(a) : Inf(Lb),
	Sup(a) = (Sup(a)>Sup(Lb)) ? Sup(a) : Sup(Lb);
    return a;
}

inline lx_interval operator | (const lx_real &a, const lx_interval &b) throw() 
{
    return lx_interval( (a<Inf(b)) ? a : Inf(b),
                       (a>Sup(b)) ? a : Sup(b) );
}

inline lx_interval operator | (const real &a, const lx_interval &b) throw() 
{
    lx_real La(a);
    return lx_interval( (La<Inf(b)) ? La : Inf(b),
                       (La>Sup(b)) ? La : Sup(b) );
}

inline lx_interval operator | (const l_real &a, const lx_interval &b) throw() 
{
    lx_real La(0,a);
    return lx_interval( (La<Inf(b)) ? La : Inf(b),
                       (La>Sup(b)) ? La : Sup(b) );
}

inline lx_interval operator | (const lx_interval &a, const lx_real &b) throw() 
{
    return lx_interval( (Inf(a)<b) ? Inf(a) : b,
                       (Sup(a)>b) ? Sup(a) : b );
}

inline lx_interval operator | (const lx_interval &a, const real &b) throw() 
{
    lx_real Lb(b);
    return lx_interval( (Inf(a)<Lb) ? Inf(a) : Lb,
                       (Sup(a)>Lb) ? Sup(a) : Lb );
}

inline lx_interval operator | (const lx_interval &a, const l_real &b) throw() 
{
    lx_real Lb(0.0,b);
    return lx_interval( (Inf(a)<Lb) ? Inf(a) : Lb,
                       (Sup(a)>Lb) ? Sup(a) : Lb );
}

inline lx_interval & operator |= (lx_interval &a, const real &b) throw() 
{
    lx_real Lb(b);
    Inf(a) = (Inf(a)<Lb) ? Inf(a) : Lb, Sup(a) = (Sup(a)>Lb) ? Sup(a) : Lb;
    return a;
}

inline lx_interval & operator |= (lx_interval &a, const l_real &b) throw() 
{
    lx_real Lb(0.0,b);
    Inf(a) = (Inf(a)<Lb) ? Inf(a) : Lb, Sup(a) = (Sup(a)>Lb) ? Sup(a) : Lb;
    return a;
}


inline lx_interval & operator |= (lx_interval &a, const lx_real &b) throw() 
{
    Inf(a) = (Inf(a)<b) ? Inf(a) : b, Sup(a) = (Sup(a)>b) ? Sup(a) : b;
    return a;
}

inline lx_interval operator |(const lx_real &a, const lx_real &b) throw()
{
    if(a>b) return lx_interval(b,a);
    else    return lx_interval(a,b);
}

// --------------------------- Intersection -----------------------------

inline lx_interval operator & (const lx_interval &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    return lx_interval( (Inf(a)>Inf(b)) ? Inf(a) : Inf(b),
                       (Sup(a)<Sup(b)) ? Sup(a) : Sup(b));
}

inline lx_interval operator & (const lx_interval &a, const l_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_interval Lb(0.0,b);
    return lx_interval( (Inf(a)>Inf(Lb)) ? Inf(a) : Inf(Lb),
                       (Sup(a)<Sup(Lb)) ? Sup(a) : Sup(Lb));
}

inline lx_interval & operator &= (lx_interval &a, const l_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_interval Lb(0.0,b);
    Inf(a) = (Inf(a)>Inf(Lb)) ? Inf(a) : Inf(Lb),
	Sup(a) = (Sup(a)<Sup(Lb)) ? Sup(a) : Sup(Lb);
    if (Inf(a)>Sup(a))
	cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("lx_interval & operator &=(lx_interval &a,const l_interval &b)"));
    return a;
}

inline lx_interval operator & (const l_interval &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_interval La(0.0,a);
    return lx_interval( (Inf(La)>Inf(b)) ? Inf(La) : Inf(b),
                       (Sup(La)<Sup(b)) ? Sup(La) : Sup(b));
}

inline lx_interval operator & (const lx_interval &a, const interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_interval Lb(0.0,l_interval(b));
    return lx_interval( (Inf(a)>Inf(Lb)) ? Inf(a) : Inf(Lb),
                       (Sup(a)<Sup(Lb)) ? Sup(a) : Sup(Lb));
}

inline lx_interval & operator &= (lx_interval &a, const interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_interval Lb(0.0,l_interval(b));
    Inf(a) = (Inf(a)>Inf(Lb)) ? Inf(a) : Inf(Lb),
	Sup(a) = (Sup(a)<Sup(Lb)) ? Sup(a) : Sup(Lb);
    if (Inf(a)>Sup(a))
	cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("lx_interval & operator &=(lx_interval &a,const interval &b)"));
    return a;
}

inline lx_interval operator & (const interval &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_interval La(0.0,l_interval(a));
    return lx_interval( (Inf(La)>Inf(b)) ? Inf(La) : Inf(b),
                       (Sup(La)<Sup(b)) ? Sup(La) : Sup(b));
}

inline lx_interval & operator &= (lx_interval &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    Inf(a)=(Inf(a)>Inf(b))?Inf(a):Inf(b),Sup(a)=(Sup(a)<Sup(b))?Sup(a):Sup(b);
    if (Inf(a)>Sup(a))
	cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("lx_interval & operator &=(lx_interval &a,const lx_interval &b)"));
   return a;
}

inline lx_interval operator & (const lx_interval &a, const lx_real &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
   return lx_interval( (Inf(a)>b) ? Inf(a) : b,
                      (Sup(a)<b) ? Sup(a) : b );
}

inline lx_interval operator & (const lx_interval &a, const real &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_real Lb(b);
    return lx_interval( (Inf(a)>Lb) ? Inf(a) : Lb,
		       (Sup(a)<Lb) ? Sup(a) : Lb );
}

inline lx_interval operator & (const lx_interval &a, const l_real &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_real Lb(0.0,b);
    return lx_interval( (Inf(a)>Lb) ? Inf(a) : Lb,
		       (Sup(a)<Lb) ? Sup(a) : Lb );
}

inline lx_interval operator & (const lx_real &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    return lx_interval( (a>Inf(b)) ? a : Inf(b),
                       (a<Sup(b)) ? a : Sup(b) );
}

inline lx_interval operator & (const real &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_real La(a);
    return lx_interval( (La>Inf(b)) ? La : Inf(b),
                       (La<Sup(b)) ? La : Sup(b) );
}

inline lx_interval operator & (const l_real &a, const lx_interval &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_real La(0.0,a);
    return lx_interval( (La>Inf(b)) ? La : Inf(b),
                       (La<Sup(b)) ? La : Sup(b) );
}

inline lx_interval & operator &= (lx_interval &a,const lx_real &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    Inf(a) = (Inf(a)>b) ? Inf(a) : b, Sup(a) = (Sup(a)<b) ? Sup(a) : b;
    if(Inf(a)>Sup(a))
	cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("lx_interval & operator &=(lx_interval &a,const lx_real &b)"));
    return a;
}

inline lx_interval & operator &= (lx_interval &a, const real &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_real Lb(b);
    Inf(a) = (Inf(a)>Lb) ? Inf(a) : Lb, Sup(a) = (Sup(a)<Lb) ? Sup(a) : Lb;
    if(Inf(a)>Sup(a))
	cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("lx_interval & operator &=(lx_interval &a,const real &b)"));
    return a;
}

inline lx_interval & operator &= (lx_interval &a, const l_real &b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL) 
{
    lx_real Lb(0.0,b);
    Inf(a) = (Inf(a)>Lb) ? Inf(a) : Lb, Sup(a) = (Sup(a)<Lb) ? Sup(a) : Lb;
    if(Inf(a)>Sup(a))
	cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("lx_interval & operator &=(lx_interval &a,const l_real &b)"));
    return a;
}

// ------------------------- SetInf, SetSup -----------------------------

inline lx_interval & SetInf(lx_interval& a, const lx_real& b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
    return a = lx_interval(b,Sup(a));
}

inline lx_interval & SetInf(lx_interval& a, const l_real& b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
    return a = lx_interval(lx_real(0.0,b),Sup(a));
}

inline lx_interval & SetInf(lx_interval& a, const real& b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
    return a = lx_interval(lx_real(0.0,l_real(b)),Sup(a));
}

inline lx_interval & SetSup(lx_interval& a, const lx_real& b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
    return a = lx_interval(Inf(a),b);
}

inline lx_interval & SetSup(lx_interval& a, const l_real& b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
    return a = lx_interval(Inf(a),lx_real(0.0,b));
}

inline lx_interval & SetSup(lx_interval& a, const real& b) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
    return a = lx_interval(Inf(a),lx_real(0.0,l_real(b)));
}


} // end namespace cxsc
