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

/* CVS $Id: lx_cinterval.inl,v 1.9 2014/01/30 17:23:47 cxsc Exp $ */

/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

namespace cxsc {
	
// --------------------------------------------------------------------------
// ------ Inline functions and operators related to type lx_cinterval -------
// --------------------------------------------------------------------------

// --------------------------------------------------------------------------
// -------------------------- Constructors ----------------------------------
// --------------------------------------------------------------------------

inline lx_cinterval::lx_cinterval(const lx_interval & a, 
                                const lx_interval & b) throw()
      : re(a), im(b) { }

inline lx_cinterval::lx_cinterval(const l_interval & a, 
                                const l_interval & b) throw()
      : re(a), im(b) { }

inline lx_cinterval::lx_cinterval(const interval & a, 
                                const interval & b) throw()
      : re(a), im(b) { }

inline lx_cinterval::lx_cinterval(const l_real & a, 
                                const l_real & b) throw()
      : re(a), im(b) { }

inline lx_cinterval::lx_cinterval(const lx_real & a, 
                                const lx_real & b) throw()
      : re(a), im(b) { } 

inline lx_cinterval::lx_cinterval(const real & a, 
                                const real & b) throw()
      : re(a), im(b) { }

inline lx_cinterval::lx_cinterval(const l_cinterval & a) throw()
      : re(Re(a)), im(Im(a)) { }

inline lx_cinterval::lx_cinterval(const complex & a) throw()
      : re(Re(a)), im(Im(a)) { }

inline lx_cinterval::lx_cinterval(const l_complex & a) throw()
      : re(Re(a)), im(Im(a)) { }

inline lx_cinterval::lx_cinterval(const lx_complex & a) throw()
      : re(Re(a)), im(Im(a)) { }

inline lx_cinterval::lx_cinterval(const lx_complex & a, const lx_complex & b)
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
    : re(Re(a),Re(b)),
      im(Im(a),Im(b))
{
   if(Inf(re)>Sup(re) || Inf(im)>Sup(im))
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("lx_cinterval::lx_cinterval(const lx_complex & a,const lx_complex & b)"));
}

inline lx_cinterval::lx_cinterval(const l_complex & a, const l_complex & b)
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
    : re(Re(a),Re(b)),
      im(Im(a),Im(b))
{
   if(Inf(re)>Sup(re) || Inf(im)>Sup(im))
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("lx_cinterval::lx_cinterval(const l_complex & a,const l_complex & b)"));
}

inline lx_cinterval::lx_cinterval(const complex & a, const complex & b)
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
    : re(Re(a),Re(b)),
      im(Im(a),Im(b))
{
   if(Inf(re)>Sup(re) || Inf(im)>Sup(im))
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("lx_cinterval::lx_cinterval(const complex & a,const complex & b)"));
}

inline lx_cinterval::lx_cinterval(const cinterval & a) throw()
      : re(Re(a)), im(Im(a)) { }

inline lx_cinterval::lx_cinterval(const real& na, const l_interval &la, 
                                  const real& nb, const l_interval &lb) 
                                  throw() : re(na,la), im(nb,lb) { }

inline lx_cinterval::lx_cinterval(const real &n, const l_interval &a, 
                                         const lx_interval &b) 
                                         throw() : re(n,a), im(b) { }

inline lx_cinterval::lx_cinterval(const lx_interval &a,
                                  const real &n, const l_interval &b) 
                                        throw() : re(a), im(n,b) { }

inline lx_cinterval::lx_cinterval(const real &nr, const string &sr, 
                                  const real &ni, const string &si) 
          throw() : re(lx_interval(nr,sr)), im(lx_interval(ni,si)) { }

inline lx_cinterval::lx_cinterval(const lx_interval & a) throw()
                                                  : re(a), im(0) { }

inline lx_cinterval::lx_cinterval(const l_interval & a) throw()
                                                  : re(a), im(0) { }

inline lx_cinterval::lx_cinterval(const interval & a) throw()
                                                  : re(a), im(0) { }

inline lx_cinterval::lx_cinterval(const lx_real & a) throw()
                                                  : re(a), im(0) { }

inline lx_cinterval::lx_cinterval(const l_real & a) throw()
                                                  : re(a), im(0) { }

inline lx_cinterval::lx_cinterval(const real & a) throw()
                                                  : re(a), im(0) { }

inline lx_cinterval::lx_cinterval(const real& n, const l_interval &a) throw()
                                                  : re(n,a), im(0) { }



// -----------------------------------------------------------------------
// ------------------------------ Assignments ----------------------------
// -----------------------------------------------------------------------

inline lx_cinterval & lx_cinterval::operator = (const lx_cinterval & a) throw()
{
    re = a.re;
    im = a.im;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const l_cinterval & a) throw()
{
    re = Re(a);
    im = Im(a);
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const cinterval & a) throw()
{
    re = Re(a);
    im = Im(a);
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const lx_interval & a) throw()
{
    re = a;
    im = 0.0;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const l_interval & a) throw()
{
    re = a;
    im = 0.0;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const interval & a) throw()
{
    re = a;
    im = 0.0;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const lx_real & a) throw()
{
    re = a;
    im = 0.0;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const l_real & a) throw()
{
    re = a;
    im = 0.0;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const real & a) throw()
{
    re = a;
    im = 0.0;
    return *this;
} 

inline lx_cinterval & lx_cinterval::operator = (const lx_complex & a ) throw()
{
    re = Re(a);
    im = Im(a);
    return *this;
}

inline lx_cinterval & lx_cinterval::operator = (const l_complex & a ) throw()
{
    re = Re(a);
    im = Im(a);
    return *this;
}

inline lx_cinterval & lx_cinterval::operator = (const complex & a ) throw()
{
    re = Re(a);
    im = Im(a);
    return *this;
}



// -----------------------------------------------------------------------
// ---------------------------- Functions --------------------------------
// -----------------------------------------------------------------------

inline lx_complex Inf(const lx_cinterval &a) throw()
{ return lx_complex(Inf(a.re),Inf(a.im)); }

inline lx_complex Sup(const lx_cinterval &a) throw()
{ return lx_complex(Sup(a.re),Sup(a.im)); }

inline lx_interval Re(const lx_cinterval &a) throw() 
{ return a.re; }

inline lx_interval Im(const lx_cinterval &a) throw() 
{ return a.im; }


inline lx_cinterval & SetRe(lx_cinterval &a, const lx_interval &b) 
{ a.re=b; return a; }
inline lx_cinterval & SetRe(lx_cinterval &a, const l_interval &b) 
{ a.re=b; return a; }
inline lx_cinterval & SetRe(lx_cinterval &a, const interval &b) 
{ a.re=b; return a; }
inline lx_cinterval & SetRe(lx_cinterval &a, const lx_real &b) 
{ a.re=b; return a; }
inline lx_cinterval & SetRe(lx_cinterval &a, const l_real &b) 
{ a.re=b; return a; }
inline lx_cinterval & SetRe(lx_cinterval &a, const real &b) 
{ a.re=b; return a; }

inline lx_cinterval & SetIm(lx_cinterval &a, const lx_interval &b) 
{ a.im=b; return a; }
inline lx_cinterval & SetIm(lx_cinterval &a, const l_interval &b) 
{ a.im=b; return a; }
inline lx_cinterval & SetIm(lx_cinterval &a, const interval &b) 
{ a.im=b; return a; }
inline lx_cinterval & SetIm(lx_cinterval &a, const lx_real &b) 
{ a.im=b; return a; }
inline lx_cinterval & SetIm(lx_cinterval &a, const l_real &b) 
{ a.im=b; return a; }
inline lx_cinterval & SetIm(lx_cinterval &a, const real &b) 
{ a.im=b; return a; }


inline lx_real InfRe(const lx_cinterval &a) throw() 
{ return Inf(a.re); }
inline lx_real InfIm(const lx_cinterval &a) throw()
{ return Inf(a.im); }
inline lx_real SupRe(const lx_cinterval &a) throw()
{ return Sup(a.re); }
inline lx_real SupIm(const lx_cinterval &a) throw()
{ return Sup(a.im); }

inline lx_complex mid(const lx_cinterval &a) throw()
{ return lx_complex(mid(a.re),mid(a.im)); }

inline lx_complex diam(const lx_cinterval &a) throw()
{ return lx_complex(diam(a.re),diam(a.im)); }

inline real expo_Re(const lx_cinterval &a) throw()
{ return expo(a.re); }

inline real expo_Im(const lx_cinterval &a) throw() 
{ return expo(a.im); }

inline l_interval li_part_Re(const lx_cinterval &a) throw()
{ return li_part(a.re); }

inline l_interval li_part_Im(const lx_cinterval &a) throw()
{ return li_part(a.im); }

inline lx_cinterval adjust(const lx_cinterval &a) throw()
{ return lx_cinterval(adjust(a.re),adjust(a.im)); }

inline lx_cinterval conj(const lx_cinterval &a) throw()
{ return lx_cinterval(a.re,-a.im); }

inline void times2pown(lx_cinterval& x, const real &n) throw()
{ 
    lx_interval a(x.re),b(x.im);
    times2pown(a,n);
    times2pown(b,n);
    x = lx_cinterval(a,b);	
}

inline lx_interval abs(const lx_cinterval &a) throw()
{ 
	return sqrtx2y2(a.re,a.im); 
}

// -----------------------------------------------------------------------
// ------------------------ Monadic Operators ----------------------------
// -----------------------------------------------------------------------

inline lx_cinterval operator-(const lx_cinterval & a) throw()
{  return lx_cinterval(-a.re,-a.im); }

inline lx_cinterval operator+(const lx_cinterval & a) throw()
{  return a; }


// -----------------------------------------------------------------------
// ----------------------- Arithmetic Operators --------------------------
// -----------------------------------------------------------------------

inline lx_cinterval operator + (const lx_cinterval &a, const lx_cinterval &b) 
       throw()
       { return lx_cinterval(a.re + b.re, a.im + b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const l_cinterval &b) 
       throw()
       { return lx_cinterval(a.re + Re(b), a.im + Im(b)); }

inline lx_cinterval operator + (const l_cinterval &a, const lx_cinterval &b) 
       throw()
       { return lx_cinterval(Re(a) + b.re, Im(a) + b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const cinterval &b) 
       throw()
       { return lx_cinterval(a.re + Re(b), a.im + Im(b)); }
 
inline lx_cinterval operator + (const cinterval &a, const lx_cinterval &b) 
       throw()
       { return lx_cinterval(Re(a) + b.re, Im(a) + b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const lx_interval &b) 
       throw()
       { return lx_cinterval(a.re + b, a.im); }

inline lx_cinterval operator + (const lx_interval &a, const lx_cinterval &b) 
       throw()
       { return lx_cinterval(a + b.re, b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const l_interval &b) 
       throw()
       { return lx_cinterval(a.re + b, a.im); }

inline lx_cinterval operator + (const l_interval &a, const lx_cinterval &b) 
       throw()
       { return lx_cinterval(a + b.re, b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const lx_real &b) throw()
       { return lx_cinterval(a.re + b, a.im); }

inline lx_cinterval operator + (const lx_real &a, const lx_cinterval &b) throw()
       { return lx_cinterval(a + b.re, b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const l_real &b) throw()
       { return lx_cinterval(a.re + b, a.im); }

inline lx_cinterval operator + (const l_real &a, const lx_cinterval &b) throw()
       { return lx_cinterval(a + b.re, b.im); } 

inline lx_cinterval operator + (const lx_cinterval &a, const real &b) throw()
       { return lx_cinterval(a.re + b, a.im); }

inline lx_cinterval operator + (const real &a, const lx_cinterval &b) throw()
       { return lx_cinterval(a + b.re, b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const complex &b) throw()
       { return lx_cinterval(a.re + Re(b), a.im + Im(b)); }
inline lx_cinterval operator + (const complex &a, const lx_cinterval &b) throw()
       { return lx_cinterval(Re(a) + b.re, Im(a) + b.im); } 

inline lx_cinterval operator + (const lx_cinterval &a, const l_complex &b) 
    throw() { return lx_cinterval(a.re + Re(b), a.im + Im(b)); }
inline lx_cinterval operator + (const l_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(Re(a) + b.re, Im(a) + b.im); }

inline lx_cinterval operator + (const lx_cinterval &a, const lx_complex &b) 
    throw() { return lx_cinterval(a.re + Re(b), a.im + Im(b)); }
inline lx_cinterval operator + (const lx_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(Re(a) + b.re, Im(a) + b.im); }


inline lx_cinterval & operator +=(lx_cinterval &a, const lx_cinterval &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const lx_interval &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const l_interval &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const l_cinterval &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const l_real &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const lx_real &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const real &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const interval &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const cinterval &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const complex &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const l_complex &b) throw()
{  return a = a+b; }

inline lx_cinterval & operator +=(lx_cinterval &a, const lx_complex &b) throw()
{  return a = a+b; }


inline lx_cinterval operator - (const lx_cinterval &a, const lx_cinterval &b) 
       throw() { return lx_cinterval(a.re - b.re, a.im - b.im); }

inline lx_cinterval operator - (const lx_cinterval &a, const l_cinterval &b) 
       throw() { return lx_cinterval(a.re - Re(b), a.im - Im(b)); }

inline lx_cinterval operator - (const l_cinterval &a, const lx_cinterval &b) 
       throw() { return lx_cinterval(Re(a) - b.re, Im(a) - b.im); }

inline lx_cinterval operator - (const lx_cinterval &a, const cinterval &b) 
       throw() { return lx_cinterval(a.re - Re(b), a.im - Im(b)); }
 
inline lx_cinterval operator - (const cinterval &a, const lx_cinterval &b) 
       throw() { return lx_cinterval(Re(a) - b.re, Im(a) - b.im); }

inline lx_cinterval operator - (const lx_cinterval &a, const lx_interval &b) 
       throw() { return lx_cinterval(a.re - b, a.im); }

inline lx_cinterval operator - (const lx_interval &a, const lx_cinterval &b) 
       throw() { return lx_cinterval(a - b.re, -b.im); }

inline lx_cinterval operator - (const lx_cinterval &a, const l_interval &b) 
       throw() { return lx_cinterval(a.re - b, a.im); }

inline lx_cinterval operator - (const l_interval &a, const lx_cinterval &b) 
       throw() { return lx_cinterval(a - b.re, -b.im); }

inline lx_cinterval operator - (const lx_cinterval &a, const lx_real &b) throw()
{ return lx_cinterval(a.re - b, a.im); }

inline lx_cinterval operator - (const lx_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a - b.re, -b.im); }

inline lx_cinterval operator - (const lx_cinterval &a, const l_real &b) throw()
{ return lx_cinterval(a.re - b, a.im); }

inline lx_cinterval operator - (const l_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a - b.re, -b.im); } 

inline lx_cinterval operator - (const lx_cinterval &a, const real &b) throw()
{ return lx_cinterval(a.re - b, a.im); }

inline lx_cinterval operator - (const real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a - b.re, -b.im); } 

inline lx_cinterval operator - (const lx_cinterval &a, const complex &b) throw()
       { return lx_cinterval(a.re - Re(b), a.im - Im(b)); }
inline lx_cinterval operator - (const complex &a, const lx_cinterval &b) throw()
       { return lx_cinterval(Re(a) - b.re, Im(a) - b.im); }
inline lx_cinterval operator - (const lx_cinterval &a, const l_complex &b) 
    throw() { return lx_cinterval(a.re - Re(b), a.im - Im(b)); }
inline lx_cinterval operator - (const l_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(Re(a) - b.re, Im(a) - b.im); }
inline lx_cinterval operator - (const lx_cinterval &a, const lx_complex &b) 
    throw() { return lx_cinterval(a.re - Re(b), a.im - Im(b)); }
inline lx_cinterval operator - (const lx_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(Re(a) - b.re, Im(a) - b.im); }


inline lx_cinterval & operator -=(lx_cinterval &a, const lx_cinterval &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const lx_interval &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const l_interval &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const l_cinterval &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const l_real &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const lx_real &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const real &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const interval &b) throw()
{  return a = a-b; }

inline lx_cinterval & operator -=(lx_cinterval &a, const cinterval &b) throw()
{  return a = a-b; }
inline lx_cinterval & operator -=(lx_cinterval &a, const complex &b) throw()
{  return a = a-b; }
inline lx_cinterval & operator -=(lx_cinterval &a, const l_complex &b) throw()
{  return a = a-b; }
inline lx_cinterval & operator -=(lx_cinterval &a, const lx_complex &b) throw()
{  return a = a-b; }


inline lx_cinterval operator * (const lx_cinterval &a,const lx_cinterval &b) 
throw()
{ return lx_cinterval(a.re*b.re - a.im*b.im, a.im*b.re + a.re*b.im); }

inline lx_cinterval operator * (const lx_cinterval &a, const l_cinterval &b) 
throw()
{ return a * lx_cinterval(b); }

inline lx_cinterval operator * (const l_cinterval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a) * b; }

inline lx_cinterval operator * (const lx_cinterval &a, const cinterval &b) 
throw()
{ return a * lx_cinterval(b); }
 
inline lx_cinterval operator * (const cinterval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a) * b; }


inline lx_cinterval operator * (const lx_cinterval &a, const lx_interval &b) 
throw()
{ return lx_cinterval(a.re*b, a.im*b); }

inline lx_cinterval operator * (const lx_interval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a*b.re, a*b.im); }

inline lx_cinterval operator * (const lx_cinterval &a, const l_interval &b) 
throw()
{ return lx_cinterval(a.re*b, a.im*b); }

inline lx_cinterval operator * (const l_interval &a, const lx_cinterval &b) 
throw() 
{ return lx_cinterval(a*b.re, a*b.im); }

inline lx_cinterval operator * (const lx_cinterval &a, const lx_real &b) throw()
{ return lx_cinterval(a.re*b, a.im*b); }

inline lx_cinterval operator * (const lx_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a*b.re, a*b.im); }

inline lx_cinterval operator * (const lx_cinterval &a, const l_real &b) throw()
{ return lx_cinterval(a.re*b, a.im*b); }

inline lx_cinterval operator * (const l_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a*b.re, a*b.im); } 

inline lx_cinterval operator * (const lx_cinterval &a, const real &b) throw()
{ return lx_cinterval(a.re*b, a.im*b); }

inline lx_cinterval operator * (const real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a*b.re, a*b.im); } 

inline lx_cinterval operator * (const lx_cinterval &a, const complex &b) throw()
{ return a * lx_cinterval(b); }
inline lx_cinterval operator * (const complex &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) * b; }
inline lx_cinterval operator * (const lx_cinterval &a, const l_complex &b) 
    throw() { return a * lx_cinterval(b); }
inline lx_cinterval operator * (const l_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(a) * b; } 
inline lx_cinterval operator * (const lx_cinterval &a, const lx_complex &b) 
    throw() { return a * lx_cinterval(b); }
inline lx_cinterval operator * (const lx_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(a) * b; } 


inline lx_cinterval & operator *=(lx_cinterval &a, const lx_cinterval &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const lx_interval &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const l_interval &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const l_cinterval &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const l_real &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const lx_real &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const real &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const interval &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const cinterval &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const complex &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const l_complex &b) throw()
{  return a = a*b; }

inline lx_cinterval & operator *=(lx_cinterval &a, const lx_complex &b) throw()
{  return a = a*b; }


inline lx_cinterval operator / (const lx_cinterval &a, const lx_cinterval &b) 
throw()
{
    lx_interval Ne(sqr(b.re) + sqr(b.im));

    return lx_cinterval( (a.re*b.re + a.im*b.im)/Ne, 
                        (a.im*b.re - a.re*b.im)/Ne );
}

inline lx_cinterval operator / (const lx_cinterval &a, const l_cinterval &b) 
throw()
{ return a / lx_cinterval(b); }

inline lx_cinterval operator / (const l_cinterval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const cinterval &b) 
throw()
{ return a / lx_cinterval(b); }
 
inline lx_cinterval operator / (const cinterval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const lx_interval &b) 
throw()
{ return lx_cinterval(a.re/b, a.im/b); }

inline lx_cinterval operator / (const lx_interval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const l_interval &b) 
throw()
{ return lx_cinterval(a.re/b, a.im/b); }

inline lx_cinterval operator / (const l_interval &a, const lx_cinterval &b) 
throw()
{ return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const l_real &b) throw()
{ return lx_cinterval(a.re/b, a.im/b); }

inline lx_cinterval operator / (const l_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const lx_real &b) throw()
{ return lx_cinterval(a.re/b, a.im/b); }

inline lx_cinterval operator / (const lx_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) / b; } 

inline lx_cinterval operator / (const lx_cinterval &a, const real &b) throw()
{ return lx_cinterval(a.re/b, a.im/b); }

inline lx_cinterval operator / (const real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) / b; } 

inline lx_cinterval operator / (const lx_cinterval &a, const complex &b) throw()
{ return a / lx_cinterval(b); }
inline lx_cinterval operator / (const complex &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const l_complex &b) 
    throw() { return a / lx_cinterval(b); }
inline lx_cinterval operator / (const l_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(a) / b; }

inline lx_cinterval operator / (const lx_cinterval &a, const lx_complex &b) 
    throw() { return a / lx_cinterval(b); }
inline lx_cinterval operator / (const lx_complex &a, const lx_cinterval &b) 
    throw() { return lx_cinterval(a) / b; } 


inline lx_cinterval & operator /=(lx_cinterval &a, const lx_cinterval &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const lx_interval &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const l_interval &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const l_cinterval &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const l_real &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const lx_real &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const real &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const interval &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const cinterval &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const complex &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const l_complex &b) throw()
{  return a = a/b; }

inline lx_cinterval & operator /=(lx_cinterval &a, const lx_complex &b) throw()
{  return a = a/b; }



inline bool operator ! (const lx_cinterval & a) throw()
{ return !a.re && !a.im; }


inline bool operator == (const lx_cinterval &a, const lx_cinterval &b) throw()
{ return a.re == b.re && a.im == b.im; }

inline bool operator == (const lx_cinterval &a, const l_cinterval &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const l_cinterval &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const lx_interval &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const lx_interval &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const l_interval &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const l_interval &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const l_real &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const l_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const lx_real &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const lx_real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const real &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const real &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const interval &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const interval &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const cinterval &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const cinterval &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const complex &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const complex &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const l_complex &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const l_complex &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }

inline bool operator == (const lx_cinterval &a, const lx_complex &b) throw()
{ return a == lx_cinterval(b); }
inline bool operator == (const lx_complex &a, const lx_cinterval &b) throw()
{ return lx_cinterval(a) == b; }


inline bool operator != (const lx_cinterval &a, const lx_cinterval &b) throw()
{ return a.re != b.re || a.im != b.im; }

inline bool operator != (const lx_cinterval &a, const l_cinterval &b) throw()
{ return !(a == b); }
inline bool operator != (const l_cinterval &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const lx_interval &b) throw()
{ return !(a == b); }
inline bool operator != (const lx_interval &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const l_interval &b) throw()
{ return !(a == b); }
inline bool operator != (const l_interval &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const l_real &b) throw()
{ return !(a == b); }
inline bool operator != (const l_real &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const lx_real &b) throw()
{ return !(a == b); }
inline bool operator != (const lx_real &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const real &b) throw()
{ return !(a == b); }
inline bool operator != (const real &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const interval &b) throw()
{ return !(a == b); }
inline bool operator != (const interval &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const cinterval &b) throw()
{ return !(a == b); }
inline bool operator != (const cinterval &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const complex &b) throw()
{ return !(a == b); }
inline bool operator != (const complex &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const l_complex &b) throw()
{ return !(a == b); }
inline bool operator != (const l_complex &a, const lx_cinterval &b) throw()
{ return !(a == b); }

inline bool operator != (const lx_cinterval &a, const lx_complex &b) throw()
{ return !(a == b); }
inline bool operator != (const lx_complex &a, const lx_cinterval &b) throw()
{ return !(a == b); }



// ------------------------- Set Operators -----------------------------------

inline bool operator < (const lx_cinterval & a, const lx_cinterval & b) throw()
{
   if (Inf(a.re) <= Inf(b.re) || Sup(a.re) >= Sup(b.re)) 
      return false;
   if (Inf(a.im) <= Inf(b.im) || Sup(a.im) >= Sup(b.im)) 
      return false;
      
   return true; 
}

inline bool operator > (const lx_cinterval & a, const lx_cinterval & b) throw() 
{ return b < a; }

inline bool operator <= (const lx_cinterval & a, const lx_cinterval & b) throw()
{
   if (Inf(a.re) < Inf(b.re) || Sup(a.re) > Sup(b.re)) 
      return false;
   if (Inf(a.im) < Inf(b.im) || Sup(a.im) > Sup(b.im)) 
      return false;
      
   return true; 
}

inline bool operator >= (const lx_cinterval & a, const lx_cinterval & b) throw()
{ return b <= a; }


inline bool operator  <(const l_cinterval & a, const lx_cinterval & b) throw() 
                                             { return lx_cinterval(a) < b; }
inline bool operator  >(const l_cinterval & a, const lx_cinterval & b) throw()
                                             { return lx_cinterval(a) > b; }
inline bool operator <=(const l_cinterval & a, const lx_cinterval & b) throw()
                                            { return lx_cinterval(a) <= b; }
inline bool operator >=(const l_cinterval & a, const lx_cinterval & b) throw()
                                            { return lx_cinterval(a) >= b; }

inline bool operator  <(const lx_cinterval & a, const l_cinterval & b) throw()
                                             { return a < lx_cinterval(b); }
inline bool operator  >(const lx_cinterval & a, const l_cinterval & b) throw() 
                                             { return a > lx_cinterval(b); }
inline bool operator <=(const lx_cinterval & a, const l_cinterval & b) throw() 
                                            { return a <= lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const l_cinterval & b) throw()
                                            { return a >= lx_cinterval(b); }

inline bool operator  <(const cinterval & a, const lx_cinterval & b) throw() 
                                             { return lx_cinterval(a) < b; }
inline bool operator  >(const cinterval & a, const lx_cinterval & b) throw()
                                             { return lx_cinterval(a) > b; }
inline bool operator <=(const cinterval & a, const lx_cinterval & b) throw()
                                            { return lx_cinterval(a) <= b; }
inline bool operator >=(const cinterval & a, const lx_cinterval & b) throw()
                                            { return lx_cinterval(a) >= b; }

inline bool operator  <(const lx_cinterval & a, const cinterval & b) throw()
                                             { return a < lx_cinterval(b); }
inline bool operator  >(const lx_cinterval & a, const cinterval & b) throw() 
                                             { return a > lx_cinterval(b); }
inline bool operator <=(const lx_cinterval & a, const cinterval & b) throw() 
                                            { return a <= lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const cinterval & b) throw()
                                            { return a >= lx_cinterval(b); }

inline bool operator  <(const lx_interval & a, const lx_cinterval & b) throw() 
                                              { return lx_cinterval(a) < b; }
inline bool operator  >(const lx_interval & a, const lx_cinterval & b) throw()
                                              { return lx_cinterval(a) > b; }
inline bool operator <=(const lx_interval & a, const lx_cinterval & b) throw()
                                             { return lx_cinterval(a) <= b; }
inline bool operator >=(const lx_interval & a, const lx_cinterval & b) throw()
                                             { return lx_cinterval(a) >= b; }

inline bool operator  <(const lx_cinterval & a, const lx_interval & b) throw()
                                             { return a < lx_cinterval(b); }
inline bool operator  >(const lx_cinterval & a, const lx_interval & b) throw() 
                                             { return a > lx_cinterval(b); }
inline bool operator <=(const lx_cinterval & a, const lx_interval & b) throw() 
                                            { return a <= lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const lx_interval & b) throw()
                                            { return a >= lx_cinterval(b); }

inline bool operator  <(const l_interval & a, const lx_cinterval & b) throw() 
                                              { return lx_cinterval(a) < b; }
inline bool operator  >(const l_interval & a, const lx_cinterval & b) throw()
                                              { return lx_cinterval(a) > b; }
inline bool operator <=(const l_interval & a, const lx_cinterval & b) throw()
                                             { return lx_cinterval(a) <= b; }
inline bool operator >=(const l_interval & a, const lx_cinterval & b) throw()
                                             { return lx_cinterval(a) >= b; }

inline bool operator  <(const lx_cinterval & a, const l_interval & b) throw()
                                             { return a < lx_cinterval(b); }
inline bool operator  >(const lx_cinterval & a, const l_interval & b) throw() 
                                             { return a > lx_cinterval(b); }
inline bool operator <=(const lx_cinterval & a, const l_interval & b) throw() 
                                            { return a <= lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const l_interval & b) throw()
                                            { return a >= lx_cinterval(b); }

inline bool operator  <(const interval & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator  >(const interval & a, const lx_cinterval & b) throw()
                                            { return lx_cinterval(a) > b; }
inline bool operator <=(const interval & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator >=(const interval & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) >= b; }

inline bool operator  <(const lx_cinterval & a, const interval & b) throw()
                                            { return a < lx_cinterval(b); }
inline bool operator  >(const lx_cinterval & a, const interval & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator <=(const lx_cinterval & a, const interval & b) throw() 
                                           { return a <= lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const interval & b) throw()
                                           { return a >= lx_cinterval(b); }

inline bool operator  <(const lx_real & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator <=(const lx_real & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator  >(const lx_cinterval & a, const lx_real & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const lx_real & b) throw()
                                           { return a >= lx_cinterval(b); }

inline bool operator  <(const l_real & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator <=(const l_real & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator  >(const lx_cinterval & a, const l_real & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const l_real & b) throw()
                                           { return a >= lx_cinterval(b); }

inline bool operator  <(const real & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator <=(const real & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator  >(const lx_cinterval & a, const real & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const real & b) throw()
                                           { return a >= lx_cinterval(b); }

inline bool operator  <(const complex & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator <=(const complex & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator  >(const lx_cinterval & a, const complex & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const complex & b) throw()
                                           { return a >= lx_cinterval(b); }

inline bool operator  <(const l_complex & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator <=(const l_complex & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator  >(const lx_cinterval & a, const l_complex & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const l_complex & b) throw()
                                           { return a >= lx_cinterval(b); }

inline bool operator  <(const lx_complex & a, const lx_cinterval & b) throw() 
                                            { return lx_cinterval(a) < b; }
inline bool operator <=(const lx_complex & a, const lx_cinterval & b) throw()
                                           { return lx_cinterval(a) <= b; }
inline bool operator  >(const lx_cinterval & a, const lx_complex & b) throw() 
                                            { return a > lx_cinterval(b); }
inline bool operator >=(const lx_cinterval & a, const lx_complex & b) throw()
                                           { return a >= lx_cinterval(b); }

// ------------------------- Intersection ------------------------------------

inline lx_cinterval operator & (const lx_cinterval& a, 
			       const lx_cinterval& b) throw()
{
    lx_cinterval tmp = a;
    SetInf(tmp.re, max(Inf(a.re),Inf(b.re)));
    SetInf(tmp.im, max(Inf(a.im),Inf(b.im)));
    SetSup(tmp.re, min(Sup(a.re),Sup(b.re)));
    SetSup(tmp.im, min(Sup(a.im),Sup(b.im)));
    if (Inf(tmp.re) > Sup(tmp.re) || Inf(tmp.im) > Sup(tmp.im))
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
      "lx_cinterval operator &(const lx_cinterval& a,const lx_cinterval& b)")); 
    return tmp;
}

inline lx_cinterval & operator &= (lx_cinterval& a, const lx_cinterval& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL) { return a = a&b; }

inline lx_cinterval operator & (const lx_cinterval& a, const lx_real& b) 
    throw() { return a & lx_cinterval(b,lx_real(0.0)); }
inline lx_cinterval operator & (const lx_real& a, const lx_cinterval& b ) 
    throw() { return lx_cinterval(a,lx_real(0.0)) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const lx_real& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const l_real& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const l_real& a, const lx_cinterval& b ) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const l_real& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const real& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const real& a, const lx_cinterval& b ) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const real& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const l_cinterval& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const l_cinterval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const l_cinterval& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const cinterval& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const cinterval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const cinterval& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const lx_interval& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_interval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const lx_interval& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const l_interval& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const l_interval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const l_interval& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const interval& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const interval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const interval& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const lx_complex& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_complex& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const lx_complex& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const l_complex& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const l_complex& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const l_complex& b) 
    throw() { return a = a & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_cinterval& a, const complex& b) 
    throw() { return a & lx_cinterval(b); }
inline lx_cinterval operator & (const complex& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) & b; }
inline lx_cinterval & operator &= (lx_cinterval& a, const complex& b) 
    throw() { return a = a & lx_cinterval(b); }

inline lx_cinterval operator & (const lx_interval& a, const complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const complex& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_interval& a, const l_complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const l_complex& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_interval& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const lx_complex& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_real& a, const cinterval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const cinterval& a, const lx_real& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_real& a, const l_cinterval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const l_cinterval& a, const lx_real& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_interval& a, const cinterval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const cinterval& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const lx_interval& a, const l_cinterval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const l_cinterval& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }

inline lx_cinterval operator & (const l_interval& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const lx_complex& a, const l_interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const l_cinterval& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const lx_complex& a, const l_cinterval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }

inline lx_cinterval operator & (const interval& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const lx_complex& a, const interval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & (const cinterval& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }
inline lx_cinterval operator & ( const lx_complex& a, const cinterval& b) 
    throw() { return lx_cinterval(a) & lx_cinterval(b); }


// -------------------------- Convex Hull ------------------------------------

inline lx_cinterval operator | (const lx_cinterval& a,
			        const lx_cinterval& b) throw()
{
   lx_cinterval tmp = a;
   SetInf(tmp.re, min(Inf(a.re), Inf(b.re)));
   SetInf(tmp.im, min(Inf(a.im), Inf(b.im)));
   SetSup(tmp.re, max(Sup(a.re), Sup(b.re)));
   SetSup(tmp.im, max(Sup(a.im), Sup(b.im)));
   return tmp;
}

inline lx_cinterval & operator |= (lx_cinterval& a, const lx_cinterval& b) 
    throw() { return a = a|b; }
inline lx_cinterval operator | (const lx_cinterval& a, const lx_real& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_real& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const lx_real& b) 
    throw() { return a = a|lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const l_real& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const l_real& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const l_real& b) 
    throw() { return a = a|lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const real& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const real& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const real& b) 
    throw() { return a = a|lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const l_cinterval& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const l_cinterval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const l_cinterval& b) 
    throw() { return a = a|lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const cinterval& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const cinterval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const cinterval& b) 
    throw() { return a = a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const lx_interval& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_interval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const lx_interval& b) 
    throw() { return a = a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const l_interval& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const l_interval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const l_interval& b) 
    throw() { return a = a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const interval& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const interval& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const interval& b) 
    throw() { return a = a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const lx_complex& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_complex& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const lx_complex& b) 
    throw() { return a = a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const l_complex& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const l_complex& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const l_complex& b) 
    throw() { return a = a | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_cinterval& a, const complex& b) 
    throw() { return a | lx_cinterval(b); }
inline lx_cinterval operator | (const complex& a, const lx_cinterval& b) 
    throw() { return lx_cinterval(a) | b; }
inline lx_cinterval & operator |= (lx_cinterval& a, const complex& b) 
    throw() { return a = a | lx_cinterval(b); }

inline lx_cinterval operator | (const lx_interval& a, const complex& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | ( const complex& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_interval& a, const l_complex& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | ( const l_complex& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_interval& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | ( const lx_complex& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_real& a, const cinterval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const cinterval& a, const lx_real& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_real& a, const l_cinterval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const l_cinterval& a, const lx_real& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_interval& a, const cinterval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const cinterval& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_interval& a, const l_cinterval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const l_cinterval& a, const lx_interval& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }


inline lx_cinterval operator | (const lx_real& a, const complex& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const complex& a, const lx_real& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_real& a, const l_complex& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const l_complex& a, const lx_real& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_real& a, const lx_complex& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }
inline lx_cinterval operator | (const lx_complex& a, const lx_real& b) 
    throw() { return lx_cinterval(a) | lx_cinterval(b); }


// ------------------------- Others --------------------------------------

inline lx_cinterval & SetInf(lx_cinterval& a, const lx_complex& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Inf(a.re) = Re(b);
    Inf(a.im) = Im(b);

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("lx_cinterval & SetInf(lx_cinterval& a, const lx_complex& b)"));
   
    return a;
}

inline lx_cinterval & SetInf(lx_cinterval& a, const l_complex& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Inf(a.re) = Re(b);
    Inf(a.im) = Im(b);

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("lx_cinterval & SetInf(lx_cinterval& a, const l_complex& b)"));
   
    return a;
}

inline lx_cinterval & SetInf(lx_cinterval& a, const complex& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Inf(a.re) = Re(b);
    Inf(a.im) = Im(b);

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("lx_cinterval & SetInf(lx_cinterval& a, const complex& b)"));
   
    return a;
}

inline lx_cinterval & SetInf(lx_cinterval& a, const lx_real & b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Inf(a.re)=b;
    Inf(a.im)=0.0;

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline lx_cinterval & SetInf(lx_cinterval& a,const lx_real & b)"));
   
   return a;
}

inline lx_cinterval & SetInf(lx_cinterval& a, const l_real & b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Inf(a.re)=b;
    Inf(a.im)=0.0;

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetInf(lx_cinterval& a,const l_real & b)"));
   
   return a;
}

inline lx_cinterval & SetInf(lx_cinterval& a, const real & b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Inf(a.re)=b;
    Inf(a.im)=0.0;

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetInf(lx_cinterval& a,const real & b)"));
   
   return a;
}

inline lx_cinterval & SetSup(lx_cinterval& a, const lx_complex& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Sup(a.re) = Re(b);
    Sup(a.im) = Im(b);

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetSup(lx_cinterval& a, const lx_complex& b)"));
   
    return a;
}


inline lx_cinterval & SetSup(lx_cinterval& a, const l_complex& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Sup(a.re) = Re(b);
    Sup(a.im) = Im(b);

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetSup(lx_cinterval& a, const l_complex& b)"));
   
    return a;
}

inline lx_cinterval & SetSup(lx_cinterval& a, const complex& b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Sup(a.re) = Re(b);
    Sup(a.im) = Im(b);

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetSup(lx_cinterval& a, const complex& b)"));
   
    return a;
}

inline lx_cinterval & SetSup(lx_cinterval& a, const lx_real & b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Sup(a.re) = b;
    Sup(a.im) = 0.0;

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetSup(lx_cinterval& a, const lx_real& b)"));
   
    return a;
}

inline lx_cinterval & SetSup(lx_cinterval& a, const l_real & b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Sup(a.re) = b;
    Sup(a.im) = 0.0;

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetSup(lx_cinterval& a, const l_real& b)"));
   
    return a;
}

inline lx_cinterval & SetSup(lx_cinterval& a, const real & b) 
    throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
    Sup(a.re) = b;
    Sup(a.im) = 0.0;

    if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
	cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL(
		"inline lx_cinterval & SetSup(lx_cinterval& a, const real& b)"));
   
    return a;
}

inline bool IsEmpty(const lx_cinterval& a) throw()
{ return (IsEmpty(a.re) || IsEmpty(a.im)); }


// -----------------------------------------------------------------------
// ----------------------- Output ----------------------------------------
// -----------------------------------------------------------------------

inline std::ostream& operator << (std::ostream& s,const lx_cinterval& a) 
    throw()
// A value a of type lx_cinterval is written to the output channel.
{     
    s << '('          
    << a.re << ", "   
    << a.im       
    << ')';
    return s;
}

inline std::string & operator << (std::string &s,const lx_cinterval& a) 
    throw()
// The value of a variable a of type lx_cinterval is copied to a string s.
// s has the form:  {ex,li}
{  
    string str;
    s+='(';
    s << a.re; 
    s+=',';
    s << a.im; 
    s+=')';
    return s;
}

} // end namespace cxsc
