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

/* CVS $Id: cinterval.inl,v 1.20 2014/01/30 17:23:44 cxsc Exp $ */

namespace cxsc {

// Inlined functions for cinterval.

// ---- implicit constructors  ------------------------------
      
inline cinterval::cinterval(const interval & a,const interval & b)
      : re(a), im(b)
{
}

inline cinterval::cinterval(const complex & a,const complex & b) 
   : re(Re(a),Re(b)),
     im(Im(a),Im(b))
{
   if(Inf(re)>Sup(re) || Inf(im)>Sup(im))
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("cinterval::cinterval(const complex & a,const complex & b)"));
}

// ---- explicit constructors  ------------------------------

inline cinterval::cinterval(const real     & a) : re(a,a), im(0,0) {}
inline cinterval::cinterval(const interval & a) : re(a), im(0,0) {}
inline cinterval::cinterval(const complex  & a) : re(Re(a),Re(a)),im(Im(a),Im(a)) {} 
      
// ---- assignments -----------------------------------------

inline cinterval & cinterval::operator =(const real & a)
{
   re=a,im=0.0;
   return *this;   
}

inline cinterval & cinterval::operator =(const interval & a)
{
   re=a,im=0.0;
   return *this;
}

inline cinterval & cinterval::operator =(const complex & a)
{
   re=Re(a),im=Im(a);
   return *this;
}

inline cinterval & cinterval::operator =(const cinterval & a)
{
   re=a.re;
   im=a.im;
   return *this;
}

inline cinterval & cinterval::operator =(const dotprecision & a)
{
   return *this=cinterval(a);
}

inline cinterval & cinterval::operator =(const idotprecision & a)
{
   return *this=cinterval(a);
}

inline cinterval & cinterval::operator =(const cdotprecision & a)
{
   return *this=cinterval(a);
}

inline cinterval & cinterval::operator =(const cidotprecision & a)
{
   return *this=cinterval(a);
}

// ---- compatiblility typecasts ----------------------------

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const real & a)
*/
inline cinterval _cinterval(const real & a)
{
   return cinterval(interval(a,a), interval(0.0,0.0));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a)
*/
inline cinterval _cinterval(const complex & a)
{
   return cinterval(a,a);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const interval & a)
*/
inline cinterval _cinterval(const interval & a)
{
   return cinterval(a,_interval(0.0,0.0));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const dotprecision &)
*/
inline cinterval _cinterval(const dotprecision & a) { return cinterval(a); }

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const cdotprecision &)
*/
inline cinterval _cinterval(const cdotprecision & a) { return cinterval(a); }

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const idotprecision &)
*/
inline cinterval _cinterval(const idotprecision & a) { return cinterval(a); }

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const cidotprecision &)
*/
inline cinterval _cinterval(const cidotprecision & a) { return cinterval(a); }

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a,const complex & b)
*/
inline cinterval _cinterval(const complex & a,const complex & b)
{
   return cinterval(interval(Re(a),Re(b)),interval(Im(a),Im(b)));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a,const complex & b)
*/
inline cinterval _cinterval(const real & a,const complex & b)
{
   return cinterval(complex(a),b);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a,const complex & b)
*/
inline cinterval _cinterval(const complex & a,const real & b)
{
   return cinterval(a,complex(b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const interval & a,const interval &b)
*/
inline cinterval _cinterval(const interval & a,const interval & b)
{
   return cinterval(a,b);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const interval & a,const interval &b)
*/
inline cinterval _cinterval(const real & a,const interval & b)
{
   return cinterval(interval(a),b);   
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const interval & a,const interval &b)
*/
inline cinterval _cinterval(const interval & a,const real & b)
{
   return cinterval(a,interval(b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a,const complex & b)
*/
inline cinterval _unchecked_cinterval(const complex & a,const complex & b)
{
   cinterval tmp;
   UncheckedSetInf(tmp,a);
   UncheckedSetSup(tmp,b);
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a,const complex & b)
*/
inline cinterval _unchecked_cinterval(const real & a,const complex & b)
{
   cinterval tmp;
   UncheckedSetInf(tmp,_complex(a));
   UncheckedSetSup(tmp,b);
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cinterval::cinterval(const complex & a,const complex & b)
*/
inline cinterval _unchecked_cinterval(const complex & a,const real & b)
{
   cinterval tmp;
   UncheckedSetInf(tmp,a);
   UncheckedSetSup(tmp,_complex(b));
   return tmp;
}

// ---- Std.Operators ---------------------------------------

inline cinterval operator -(const cinterval & a)
{
   return cinterval(-a.re,-a.im);
}

inline cinterval operator +(const cinterval & a)
{
   return a;
}

inline cinterval operator +(const cinterval & a,const cinterval & b)
{
   return cinterval(a.re+b.re,a.im+b.im);
}

inline cinterval operator -(const cinterval & a,const cinterval & b)
{
   return cinterval(a.re-b.re,a.im-b.im);
}

inline cinterval operator &(const cinterval & a,const cinterval & b)
{
   cinterval tmp = a;
   SetInf(tmp.re, max(Inf(a.re), Inf(b.re)));
   SetInf(tmp.im, max(Inf(a.im), Inf(b.im)));
   SetSup(tmp.re, min(Sup(a.re), Sup(b.re)));
   SetSup(tmp.im, min(Sup(a.im), Sup(b.im)));
   if (Inf(tmp.re) > Sup(tmp.re) || Inf(tmp.im) > Sup(tmp.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline cinterval operator &(const cinterval & a,const cinterval & b)"));
  return tmp;
}

inline cinterval operator |(const cinterval & a,const cinterval & b)
{
   cinterval tmp = a;
   SetInf(tmp.re, min(Inf(a.re), Inf(b.re)));
   SetInf(tmp.im, min(Inf(a.im), Inf(b.im)));
   SetSup(tmp.re, max(Sup(a.re), Sup(b.re)));
   SetSup(tmp.im, max(Sup(a.im), Sup(b.im)));
   return tmp;
}

inline cinterval & operator +=(cinterval & a, const cinterval & b) { return a=a+b; }
inline cinterval & operator -=(cinterval & a, const cinterval & b) { return a=a-b; }
inline cinterval & operator *=(cinterval & a, const cinterval & b) { return a=a*b; }
inline cinterval & operator /=(cinterval & a, const cinterval & b) { return a=a/b; }
inline cinterval & operator |=(cinterval & a, const cinterval & b) { return a=a|b; }
inline cinterval & operator &=(cinterval & a, const cinterval & b) { return a=a&b; }

// CI-R

inline cinterval operator +(const cinterval & a,const real & b) { return a+_cinterval(b); }
inline cinterval operator +(const real & a,const cinterval & b) { return _cinterval(a)+b; }
inline cinterval operator -(const cinterval & a,const real & b) { return a-_cinterval(b); }
inline cinterval operator -(const real & a,const cinterval & b) { return _cinterval(a)-b; }
inline cinterval operator *(const cinterval & a,const real & b) { return a*_cinterval(b); }
inline cinterval operator *(const real & a,const cinterval & b)
{ // return _cinterval(a)*b;
     return cinterval(b.re*a, b.im*a); // Blomquist 07.11.02;
}
inline cinterval operator /(const cinterval & a,const real & b)
{ // return a/_cinterval(b);
     return cinterval(a.re/b, a.im/b); // Blomquist 07.11.02;
}
inline cinterval operator /(const real & a,const cinterval & b) { return _cinterval(a)/b; }
inline cinterval operator |(const cinterval & a,const real & b) { return a|_cinterval(b); }
inline cinterval operator |(const real & a,const cinterval & b) { return _cinterval(a)|b; }
inline cinterval operator &(const cinterval & a,const real & b) { return a&_cinterval(b); }
inline cinterval operator &(const real & a,const cinterval & b) { return _cinterval(a)&b; }

inline cinterval & operator +=(cinterval & a, const real & b) { return a=a+_cinterval(b); }
inline cinterval & operator -=(cinterval & a, const real & b) { return a=a-_cinterval(b); }
inline cinterval & operator *=(cinterval & a, const real & b)
{ // return a=a*_cinterval(b);
     return a = a * b;  // Blomquist 07.11.02;
}
inline cinterval & operator /=(cinterval & a, const real & b)
{ // return a=a/_cinterval(b);
     return a = a / b;  // Blomquist 07.11.02;
}
inline cinterval & operator |=(cinterval & a, const real & b) { return a=a|_cinterval(b); }
inline cinterval & operator &=(cinterval & a, const real & b) { return a=a&_cinterval(b); }

// CI-C

inline cinterval operator +(const cinterval & a,const complex & b) { return a+_cinterval(b); }
inline cinterval operator +(const complex & a,const cinterval & b) { return _cinterval(a)+b; }
inline cinterval operator -(const cinterval & a,const complex & b) { return a-_cinterval(b); }
inline cinterval operator -(const complex & a,const cinterval & b) { return _cinterval(a)-b; }
inline cinterval operator *(const cinterval & a,const complex & b) { return a*_cinterval(b); }
inline cinterval operator *(const complex & a,const cinterval & b) { return _cinterval(a)*b; }
inline cinterval operator /(const cinterval & a,const complex & b) { return a/_cinterval(b); }
inline cinterval operator /(const complex & a,const cinterval & b) { return _cinterval(a)/b; }
inline cinterval operator |(const cinterval & a,const complex & b) { return a|_cinterval(b); }
inline cinterval operator |(const complex & a,const cinterval & b) { return _cinterval(a)|b; }
inline cinterval operator &(const cinterval & a,const complex & b) { return a&_cinterval(b); }
inline cinterval operator &(const complex & a,const cinterval & b) { return _cinterval(a)&b; }

inline cinterval & operator +=(cinterval & a, const complex & b) { return a=a+_cinterval(b); }
inline cinterval & operator -=(cinterval & a, const complex & b) { return a=a-_cinterval(b); }
inline cinterval & operator *=(cinterval & a, const complex & b) { return a=a*_cinterval(b); }
inline cinterval & operator /=(cinterval & a, const complex & b) { return a=a/_cinterval(b); }
inline cinterval & operator |=(cinterval & a, const complex & b) { return a=a|_cinterval(b); }
inline cinterval & operator &=(cinterval & a, const complex & b) { return a=a&_cinterval(b); }

// CI-I

inline cinterval operator +(const cinterval & a,const interval & b) { return a+_cinterval(b); }
inline cinterval operator +(const interval & a,const cinterval & b) { return _cinterval(a)+b; }
inline cinterval operator -(const cinterval & a,const interval & b) { return a-_cinterval(b); }
inline cinterval operator -(const interval & a,const cinterval & b) { return _cinterval(a)-b; }
inline cinterval operator *(const cinterval & a,const interval & b)
{ // return a*_cinterval(b);
     return cinterval(a.re*b,a.im*b);  // Blomquist, 07.11.02;
}
inline cinterval operator *(const interval & a,const cinterval & b)
{ // return _cinterval(a)*b;
     return cinterval(b.re*a,b.im*a);
}
inline cinterval operator /(const cinterval & a,const interval & b)
{ // return a/_cinterval(b);
     return cinterval(a.re/b,a.im/b);
}
inline cinterval operator /(const interval & a,const cinterval & b) { return _cinterval(a)/b; }
inline cinterval operator |(const cinterval & a,const interval & b) { return a|_cinterval(b); }
inline cinterval operator |(const interval & a,const cinterval & b) { return _cinterval(a)|b; }
inline cinterval operator &(const cinterval & a,const interval & b) { return a&_cinterval(b); }
inline cinterval operator &(const interval & a,const cinterval & b) { return _cinterval(a)&b; }

inline cinterval & operator +=(cinterval & a, const interval & b) { return a=a+_cinterval(b); }
inline cinterval & operator -=(cinterval & a, const interval & b) { return a=a-_cinterval(b); }
inline cinterval & operator *=(cinterval & a, const interval & b) { return a=a*_cinterval(b); }
inline cinterval & operator /=(cinterval & a, const interval & b) { return a=a/_cinterval(b); }
inline cinterval & operator |=(cinterval & a, const interval & b) { return a=a|_cinterval(b); }
inline cinterval & operator &=(cinterval & a, const interval & b) { return a=a&_cinterval(b); }

// C-R

inline cinterval operator |(const complex & a,const real & b) { return _cinterval(a)|_cinterval(b); }
inline cinterval operator |(const real & a,const complex & b) { return _cinterval(a)|_cinterval(b); }

// C-I

inline cinterval operator +(const complex & a,const interval & b) { return _cinterval(a)+_cinterval(b); }
inline cinterval operator +(const interval & a,const complex & b) { return _cinterval(a)+_cinterval(b); }
inline cinterval operator -(const complex & a,const interval & b) { return _cinterval(a)-_cinterval(b); }
inline cinterval operator -(const interval & a,const complex & b) { return _cinterval(a)-_cinterval(b); }
inline cinterval operator *(const complex & a,const interval & b)
{ // return _cinterval(a)*_cinterval(b);
     return _cinterval(a)*b;  // Blomquist, 07.11.02;
}
inline cinterval operator *(const interval & a,const complex & b)
{ // return _cinterval(a)*_cinterval(b);
     return _cinterval(b) * a;  // Blomquist, 07.11.02;
}
inline cinterval operator /(const complex & a,const interval & b)
{ // return _cinterval(a)/_cinterval(b);
     return _cinterval(a) / b;  // Blomquist, 07.11.02;
}
inline cinterval operator /(const interval & a,const complex & b) { return _cinterval(a)/_cinterval(b); }
inline cinterval operator |(const complex & a,const interval & b) { return _cinterval(a)|_cinterval(b); }
inline cinterval operator |(const interval & a,const complex & b) { return _cinterval(a)|_cinterval(b); }
inline cinterval operator &(const complex & a,const interval & b) { return _cinterval(a)&_cinterval(b); }
inline cinterval operator &(const interval & a,const complex & b) { return _cinterval(a)&_cinterval(b); }
      
// C-C

inline cinterval operator |(const complex & a,const complex & b) { return _cinterval(a)|_cinterval(b); }

// ---- Comp.Operat.  ---------------------------------------
inline bool operator!  (const cinterval & a) 
{
   return !a.re && !a.im;
}

inline bool operator== (const cinterval & a, const cinterval & b)
{
   return a.re==b.re && a.im==b.im;
}

inline bool operator!= (const cinterval & a, const cinterval & b)
{
   return a.re!=b.re || a.im!=b.im;
}

// CI-R

inline bool operator== (const cinterval & a, const real & b) { return a==_cinterval(b); }
inline bool operator== (const real & a, const cinterval & b) { return _cinterval(a)==b; }
inline bool operator!= (const cinterval & a, const real & b) { return a!=_cinterval(b); }
inline bool operator!= (const real & a, const cinterval & b) { return _cinterval(a)!=b; }

// CI-C

inline bool operator== (const cinterval & a, const complex & b) { return a==_cinterval(b); }
inline bool operator== (const complex & a, const cinterval & b) { return _cinterval(a)==b; }
inline bool operator!= (const cinterval & a, const complex & b) { return a!=_cinterval(b); }
inline bool operator!= (const complex & a, const cinterval & b) { return _cinterval(a)!=b; }

// CI-I

inline bool operator== (const cinterval & a, const interval & b) { return a==_cinterval(b); }
inline bool operator== (const interval & a, const cinterval & b) { return _cinterval(a)==b; }
inline bool operator!= (const cinterval & a, const interval & b) { return a!=_cinterval(b); }
inline bool operator!= (const interval & a, const cinterval & b) { return _cinterval(a)!=b; }

// ---- Set Operators ----
inline bool operator  <(const cinterval & a,const cinterval & b)
{
   if (Inf(a.re) <= Inf(b.re) || Sup(a.re) >= Sup(b.re)) 
      return false;
   if (Inf(a.im) <= Inf(b.im) || Sup(a.im) >= Sup(b.im)) 
      return false;
      
   return true; 
}

inline bool operator  >(const cinterval & a,const cinterval & b) { return b<a; }

inline bool operator <=(const cinterval & a,const cinterval & b)
{
   if (Inf(a.re) < Inf(b.re) || Sup(a.re) > Sup(b.re)) 
      return false;
   if (Inf(a.im) < Inf(b.im) || Sup(a.im) > Sup(b.im)) 
      return false;
      
   return true; 
}

inline bool operator >=(const cinterval & a,const cinterval & b) { return b<=a; }

// CI-R

inline bool operator  <(const real & a,const cinterval & b) { return _cinterval(a)<b; }
inline bool operator  >(const real & a,const cinterval & b) { return _cinterval(a)>b; }
inline bool operator <=(const real & a,const cinterval & b) { return _cinterval(a)<=b; }
inline bool operator >=(const real & a,const cinterval & b) { return _cinterval(a)>=b; }

inline bool operator  <(const cinterval & a,const real & b) { return a<_cinterval(b); }
inline bool operator  >(const cinterval & a,const real & b) { return a>_cinterval(b); }
inline bool operator <=(const cinterval & a,const real & b) { return a<=_cinterval(b); }
inline bool operator >=(const cinterval & a,const real & b) { return a>=_cinterval(b); }

// CI-C

inline bool operator  <(const complex & a,const cinterval & b) { return _cinterval(a)<b; }
inline bool operator  >(const complex & a,const cinterval & b) { return _cinterval(a)>b; }
inline bool operator <=(const complex & a,const cinterval & b) { return _cinterval(a)<=b; }
inline bool operator >=(const complex & a,const cinterval & b) { return _cinterval(a)>=b; }

inline bool operator  <(const cinterval & a,const complex & b) { return a<_cinterval(b); }
inline bool operator  >(const cinterval & a,const complex & b) { return a>_cinterval(b); }
inline bool operator <=(const cinterval & a,const complex & b) { return a<=_cinterval(b); }
inline bool operator >=(const cinterval & a,const complex & b) { return a>=_cinterval(b); }

// CI-C

inline bool operator  <(const interval & a,const cinterval & b) { return _cinterval(a)<b; }
inline bool operator  >(const interval & a,const cinterval & b) { return _cinterval(a)>b; }
inline bool operator <=(const interval & a,const cinterval & b) { return _cinterval(a)<=b; }
inline bool operator >=(const interval & a,const cinterval & b) { return _cinterval(a)>=b; }

inline bool operator  <(const cinterval & a,const interval & b) { return a<_cinterval(b); }
inline bool operator  >(const cinterval & a,const interval & b) { return a>_cinterval(b); }
inline bool operator <=(const cinterval & a,const interval & b) { return a<=_cinterval(b); }
inline bool operator >=(const cinterval & a,const interval & b) { return a>=_cinterval(b); }

// ---- Others   -------------------------------------------
inline complex    Inf(const cinterval & a) { return _complex(Inf(a.re),Inf(a.im)); }
inline complex    Sup(const cinterval & a) { return _complex(Sup(a.re),Sup(a.im)); }

inline cinterval & SetInf(cinterval & a,const complex & b)
{
   Inf(a.re)=Re(b);
   Inf(a.im)=Im(b);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline cinterval & SetInf(cinterval & a,const complex & b)"));
   
   return a;
}

inline cinterval & SetSup(cinterval & a,const complex & b)
{
   Sup(a.re)=Re(b);
   Sup(a.im)=Im(b);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline cinterval & SetSup(cinterval & a,const complex & b)"));
   
   return a;
}

inline cinterval & SetInf(cinterval & a,const real & b)
{
   Inf(a.re)=b;
   Inf(a.im)=0.0;

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline cinterval & SetInf(cinterval & a,const real & b)"));
   
   return a;
}

inline cinterval & SetSup(cinterval & a,const real & b)
{
   Sup(a.re)=b;
   Sup(a.im)=0.0;

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline cinterval & SetSup(cinterval & a,const real & b)"));
   
   return a;
}

inline cinterval & UncheckedSetInf(cinterval & a,const complex & b)
{
   Inf(a.re)=Re(b);
   Inf(a.im)=Im(b);
   return a;
}

inline cinterval & UncheckedSetInf(cinterval & a,const real & b)
{
   Inf(a.re)=b;
   Inf(a.im)=0.0;
   return a;
}

inline cinterval & UncheckedSetSup(cinterval & a,const complex & b)
{
   Sup(a.re)=Re(b);
   Sup(a.im)=Im(b);
   return a;
}

inline cinterval & UncheckedSetSup(cinterval & a,const real & b)
{
   Sup(a.re)=b;
   Sup(a.im)=0.0;
   return a;
}

inline interval & Re(cinterval & a)       { return a.re; }
inline interval   Re(const cinterval & a) { return a.re; }
inline interval & Im(cinterval & a)       { return a.im; }
inline interval   Im(const cinterval & a) { return a.im; }
      
inline cinterval & SetRe(cinterval & a,const interval & b) { a.re=b; return a; }
inline cinterval & SetIm(cinterval & a,const interval & b) { a.im=b; return a; } 
inline cinterval & SetRe(cinterval & a,const real     & b) { a.re=b; return a; }
inline cinterval & SetIm(cinterval & a,const real     & b) { a.im=b; return a; } 

inline real InfRe(const cinterval &a) { return Inf(a.re); }
inline real InfIm(const cinterval &a) { return Inf(a.im); }
inline real SupRe(const cinterval &a) { return Sup(a.re); }
inline real SupIm(const cinterval &a) { return Sup(a.im); }
      
inline real & InfRe(cinterval &a) { return Inf(a.re); }
inline real & InfIm(cinterval &a) { return Inf(a.im); }
inline real & SupRe(cinterval &a) { return Sup(a.re); }
inline real & SupIm(cinterval &a) { return Sup(a.im); }










inline cinterval conj(const cinterval & a) { return cinterval(a.re,-a.im); }

inline complex mid(const cinterval &a) { return complex(mid(a.re),mid(a.im)); }
inline complex diam(const cinterval &a){ return complex(diam(a.re),diam(a.im)); }

cinterval mult_operator(const cinterval & a,const cinterval & b);
cinterval div_operator(const cinterval & a,const cinterval & b);

inline cinterval operator *(const cinterval & a,const cinterval & b)
{
#ifdef CXSC_FAST_COMPLEX_OPERATIONS
  return cinterval(Re(a)*Re(b)-Im(a)*Im(b), Re(a)*Im(b)+Im(a)*Re(b));
#else
  return mult_operator(a,b);
#endif
} 

inline cinterval operator / (const cinterval & a, const cinterval & b)
{
    if (0.0 <= Re(b) && 0.0 <= Im(b) ) {
      cxscthrow(DIV_BY_ZERO("cinterval operator / (const cinterval&, const cinterval&)"));
      return a; // dummy result
    }      
#ifdef CXSC_FAST_COMPLEX_OPERATIONS
    return a * (1.0 / b);
#else
    return div_operator(a,b);
#endif
}


} // namespace cxsc

