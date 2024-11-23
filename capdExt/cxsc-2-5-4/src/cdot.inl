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

/* CVS $Id: cdot.inl,v 1.31 2014/01/30 17:23:43 cxsc Exp $ */

namespace cxsc {

// ---- Konstruktoren ----

inline cdotprecision::cdotprecision(const dotprecision &a,const dotprecision &b)
                        : re(a), im(b), k(0) { re.set_k(0); im.set_k(0); }

inline cdotprecision::cdotprecision(const real &a,const real &b)
                        : re(a), im(b), k(0) {}

inline cdotprecision::cdotprecision(const cdotprecision &a)
                        : re(a.re), im(a.im), k(a.k) {}

inline cdotprecision::cdotprecision(const l_real &a, const l_real &b)
                        : re(a), im(b), k(0) {}

inline cdotprecision::cdotprecision(const l_complex &a)
                        : re(Re(_cdotprecision(a))), im(Im(_cdotprecision(a))), k(0){}

inline cdotprecision::cdotprecision(const real         &r)   
                        : re(r), im(0), k(0) {}
			
inline cdotprecision::cdotprecision(const complex      &c)
                        : re(Re(c)),im(Im(c)), k(0) {}
			
inline cdotprecision::cdotprecision(const dotprecision &r)
                        : re(r), im(0), k(0) { re.set_k(0); im.set_k(0); }
			
inline cdotprecision::cdotprecision(const l_real &r)
                        : re(r), im(0), k(0) {}




inline cdotprecision& cdotprecision::operator= (const real & a)
{
   re=a;im=0; return *this;
}

inline cdotprecision& cdotprecision::operator= (const complex & a) 
{
   re=Re(a),im=Im(a); return *this;
}
inline cdotprecision& cdotprecision::operator= (const dotprecision & a)
{
   re=a;im=0; return *this;
}
inline cdotprecision& cdotprecision::operator= (const cdotprecision & a)
{
   re=a.re,im=a.im; return *this;
}
inline cdotprecision& cdotprecision::operator= (const l_real & a) 
{
   re=a;im=0; return *this;
}


// ---- Typwandlungen ----

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const dotprecision &r)
*/
inline cdotprecision _cdotprecision(const dotprecision& a) 
{
  return cdotprecision (a);
}
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const real &r)
*/
inline cdotprecision _cdotprecision(const real& a) 
{
  return cdotprecision (a);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const l_real &r)
*/
inline cdotprecision _cdotprecision(const l_real& a) 
{ 
  return cdotprecision (a);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const complex &c)
*/
inline cdotprecision _cdotprecision(const complex& a) 
{
  return cdotprecision (a);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const dotprecision&, const dotprecision&)
*/
inline cdotprecision _cdotprecision(const dotprecision& a, const dotprecision& b) 
{
  return cdotprecision (a,b);
}
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const real &, const real &)
*/
inline cdotprecision _cdotprecision(const real& a, const real& b) 
{
  return cdotprecision (a,b);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const l_real &, const l_real &)
*/
inline cdotprecision _cdotprecision(const l_real& a, const l_real& b)
{
  return cdotprecision (a,b);
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cdotprecision::cdotprecision(const l_complex &)
*/
inline cdotprecision _cdotprecision(const l_complex& lc)
{
   return cdotprecision(lc);
}

// ---- Standardfunkt ---- (arithmetische Operatoren)

inline cdotprecision operator-(const cdotprecision &a) throw() 
{ 
  return cdotprecision (-a.re, -a.im); 
}
inline cdotprecision operator+(const cdotprecision &a) throw() 
{ 
  return a; 
}

inline cdotprecision operator+(const cdotprecision &a,const cdotprecision &b) throw() 
{ 
  return cdotprecision(a.re+b.re,a.im+b.im); 
}
inline cdotprecision operator-(const cdotprecision &a,const cdotprecision &b) throw() 
{ 
  return cdotprecision(a.re-b.re,a.im-b.im); 
}

inline cdotprecision operator +(const cdotprecision &a,const complex &b) throw() { return cdotprecision(a.re+Re(b),a.im+Im(b)); }
inline cdotprecision operator +(const complex &b,const cdotprecision &a) throw() { return cdotprecision(a.re+Re(b),a.im+Im(b)); } 
inline cdotprecision operator -(const cdotprecision &a,const complex &b) throw() { return cdotprecision(a.re-Re(b),a.im-Im(b)); } 
inline cdotprecision operator -(const complex &a,const cdotprecision &b) throw() { return cdotprecision(Re(a)-b.re,Im(a)-b.im); }

inline cdotprecision operator +(const cdotprecision &a,const dotprecision &b) throw() { return cdotprecision(a.re+b,a.im); }
inline cdotprecision operator +(const dotprecision &b,const cdotprecision &a) throw() { return cdotprecision(a.re+b,a.im); } 
inline cdotprecision operator -(const cdotprecision &a,const dotprecision &b) throw() { return cdotprecision(a.re-b,a.im); } 
inline cdotprecision operator -(const dotprecision &a,const cdotprecision &b) throw() { return cdotprecision(a-b.re,-b.im); }

inline cdotprecision operator +(const cdotprecision &a,const real &b) throw() { return cdotprecision(a.re+b,a.im); }
inline cdotprecision operator +(const real &b,const cdotprecision &a) throw() { return cdotprecision(a.re+b,a.im); } 
inline cdotprecision operator -(const cdotprecision &a,const real &b) throw() { return cdotprecision(a.re-b,a.im); } 
inline cdotprecision operator -(const real &a,const cdotprecision &b) throw() { return cdotprecision(a-b.re,b.im); }

inline cdotprecision operator +(const cdotprecision &a, const l_real &b) throw() 
{ 
  return cdotprecision(a.re+b,a.im);
}
inline cdotprecision operator +(const l_real &b, const cdotprecision &a) throw() 
{ 
  return cdotprecision(a.re+b,a.im);
}
inline cdotprecision operator -(const cdotprecision &a, const l_real &b) throw() 
{ 
  return cdotprecision(a.re-b,a.im);
}
inline cdotprecision operator -(const l_real &b, const cdotprecision &a) throw() 
{ 
  return cdotprecision(b-a.re,-a.im); 
}

inline cdotprecision & operator +=(cdotprecision &a,const cdotprecision &b) throw() { a.re+=b.re;a.im+=b.im; return a;}
inline cdotprecision & operator +=(cdotprecision &a,const dotprecision &b) throw() { a.re+=b; return a;}
inline cdotprecision & operator -=(cdotprecision &a,const cdotprecision &b) throw() { a.re-=b.re;a.im-=b.im; return a;}
inline cdotprecision & operator -=(cdotprecision &a,const dotprecision &b) throw() { a.re-=b; return a;}
inline cdotprecision & operator +=(cdotprecision &a,const complex &b) throw() 
{
   a.re+=Re(b);
   a.im+=Im(b);
   return a;
}
inline cdotprecision & operator -=(cdotprecision &a,const complex &b) throw()
{
   a.re-=Re(b);
   a.im-=Im(b);
   return a;
}
inline cdotprecision & operator +=(cdotprecision &a,const real &b) throw()
{
   a.re+=b;
   return a;
}
inline cdotprecision & operator -=(cdotprecision &a,const real &b) throw()
{
   a.re-=b;
   return a;
}
inline cdotprecision & operator +=(cdotprecision &a,const l_real &b) throw()
{ // Blomquist 17.09.02.
   a.re+=b;
   return a;
}
inline cdotprecision & operator -=(cdotprecision &a,const l_real &b) throw()
{ // Blomquist 17.09.02.
   a.re-=b;
   return a;
}
inline cdotprecision & operator += (cdotprecision &cd, const l_complex &lc) throw() 
{
   cd = cd + cdotprecision(lc); return cd;
}
inline cdotprecision & operator -= (cdotprecision &cd, const l_complex &lc) throw() 
{
   cd = cd - cdotprecision(lc); return cd;
}
      
// --- Vergleichsoperationen ----
inline bool operator ==(const cdotprecision &a,const cdotprecision &b) throw() {   return(a.re==b.re && a.im==b.im); }
inline bool operator !=(const cdotprecision &a,const cdotprecision &b) throw() {   return(a.re!=b.re || a.im!=b.im); }
inline bool operator ==(const dotprecision &r,const cdotprecision &a)     throw() {   return(r==a.re && !a.im); }
inline bool operator !=(const dotprecision &r,const cdotprecision &a)     throw() {   return(r!=a.re || !!a.im); }
inline bool operator ==(const cdotprecision &a,const dotprecision &r)     throw() {   return(r==a.re && !a.im); }
inline bool operator !=(const cdotprecision &a,const dotprecision &r)     throw() {   return(r!=a.re || !!a.im); }
inline bool operator ==(const complex &c,const cdotprecision &a)     throw() {   return(Re(c)==a.re && Im(c)==a.im); }
inline bool operator !=(const complex &c,const cdotprecision &a)     throw() {   return(Re(c)!=a.re || Im(c)!=a.im); }
inline bool operator ==(const cdotprecision &a,const complex &c)     throw() {   return(Re(c)==a.re && Im(c)==a.im); }
inline bool operator !=(const cdotprecision &a,const complex &c)     throw() {   return(Re(c)!=a.re || Im(c)!=a.im); }
inline bool operator ==(const real &c,const cdotprecision &a)     throw() {   return(c==a.re && !a.im); }
inline bool operator !=(const real &c,const cdotprecision &a)     throw() {   return(c!=a.re || !!a.im); }
inline bool operator ==(const cdotprecision &a,const real &c)     throw() {   return(c==a.re && !a.im); }
inline bool operator !=(const cdotprecision &a,const real &c)     throw() {   return(c!=a.re || !!a.im); }

inline bool operator ==(const l_real &c,const cdotprecision &a) throw() 
{ 
  return(c==a.re && !a.im);
}

inline bool operator !=(const l_real &c,const cdotprecision &a) throw()
{ 
  return(c!=a.re || !!a.im);
}

inline bool operator ==(const cdotprecision &a,const l_real &c) throw()
{ 
  return(c==a.re && !a.im);
}

inline bool operator !=(const cdotprecision &a,const l_real &c) throw()
{ 
  return(c!=a.re || !!a.im);
}

// ----- Funktionen -----
inline dotprecision & Re(cdotprecision& a)
{ 
  return a.re; 
}

inline dotprecision & Im(cdotprecision& a) throw() 
{ 
  return a.im;
}

inline const dotprecision & Re(const cdotprecision& a)
{ 
  return a.re; 
}

inline const dotprecision & Im(const cdotprecision& a) throw()
{ 
  return a.im; 
}

inline cdotprecision& SetRe (cdotprecision& a, const dotprecision& b)  throw()
{                         // ggf. exception
   a.re=b; 
   return a;
}

inline cdotprecision& SetIm (cdotprecision& a, const dotprecision& b) throw()
{
   a.im=b;
   return a;
}

inline cdotprecision conj(const cdotprecision& a) throw()
{
   return cdotprecision(a.re,-a.im);
}

inline bool operator !(const cdotprecision &a) throw() 
{ 
   return !a.re && !a.im;
}  

inline void accumulate  (cdotprecision & a, const complex & b, const real & c) throw()
{ 
  accumulate(a,b,complex(c));
}

inline void accumulate  (cdotprecision & a, const real & b, const complex & c) throw()
{ 
  accumulate(a,complex(b),c);
}

inline void accumulate  (cdotprecision & a, const real & b, const real & c) throw()
{ 
  accumulate(a,complex(b),complex(c)); 
}

} // namespace cxsc

