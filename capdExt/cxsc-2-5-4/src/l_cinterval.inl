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

/* CVS $Id: l_cinterval.inl,v 1.14 2014/01/30 17:23:46 cxsc Exp $ */

namespace cxsc {


// Inlined functions for l_cinterval.

// ---- implicit constructors  ------------------------------
      
inline l_cinterval::l_cinterval(const interval & a,const interval & b) throw()
      : re(a), im(b)
{
}

inline l_cinterval::l_cinterval(const l_interval & a, 
                                const l_interval & b) throw()
      : re(a), im(b)
{
}

inline l_cinterval::l_cinterval(const complex & a, const complex & b)  
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
   : re(Re(a),Re(b)),
     im(Im(a),Im(b))
{
   if(Inf(re)>Sup(re) || Inf(im)>Sup(im))
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("l_cinterval::l_cinterval(const complex & a,const complex & b)"));
}

inline l_cinterval::l_cinterval(const l_complex & a, const l_complex & b)  
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
   : re(Re(a),Re(b)),
     im(Im(a),Im(b))
{
   if(Inf(re)>Sup(re) || Inf(im)>Sup(im))
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("l_cinterval::l_cinterval(const complex & a,const complex & b)"));
}

// ---- explicit constructors  ------------------------------

inline l_cinterval::l_cinterval(const real & a) throw() : re(a,a), im(0,0) {}
inline l_cinterval::l_cinterval(const l_real & a) throw() : re(a,a), im(0,0) {}


inline l_cinterval::l_cinterval(const interval & a) 
                               throw() : re(a), im(0,0) {}
inline l_cinterval::l_cinterval(const l_interval & a) 
                               throw() : re(a), im(0,0) {}
inline l_cinterval::l_cinterval(const complex & a) 
                           throw() : re(Re(a),Re(a)),im(Im(a),Im(a)) {}
inline l_cinterval::l_cinterval(const l_complex & a) 
                           throw() : re(Re(a),Re(a)),im(Im(a),Im(a)) {} 
inline l_cinterval::l_cinterval(const cinterval & a) 
       throw() : re(Inf(Re(a)),Sup(Re(a))),im(Inf(Im(a)),Sup(Im(a))) {} 

     
// ---- assignments -----------------------------------------

inline l_cinterval & l_cinterval::operator =(const real & a) throw()
{
   re=a,im=0.0;
   return *this;   
}

inline l_cinterval & l_cinterval::operator =(const l_real & a) throw()
{
   re=a,im=0.0;
   return *this;   
}


inline l_cinterval & l_cinterval::operator =(const interval & a) throw()
{
   re=a,im=0.0;
   return *this;
}

inline l_cinterval & l_cinterval::operator =(const l_interval & a) throw()
{
   re=a,im=0.0;
   return *this;
}

inline l_cinterval & l_cinterval::operator =(const complex & a) throw()
{
   re=Re(a),im=Im(a);
   return *this;
}

inline l_cinterval & l_cinterval::operator =(const l_complex & a) throw()
{
   re=Re(a),im=Im(a);
   return *this;
}

inline l_cinterval & l_cinterval::operator =(const cinterval & a) throw()
{
     re = Re(a);	
     im = Im(a);	
     return *this;
}

inline l_cinterval & l_cinterval::operator =(const l_cinterval & a) throw()
{
     re = a.re;	
     im = a.im;	
     return *this;
}
 
inline l_cinterval & l_cinterval::operator =(const dotprecision & a) throw()
{
   return *this = l_cinterval(a);
}

inline l_cinterval & l_cinterval::operator =(const idotprecision & a) throw()
{
   return *this = l_cinterval(a);
}

inline l_cinterval & l_cinterval::operator =(const cdotprecision & a) throw()
{
   return *this = l_cinterval(a);
}

inline l_cinterval & l_cinterval::operator =(const cidotprecision & a) throw()
{
   return *this = l_cinterval(a);
}


// ---- Std.Operators ---------------------------------------

inline l_cinterval operator -(const l_cinterval & a) throw ()
{
   return l_cinterval(-a.re,-a.im);
}

inline l_cinterval operator +(const l_cinterval & a) throw ()
{
   return a;
}

inline bool operator! (const l_cinterval & a)  throw()
{
   return !a.re && !a.im;
}

inline l_cinterval operator +(const l_cinterval & a, 
                              const l_cinterval & b) throw()
{
   return l_cinterval(a.re + b.re, a.im + b.im);
}

inline l_cinterval operator -(const l_cinterval & a, 
                              const l_cinterval & b) throw()
{
   return l_cinterval(a.re - b.re, a.im - b.im);
}

inline l_cinterval operator &(const l_cinterval & a, const l_cinterval & b) 
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   l_cinterval tmp = a;
   l_real x,y;
   y = Inf(a.re);   x = Inf(b.re);
   if (x>y) y = x;  // y = max(Inf(a.re), Inf(b.re))
   SetInf(tmp.re, y);
   y = Inf(a.im);   x = Inf(b.im);
   if (x>y) y = x;  // y = max(Inf(a.im), Inf(b.im))
   SetInf(tmp.im, y);
   y = Sup(a.re);   x = Sup(b.re);
   if (x<y) y = x;  // y = min(Sup(a.re), Sup(b.re)) 
   SetSup(tmp.re, y);
   y = Sup(a.im);   x = Sup(b.im);
   if (x<y) y = x;  // y = min(Sup(a.im), Sup(b.im))
   SetSup(tmp.im, y);
   if (Inf(tmp.re) > Sup(tmp.re) || Inf(tmp.im) > Sup(tmp.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval operator &(const l_cinterval & a,const l_cinterval & b)"));
  return tmp;
}

inline l_cinterval operator |(const l_cinterval & a, const l_cinterval & b) 
                                                                    throw()
{
   l_cinterval tmp = a;
   l_real x,y;
   y = Inf(a.re);   x = Inf(b.re);
   if (x<y) y = x; // y = min(Inf(a.re), Inf(b.re))
   SetInf(tmp.re, y);
   y = Inf(a.im);   x = Inf(b.im);
   if (x<y) y = x; // y = min(Inf(a.im), Inf(b.im))
   SetInf(tmp.im, y);
   y = Sup(a.re);   x = Sup(b.re);
   if (x>y) y = x; // y = max(Sup(a.re), Sup(b.re)) 
   SetSup(tmp.re, y);
   y = Sup(a.im);   x = Sup(b.im);
   if (x>y) y = x; // y = max(Sup(a.im), Sup(b.im))
   SetSup(tmp.im, y);
   return tmp;
}

inline l_cinterval & operator +=(l_cinterval & a, const l_cinterval & b) 
                                             throw() { return a = a+b; }
inline l_cinterval & operator -=(l_cinterval & a, const l_cinterval & b) 
                                             throw() { return a = a-b; }
inline l_cinterval & operator *=(l_cinterval & a, const l_cinterval & b) 
                                             throw() { return a = a*b; }
inline l_cinterval & operator /=(l_cinterval & a, const l_cinterval & b) 
                                             throw() { return a = a/b; }
inline l_cinterval & operator |=(l_cinterval & a, const l_cinterval & b) 
                                             throw() { return a = a|b; }
inline l_cinterval & operator &=(l_cinterval & a, const l_cinterval & b) 
               throw(ERROR_CINTERVAL_EMPTY_INTERVAL) { return a = a&b; }

// LCI <--> R

inline l_cinterval operator +(const l_cinterval & a, const real & b) throw() 
                                              { return a + l_cinterval(b); }
inline l_cinterval operator +(const real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a, const real & b) throw() 
                                              { return a - l_cinterval(b); }
inline l_cinterval operator -(const real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a, const real & b) throw() 
                                     { return l_cinterval(a.re*b, a.im*b); }
inline l_cinterval operator *(const real & a, const l_cinterval & b) throw()
                                     { return l_cinterval(b.re*a, b.im*a); }
inline l_cinterval operator /(const l_cinterval & a, const real & b) throw()
                                     { return l_cinterval(a.re/b, a.im/b); }
inline l_cinterval operator /(const real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a, const real & b) throw() 
                                                { return a|l_cinterval(b); }
inline l_cinterval operator |(const real & a, const l_cinterval & b) throw() 
                                                { return l_cinterval(a)|b; }
inline l_cinterval operator &(const l_cinterval & a, const real & b) throw() 
                                              { return a & l_cinterval(b); }
inline l_cinterval operator &(const real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const real & b) throw() 
                                            { return a = a+l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const real & b) throw() 
                                           { return a = a -l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const real & b) throw()
                                                       { return a = a * b; }
inline l_cinterval & operator /=(l_cinterval & a, const real & b) throw()
                                                       { return a = a / b; }
inline l_cinterval & operator |=(l_cinterval & a, const real & b) throw() 
                                            { return a = a|l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const real & b) throw() 
                                            { return a = a&l_cinterval(b); }

// LCI <--> LR

inline l_cinterval operator +(const l_cinterval & a, const l_real & b) throw() 
                                              { return a + l_cinterval(b); }
inline l_cinterval operator +(const l_real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a, const l_real & b) throw() 
                                              { return a - l_cinterval(b); }
inline l_cinterval operator -(const l_real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a, const l_real & b) throw() 
                                     { return l_cinterval(a.re*b, a.im*b); }
inline l_cinterval operator *(const l_real & a, const l_cinterval & b) throw()
                                     { return l_cinterval(b.re*a, b.im*a); }
inline l_cinterval operator /(const l_cinterval & a, const l_real & b) throw()
                                     { return l_cinterval(a.re/b, a.im/b); }
inline l_cinterval operator /(const l_real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a, const l_real & b) throw() 
                                                { return a|l_cinterval(b); }
inline l_cinterval operator |(const l_real & a, const l_cinterval & b) throw() 
                                                { return l_cinterval(a)|b; }
inline l_cinterval operator &(const l_cinterval & a, const l_real & b) throw() 
                                              { return a & l_cinterval(b); }
inline l_cinterval operator &(const l_real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const l_real & b) throw() 
                                            { return a = a+l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const l_real & b) throw() 
                                           { return a = a -l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const l_real & b) throw()
                                                       { return a = a * b; }
inline l_cinterval & operator /=(l_cinterval & a, const l_real & b) throw()
                                                       { return a = a / b; }
inline l_cinterval & operator |=(l_cinterval & a, const l_real & b) throw() 
                                            { return a = a|l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const l_real & b) throw() 
                                            { return a = a&l_cinterval(b); }



// LCI <--> C

inline l_cinterval operator +(const l_cinterval & a, const complex & b) 
                                      throw() { return a + l_cinterval(b); }
inline l_cinterval operator +(const complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a, const complex & b) 
                                      throw() { return a - l_cinterval(b); }
inline l_cinterval operator -(const complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a, const complex & b) 
                                      throw() { return a * l_cinterval(b); }
inline l_cinterval operator *(const complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) * b; }
inline l_cinterval operator /(const l_cinterval & a, const complex & b) 
                                      throw() { return a / l_cinterval(b); }
inline l_cinterval operator /(const complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a, const complex & b) 
                                      throw() { return a | l_cinterval(b); }
inline l_cinterval operator |(const complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) | b; }
inline l_cinterval operator &(const l_cinterval & a, const complex & b) 
                                      throw() { return a & l_cinterval(b); }
inline l_cinterval operator &(const complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const complex & b) throw() 
                                          { return a = a + l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const complex & b) throw() 
                                          { return a = a - l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const complex & b) throw() 
                                          { return a = a * l_cinterval(b); }
inline l_cinterval & operator /=(l_cinterval & a, const complex & b) throw() 
                                          { return a = a / l_cinterval(b); }
inline l_cinterval & operator |=(l_cinterval & a, const complex & b) throw() 
                                          { return a = a | l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const complex & b) throw() 
                                          { return a = a & l_cinterval(b); }


// LCI <--> LC

inline l_cinterval operator +(const l_cinterval & a, const l_complex & b) 
                                      throw() { return a + l_cinterval(b); }
inline l_cinterval operator +(const l_complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a, const l_complex & b) 
                                      throw() { return a - l_cinterval(b); }
inline l_cinterval operator -(const l_complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a, const l_complex & b) 
                                      throw() { return a * l_cinterval(b); }
inline l_cinterval operator *(const l_complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) * b; }
inline l_cinterval operator /(const l_cinterval & a, const l_complex & b) 
                                      throw() { return a / l_cinterval(b); }
inline l_cinterval operator /(const l_complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a, const l_complex & b) 
                                      throw() { return a | l_cinterval(b); }
inline l_cinterval operator |(const l_complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) | b; }
inline l_cinterval operator &(const l_cinterval & a, const l_complex & b) 
                                      throw() { return a & l_cinterval(b); }
inline l_cinterval operator &(const l_complex & a, const l_cinterval & b) 
                                      throw() { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const l_complex & b) throw() 
                                          { return a = a + l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const l_complex & b) throw() 
                                          { return a = a - l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const l_complex & b) throw() 
                                          { return a = a * l_cinterval(b); }
inline l_cinterval & operator /=(l_cinterval & a, const l_complex & b) throw() 
                                          { return a = a / l_cinterval(b); }
inline l_cinterval & operator |=(l_cinterval & a, const l_complex & b) throw() 
                                          { return a = a | l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const l_complex & b) throw() 
                                          { return a = a & l_cinterval(b); }



// LCI <--> I

inline l_cinterval operator +(const l_cinterval & a, 
                              const interval & b) throw()
                                            { return a + l_cinterval(b); } 
inline l_cinterval operator +(const interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a,
                              const interval & b) throw() 
                                            { return a - l_cinterval(b); }
inline l_cinterval operator -(const interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a,
                              const interval & b) throw()
                                   { return l_cinterval(a.re*b, a.im*b); }
inline l_cinterval operator *(const interval & a,
                              const l_cinterval & b) throw()
                                   { return l_cinterval(b.re*a, b.im*a); }
inline l_cinterval operator /(const l_cinterval & a,
                              const interval & b) throw()
                                   { return l_cinterval(a.re/b, a.im/b); }
inline l_cinterval operator /(const interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a,
                              const interval & b) throw() 
                                            { return a | l_cinterval(b); }
inline l_cinterval operator |(const interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) | b; }
inline l_cinterval operator &(const l_cinterval & a,
                              const interval & b) throw() 
                                            { return a & l_cinterval(b); }
inline l_cinterval operator &(const interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const interval & b) 
                                throw() { return a = a + l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const interval & b) 
                                throw() { return a = a - l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const interval & b) 
                                throw() { return a = a * l_cinterval(b); }
inline l_cinterval & operator /=(l_cinterval & a, const interval & b) 
                                throw() { return a = a / l_cinterval(b); }
inline l_cinterval & operator |=(l_cinterval & a, const interval & b) 
                                throw() { return a = a | l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const interval & b) 
                                throw() { return a = a & l_cinterval(b); }

// LCI <--> LI

inline l_cinterval operator +(const l_cinterval & a, 
                              const l_interval & b) throw()
                                            { return a + l_cinterval(b); } 
inline l_cinterval operator +(const l_interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a,
                              const l_interval & b) throw() 
                                            { return a - l_cinterval(b); }
inline l_cinterval operator -(const l_interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a,
                              const l_interval & b) throw()
                                   { return l_cinterval(a.re*b, a.im*b); }
inline l_cinterval operator *(const l_interval & a,
                              const l_cinterval & b) throw()
                                   { return l_cinterval(b.re*a, b.im*a); }
inline l_cinterval operator /(const l_cinterval & a,
                              const l_interval & b) throw()
                                   { return l_cinterval(a.re/b, a.im/b); }
inline l_cinterval operator /(const l_interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a,
                              const l_interval & b) throw() 
                                            { return a | l_cinterval(b); }
inline l_cinterval operator |(const l_interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) | b; }
inline l_cinterval operator &(const l_cinterval & a,
                              const l_interval & b) throw() 
                                            { return a & l_cinterval(b); }
inline l_cinterval operator &(const l_interval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const l_interval & b) 
                                throw() { return a = a + l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const l_interval & b) 
                                throw() { return a = a - l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const l_interval & b) 
                                throw() { return a = a * l_cinterval(b); }
inline l_cinterval & operator /=(l_cinterval & a, const l_interval & b) 
                                throw() { return a = a / l_cinterval(b); }
inline l_cinterval & operator |=(l_cinterval & a, const l_interval & b) 
                                throw() { return a = a | l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const l_interval & b) 
                                throw() { return a = a & l_cinterval(b); }

// LCI <--> CI

inline l_cinterval operator +(const l_cinterval & a, 
                              const cinterval & b) throw()
                                            { return a + l_cinterval(b); } 
inline l_cinterval operator +(const cinterval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) + b; }
inline l_cinterval operator -(const l_cinterval & a,
                              const cinterval & b) throw() 
                                            { return a - l_cinterval(b); }
inline l_cinterval operator -(const cinterval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) - b; }
inline l_cinterval operator *(const l_cinterval & a,
                              const cinterval & b) throw()
                                            { return a * l_cinterval(b); }
inline l_cinterval operator *(const cinterval & a,
                              const l_cinterval & b) throw()
                                            { return l_cinterval(a) * b; }
inline l_cinterval operator /(const l_cinterval & a,
                              const cinterval & b) throw()
                                            { return a / l_cinterval(b); }
inline l_cinterval operator /(const cinterval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) / b; }
inline l_cinterval operator |(const l_cinterval & a,
                              const cinterval & b) throw() 
                                            { return a | l_cinterval(b); }
inline l_cinterval operator |(const cinterval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) | b; }
inline l_cinterval operator &(const l_cinterval & a,
                              const cinterval & b) throw() 
                                            { return a & l_cinterval(b); }
inline l_cinterval operator &(const cinterval & a,
                              const l_cinterval & b) throw() 
                                            { return l_cinterval(a) & b; }

inline l_cinterval & operator +=(l_cinterval & a, const cinterval & b) 
                                throw() { return a = a + l_cinterval(b); }
inline l_cinterval & operator -=(l_cinterval & a, const cinterval & b) 
                                throw() { return a = a - l_cinterval(b); }
inline l_cinterval & operator *=(l_cinterval & a, const cinterval & b) 
                                throw() { return a = a * l_cinterval(b); }
inline l_cinterval & operator /=(l_cinterval & a, const cinterval & b) 
                                throw() { return a = a / l_cinterval(b); }
inline l_cinterval & operator |=(l_cinterval & a, const cinterval & b) 
                                throw() { return a = a | l_cinterval(b); }
inline l_cinterval & operator &=(l_cinterval & a, const cinterval & b) 
                                throw() { return a = a & l_cinterval(b); }

// C-R

inline l_cinterval operator |(const l_complex & a, const real & b) throw() 
                               { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const real & a, const l_complex & b) throw() 
                               { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const complex & a, const l_real & b) throw() 
                               { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_real & a, const complex & b) throw() 
                               { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_complex & a, const l_real & b) throw() 
                                 { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_real & a, const l_complex & b) throw() 
                                 { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const cinterval & a, const l_real & b) throw() 
                                 { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_real & a, const cinterval & b) throw() 
                                 { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const cinterval & a, const l_complex & b) 
                    throw() { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_complex & a, const cinterval & b) 
                    throw() { return l_cinterval(a) | l_cinterval(b); }

// LC <--> I

inline l_cinterval operator +(const l_complex & a, const interval & b) throw() 
                                   { return l_cinterval(a) + l_cinterval(b); }
inline l_cinterval operator +(const interval & a, const l_complex & b) throw() 
                                   { return l_cinterval(a) + l_cinterval(b); }
inline l_cinterval operator -(const l_complex & a, const interval & b) throw() 
                                   { return l_cinterval(a) - l_cinterval(b); }
inline l_cinterval operator -(const interval & a, const l_complex & b) throw() 
                                   { return l_cinterval(a) - l_cinterval(b); }
inline l_cinterval operator *(const l_complex & a, const interval & b) throw()
                                                { return l_cinterval(a) * b; }
inline l_cinterval operator *(const interval & a, const l_complex & b) throw()
                                                { return l_cinterval(b) * a; }
inline l_cinterval operator /(const l_complex & a, const interval & b) throw()
                                                { return l_cinterval(a) / b; }
inline l_cinterval operator /(const interval & a, const l_complex & b) throw()
                                   { return l_cinterval(a) / l_cinterval(b); }
inline l_cinterval operator |(const l_complex & a, const interval & b) throw()
                                   { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const interval & a, const l_complex & b) throw()
                                   { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator &(const l_complex & a, const interval & b) throw()
                                   { return l_cinterval(a) & l_cinterval(b); }
inline l_cinterval operator &(const interval & a, const l_complex & b) throw()
                                   { return l_cinterval(a) & l_cinterval(b); }

// C <--> LI

inline l_cinterval operator +(const complex & a, const l_interval & b) throw() 
                                   { return l_cinterval(a) + l_cinterval(b); }
inline l_cinterval operator +(const l_interval & a, const complex & b) throw() 
                                   { return l_cinterval(a) + l_cinterval(b); }
inline l_cinterval operator -(const complex & a, const l_interval & b) throw() 
                                   { return l_cinterval(a) - l_cinterval(b); }
inline l_cinterval operator -(const l_interval & a, const complex & b) throw() 
                                   { return l_cinterval(a) - l_cinterval(b); }
inline l_cinterval operator *(const complex & a, const l_interval & b) throw()
                                                { return l_cinterval(a) * b; }
inline l_cinterval operator *(const l_interval & a, const complex & b) throw()
                                                { return l_cinterval(b) * a; }
inline l_cinterval operator /(const complex & a, const l_interval & b) throw()
                                                { return l_cinterval(a) / b; }
inline l_cinterval operator /(const l_interval & a, const complex & b) throw()
                                   { return l_cinterval(a) / l_cinterval(b); }
inline l_cinterval operator |(const complex & a, const l_interval & b) throw()
                                   { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_interval & a, const complex & b) throw()
                                   { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator &(const complex & a, const l_interval & b) throw()
                                   { return l_cinterval(a) & l_cinterval(b); }
inline l_cinterval operator &(const l_interval & a, const complex & b) throw()
                                   { return l_cinterval(a) & l_cinterval(b); }

// LC <--> LI

inline l_cinterval operator +(const l_complex & a, const l_interval & b) 
                           throw() { return l_cinterval(a) + l_cinterval(b); }
inline l_cinterval operator +(const l_interval & a, const l_complex & b) 
                           throw() { return l_cinterval(a) + l_cinterval(b); }
inline l_cinterval operator -(const l_complex & a, const l_interval & b) 
                           throw() { return l_cinterval(a) - l_cinterval(b); }
inline l_cinterval operator -(const l_interval & a, const l_complex & b) 
                           throw() { return l_cinterval(a) - l_cinterval(b); }
inline l_cinterval operator *(const l_complex & a, const l_interval & b) 
                                        throw() { return l_cinterval(a) * b; }
inline l_cinterval operator *(const l_interval & a, const l_complex & b) 
                                        throw() { return l_cinterval(b) * a; }
inline l_cinterval operator /(const l_complex & a, const l_interval & b) 
                                        throw() { return l_cinterval(a) / b; }
inline l_cinterval operator /(const l_interval & a, const l_complex & b) 
                           throw() { return l_cinterval(a) / l_cinterval(b); }
inline l_cinterval operator |(const l_complex & a, const l_interval & b) 
                           throw() { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_interval & a, const l_complex & b) 
                           throw() { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator &(const l_complex & a, const l_interval & b) 
                           throw() { return l_cinterval(a) & l_cinterval(b); }
inline l_cinterval operator &(const l_interval & a, const l_complex & b) 
                           throw() { return l_cinterval(a) & l_cinterval(b); }


// LC <--> C

inline l_cinterval operator |(const l_complex & a, const complex & b) 
                           throw() { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const complex & a, const l_complex & b) 
                           throw() { return l_cinterval(a) | l_cinterval(b); }
inline l_cinterval operator |(const l_complex & a, const l_complex & b) 
                           throw() { return l_cinterval(a) | l_cinterval(b); }


// ---- Comp.Operat.  ---------------------------------------

inline bool operator== (const l_cinterval & a, const l_cinterval & b) throw()
{
   return a.re==b.re && a.im==b.im;
}
inline bool operator!= (const l_cinterval & a, const l_cinterval & b) throw()
{
   return a.re!=b.re || a.im!=b.im;
}


// LCI-R

inline bool operator== (const l_cinterval & a, const real & b) throw() 
                                              { return a == l_cinterval(b); }
inline bool operator== (const real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const real & b) throw() 
                                              { return a != l_cinterval(b); }
inline bool operator!= (const real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) != b; }

// LCI-LR

inline bool operator== (const l_cinterval & a, const l_real & b) throw() 
                                              { return a == l_cinterval(b); }
inline bool operator== (const l_real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const l_real & b) throw() 
                                              { return a != l_cinterval(b); }
inline bool operator!= (const l_real & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) != b; }

// LCI <--> I

inline bool operator== (const l_cinterval & a, const interval & b) throw() 
                                              { return a == l_cinterval(b); }
inline bool operator== (const interval & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const interval & b) throw() 
                                              { return a != l_cinterval(b); }
inline bool operator!= (const interval & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) != b; }

// LCI <--> LI

inline bool operator== (const l_cinterval & a, const l_interval & b) throw() 
                                             { return a == l_cinterval(b); }
inline bool operator== (const l_interval & a, const l_cinterval & b) throw() 
                                             { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const l_interval & b) throw() 
                                             { return a != l_cinterval(b); }
inline bool operator!= (const l_interval & a, const l_cinterval & b) throw() 
                                             { return l_cinterval(a) != b; }

// LCI <--> C

inline bool operator== (const l_cinterval & a, const complex & b) throw() 
                                            { return a == l_cinterval(b); }
inline bool operator== (const complex & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const complex & b) throw() 
                                            { return a != l_cinterval(b); }
inline bool operator!= (const complex & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) != b; }

// LCI <--> LC

inline bool operator== (const l_cinterval & a, const l_complex & b) throw() 
                                            { return a == l_cinterval(b); }
inline bool operator== (const l_complex & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const l_complex & b) throw() 
                                            { return a != l_cinterval(b); }
inline bool operator!= (const l_complex & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) != b; }

// LCI <--> CI

inline bool operator== (const l_cinterval & a, const cinterval & b) throw() 
                                            { return a == l_cinterval(b); }
inline bool operator== (const cinterval & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) == b; }
inline bool operator!= (const l_cinterval & a, const cinterval & b) throw() 
                                            { return a != l_cinterval(b); }
inline bool operator!= (const cinterval & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) != b; }


// ---- Set Operators ----
inline bool operator  <(const l_cinterval & a, const l_cinterval & b) throw()
{
   if (Inf(a.re) <= Inf(b.re) || Sup(a.re) >= Sup(b.re)) 
      return false;
   if (Inf(a.im) <= Inf(b.im) || Sup(a.im) >= Sup(b.im)) 
      return false;
      
   return true; 
}

inline bool operator  >(const l_cinterval & a, const l_cinterval & b) throw() 
                                                            { return b < a; }

inline bool operator <=(const l_cinterval & a, const l_cinterval & b) throw()
{
   if (Inf(a.re) < Inf(b.re) || Sup(a.re) > Sup(b.re)) 
      return false;
   if (Inf(a.im) < Inf(b.im) || Sup(a.im) > Sup(b.im)) 
      return false;
      
   return true; 
}

inline bool operator >= (const l_cinterval & a, const l_cinterval & b) throw()
                                                            { return b <= a; }

// lCI <--> R

inline bool operator  <(const real & a, const l_cinterval & b) throw() 
                                        { return l_cinterval(a) < b; }
inline bool operator  >(const real & a, const l_cinterval & b) throw() 
                                        { return l_cinterval(a) > b; }
inline bool operator <=(const real & a, const l_cinterval & b) throw() 
                                       { return l_cinterval(a) <= b; }
inline bool operator >=(const real & a, const l_cinterval & b) throw() 
                                       { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const real & b) throw() 
                                        { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const real & b) throw() 
                                        { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const real & b) throw() 
                                       { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const real & b) throw() 
                                       { return a >= l_cinterval(b); }

// lCI <--> LR

inline bool operator  <(const l_real & a, const l_cinterval & b) throw() 
                                          { return l_cinterval(a) < b; }
inline bool operator  >(const l_real & a, const l_cinterval & b) throw() 
                                          { return l_cinterval(a) > b; }
inline bool operator <=(const l_real & a, const l_cinterval & b) throw() 
                                         { return l_cinterval(a) <= b; }
inline bool operator >=(const l_real & a, const l_cinterval & b) throw() 
                                         { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const l_real & b) throw() 
                                          { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const l_real & b) throw() 
                                          { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const l_real & b) throw() 
                                         { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const l_real & b) throw() 
                                         { return a >= l_cinterval(b); }

// LCI <--> I

inline bool operator  <(const interval & a, const l_cinterval & b) throw() 
                                            { return l_cinterval(a) < b; }
inline bool operator  >(const interval & a, const l_cinterval & b) throw()
                                            { return l_cinterval(a) > b; }
inline bool operator <=(const interval & a, const l_cinterval & b) throw()
                                           { return l_cinterval(a) <= b; }
inline bool operator >=(const interval & a, const l_cinterval & b) throw() 
                                           { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const interval & b) throw() 
                                            { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const interval & b) throw() 
                                            { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const interval & b) throw() 
                                           { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const interval & b) throw() 
                                           { return a >= l_cinterval(b); }

// LCI <--> LI

inline bool operator  <(const l_interval & a, const l_cinterval & b) throw() 
                                              { return l_cinterval(a) < b; }
inline bool operator  >(const l_interval & a, const l_cinterval & b) throw()
                                              { return l_cinterval(a) > b; }
inline bool operator <=(const l_interval & a, const l_cinterval & b) throw()
                                             { return l_cinterval(a) <= b; }
inline bool operator >=(const l_interval & a, const l_cinterval & b) throw() 
                                             { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const l_interval & b) throw() 
                                              { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const l_interval & b) throw() 
                                              { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const l_interval & b) throw() 
                                             { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const l_interval & b) throw() 
                                             { return a >= l_cinterval(b); }

// LCI <--> C

inline bool operator  <(const complex & a, const l_cinterval & b) throw() 
                                           { return l_cinterval(a) < b; }
inline bool operator  >(const complex & a, const l_cinterval & b) throw()
                                           { return l_cinterval(a) > b; }
inline bool operator <=(const complex & a, const l_cinterval & b) throw()
                                          { return l_cinterval(a) <= b; }
inline bool operator >=(const complex & a, const l_cinterval & b) throw()
                                          { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const complex & b) throw()
                                           { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const complex & b) throw() 
                                           { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const complex & b) throw() 
                                          { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const complex & b) throw()
                                          { return a >= l_cinterval(b); }

// LCI <--> LC

inline bool operator  <(const l_complex & a, const l_cinterval & b) throw() 
                                             { return l_cinterval(a) < b; }
inline bool operator  >(const l_complex & a, const l_cinterval & b) throw()
                                             { return l_cinterval(a) > b; }
inline bool operator <=(const l_complex & a, const l_cinterval & b) throw()
                                            { return l_cinterval(a) <= b; }
inline bool operator >=(const l_complex & a, const l_cinterval & b) throw()
                                            { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const l_complex & b) throw()
                                             { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const l_complex & b) throw() 
                                             { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const l_complex & b) throw() 
                                            { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const l_complex & b) throw()
                                            { return a >= l_cinterval(b); }

// LCI <--> CI

inline bool operator  <(const cinterval & a, const l_cinterval & b) throw() 
                                             { return l_cinterval(a) < b; }
inline bool operator  >(const cinterval & a, const l_cinterval & b) throw()
                                             { return l_cinterval(a) > b; }
inline bool operator <=(const cinterval & a, const l_cinterval & b) throw()
                                            { return l_cinterval(a) <= b; }
inline bool operator >=(const cinterval & a, const l_cinterval & b) throw()
                                            { return l_cinterval(a) >= b; }

inline bool operator  <(const l_cinterval & a, const cinterval & b) throw()
                                             { return a < l_cinterval(b); }
inline bool operator  >(const l_cinterval & a, const cinterval & b) throw() 
                                             { return a > l_cinterval(b); }
inline bool operator <=(const l_cinterval & a, const cinterval & b) throw() 
                                            { return a <= l_cinterval(b); }
inline bool operator >=(const l_cinterval & a, const cinterval & b) throw()
                                            { return a >= l_cinterval(b); }


// ---- Others   -------------------------------------------

inline l_complex    Inf(const l_cinterval & a) throw() 
                                 { return l_complex(Inf(a.re),Inf(a.im)); }
inline l_complex    Sup(const l_cinterval & a) throw() 
                                 { return l_complex(Sup(a.re),Sup(a.im)); }

inline l_cinterval & SetInf(l_cinterval & a, const complex & b) 
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Inf(a.re) = Re(b);
   Inf(a.im) = Im(b);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetInf(l_cinterval & a, const complex & b)"));
   
   return a;
}

inline l_cinterval & SetSup(l_cinterval & a, const complex & b) 
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Sup(a.re)=Re(b);
   Sup(a.im)=Im(b);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetSup(l_cinterval & a, const complex & b)"));
   
   return a;
}

inline l_cinterval & SetInf(l_cinterval & a, const l_complex & b) 
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Inf(a.re) = Re(b);
   Inf(a.im) = Im(b);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetInf(l_cinterval & a, const l_complex & b)"));
   
   return a;
}

inline l_cinterval & SetSup(l_cinterval & a, const l_complex & b) 
                                      throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Sup(a.re)=Re(b);
   Sup(a.im)=Im(b);

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetSup(l_cinterval & a, const l_complex & b)"));
   
   return a;
}

inline l_cinterval & SetInf(l_cinterval & a, const real & b) 
                         throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Inf(a.re)=b;
   Inf(a.im)=0.0;

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetSup(l_cinterval & a, const real & b)"));
   
   return a;
}

inline l_cinterval & SetSup(l_cinterval & a, const real & b) 
                         throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Sup(a.re)=b;
   Sup(a.im)=0.0;

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetSup(l_cinterval & a, const real & b)"));
   
   return a;
}

inline l_cinterval & SetInf(l_cinterval & a, const l_real & b) 
                         throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Inf(a.re)=b;
   Inf(a.im)=0.0;

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetSup(l_cinterval & a, const l_real & b)"));
   
   return a;
}

inline l_cinterval & SetSup(l_cinterval & a, const l_real & b) 
                         throw(ERROR_CINTERVAL_EMPTY_INTERVAL)
{
   Sup(a.re)=b;
   Sup(a.im)=0.0;

   if (Inf(a.re) > Sup(a.re) || Inf(a.im) > Sup(a.im)) 
      cxscthrow(ERROR_CINTERVAL_EMPTY_INTERVAL("inline l_cinterval & SetSup(l_cinterval & a, const l_real & b)"));
   
   return a;
}


inline l_cinterval & UncheckedSetInf(l_cinterval & a, const complex & b)
                                                                   throw()
{
   Inf(a.re)=Re(b);
   Inf(a.im)=Im(b);
   return a;
}

inline l_cinterval & UncheckedSetInf(l_cinterval & a, const real & b) 
                                                                throw()
{
   Inf(a.re)=b;
   Inf(a.im)=0.0;
   return a;
}

inline l_cinterval & UncheckedSetSup(l_cinterval & a, const complex & b) 
                                                                   throw()
{
   Sup(a.re)=Re(b);
   Sup(a.im)=Im(b);
   return a;
}

inline l_cinterval & UncheckedSetSup(l_cinterval & a, const real & b) 
                                                                throw()
{
   Sup(a.re)=b;
   Sup(a.im)=0.0;
   return a;
}

inline l_cinterval & UncheckedSetInf(l_cinterval & a, const l_complex & b)
                                                                   throw()
{
   Inf(a.re)=Re(b);
   Inf(a.im)=Im(b);
   return a;
}

inline l_cinterval & UncheckedSetInf(l_cinterval & a, const l_real & b) 
                                                                throw()
{
   Inf(a.re)=b;
   Inf(a.im)=0.0;
   return a;
}

inline l_cinterval & UncheckedSetSup(l_cinterval & a, const l_complex & b) 
                                                                   throw()
{
   Sup(a.re)=Re(b);
   Sup(a.im)=Im(b);
   return a;
}

inline l_cinterval & UncheckedSetSup(l_cinterval & a, const l_real & b) 
                                                                throw()
{
   Sup(a.re)=b;
   Sup(a.im)=0.0;
   return a;
}


inline l_cinterval conj(const l_cinterval & a) throw() 
{ return l_cinterval(a.re,-a.im); }

inline l_complex mid(const l_cinterval &a) throw() 
{ return l_complex(mid(a.re), mid(a.im)); }

inline l_complex diam(const l_cinterval &a) throw()
{ return l_complex(diam(a.re),diam(a.im)); }

inline l_cinterval adjust(const l_cinterval & a) throw() 
{ 
  return l_cinterval(adjust(Re(a)),adjust(Im(a)));	
	
// return l_cinterval(a.re,-a.im); 
}

inline void times2pown(l_cinterval& x, const int& n) throw() 
// Blomquist, 08.03.07;
{
    if ( n<-1074 || n>1023 ) 
    { std::cerr << "Error in:  " 
           << "times2pown(l_cinterval& x, const int& n): " << std::endl
           << " -1074 <= n <= +1023 not fulfilled" << std::endl; exit(0); 
    }
    l_interval u(Re(x)),v(Im(x));
    times2pown(u,n);
    times2pown(v,n);
    x = l_cinterval(u,v);
}

inline void Times2pown(l_cinterval& x, const int& n) throw() 
// Blomquist, 28.03.07;
{

    l_interval u(Re(x)),v(Im(x));
    Times2pown(u,n);
    Times2pown(v,n);
    x = l_cinterval(u,v);
}


} // namespace cxsc




