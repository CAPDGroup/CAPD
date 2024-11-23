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

/* CVS $Id: cdot.hpp,v 1.35 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_CDOT_HPP_INCLUDED
#define _CXSC_CDOT_HPP_INCLUDED

#include <iostream>
#include <string>
#include "dot.hpp"
#include "complex.hpp"
#include "l_real.hpp" 
#include "l_complex.hpp"
#include "cinterval.hpp"

namespace cxsc {

class l_complex;  // forward declaration

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- class cdotprecision (declaration)                                 ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! The Data Type cdotprecision
/*!
The data types dotprecision, idotprecision, cdotprecision and cidotprecision
are based on the scalar data types real, interval, complex, and cinterval, 
respectively. Variables of these data types permit the exact representation of
products of two arbitrary numbers of the corresponding scalar base type and 
the exact summation of an arbitrary number of such products in a dotprecision
accumulator, i.e. in a fixed-point format of suitable size.

Since C-XSC Version 2.3.0, the precision for dot products computed with the dotprecision types can be choosen by the user. The default is precision k=0, which means maximum precision (the behaviour of older C-XSC Versions). For k=1, pure floating point operations are used and an error bound is computed using manipulation of the rounding mode of the processor. For k>=2, the so called DotK algorithm is used, simulating higher precision computations and also computing an error bound. When calling the rnd() function, the error bound will be added to the result interval. The resulting intervals will be wider for lower k, but computations will be significantly faster than with maximum precision (k=0).

\sa cxsc::dotprecision
*/
class cdotprecision
{  
   private:
      // ---- Datenelemente ---------------------------------------
      dotprecision re,im;
      int k;

   public:
      // ---- Constructors  ---------------------------------------
      //! Constructor of class cdotprecision
      cdotprecision() : k(0) {}
      //! Constructor of class cdotprecision
      inline cdotprecision(const dotprecision&, const dotprecision&);
      //! Constructor of class cdotprecision
      inline cdotprecision(const real &, const real &);
      //! Constructor of class cdotprecision
      inline cdotprecision(const cdotprecision&);
      //! Constructor of class cdotprecision
      inline cdotprecision(const l_real &, const l_real &);
      //! Constructor of class cdotprecision
      inline cdotprecision(const l_complex &);
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const real         &r);
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const complex      &c);
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const dotprecision &r);
      //! Constructor of class cdotprecision
      explicit inline cdotprecision(const l_real &r);

      //! Get currently set precision for computation of dot products
      inline int get_k() const { return k; }
      //! Set precision for computation of dot products
      inline void set_k(unsigned int i) { k=i; re.set_k(i); im.set_k(i); }
      //! Get currently set precision for computation of dot products
      inline int get_dotprec() const { return k; }
      //! Set precision for computation of dot products
      inline void set_dotprec(unsigned int i) { k=i; re.set_k(i); im.set_k(i); }

      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const real & a);
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const complex & a);
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const dotprecision & a);
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const cdotprecision & a);
      //! Implementation of standard assigning operator
      inline cdotprecision& operator= (const l_real & a);

      // ---- Destruktor    ----
      // ~cdotprecision() {} unnoetig

      // ---- Typwandlungen ----
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const dotprecision&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const real&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const l_real&); 
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const complex&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const dotprecision&, const dotprecision&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const real&,const real&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const l_real&, const l_real&); 
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cdotprecision _cdotprecision(const l_complex&);

      // ---- Ausgabefunkt. ---------------------------------------
      //! Implementation of standard input method
      friend std::istream& operator >> (std::istream& s, cdotprecision& a)       throw();
      //! Implementation of standard output method
      friend std::ostream& operator << (std::ostream& s, const cdotprecision& a) throw();
      //! Implementation of standard input method
      friend std::string&  operator >> (std::string&  s, cdotprecision& a)       throw();
      //! Implementation of standard output method
      friend std::string&  operator << (std::string&  s, const cdotprecision& a) throw();
      //! Implementation of standard input method
      friend void          operator >> (const std::string &s,cdotprecision& a)   throw();
      //! Implementation of standard input method
      friend void          operator >> (const char *s       ,cdotprecision& a)   throw();

      // ---- Standardfunkt ---- (arithmetische Operatoren)
      //! Implementation of standard algebraic negative sign operation
      friend inline cdotprecision operator -(const cdotprecision &) throw();
      //! Implementation of standard algebraic positive sign operation
      friend inline cdotprecision operator +(const cdotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const cdotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const complex &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const complex &,const cdotprecision &) throw();
      
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const dotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const dotprecision &,const cdotprecision &) throw();
      
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const real &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const real &,const cdotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cdotprecision operator +(const l_real &, const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cdotprecision operator -(const l_real &, const cdotprecision &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const complex &) throw();
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const real &) throw();
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const l_real &) throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cdotprecision & operator +=(cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cdotprecision & operator -=(cdotprecision &,const dotprecision &) throw();

      // ---- Vergleichsop. ----

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const dotprecision &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const dotprecision &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const dotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const dotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const complex &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const complex &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const complex &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const complex &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const real &,const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const real &,const cdotprecision &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &,const real &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &,const real &) throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard equality operation
      friend inline bool operator ==(const l_real &, const cdotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cdotprecision &, const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const l_real &, const cdotprecision &) throw();

      //! Implementation of standard negation operation
      friend inline bool operator !(const cdotprecision&) throw();
      
      // ---- Funktionen    ----
      
      //! Returns the real part of the complex dotprecision value
      friend inline dotprecision & Re(cdotprecision& a);
      //! Returns the imaginary part of the complex dotprecision value
      friend inline dotprecision & Im(cdotprecision& a)  throw();
      //! Returns the real part of the complex dotprecision value
      friend inline const dotprecision & Re(const cdotprecision& a);
      //! Returns the imaginary part of the complex dotprecision value
      friend inline const dotprecision & Im(const cdotprecision& a)  throw();
      //! Sets the real part of a complex dotprecision value
      friend inline cdotprecision& SetRe (cdotprecision& a, const dotprecision& b)  throw();
      //! Sets the imaginary part of a complex dotprecision value
      friend inline cdotprecision& SetIm (cdotprecision& a, const dotprecision& b) throw();
      
      //! Returns the conjugated complex dotprecision value
      friend inline cdotprecision conj(const cdotprecision &a) throw();

      //! Converting the exact complex dotprecision value with one rounding into a complex value
      friend void rnd(const cdotprecision &,complex &,rndtype ) throw();
      //! Converting the exact complex dotprecision value with one rounding into the nearest lower and upper complex value
      friend void rnd(const cdotprecision &,complex &,complex &) throw();
      //! Converting the exact complex dotprecision value into an interval with the nearest lower and upper bound
      friend void rnd (const cdotprecision&, cinterval&) throw();
      //! Converting the exact complex dotprecision value with one rounding into a complex value
      friend complex rnd(const cdotprecision &,rndtype) throw();
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend        void accumulate  (cdotprecision&, const complex&, const complex&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cdotprecision&, const complex&, const real&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cdotprecision&, const real&, const complex&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cdotprecision&, const real&, const real&) throw();
};

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- friend functions of class cdotprecision (not inline)              ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

inline cdotprecision _cdotprecision(const l_complex&);

std::istream& operator >> (std::istream& s, cdotprecision& a)       throw();
std::ostream& operator << (std::ostream& s, const cdotprecision& a) throw();
std::string&  operator >> (std::string&  s, cdotprecision& a)       throw();
std::string&  operator << (std::string&  s, const cdotprecision& a) throw();
void          operator >> (const std::string &s,cdotprecision& a)   throw();
void          operator >> (const char *s       ,cdotprecision& a)   throw();

void rnd(const cdotprecision &,complex &,rndtype = RND_NEXT) throw();
void rnd(const cdotprecision &,complex &,complex &) throw();
void rnd(const cdotprecision&, cinterval&) throw();
complex rnd(const cdotprecision &,rndtype = RND_NEXT) throw();

void accumulate  (cdotprecision&, const complex&, const complex&) throw();


// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- global functions associated with class cdotprecision              ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! Implementation of standard algebraic addition and allocation operation
inline cdotprecision & operator += (cdotprecision &cd, const l_complex &lc) throw(); 
inline cdotprecision & operator -= (cdotprecision &cd, const l_complex &lc) throw(); 

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ----  global CDotprecision Akku's                                      ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//#define MAXCDOTAKKU     (MAXDOTAKKU / 2)
//extern cdotprecision cdotakku[MAXCDOTAKKU];
  
// ---------------------------------------------------------------------------

} // namespace cxsc 


// ---------------------------------------------------------------------------
// ----                                                                   ----
// ----  include definitions of inline functions                          ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

#include "cdot.inl"

#endif // _CXSC_CDOT_HPP_INCLUDED

