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
/* CVS $Id: l_complex.hpp,v 1.19 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_L_COMPLEX_HPP_INCLUDED
#define _CXSC_L_COMPLEX_HPP_INCLUDED

#include <iostream>
#include <string>

#include "except.hpp"
#include "complex.hpp"
#include "l_real.hpp"
#include "l_rmath.hpp"
#include "cdot.hpp"

namespace cxsc {

//! The Multiple-Precision Data Type l_complex
/*!
The multiple-precision data type l_complex is a variant of the scalar type complex, which provides support for longer numbers, thus increasing the accuracy of the data type.

\sa complex
*/
class l_complex
{
   private:
    // ------------- Datenelemente -------------------------------------------
    l_real re, im;

   public:
    //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
    l_complex _l_complex(const cdotprecision &);
    // ------------- Constructors --------------------------------------------
    //! Constructor of class l_complex
    l_complex(void)  throw() {}
    //! Constructor of class l_complex
    l_complex(const l_real& a, const l_real& b) throw() : re(a), im(b) { }
    //! Constructor of class l_complex
    l_complex(const real& a, const real& b) throw() : re(a), im(b) { }

    //! Implementation of standard assigning operator
    l_complex & operator = (const l_real & lr) throw() 
        { re = lr; im = 0; return *this; }
    //! Implementation of standard assigning operator
    l_complex & operator = (const real & r) throw() 
      { re = r; im = 0; return *this; } 
    //! Implementation of standard assigning operator
    l_complex & operator = (const complex & c) throw() 
        { re = Re(c); im = Im(c); return *this; } 
    //! Implementation of standard assigning operator
    l_complex & operator = (const dotprecision & d) throw() 
        { re = d; im = 0.0; return *this; }
    //! Implementation of standard assigning operator
    l_complex & operator = (const cdotprecision & cd) throw() 
        { re = Re(_l_complex(cd)); im = Im(_l_complex(cd)); return *this; }
	 
	 //! Implementation of standard assigning operator
	 l_complex & operator = (const lx_complex &) throw();

    // ------------- Type-Casts ----------------------------------------------
    //! Constructor of class l_complex
    explicit inline l_complex(const l_real  &r) throw() : re(r), im(0.0) { }
    //! Constructor of class l_complex
    explicit inline l_complex(const   real  &r) throw() : re(r), im(0.0) { }
    //! Constructor of class l_complex
    explicit inline l_complex(const complex &r) throw() : re(Re(r)), im(Im(r))
                                                                         { }
    //! Constructor of class l_complex
    explicit inline l_complex(const dotprecision &d) throw()
                  : re(d), im(0.0) { }
    //! Constructor of class l_complex
    explicit inline l_complex(const cdotprecision &cd) throw()
	          : re(Re(_l_complex(cd))), im(Im(_l_complex(cd))) { }

//    friend inline l_complex _l_complex(const l_real &a) throw()
//	{ return l_complex(a); }
//    friend inline l_complex _l_complex(const l_real &a, const l_real &b)
//        throw()   { return l_complex(a,b); }
//    friend inline l_complex _l_complex(const real &a) throw()
//	{ return l_complex(a); }
//    friend inline l_complex _l_complex(const real &a, const real &b)
//        throw()   { return l_complex(a,b); }
//    friend inline l_complex _l_complex(const complex &c)
//        throw()   { return l_complex(c); }
//    friend inline l_complex _l_complex(const dotprecision &d)
//        throw()   { return l_complex(d); }
    /*!
    \deprecated use standard contructors for typecasting

    \sa cxsc::l_complex::l_complex(const cdotprecision &cd)
    */
    friend inline l_complex _l_complex(const cdotprecision &cd)
        throw()   { return l_complex(cd); }

    // ----------------------------------------------------------------------
    //! Returns the precision of the long datatype value
    friend int StagPrec(const l_complex&) throw();

// ------------- Arithmetic Operators ---------------------------------------
    // ----------------------------------------------------------------------
    // ------------- Unary Operators ----------------------------------------

    //! Implementation of standard algebraic negative sign operation
    friend l_complex operator-(const l_complex&);
    //! Implementation of standard algebraic positive sign operation
    friend l_complex operator+(const l_complex&);

    // ------------- Binary Operators ---------------------------------------
    // ----------------- l_complex +/- l_complex ----------------------------
    //! Implementation of standard algebraic addition operation
    friend inline l_complex operator +
                   (const l_complex &, const l_complex &) throw();

    //! Implementation of standard algebraic subtraction operation
    friend inline l_complex operator -
                   (const l_complex &, const l_complex &) throw();
    
    // ----------------- l_complex + complex --------------------------------
    //! Implementation of standard algebraic positive sign operation
    inline l_complex operator + (const complex &) const throw();

    //! Implementation of standard algebraic addition operation
    friend inline l_complex operator +
                   (const complex &, const l_complex &) throw();

    // ---------------- l_complex + real ------------------------------------
    //! Implementation of standard algebraic positive sign operation
    inline l_complex operator + (const real &) const throw();

    //! Implementation of standard algebraic addition operation
    friend inline l_complex operator +
                   (const real &, const l_complex &) throw();

    // ---------------- l_complex + l_real ----------------------------------
    //! Implementation of standard algebraic positive sign operation
    inline l_complex operator + (const l_real &) const throw();

    //! Implementation of standard algebraic addition operation
    friend inline l_complex operator +
                   (const l_real &, const l_complex &) throw();

    // ---------------- l_complex - l_real ----------------------------------
    //! Implementation of standard algebraic negative sign operation
    inline l_complex operator - (const l_real &) const throw();

    //! Implementation of standard algebraic subtraction operation
    friend inline l_complex operator -
                   (const l_real &, const l_complex &) throw();

    // ----------------- l_complex - complex -------------------------------
    //! Implementation of standard algebraic negative sign operation
    inline l_complex operator - (const complex &) const throw();

    //! Implementation of standard algebraic subtraction operation
    friend inline l_complex operator -
                   (const complex &, const l_complex &) throw();

    // ---------------- l_complex - real ------------------------------------
    //! Implementation of standard algebraic negative sign operation
    inline l_complex operator - (const real &) const throw();

    //! Implementation of standard algebraic subtraction operation
    friend inline l_complex operator -
                   (const real &, const l_complex &) throw();

    // ---------------- l_complex + cdotprecision ---------------------------
    //! Implementation of standard algebraic addition operation
    friend cdotprecision operator+
	   (const l_complex &, const cdotprecision &) throw();

    //! Implementation of standard algebraic addition operation
    friend cdotprecision operator+ 
      (const cdotprecision &, const l_complex &) throw();

    // ---------------- l_complex - cdotprecision ---------------------------
    //! Implementation of standard algebraic subtraction operation
    friend cdotprecision operator- 
      (const l_complex &, const cdotprecision &) throw();

    //! Implementation of standard algebraic subtraction operation
    friend cdotprecision operator- 
      (const cdotprecision &, const l_complex &) throw(); 

    // ---------------- l_complex + dotprecision ----------------------------
    //! Implementation of standard algebraic addition operation
    friend cdotprecision operator+ 
      (const l_complex &, const dotprecision &) throw();

    //! Implementation of standard algebraic addition operation
    friend cdotprecision operator+ 
      (const dotprecision &, const l_complex &) throw();

    // ---------------- l_complex - dotprecision ----------------------------
    //! Implementation of standard algebraic subtraction operation
    friend cdotprecision operator- 
      (const l_complex &, const dotprecision &) throw();

    //! Implementation of standard algebraic subtraction operation
    friend cdotprecision operator- 
      (const dotprecision &, const l_complex &) throw();


// ------------- Multiplication ---------------------------------------------
// ------------------ l_complex * l_complex ---------------------------------

//! Implementation of standard algebraic multiplication operation
friend l_complex operator * (const l_complex& a, const l_complex& b)
    throw();

// ------------------ l_complex * complex -----------------------------------
//! Implementation of standard algebraic multiplication operation
friend l_complex operator * (const l_complex& a, const complex& b)
    throw();

//! Implementation of standard algebraic multiplication operation
friend l_complex operator * ( const complex& b, const l_complex& a )
    throw();

// ------------------ l_complex * real ---------------------------------------
    //! Implementation of standard algebraic multiplication operation
    inline l_complex operator * (const real &) const throw();

    //! Implementation of standard algebraic multiplication operation
    friend inline l_complex operator *
                   (const real &, const l_complex &) throw();

// ------------------ l_complex * l_real -------------------------------------
    //! Implementation of standard algebraic multiplication operation
    inline l_complex operator * (const l_real &) const throw();

    //! Implementation of standard algebraic multiplication operation
    friend inline l_complex operator *
                   (const l_real &, const l_complex &) throw();


// ----------------- Others --------------------------------------------------
    //! Returns the real part of the complex value
    friend l_real & Re(l_complex& a); // { return a.re; }
    //! Returns the real part of the complex value
    friend l_real   Re(const l_complex& a); // { return a.re; }
    //! Returns the imaginary part of the complex value
    friend l_real & Im(l_complex& a); // { return a.im; }
    //! Returns the imaginary part of the complex value
    friend l_real   Im(const l_complex& a); // { return a.im; }

    //! Returns the conjugated complex value
    friend inline l_complex conj(const l_complex&) throw(); // conjugated value
    //! Sets the real part of a complex value
    friend l_complex & SetRe(l_complex & a,const l_real & b); 
	// { a.re=b; return a; } // The real part of a is substituted by b.
    // SetRe(lc,lr); --> Re(lc)=lr;
    // lc1 = SetRe(lc,lr); --> Re(lc)=lr; and lc1 = lc;
    //! Sets the imaginary part of a complex value
    friend l_complex & SetIm(l_complex & a,const l_real & b); 
        // { a.im=b; return a; } // See SetRe(...);


// ----- accumulate(cdotprecision,l_complex,l_complex|complex|real|l_real ----
    //! The accurate scalar product of the last two arguments added to the value of the first argument
    friend void accumulate(cdotprecision&, const l_complex&, 
                                           const l_complex&) throw();

    //! The accurate scalar product of the last two arguments added to the value of the first argument
    friend void accumulate(cdotprecision&, const l_complex&, 
                                           const complex&) throw();

    //! The accurate scalar product of the last two arguments added to the value of the first argument
    friend void accumulate(cdotprecision&, const l_complex&, 
                                           const real&) throw();

    //! The accurate scalar product of the last two arguments added to the value of the first argument
    friend void accumulate(cdotprecision&, const l_complex&, 
                                           const l_real&) throw();

// ---------------- cdotprecision +(-)= l_complex ---------------------------- 
//! Implementation of standard algebraic addition and allocation operation
friend inline cdotprecision & operator += (cdotprecision &cd, 
                                           const l_complex &lc) throw(); 
//! Implementation of standard algebraic subtraction and allocation operation
friend inline cdotprecision & operator -= (cdotprecision &cd, 
                                           const l_complex &lc) throw(); 

// ---------------- l_complex +(-)= l_complex|complex|real|l_real ------------
//! Implementation of standard algebraic addition and allocation operation
friend inline l_complex & operator += (l_complex &,const l_complex &) throw(); 
//! Implementation of standard algebraic subtraction and allocation operation
friend inline l_complex & operator -= (l_complex &,const l_complex &) throw(); 
//! Implementation of standard algebraic addition and allocation operation
friend inline l_complex & operator += (l_complex &,const complex &) throw(); 
//! Implementation of standard algebraic subtraction and allocation operation
friend inline l_complex & operator -= (l_complex &,const complex &) throw();
//! Implementation of standard algebraic addition and allocation operation
friend inline l_complex & operator += (l_complex &,const real &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
friend inline l_complex & operator -= (l_complex &,const real &) throw();
//! Implementation of standard algebraic addition and allocation operation
friend inline l_complex & operator += (l_complex &,const l_real &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
friend inline l_complex & operator -= (l_complex &,const l_real &) throw(); 

// ---------------- l_complex *= l_complex|complex|real|l_real ---------------
//! Implementation of standard algebraic multiplication and allocation operation
friend inline l_complex & operator *= (l_complex &,const l_complex &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
friend inline l_complex & operator *= (l_complex &,const complex &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
friend inline l_complex & operator *= (l_complex &,const real &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
friend inline l_complex & operator *= (l_complex &,const l_real &) throw();

// ---------------- Compare Operators ----------------------------------------
friend inline bool operator! (const l_complex &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const l_complex &, const l_complex &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const l_complex &, const l_complex &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const l_complex &, const complex &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const complex &, const l_complex &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const l_complex &, const complex &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const complex &, const l_complex &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const l_complex &, const real &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const real &, const l_complex &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const l_complex &, const real &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const real &, const l_complex &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const l_complex &, const l_real &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const l_real &, const l_complex &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const l_complex &, const l_real &) throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const l_real &, const l_complex &) throw();
//! Implementation of standard equality operation
friend inline bool operator== (const l_complex &, const dotprecision &)
                                                                  throw();
//! Implementation of standard equality operation
friend inline bool operator== (const dotprecision &, const l_complex &)
                                                                  throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const l_complex &, const dotprecision &)
                                                                  throw();
//! Implementation of standard negated equality operation
friend inline bool operator!= (const dotprecision &, const l_complex &)
                                                                  throw();

//! Implementation of standard equality operation
friend inline bool operator ==(const cdotprecision &, const l_complex &)
    throw(); // {l_complex.inl}
//! Implementation of standard equality operation
friend inline bool operator ==(const l_complex &, const cdotprecision &)
    throw(); // {l_complex.inl}
//! Implementation of standard negated equality operation
friend inline bool operator !=(const cdotprecision &, const l_complex &)
    throw(); // {l_complex.inl}
//! Implementation of standard negated equality operation
friend inline bool operator !=(const l_complex &, const cdotprecision &)
    throw(); // {l_complex.inl}

// -------------- Division: Directed rounding, Blomquist 28.11.02 ------------

    //! Division of two real values and rounding to the nearest value
    friend l_complex divn (const l_complex &, const l_complex &);
    //! Division of two real values and rounding the result downwards
    friend l_complex divd (const l_complex &, const l_complex &);
    //! Division of two real values and rounding the result upwards
    friend l_complex divu (const l_complex &, const l_complex &);

// -------------- Division:  Blomquist 28.11.02 ------------------------------
    //! Implementation of standard algebraic division operation
    friend l_complex operator / (const l_complex &,const l_complex &) throw();

    //! Implementation of standard algebraic division operation
    friend inline l_complex operator / (const l_complex & a,
                                        const complex & b) throw();
    //! Implementation of standard algebraic division operation
    friend inline l_complex operator / (const l_complex & a,
                                        const l_real & b) throw();
    //! Implementation of standard algebraic division operation
    friend inline l_complex operator / (const l_complex & a,
                                        const real & b) throw();

    //! Implementation of standard algebraic division operation
    friend inline l_complex operator / (const complex& a,
                                        const l_complex& b) throw();
    //! Implementation of standard algebraic division operation
    friend inline l_complex operator / (const real& a,
                                        const l_complex& b) throw();
    //! Implementation of standard algebraic division operation
    friend inline l_complex operator / (const l_real& a,
                                        const l_complex& b) throw();

    //! Implementation of standard algebraic division and allocation operation
    friend inline l_complex& operator /=(l_complex&,const l_complex&) throw();
    //! Implementation of standard algebraic division and allocation operation
    friend inline l_complex& operator /=(l_complex&,const complex&) throw();
    //! Implementation of standard algebraic division and allocation operation
    friend inline l_complex& operator /=(l_complex&,const real&) throw();
    //! Implementation of standard algebraic division and allocation operation
    friend inline l_complex& operator /=(l_complex&,const l_real&) throw();

//! The absolute value of a l_complex value
friend l_real abs2(const l_complex &a) throw(); // a.re*a.re + a.im*a.im;
//! The absolute value of a l_complex value
friend l_real abs (const l_complex &z) throw(); 

// ----------------------- Output --------------------------------------------

//! Implementation of standard output method
friend std::ostream& operator << (std::ostream& s,const l_complex& z ) throw()
// A complex number z of type l_complex is written to the output channel.
{     
    s << '('          
    << z.re << ", "   
    << z.im       
    << ')';
    return s;
}

//! Implementation of standard output method
friend std::string & operator << (std::string &s, const l_complex& a) throw()
// The value of a variable a of type l_complex is copied to a string s.
// s has the form:  (Re(a),Im(a))
{  
    s+='(';
    s << a.re;
    s+=", ";
    s << a.im; 
    s+=')';
    return s;
}

// ----------------------- Input ---------------------------------------------

//! Implementation of standard input method
friend std::istream & operator >> (std::istream &s, l_complex &a) throw()
// An input of a complex number z of the form (Re(z),Im(z)) is copied to
// the variable a of type l_complex.
{  
   char c;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, '(');
   if (inpdotflag) 
      s.putback(c);

   s >> a.re;

   skipeolnflag = inpdotflag = true;
   c = skipwhitespacessinglechar (s, ',');
   if (inpdotflag) s.putback(c);

   s >> a.im >> RestoreOpt;

   if (!waseolnflag) 
   {
      skipeolnflag = false, inpdotflag = true;
      c = skipwhitespaces (s);
      if (inpdotflag && c != ')') 
         s.putback(c);
   }
   return s;
}

//! Implementation of standard input method
friend std::string & operator >> (std::string &s, l_complex &a) throw()
// A complex number z of the form (Re(z),Im(z)), represented in a string s
// is copied to a of type l_complex.
{   
   s = skipwhitespacessinglechar (s, '(');
   s >> SaveOpt >> a.re;
   s = skipwhitespacessinglechar (s, ',');
   s >> a.im >> RestoreOpt;
   s = skipwhitespaces (s);

   if (s[0] == ')') 
      s.erase(0,1);
   return s;
}

}; // end of class l_complex

inline l_complex _l_complex(const l_real &a) throw()
	{ return l_complex(a); }
inline l_complex _l_complex(const l_real &a, const l_real &b)
        throw()   { return l_complex(a,b); }
inline l_complex _l_complex(const real &a) throw()
	{ return l_complex(a); }
inline l_complex _l_complex(const real &a, const real &b)
        throw()   { return l_complex(a,b); }
inline l_complex _l_complex(const complex &c)
        throw()   { return l_complex(c); }
inline l_complex _l_complex(const dotprecision &d)
        throw()   { return l_complex(d); }
inline l_complex conj(const l_complex&) throw();
//inline l_complex _l_complex(const cdotprecision &cd)
//        throw()   { return l_complex(cd); }

l_real & Re(l_complex& a);
l_real   Re(const l_complex& a);
l_real & Im(l_complex& a);
l_real   Im(const l_complex& a);

l_complex & SetRe(l_complex & a,const l_real & b);
	//{ a.re=b; return a; } // The real part of a is substituted by b.
    // SetRe(lc,lr); --> Re(lc)=lr;
    // lc1 = SetRe(lc,lr); --> Re(lc)=lr; and lc1 = lc;
l_complex & SetIm(l_complex & a,const l_real & b); 
//        { a.im=b; return a; } // See SetRe(...);

l_complex _l_complex(const cdotprecision &) throw();

}  // end namespace cxsc

#include "l_complex.inl"

#endif // _CXSC_L_COMPLEX_HPP_INCLUDED
