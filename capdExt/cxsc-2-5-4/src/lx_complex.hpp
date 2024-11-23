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

/* CVS $Id: lx_complex.hpp,v 1.9 2014/01/30 17:23:47 cxsc Exp $ */


/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

/*
**  Implementation of the classes
**
**  lx_complex    with all tools and elementary functions for complex
**                point and interval aruments
**
*/

#ifndef _CXSC_LX_COMPLEX_HPP_INCLUDED
#define _CXSC_LX_COMPLEX_HPP_INCLUDED

#include <iostream>
#include <except.hpp>
#include <l_cinterval.hpp>
#include <l_complex.hpp>
#include "lx_interval.hpp"

namespace cxsc {
// --------------------------------------------------------------------------
//      Class lx_complex
// --------------------------------------------------------------------------
	 
	class lx_complex
	{
		private:
    // ------------- Data Elements -------------------------------------------
			lx_real re;
			lx_real im;

		public:
    // ------------- Constructors --------------------------------------------

	//! Constructor of class lx_complex
			lx_complex(void)  throw() {}
	//! Constructor of class lx_complex
			lx_complex(const real& nr, const l_real &ar, const real& ni, const l_real &ai)
					throw() : re(lx_real(nr,ar)), im(lx_real(ni,ai))  { }
	//! Constructor of class lx_complex
			lx_complex(const real& n, const real &a) 
					throw() : re(lx_real(n,a)), im(lx_real(0.0)) { }
	//! Constructor of class lx_complex
			lx_complex(const lx_real &a, const lx_real &b)
					throw() : re(a), im(b) { }
	//! Constructor of class lx_complex
			lx_complex(const lx_real &a) throw() : re(a), im(lx_real(0.0)) { }
	//! Constructor of class lx_complex
			explicit lx_complex(const l_real &a, const l_real &b) 
					throw() : re(a), im(b) { }
	//! Constructor of class lx_complex
			explicit lx_complex(const l_real &a) 
					throw() : re(a), im(lx_real(0.0)) { }
	//! Constructor of class lx_complex
			explicit lx_complex(const real &a) 
					throw() : re(a), im(lx_real(0.0)) { }
	//! Constructor of class lx_complex
			explicit lx_complex(const complex &a) 
					throw() : re(Re(a)), im(Im(a)) { }
	//! Constructor of class lx_complex
			explicit lx_complex(const l_complex &a) 
					throw() : re(lx_real(0,Re(a))), im(lx_real(0,Im(a))) { }
	//! Constructor of class lx_complex
			lx_complex(const real& pr, const string& sr,
						  const real& pi, const string& si)
					throw() : re(lx_real(pr,sr)), im(lx_real(pi,si)) { }

// ------------- Assignments ---------------------------------------------
//! Implementation of standard assigning operator
			lx_complex & operator = (const lx_real & lr) throw() 
			{ re = lr; im = 0; return *this; }
//! Implementation of standard assigning operator
			lx_complex & operator = (const l_real & lr) throw() 
			{ re = lr; im = 0; return *this; }
//! Implementation of standard assigning operator
			lx_complex & operator = (const real & lr) throw() 
			{ re = lr; im = 0; return *this; }
//! Implementation of standard assigning operator
			lx_complex & operator = (const l_complex & c) throw() 
			{ re = Re(c); im = Im(c); return *this; }
//! Implementation of standard assigning operator
			lx_complex & operator = (const complex & c) throw() 
			{ re = Re(c); im = Im(c); return *this; } 
    
    // ------------- Functions -----------------------------------------------

//! Returns the real part of a complex value
			friend inline lx_real Re(const lx_complex&);
//! Returns the imaginary part of a complex value
			friend inline lx_real Im(const lx_complex&);
//! Returns the precision of the real part
			friend inline int StagPrec(const lx_complex&) throw();
//! Returns the exponent of base 2 of the real part
			friend inline real expoRe(const lx_complex&) throw();
//! Returns the exponent of base 2 of the imaginary part
			friend inline real expoIm(const lx_complex&) throw();
//! Returns the lr-part of the real part of a complex value
			friend inline l_real lr_partRe(const lx_complex&) throw();
//! Returns the lr-part of the imaginary part of a complex value
			friend inline l_real lr_partIm(const lx_complex&) throw();

//! Sets a new real part of a complex value
			friend inline lx_complex & SetRe(lx_complex&, const lx_real&);
//! Sets a new real part of a complex value
			friend inline lx_complex & SetRe(lx_complex&, const l_real&);
//! Sets a new real part of a complex value
			friend inline lx_complex & SetRe(lx_complex&, const real&);

//! Sets a new imaginary part of a complex value
			friend inline lx_complex & SetIm(lx_complex&, const lx_real&);
//! Sets a new imaginary part of a complex value
			friend inline lx_complex & SetIm(lx_complex&, const l_real&);
//! Sets a new imaginary part of a complex value
			friend inline lx_complex & SetIm(lx_complex&, const real&);
	 
//! Returns the conjugated complex value
			friend inline lx_complex conj(const lx_complex& a) throw();
//! Returns the absolute value of a complex number
			friend lx_real abs  (const lx_complex&) throw();
//! Returns the square of the absolute value of a complex number
			friend lx_real abs2 (const lx_complex&) throw();
	 
//! Implementation of standard negation operation
			friend inline bool operator ! (const lx_complex& ) throw();
	 
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_complex&, const lx_complex&) 
					throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_complex&, const l_complex&) 
					throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_complex&, const complex&) 
					throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const l_complex&, const lx_complex&) 
					throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const complex&,   const lx_complex&) 
					throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_complex&, const lx_real&) throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_complex&, const l_real&) throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_complex&, const real&)   throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const lx_real&, const lx_complex&) throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const l_real&, const lx_complex&) throw();
//! Implementation of standard equality operation
			friend inline bool operator == (const real&,   const lx_complex&) throw();

//! Implementation of standard algebraic negative sign operation
			friend inline lx_complex operator - (const lx_complex &) throw();
	 

// ----------------------- Output --------------------------------------------
	 
//! Implementation of standard output method
			friend inline std::ostream& operator << (std::ostream& s, const lx_complex& a) 
					throw(); // A value a of type lx_complex is written to the output channel.
//! Implementation of standard output method
			friend inline std::string & operator << (std::string &s, const lx_complex& a) 
					throw();
// The value of a variable a of type lx_complex is copied to a string s.
// s has the form:  ({2**()*...} , {2**()*...})

// ----------------------- Input --------------------------------------------

//! Implementation of standard input method
			friend std::string & operator >> (std::string& s, lx_complex& a) throw();

//! Implementation of standard input method
friend std::istream & operator >> (std::istream &s, lx_complex &a) throw()
// An input of a complex number z of the form ({...,...},{...,...}) is
// copied to the variable a of type lx_complex.
{
	char c;
	std::cout << "Real part: {Exponent to base 10, real number} = ?" 
			<< std::endl;
	s >> a.re;
	std::cout << "Img. part: {Exponent to base 10, real number} = ?" 
			<< std::endl;	
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

//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const lx_complex&, const lx_complex&)
					throw();
//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const lx_complex&, const lx_real&) throw();
//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const lx_real&, const lx_complex&) throw();
//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const lx_complex&, const l_real&) throw();
//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const l_real&, const lx_complex&) throw();
//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const lx_complex&, const real&) throw();
//! Implementation of standard algebraic addition operation
			friend inline lx_complex operator + (const real&, const lx_complex&) throw();

//! Implementation of standard algebraic multiplication operation
			friend inline lx_complex operator * (const lx_complex&, const lx_complex&)
					throw();
		
//! Implementation of standard algebraic division operation
			friend inline lx_complex operator / (const lx_complex&, const lx_complex&)
					throw();

}; // end of class lx_complex

// ***********************************************************************
// ---------Functions related to type lx_complex -------------------------
// ***********************************************************************

    inline std::ostream& operator << (std::ostream& s, const lx_complex&) 
			 throw();
	 inline std::string & operator << (std::string&  s, const lx_complex&) 
			 throw(); 
	 
// ------- friend functions declared in class lx_complex: -----------

	 inline lx_real Re(const lx_complex&);
	 inline lx_real Im(const lx_complex&);
	 inline int StagPrec(const lx_complex&) throw();
	 inline real expoRe(const lx_complex&) throw();
	 inline real expoIm(const lx_complex&) throw();
	 inline l_real lr_partRe(const lx_complex&) throw();
	 inline l_real lr_partIm(const lx_complex&) throw();

	 inline lx_complex & SetRe(lx_complex&, const lx_real&);
	 inline lx_complex & SetRe(lx_complex&, const l_real&);
	 inline lx_complex & SetRe(lx_complex&, const real&);

	 inline lx_complex & SetIm(lx_complex&, const lx_real&);
	 inline lx_complex & SetIm(lx_complex&, const l_real&);
	 inline lx_complex & SetIm(lx_complex&, const real&);
	 
	 inline lx_complex conj(const lx_complex&) throw();
	 
	 lx_real abs  (const lx_complex&) throw();
	 lx_real abs2 (const lx_complex&) throw();

	 inline bool operator == (const lx_complex&, const lx_complex&) throw();

	 inline bool operator == (const lx_complex&, const l_complex&) throw();
	 inline bool operator == (const lx_complex&, const complex&) throw();
	 inline bool operator == (const l_complex&, const lx_complex&) throw();
	 inline bool operator == (const complex&,   const lx_complex&) throw();

	 inline bool operator == (const lx_complex&, const lx_real&) throw();
	 inline bool operator == (const lx_complex&, const l_real&) throw();
	 inline bool operator == (const lx_complex&, const real&)   throw();
	 inline bool operator == (const lx_real&, const lx_complex&) throw();
	 inline bool operator == (const l_real&, const lx_complex&) throw();
	 inline bool operator == (const real&,   const lx_complex&) throw();

// ----------------------------------------------------------------------
	 
//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_complex&, const lx_complex&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_complex&, const l_complex&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_complex&, const complex&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const l_complex&, const lx_complex&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const complex&,   const lx_complex&) throw();

//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_complex&, const lx_real&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_complex&, const l_real&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_complex&, const real&)   throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const lx_real&, const lx_complex&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const l_real&, const lx_complex&) throw();
//! Implementation of standard negated equality operation
	 inline bool operator != (const real&,   const lx_complex&) throw();

	 inline lx_complex operator - (const lx_complex &) throw();

// -------------------------- Output ------------------------------------

	 inline std::ostream& operator << (std::ostream& s, const lx_complex& a) 
			 throw(); // A value a of type lx_complex is written to the output channel.
	 inline std::string & operator << (std::string& s, const lx_complex& a) 
			 throw();
    // The value of a variable a of type lx_complex is copied to a string s.
    // s has the form:  { ? , ? }
	 
// -------------------------- Input ------------------------------------
	 
//! Implementation of standard input method
	 std::string & operator >> (std::string& s, lx_complex& a) throw();
// Writes string s to variable a of type lx_complex;
// and returns an empty string s;
// Example:  s = "{-4000,2}" delivers a value a
// with:    10^(-4000)*2 ~ a;
	
//! Implementation of standard input method
	 void operator >> (const std::string &s, lx_complex &a) throw();
//! Implementation of standard input method
	 void operator >> (const char *s, lx_complex& a) throw();

// ---- function and operator declarations outside the class lx_complex ----
	
//! Implementation of standard algebraic positive sign operation
	 inline lx_complex operator + (const lx_complex&) throw();

//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const lx_complex&, const l_complex&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const lx_complex&, const complex&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const l_complex&, const lx_complex&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const complex&, const lx_complex&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const lx_complex&, const lx_real&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const lx_real&, const lx_complex&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const lx_complex&, const l_real&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const l_real&, const lx_complex&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const lx_complex&, const real&) throw();
//! Implementation of standard algebraic addition operation
	 inline lx_complex operator + (const real&, const lx_complex&) throw();
	 
	 //! Implementation of standard algebraic addition and allocation operation
	 inline lx_complex & operator +=(lx_complex &, const lx_complex &) throw();
	 //! Implementation of standard algebraic addition and allocation operation
	 inline lx_complex & operator +=(lx_complex &, const l_complex &) throw();
	 //! Implementation of standard algebraic addition and allocation operation
	 inline lx_complex & operator +=(lx_complex &, const complex &) throw();
	 //! Implementation of standard algebraic addition and allocation operation
	 inline lx_complex & operator +=(lx_complex &, const lx_real &) throw();
	 //! Implementation of standard algebraic addition and allocation operation
	 inline lx_complex & operator +=(lx_complex &, const l_real &) throw();
	 //! Implementation of standard algebraic addition and allocation operation
	 inline lx_complex & operator +=(lx_complex &, const real &) throw();
	 
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_complex&, const l_complex&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_complex&, const complex&) throw();
    //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const l_complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_complex&, const lx_real&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_complex&, const l_real&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_complex&, const real&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const lx_real&, const lx_complex&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const l_real&, const lx_complex&) throw();
	 //! Implementation of standard algebraic subtraction operation
	 inline lx_complex operator - (const real&, const lx_complex&) throw();
	 
	 //! Implementation of standard algebraic subtraction and allocation operation
	 inline lx_complex & operator -=(lx_complex &, const lx_complex &) throw();
	 //! Implementation of standard algebraic subtraction and allocation operation
	 inline lx_complex & operator -=(lx_complex &, const l_complex &) throw();
	 //! Implementation of standard algebraic subtraction and allocation operation
	 inline lx_complex & operator -=(lx_complex &, const complex &) throw();
	 //! Implementation of standard algebraic subtraction and allocation operation
	 inline lx_complex & operator -=(lx_complex &, const lx_real &) throw();
	 //! Implementation of standard algebraic subtraction and allocation operation
	 inline lx_complex & operator -=(lx_complex &, const l_real &) throw();
	 //! Implementation of standard algebraic subtraction and allocation operation
	 inline lx_complex & operator -=(lx_complex &, const real &) throw();
	 
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_complex&, const l_complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_complex&, const complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const l_complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const complex&,   const lx_complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_complex&, const lx_real&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_complex&, const l_real&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_complex&, const real&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const lx_real&, const lx_complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const l_real&,  const lx_complex&) throw();
	 //! Implementation of standard algebraic multiplication operation
	 inline lx_complex operator * (const real&,    const lx_complex&) throw();
	 
	 //! Implementation of standard algebraic multiplication and allocation operation
	 inline lx_complex & operator *=(lx_complex &, const lx_complex &) throw();
	 //! Implementation of standard algebraic multiplication and allocation operation
	 inline lx_complex & operator *=(lx_complex &, const l_complex &) throw();
	 //! Implementation of standard algebraic multiplication and allocation operation
	 inline lx_complex & operator *=(lx_complex &, const complex &) throw();
	 //! Implementation of standard algebraic multiplication and allocation operation
	 inline lx_complex & operator *=(lx_complex &, const lx_real &) throw();
	 //! Implementation of standard algebraic multiplication and allocation operation
	 inline lx_complex & operator *=(lx_complex &, const l_real &) throw();
	 //! Implementation of standard algebraic multiplication and allocation operation
	 inline lx_complex & operator *=(lx_complex &, const real &) throw();
	 
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_complex&, const l_complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_complex&, const complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const l_complex&, const lx_complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const complex&,   const lx_complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_complex&, const lx_real&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_complex&, const l_real&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_complex&, const real&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const lx_real&, const lx_complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const l_real&,  const lx_complex&) throw();
	 //! Implementation of standard algebraic division operation
	 inline lx_complex operator / (const real&,    const lx_complex&) throw();
	 
	 //! Implementation of standard algebraic division and allocation operation
	 inline lx_complex & operator /=(lx_complex &, const lx_complex &) throw();
	 //! Implementation of standard algebraic division and allocation operation
	 inline lx_complex & operator /=(lx_complex &, const l_complex &) throw();
	 //! Implementation of standard algebraic division and allocation operation
	 inline lx_complex & operator /=(lx_complex &, const complex &) throw();
	 //! Implementation of standard algebraic division and allocation operation
	 inline lx_complex & operator /=(lx_complex &, const lx_real &) throw();
	 //! Implementation of standard algebraic division and allocation operation
	 inline lx_complex & operator /=(lx_complex &, const l_real &) throw();
	 //! Implementation of standard algebraic division and allocation operation
	 inline lx_complex & operator /=(lx_complex &, const real &) throw();


	 
// ---------------------------------------------------------------------------
// ----- Elementary functions related to lx_complex
// ---------------------------------------------------------------------------

//! Calculates \f$ \mbox{sqr}(z) \f$
	 lx_complex sqr(const lx_complex&) throw();
//! Calculates \f$ \sqrt{z} \f$
	 lx_complex sqrt(const lx_complex&) throw();
//! Calculates \f$ \sqrt[n]{z} \f$
	 lx_complex sqrt(const lx_complex& ,int) throw();
//! Calculates \f$ \exp(z) \f$
	 lx_complex exp(const lx_complex&) throw();
//! Calculates \f$ 2^z \f$
	 lx_complex exp2(const lx_complex&) throw();
//! Calculates \f$ 10^z \f$
	 lx_complex exp10(const lx_complex&) throw();
//! Calculates \f$ \sin(z) \f$
	 lx_complex sin(const lx_complex&) throw();
//! Calculates \f$ \cos(z) \f$
	 lx_complex cos(const lx_complex&) throw();
//! Calculates \f$ \tan(z) \f$
	 lx_complex tan(const lx_complex&) throw();
//! Calculates \f$ \cot(z) \f$
	 lx_complex cot(const lx_complex&) throw();
//! Calculates \f$ \arcsin(z) \f$
	 lx_complex asin(const lx_complex&) throw();
//! Calculates \f$ \arccos(z) \f$
	 lx_complex acos(const lx_complex&) throw();
//! Calculates \f$ \arctan(z) \f$
	 lx_complex atan(const lx_complex&) throw();
//! Calculates \f$ \mbox{arccot}(z) \f$
	 lx_complex acot(const lx_complex&) throw();
//! Calculates \f$ \sinh(z) \f$
	 lx_complex sinh(const lx_complex&) throw();
//! Calculates \f$ \cosh(z) \f$
	 lx_complex cosh(const lx_complex&) throw();
//! Calculates \f$ \tanh(z) \f$
	 lx_complex tanh(const lx_complex&) throw();
//! Calculates \f$ \coth(z) \f$
	 lx_complex coth(const lx_complex&) throw();
//! Calculates \f$ \mbox{arcsinh}(z) \f$
	 lx_complex asinh(const lx_complex&) throw();
//! Calculates \f$ \mbox{arccosh}(z) \f$
	 lx_complex acosh(const lx_complex&) throw();
//! Calculates \f$ \mbox{arctanh}(z) \f$
	 lx_complex atanh(const lx_complex&) throw();
//! Calculates \f$ \mbox{arccoth}(z) \f$
	 lx_complex acoth(const lx_complex&) throw();
//! Calculates \f$ \sqrt{z} \f$ and returns all possible solutions
	 std::list<lx_complex>sqrt_all(const lx_complex&);
//! Calculates \f$ \mbox{arg}(z) \f$
	 lx_real arg(const lx_complex&) throw();
//! Calculates \f$ \mbox{arg}(z) \f$
	 lx_real Arg(const lx_complex&) throw();
//! Calculates \f$ \sqrt[n]{z} \f$ and returns all possible solutions
	 std::list<lx_complex>sqrt_all(const lx_complex&, int);
//! Calculates \f$ \ln(z) \f$
	 lx_complex ln(const lx_complex&) throw();
//! Calculates \f$ \mbox{log2}(z) \f$
	 lx_complex log2(const lx_complex&) throw();
//! Calculates \f$ \mbox{log10}(z) \f$
	 lx_complex log10(const lx_complex&) throw();
//! Calculates \f$ z^n \f$
	 lx_complex power_fast(const lx_complex&, const real&) throw();
//! Calculates \f$ z^n \f$
	 lx_complex power(const lx_complex&, const real&) throw();
//! Calculates \f$ z^y \f$
	 lx_complex pow(const lx_complex&, const lx_real&) throw();
//! Calculates \f$ z_1^{z_2} \f$
	 lx_complex pow(const lx_complex&, const lx_complex&) throw();

//! Calculates \f$ \mbox{sqrt}(1+[z]^2) \f$
	 lx_complex sqrt1px2(const lx_complex&) throw();
//! Calculates \f$ \mbox{sqrt}(1-[z]^2) \f$	 
	 lx_complex sqrt1mx2(const lx_complex&) throw();
//! Calculates \f$ \mbox{sqrt}([z]^2-1) \f$
	 lx_complex sqrtx2m1(const lx_complex&) throw();
//! Calculates \f$ \mbox{sqrt}([z]+1)-1 \f$
	 lx_complex sqrtp1m1(const lx_complex&) throw();
//! Calculates \f$ \mbox{exp}([z])-1 \f$
	 lx_complex expm1(const lx_complex&) throw();
//! Calculates \f$ \mbox{ln}(1+[z]) \f$	 
	 lx_complex lnp1(const lx_complex&) throw();

}  // end namespace cxsc

#include "lx_complex.inl"

#endif // _CXSC_LX_COMPLEX_HPP_INCLUDED
