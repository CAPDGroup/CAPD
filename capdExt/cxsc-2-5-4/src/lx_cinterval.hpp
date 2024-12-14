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

/* CVS $Id: lx_cinterval.hpp,v 1.9 2014/01/30 17:23:47 cxsc Exp $ */


/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

/*
**  Implementation of the classes
**
**  lx_cinterval  with all tools and elementary functions for complex
**                point and interval aruments
**
*/

#ifndef _CXSC_LX_CINTERVAL_HPP_INCLUDED
#define _CXSC_LX_CINTERVAL_HPP_INCLUDED

#include <iostream>
#include <except.hpp>
#include <l_cinterval.hpp>
#include <l_complex.hpp>
#include "lx_interval.hpp"
#include "lx_complex.hpp"

namespace cxsc {
	
// --------------------------------------------------------------------------
//      Class lx_cinterval
// --------------------------------------------------------------------------

class lx_cinterval
{
private:
    // ----------------- private data elements -------------------------------
    lx_interval re, im;
    // (re,im) is a complex number:  re + i*im, i = sqrt(-1).
public:
    // ------------- Constructors --------------------------------------------

    //! Constructor of class lx_cinterval
    inline lx_cinterval(void)  { }
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const lx_interval &, const lx_interval &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const l_interval &, const l_interval &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const interval &,   const interval &)  ;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const l_real &,     const l_real &)    ;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const lx_real &,     const lx_real &)    ;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const real &,       const real &)      ;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const l_cinterval &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const cinterval &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const complex &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const l_complex &);
    //! Constructor of class lx_cinterval 
    inline lx_cinterval(const lx_complex &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const lx_complex&, const lx_complex&)
	;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const l_complex&, const l_complex&)
	;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const complex&, const complex&)
	;
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const real&, const l_interval&, const real&, const l_interval&);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const real&, const l_interval&);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const real&, const l_interval&, const lx_interval&);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const lx_interval&, const real&, const l_interval&);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const real&, const string&, const real&, const string&);
    //! Constructor of class lx_cinterval
    explicit inline lx_cinterval(const lx_interval &);
    //! Constructor of class lx_cinterval
    explicit inline lx_cinterval(const l_interval &);
    //! Constructor of class lx_cinterval
    inline lx_cinterval(const interval &);
    //! Constructor of class lx_cinterval
    explicit inline lx_cinterval(const lx_real &);
    //! Constructor of class lx_cinterval
    explicit inline lx_cinterval(const l_real &);
    //! Constructor of class lx_cinterval
    explicit inline lx_cinterval(const real &);
    //! Constructor of class lx_cinterval

    // ------------- Assignments ---------------------------------------------

    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const lx_cinterval & );
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const l_cinterval & );
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const cinterval & )  ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const lx_interval & ) ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const l_interval & ) ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const interval & )   ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const lx_real & )     ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const l_real & )     ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const real & )       ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const lx_complex & )  ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const l_complex & )  ;
    //! Implementation of standard assigning operator
    inline lx_cinterval & operator = (const complex & )    ;

// ----------------------- Output --------------------------------------------

//! Implementation of standard output method
friend inline std::ostream& operator << (std::ostream& s,const lx_cinterval& a) 
   ;
// A value a of type lx_cinterval is written to the output channel.

//! Implementation of standard output method
friend inline std::string & operator << (std::string &s,const lx_cinterval& a) 
   ;
// The value of a variable a of type lx_cinterval is copied to a string s.
// s has the form:  {ex,li}


// ---------------------- Arithmetic operators ------------------------------

//! Implementation of standard algebraic negative sign operation
friend inline lx_cinterval operator -(const lx_cinterval &);

//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &,const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &,const l_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const l_cinterval &,const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const cinterval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const lx_interval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_interval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const l_interval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const l_interval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const lx_real &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const l_real &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const l_real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const real &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const complex &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const complex &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const l_complex &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const l_complex &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_cinterval &, const lx_complex &) 
   ;
//! Implementation of standard algebraic addition operation
friend inline lx_cinterval operator + (const lx_complex &, const lx_cinterval &) 
   ;


//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &,const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &,const l_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const l_cinterval &,const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const cinterval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const lx_interval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_interval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const l_interval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const l_interval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const lx_real &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const l_real &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const l_real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const real &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const complex &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const complex &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const l_complex &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const l_complex &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_cinterval &, const lx_complex &) 
   ;
//! Implementation of standard algebraic subtraction operation
friend inline lx_cinterval operator - (const lx_complex &, const lx_cinterval &) 
   ;


//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_cinterval &,const lx_cinterval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_cinterval &, const lx_interval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_interval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_cinterval &, const l_interval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const l_interval &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_cinterval &, const lx_real &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_cinterval &, const l_real &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const l_real &, const lx_cinterval &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const lx_cinterval &, const real &) 
   ;
//! Implementation of standard algebraic multiplication operation
friend inline lx_cinterval operator * (const real &, const lx_cinterval &) 
   ;

//! Implementation of standard algebraic division operation
friend inline lx_cinterval operator / (const lx_cinterval &,const lx_cinterval &) 
   ;
//! Implementation of standard algebraic division operation
friend inline lx_cinterval operator / (const lx_cinterval &, const lx_interval &) 
   ;
//! Implementation of standard algebraic division operation
friend inline lx_cinterval operator / (const lx_cinterval &, const l_interval &) 
   ;
//! Implementation of standard algebraic division operation
friend inline lx_cinterval operator / (const lx_cinterval &, const l_real &) 
   ;
//! Implementation of standard algebraic division operation
friend inline lx_cinterval operator / (const lx_cinterval &, const lx_real &) 
   ;
//! Implementation of standard algebraic division operation
friend inline lx_cinterval operator / (const lx_cinterval &, const real &) 
   ;

//! Implementation of standard negation operation
friend inline bool operator ! (const lx_cinterval&);

//! Implementation of standard equality operation
friend inline bool operator == (const lx_cinterval &, const lx_cinterval &) 
   ;

//! Implementation of negated standard equality operation
friend inline bool operator != (const lx_cinterval &, const lx_cinterval &) 
   ;


    // --------------------- Functions ---------------------------------------

    //! Returns the real interval of the complex interval
    friend inline lx_interval Re(const lx_cinterval &);
    //! Returns the imaginary interval of the complex interval
    friend inline lx_interval Im(const lx_cinterval &);

    //! Returns the infinum of a complex interval
    friend inline lx_complex Inf(const lx_cinterval &);
    //! Returns the supremum of a complex interval
    friend inline lx_complex Sup(const lx_cinterval &);

    //! Sets the real interval of the complex interval
    friend inline lx_cinterval & SetRe(lx_cinterval&, const lx_interval&);
    //! Sets the real interval of the complex interval
    friend inline lx_cinterval & SetRe(lx_cinterval&, const l_interval&);
    //! Sets the real interval of the complex interval
    friend inline lx_cinterval & SetRe(lx_cinterval&, const interval&);
    //! Sets the real interval of the complex interval
    friend inline lx_cinterval & SetRe(lx_cinterval&, const lx_real&);
    //! Sets the real interval of the complex interval
    friend inline lx_cinterval & SetRe(lx_cinterval&, const l_real&);
    //! Sets the real interval of the complex interval
    friend inline lx_cinterval & SetRe(lx_cinterval&, const real&); 

    //! Sets the imaginary interval of the complex interval
    friend inline lx_cinterval & SetIm(lx_cinterval&, const lx_interval&);
    //! Sets the imaginary interval of the complex interval
    friend inline lx_cinterval & SetIm(lx_cinterval&, const l_interval&);
    //! Sets the imaginary interval of the complex interval
    friend inline lx_cinterval & SetIm(lx_cinterval&, const interval&);
    //! Sets the imaginary interval of the complex interval
    friend inline lx_cinterval & SetIm(lx_cinterval&, const lx_real&);
    //! Sets the imaginary interval of the complex interval
    friend inline lx_cinterval & SetIm(lx_cinterval&, const l_real&);
    //! Sets the imaginary interval of the complex interval
    friend inline lx_cinterval & SetIm(lx_cinterval&, const real&);

    //! Returns the infimum of the real interval of the complex interval
    friend inline lx_real InfRe(const lx_cinterval&); 
    //! Returns the infimum of the imaginary interval of the complex interval
    friend inline lx_real InfIm(const lx_cinterval&);
    //! Returns the supremum of the real interval of the complex interval
    friend inline lx_real SupRe(const lx_cinterval&);
    //! Returns the supremum of the imaginary interval of the complex interval
    friend inline lx_real SupIm(const lx_cinterval&);

    //! Returns the rounded middle of the complex interval
    friend inline lx_complex mid(const lx_cinterval &);
    //! Returns the rounded diameter of the complex interval
    friend inline lx_complex diam(const lx_cinterval &);
    //! Returns the exponent of the real part of the complex interval
    friend inline real expo_Re(const lx_cinterval &);
    //! Returns the exponent of the imaginary part of the complex interval
    friend inline real expo_Im(const lx_cinterval &);
    //! Returns the li_part of the real part of the complex interval
    friend inline l_interval li_part_Re(const lx_cinterval &);
    //! Returns the li_part of the imaginary part of the complex interval
    friend inline l_interval li_part_Im(const lx_cinterval &);
    //! Returns the absolute value of a complex interval
    friend inline lx_interval abs(const lx_cinterval &);
    //! matches the precision of a complex interval to the actual stagprec value
    friend inline lx_cinterval adjust(const lx_cinterval &);
    //! Returns the conjugated complex interval 
    friend inline lx_cinterval conj(const lx_cinterval &);
    //! Multiplication of an interval with \f$ 2^n \f$
    friend inline void times2pown(lx_cinterval& , const real&);
    //! Returns 1 if the argument is an empty interval
    friend inline bool IsEmpty(const lx_cinterval&);

// ------------------------- Set Operators -----------------------------------

friend inline bool operator < (const lx_cinterval &, const lx_cinterval &) 
   ;
friend inline bool operator <= (const lx_cinterval &, const lx_cinterval &) 
   ;

// ------------------------- Intersection ------------------------------------

friend inline lx_cinterval operator & (const lx_cinterval& a, 
				      const lx_cinterval& b);

// -------------------------- Convex Hull ------------------------------------

friend inline lx_cinterval operator | (const lx_cinterval& a,
			       const lx_cinterval& b);

// ---------------------------- Others --------------------------------------

friend inline lx_cinterval & SetInf(lx_cinterval& a, const lx_complex& b) 
   ;
friend inline lx_cinterval & SetInf(lx_cinterval& a, const l_complex& b) 
   ;
friend inline lx_cinterval & SetInf(lx_cinterval& a, const complex& b) 
   ;
friend inline lx_cinterval & SetInf(lx_cinterval& a, const lx_real & b) 
   ;
friend inline lx_cinterval & SetInf(lx_cinterval& a, const l_real & b) 
	;
friend inline lx_cinterval & SetInf(lx_cinterval& a, const real & b) 
	;

friend inline lx_cinterval & SetSup(lx_cinterval& a, const lx_complex& b) 
   ;
friend inline lx_cinterval & SetSup(lx_cinterval& a, const l_complex& b) 
   ;
friend inline lx_cinterval & SetSup(lx_cinterval& a, const complex& b) 
   ;
friend inline lx_cinterval & SetSup(lx_cinterval& a, const lx_real & b) 
   ;
friend inline lx_cinterval & SetSup(lx_cinterval& a, const l_real & b) 
	;
friend inline lx_cinterval & SetSup(lx_cinterval& a, const real & b) 
	;

}; // end of class lx_cinterval

// ***************************************************************************
// ***************************************************************************

// ---------------------------------------------------------------------------
// ------- friend functions declared inside the class lx_cinterval ------------
// ---------------------------------------------------------------------------

    //! Implementation of standard algebraic negative sign operation
    inline lx_cinterval operator-(const lx_cinterval &);

    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &,const l_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const l_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const l_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const lx_real &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const l_real &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const l_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const real &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const complex &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const l_complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_cinterval &, const lx_complex &) 
	;
    //! Implementation of standard algebraic addition operation
    inline lx_cinterval operator + (const lx_complex &, const lx_cinterval &) 
	;

    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &,const l_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const l_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const l_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const lx_real &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const l_real &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const l_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const real &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const complex &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const l_complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_cinterval &, const lx_complex &) 
	;
    //! Implementation of standard algebraic subtraction operation
    inline lx_cinterval operator - (const lx_complex &, const lx_cinterval &) 
	;

    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const l_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const l_real &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const l_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const lx_real &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const real &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const real &, const lx_cinterval &) 
	;

    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const l_real &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const lx_real &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const real &)
	;

    //! Implementation of standard negation operation
    inline bool operator ! (const lx_cinterval&);
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const lx_cinterval &) 
	;
//! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const lx_cinterval &) 
	;

    // ---------------------- Set Operators ---------------------------------

    //! Implementation of standard less-than operation
    inline bool operator < (const lx_cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <= (const lx_cinterval &, const lx_cinterval &) 
   ;

// -------------- friend Functions declared in lx_cinterval ------------------

    //! Returns the infimum of the real and imaginary part
    inline lx_complex Inf(const lx_cinterval &);
    //! Returns the supremum of the real and imaginary part
    inline lx_complex Sup(const lx_cinterval &);
    //! Returns the real part of the complex interval

    inline lx_interval Re(const lx_cinterval &);
    //! Returns the imaginary part of the complex interval
    inline lx_interval Im(const lx_cinterval &);
    //! Returns the complex valued centre of the complex interval

    //! Sets the real interval of the complex interval
    inline lx_cinterval & SetRe(lx_cinterval&, const lx_interval&);
    //! Sets the real interval of the complex interval
    inline lx_cinterval & SetRe(lx_cinterval&, const l_interval&);
    //! Sets the real interval of the complex interval
    inline lx_cinterval & SetRe(lx_cinterval&, const interval&);
    //! Sets the real interval of the complex interval
    inline lx_cinterval & SetRe(lx_cinterval&, const lx_real&);
    //! Sets the real interval of the complex interval
    inline lx_cinterval & SetRe(lx_cinterval&, const l_real&);
    //! Sets the real interval of the complex interval
    inline lx_cinterval & SetRe(lx_cinterval&, const real&);

    //! Sets the imaginary interval of the complex interval
    inline lx_cinterval & SetIm(lx_cinterval&, const lx_interval&);
    //! Sets the imaginary interval of the complex interval
    inline lx_cinterval & SetIm(lx_cinterval&, const l_interval&);
    //! Sets the imaginary interval of the complex interval
    inline lx_cinterval & SetIm(lx_cinterval&, const interval&);
    //! Sets the imaginary interval of the complex interval
    inline lx_cinterval & SetIm(lx_cinterval&, const lx_real&);
    //! Sets the imaginary interval of the complex interval
    inline lx_cinterval & SetIm(lx_cinterval&, const l_real&);
    //! Sets the imaginary interval of the complex interval
    inline lx_cinterval & SetIm(lx_cinterval&, const real&);

    //! Returns the infimum of the real interval of the complex interval
    inline lx_real InfRe(const lx_cinterval&); 
    //! Returns the infimum of the imaginary interval of the complex interval
    inline lx_real InfIm(const lx_cinterval&);
    //! Returns the supremum of the real interval of the complex interval
    inline lx_real SupRe(const lx_cinterval&);
    //! Returns the supremum of the imaginary interval of the complex interval
    inline lx_real SupIm(const lx_cinterval&);

    //! Returns the complex middle of the complex interval
    inline lx_complex mid(const lx_cinterval &);
    //! Returns the complex valued diameter of the complex interval
    inline lx_complex diam(const lx_cinterval &);
    //! Returns the exponent n of the real part of the complex interval
    inline real expo_Re(const lx_cinterval &a);
    //! Returns the exponent n of the imaginary part of the complex interval
    inline real expo_Im(const lx_cinterval &a);
    //! Returns the l_interval of the real part of the complex interval
    inline l_interval li_part_Re(const lx_cinterval &);
    //! Returns the l_interval of the imaginary part of the complex interval
    inline l_interval li_part_Im(const lx_cinterval &);
    //! Sets the precision of a specific long datatype value
    inline lx_cinterval adjust(const lx_cinterval &);
    //! Returns the conjugated complex interval
    inline lx_cinterval conj(const lx_cinterval &);
    //! Multiplication of interval with \f$ 2^n \f$
    inline void times2pown(lx_cinterval& , const real&);
    //! Returns the absolute value of the complex interval
    inline lx_interval abs(const lx_cinterval &); 
    //! Returns the intersection of the two complex interval operands
    inline lx_cinterval operator & (const lx_cinterval& a, 
				    const lx_cinterval& b);
    //! Returns the convex hull of the two complex interval operands
    inline lx_cinterval operator | (const lx_cinterval& a,
				    const lx_cinterval& b);

    //! Returns the complex interval with the new given infimum value
    inline lx_cinterval & SetInf(lx_cinterval& a, const lx_complex& b) 
	;
    //! Returns the complex interval with the new given infimum value
    inline lx_cinterval & SetInf(lx_cinterval& a, const l_complex& b) 
	;
    //! Returns the complex interval with the new given infimum value
    inline lx_cinterval & SetInf(lx_cinterval& a, const complex& b) 
	;
    //! Returns the complex interval with the new given infimum value
    inline lx_cinterval & SetInf(lx_cinterval& a, const lx_real & b) 
	;
    //! Returns the complex interval with the new given infimum value
    inline lx_cinterval & SetInf(lx_cinterval& a, const l_real & b) 
	;
    //! Returns the complex interval with the new given infimum value
    inline lx_cinterval & SetInf(lx_cinterval& a, const real & b) 
	;

    //! Returns the complex interval with the new given supremum value
    inline lx_cinterval & SetSup(lx_cinterval& a, const lx_complex& b) 
	;
    //! Returns the complex interval with the new given supremum value
    inline lx_cinterval & SetSup(lx_cinterval& a, const l_complex& b) 
	;
    //! Returns the complex interval with the new given supremum value
    inline lx_cinterval & SetSup(lx_cinterval& a, const complex& b) 
	;
    //! Returns the complex interval with the new given supremum value
    inline lx_cinterval & SetSup(lx_cinterval& a, const lx_real & b) 
	;
    //! Returns the complex interval with the new given supremum value
    inline lx_cinterval & SetSup(lx_cinterval& a, const l_real & b) 
	;
    //! Returns the complex interval with the new given supremum value
    inline lx_cinterval & SetSup(lx_cinterval& a, const real & b) 
	;

    //! Returns 1 if the argument is an empty interval
    inline bool IsEmpty(const lx_cinterval&);

// ***************************************************************************
// ---------------------------------------------------------------------------
// -------- Functions declared only outside the class lx_cinterval ------------
// ---------------------------------------------------------------------------
// ***************************************************************************

    //! Implementation of standard algebraic positive sign operation 
    inline lx_cinterval operator+(const lx_cinterval &);

    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const lx_cinterval &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const lx_interval &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const l_interval &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const l_cinterval &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const l_real &b);
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const lx_real &b);
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const real &b);
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const interval &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const cinterval &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const complex &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const l_complex &b) 
	;
    //! Implementation of standard algebraic addition and allocation operation
    inline lx_cinterval & operator +=(lx_cinterval &a, const lx_complex &b) 
	;

    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const lx_cinterval &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const lx_interval &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const l_interval &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const l_cinterval &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const l_real &b);
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const lx_real &b);
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const real &b);
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const interval &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const cinterval &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const complex &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const l_complex &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation
    inline lx_cinterval & operator -=(lx_cinterval &a, const lx_complex &b) 
	;
    //! Implementation of standard algebraic subtraction and allocation operation

    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &,const l_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const l_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const complex &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const l_complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_cinterval &, const lx_complex &) 
	;
    //! Implementation of standard algebraic multiplication operation
    inline lx_cinterval operator * (const lx_complex &, const lx_cinterval &) 
	;

    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const lx_cinterval &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const lx_interval &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const l_interval &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const l_cinterval &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const l_real &b);
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const lx_real &b);
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const real &b);
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const interval &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const cinterval &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const complex &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const l_complex &b) 
	;
    //! Implementation of standard algebraic multiplication and allocation operation
    inline lx_cinterval & operator *=(lx_cinterval &a, const lx_complex &b) 
	;

    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &,const l_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const l_cinterval &,const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const l_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const l_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const real &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const complex &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const l_complex &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_cinterval &, const lx_complex &) 
	;
    //! Implementation of standard algebraic division operation
    inline lx_cinterval operator / (const lx_complex &, const lx_cinterval &) 
	;

    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const l_cinterval &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const l_real &);
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const lx_real &);
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const real &);
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const interval &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const complex &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard algebraic division and allocation operation
    inline lx_cinterval & operator /=(lx_cinterval &, const lx_complex &) 
	;


    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const l_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const l_cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const l_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const l_real &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const l_real &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const lx_real &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_real &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const real &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const real &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const interval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const interval &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const complex &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const complex &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const l_complex &, const lx_cinterval &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_cinterval &, const lx_complex &) 
	;
    //! Implementation of standard equality operation
    inline bool operator == (const lx_complex &, const lx_cinterval &) 
	;


    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const l_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const l_cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const lx_interval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const l_interval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const l_interval &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const l_real &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const l_real &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const lx_real &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_real &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const real &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const real &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const interval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const interval &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const cinterval &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const complex &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const complex &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const l_complex &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const l_complex &, const lx_cinterval &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_cinterval &, const lx_complex &) 
	;
    //! Implementation of standard negated equality operation
    inline bool operator != (const lx_complex &, const lx_cinterval &) 
	;


// ------------------------- Set Operators ------------------------------

    //! Implementation of standard greater-than operation
    inline bool operator > (const lx_cinterval &, const lx_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >= (const lx_cinterval &, const lx_cinterval &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const l_cinterval &, const lx_cinterval &);
    //! Implementation of standard greater-than operation 
    inline bool operator  >(const l_cinterval &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const l_cinterval &, const lx_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const l_cinterval &, const lx_cinterval &); 
    //! Implementation of standard less-than operation
    inline bool operator  <(const lx_cinterval &, const l_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const l_cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_cinterval &, const l_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const l_cinterval &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const cinterval &, const lx_cinterval &); 
    //! Implementation of standard greater-than operation
    inline bool operator  >(const cinterval &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const cinterval &, const lx_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const cinterval &, const lx_cinterval &); 
    //! Implementation of standard less-than operation
    inline bool operator  <(const lx_cinterval &, const cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_cinterval &, const cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const cinterval &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const lx_interval &, const lx_cinterval &);
    //! Implementation of standard greater-than operation 
    inline bool operator  >(const lx_interval &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_interval &, const lx_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_interval &, const lx_cinterval &);
    //! Implementation of standard less-than operation 
    inline bool operator  <(const lx_cinterval &, const lx_interval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const lx_interval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_cinterval &, const lx_interval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const lx_interval &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const l_interval &, const lx_cinterval &); 
    //! Implementation of standard greater-than operation
    inline bool operator  >(const l_interval &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const l_interval &, const lx_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const l_interval &, const lx_cinterval &);
    //! Implementation of standard less-than operation 
    inline bool operator  <(const lx_cinterval &, const l_interval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const l_interval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_cinterval &, const l_interval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const l_interval &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const interval &, const lx_cinterval &); 
    //! Implementation of standard greater-than operation
    inline bool operator  >(const interval &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const interval &, const lx_cinterval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const interval &, const lx_cinterval &); 
    //! Implementation of standard less-than operation
    inline bool operator  <(const lx_cinterval &, const interval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const interval &);
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_cinterval &, const interval &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const interval &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const lx_real &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation 
    inline bool operator <=(const lx_real &, const lx_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const lx_real &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const lx_real &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const l_real &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation 
    inline bool operator <=(const l_real &, const lx_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const l_real &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const l_real &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const real &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation 
    inline bool operator <=(const real &, const lx_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const real &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const real &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const complex &, const lx_cinterval &);
    //! Implementation of standard less-or-equal-than operation 
    inline bool operator <=(const complex &, const lx_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const complex &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const complex &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const l_complex &, const lx_cinterval &); 
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const l_complex &, const lx_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const l_complex &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const l_complex &);

    //! Implementation of standard less-than operation
    inline bool operator  <(const lx_complex &, const lx_cinterval &); 
    //! Implementation of standard less-or-equal-than operation
    inline bool operator <=(const lx_complex &, const lx_cinterval &);
    //! Implementation of standard greater-than operation
    inline bool operator  >(const lx_cinterval &, const lx_complex &);
    //! Implementation of standard greater-or-equal-than operation
    inline bool operator >=(const lx_cinterval &, const lx_complex &);

// ------------------------- Convex Hull ---------------------------------

    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const lx_cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const lx_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_real&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const lx_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const l_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_real&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const l_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const real&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const l_cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_cinterval&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const l_cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const cinterval& a, const lx_cinterval& b) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const lx_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_interval&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const lx_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const l_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_interval&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const l_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const interval&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const lx_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_complex&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const lx_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const l_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_complex&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const l_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_cinterval&, const complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const complex&, const lx_cinterval&) 
	;
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_cinterval & operator |= (lx_cinterval&, const complex&) 
	;

    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_interval&, const complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | ( const complex&, const lx_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_interval&, const lx_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | ( const lx_complex&, const lx_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_interval&, const l_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | ( const l_complex&, const lx_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_real&, const cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const cinterval&, const lx_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_real&, const l_cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_cinterval&, const lx_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_interval&, const cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const cinterval&, const lx_interval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_interval&, const l_cinterval&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_cinterval&, const lx_interval&) 
	;

    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_real&, const complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const complex&, const lx_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_real&, const l_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const l_complex&, const lx_real&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_real&, const lx_complex&) 
	;
    //! Returns the convex hull of the arguments
    inline lx_cinterval operator | (const lx_complex&, const lx_real&) 
	;

// ------------------------- Intersection ----------------------------------

    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const lx_cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const lx_real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_real&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const lx_real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const l_real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_real&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const l_real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const real&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const l_cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_cinterval&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const l_cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const cinterval& a, const lx_cinterval& b) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const lx_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_interval&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const lx_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const l_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_interval&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const l_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const interval&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_complex&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const l_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_complex&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const l_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_cinterval&, const complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const complex&, const lx_cinterval&) 
	;
    //! Allocates the intersection of the arguments to the first argument
    inline lx_cinterval & operator &= (lx_cinterval&, const complex&) 
	;

    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_interval&, const complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const complex&, const lx_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_interval&, const l_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const l_complex&, const lx_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_interval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const lx_complex&, const lx_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_real&, const cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const cinterval&, const lx_real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_real&, const l_cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_cinterval&, const lx_real&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_interval&, const cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const cinterval&, const lx_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const lx_interval&, const l_cinterval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_cinterval&, const lx_interval&) 
	;


    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_interval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const lx_complex&, const l_interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const l_cinterval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const lx_complex&, const l_cinterval&) 
	;

    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const interval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const lx_complex&, const interval&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & (const cinterval&, const lx_complex&) 
	;
    //! Returns the intersection of the arguments
    inline lx_cinterval operator & ( const lx_complex&, const cinterval&) 
	;

// ------------------------ Input --------------------------------------

    //! Implementation of standard input method
    std::string & operator >> (std::string&, lx_cinterval&);
    //! Implementation of standard input method
    void operator >> (const std::string&, lx_cinterval&);
    //! Implementation of standard input method
    void operator >> (const char *s, lx_cinterval&);
    //! Implementation of standard input method
    std::istream & operator >> (std::istream&, lx_cinterval&);
	 

// ---------------------------------------------------------------------------
// ----- Elementary functions related to lx_cinterval
// ---------------------------------------------------------------------------
	 
    //! Calculates \f$ [z]^2  \f$
	 lx_cinterval sqr(const lx_cinterval&);
    //! Calculates \f$ \sqrt{[z]}  \f$
	 lx_cinterval sqrt(const lx_cinterval&);
    //! Calculates \f$ \sqrt[n]{[z]} \f$
	 lx_cinterval sqrt(const lx_cinterval& ,int); 
    //! Calculates \f$ \exp([z]) \f$
	 lx_cinterval exp(const lx_cinterval&);
	 //! Calculates \f$ 2^{[z]} \f$
	 lx_cinterval exp2(const lx_cinterval&);
	 //! Calculates \f$ 10^{[z]} \f$
	 lx_cinterval exp10(const lx_cinterval&);
    //! Calculates \f$ \sin([z]) \f$
	 lx_cinterval sin(const lx_cinterval&);
    //! Calculates \f$ \cos([z]) \f$
	 lx_cinterval cos(const lx_cinterval&);

    //! Calculates \f$ \cosh([z]) \f$
	 lx_cinterval cosh(const lx_cinterval&);
    //! Calculates \f$ \sinh([z]) \f$
	 lx_cinterval sinh(const lx_cinterval&);

    //! Calculates \f$ \mbox{Arg}([z]) \f$
	 lx_interval Arg(const lx_cinterval&);
    //! Calculates \f$ \mbox{arg}([z]) \f$
	 lx_interval arg(const lx_cinterval&);

    //! Calculates \f$ \ln([z]) \f$
	 lx_cinterval Ln(const lx_cinterval& );
    //! Calculates \f$ \ln([z]) \f$
	 lx_cinterval ln(const lx_cinterval& );
	 
    //! Calculates \f$ \mbox{log2}([z]) \f$
	 lx_cinterval log2(const lx_cinterval& );
    //! Calculates \f$ \mbox{log10}([z]) \f$
	 lx_cinterval log10(const lx_cinterval& );

    //! Calculates \f$ [z]^n \f$
	 lx_cinterval power_fast( const lx_cinterval&, const real& );
    //! Calculates \f$ [z]^n \f$
	 lx_cinterval power( const lx_cinterval&, const real& );
    //! Calculates \f$ [z]^{[p]} \f$
	 lx_cinterval pow( const lx_cinterval& , const lx_interval& );
    //! Calculates \f$ [z]^{[w]} \f$
	 lx_cinterval pow( const lx_cinterval& , const lx_cinterval& );

    //! Calculates \f$ \tan([z]) \f$
	 lx_cinterval tan (const lx_cinterval& );
    //! Calculates \f$ \cot([z]) \f$
	 lx_cinterval cot (const lx_cinterval& );
    //! Calculates \f$ \tanh([z]) \f$
	 lx_cinterval tanh(const lx_cinterval& );
    //! Calculates \f$ \coth([z]) \f$
	 lx_cinterval coth(const lx_cinterval& );

    //! Calculates \f$ \arcsin([z]) \f$
	 lx_cinterval asin(const lx_cinterval& );
    //! Calculates \f$ \arccos([z]) \f$
	 lx_cinterval acos(const lx_cinterval& );
    //! Calculates \f$ \arctan([z]) \f$
	 lx_cinterval atan(const lx_cinterval& );
    //! Calculates \f$ \mbox{arccot}([z]) \f$
	 lx_cinterval acot(const lx_cinterval& );

    //! Calculates \f$ \mbox{arcsinh}([z]) \f$
	 lx_cinterval asinh( const lx_cinterval& );
    //! Calculates \f$ \mbox{arccosh}([z]) \f$
	 lx_cinterval acosh( const lx_cinterval& );
    //! Calculates \f$ \mbox{arctanh}([z]) \f$
	 lx_cinterval atanh( const lx_cinterval& );
    //! Calculates \f$ \mbox{arccoth}([z]) \f$
	 lx_cinterval acoth( const lx_cinterval& );
	 
    //! Calculates \f$ \mbox{sqrt}(1+[z]^2) \f$
	 lx_cinterval sqrt1px2(const lx_cinterval& z);
    //! Calculates \f$ \mbox{sqrt}(1-[z]^2) \f$	 
	 lx_cinterval sqrt1mx2(const lx_cinterval& z);
	 //! Calculates \f$ \mbox{sqrt}([z]^2-1) \f$
	 lx_cinterval sqrtx2m1(const lx_cinterval& z);
	 //! Calculates \f$ \mbox{sqrt}([z]+1)-1 \f$
	 lx_cinterval sqrtp1m1(const lx_cinterval& z);
	 //! Calculates \f$ \mbox{exp}([z])-1 \f$
	 lx_cinterval expm1(const lx_cinterval& z);
	 //! Calculates \f$ \mbox{ln}(1+[z]) \f$	 
	 lx_cinterval lnp1(const lx_cinterval& z);
	 //! Calculates \f$ \mbox{sqrt}([z]) \f$ and returns all possible solutions
	 std::list<lx_cinterval> sqrt_all( const lx_cinterval& z);
	 //! Calculates \f$ \mbox{sqrt}[n][z] \f$ and returns all possible solutions
	 std::list<lx_cinterval> sqrt_all( const lx_cinterval& z, int n );
	 //! Calculates \f$ [z]^{[y]} \f$ and returns all possible solutions
	 std::list<lx_cinterval> pow_all( const lx_cinterval& z, const lx_interval& p );
}  // end namespace cxsc

#include "lx_cinterval.inl"

#endif // _CXSC_LX_CINTERVAL_HPP_INCLUDED
