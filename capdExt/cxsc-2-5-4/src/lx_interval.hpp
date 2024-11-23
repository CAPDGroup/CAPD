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

/* CVS $Id: lx_interval.hpp,v 1.10 2014/01/30 17:23:47 cxsc Exp $ */


/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

/*
**  Implementation of the classes
**
**  lx_interval  with all tools and elementary functions for real
**               point and interval aruments
**
*/

#ifndef _CXSC_LX_INTERVAL_HPP_INCLUDED
#define _CXSC_LX_INTERVAL_HPP_INCLUDED

#include <l_imath.hpp>
#include <lx_real.hpp>
#include <sstream>
#include <cmath>
#include <iostream>

namespace cxsc {
	
class lx_interval {
	
private:
    // ------------- Data Elements -------------------------------------------
    real ex;
    l_interval li;
    // The mathematical value of an object of type lx_interval is
    // interpreted as:  2^(ex) * li;

public:
    // ------------- Constructors --------------------------------------------

    //! Constructor of class lx_interval
    lx_interval(void)  throw() {}

    //! Constructor of class lx_interval
    lx_interval(const real& n, const l_interval& a) throw() 
    {
	if ( !(Is_Integer(n)) ) 
	    cxscthrow(REAL_NOT_ALLOWED("lx_interval(const real&, const l_interval&)"));
	else
	{
	    ex = n; li = a;
	}
    }

    //! Constructor of class lx_interval
    lx_interval(const real& n, const l_real& a) throw() 
    { 
	if ( !(Is_Integer(n)) ) 
	    cxscthrow(REAL_NOT_ALLOWED("lx_interval(const real& n, const l_real& a)"));
	else
	{
	    ex = n; li = a;
	}
    }

    //! Constructor of class lx_interval
    lx_interval(const real& n, const interval& a) throw()
    { 
	if ( !(Is_Integer(n)) ) 
	    cxscthrow(REAL_NOT_ALLOWED("lx_interval(const real&, const interval&)"));
	else
	{
	    ex = n; li = a;
	}
    }

    //! Constructor of class lx_interval
    explicit lx_interval(const real& n, const real& a) throw()
    { 
		if ( !(Is_Integer(n)) ) 
	    	cxscthrow(REAL_NOT_ALLOWED("lx_interval(const real&, const real&)"));
		else
		{
	    	ex = n; li = a;
		}
    }
    //! Constructor of class lx_interval
    explicit lx_interval(const l_interval& a) throw() : ex(0), li(a) { }
    //! Constructor of class lx_interval
    explicit lx_interval(const l_real& a)     throw() : ex(0), li(a) { }
    //! Constructor of class lx_interval
    lx_interval(const l_real& a, const l_real& b)  throw() : ex(0),li(a,b) { }
    //! Constructor of class lx_interval
    explicit lx_interval(const interval& a)   throw() : ex(0), li(a)   { }
    //! Constructor of class lx_interval
    explicit lx_interval(const real& a)       throw() : ex(0), li(a)   { }
    //! Constructor of class lx_interval
    lx_interval(const lx_real&, const lx_real&) throw();
    //! Constructor of class lx_interval
    explicit lx_interval(const lx_real& a)     
                                 throw() : ex(expo(a)), li(lr_part(a)) { }
    //! Constructor of class lx_interval
    lx_interval(const real&, const string&) throw(); 

    // ------------- Assignments ---------------------------------------------

    //! Implementation of standard assigning operator
    inline lx_interval & operator = (const lx_interval & a) throw();
    //! Implementation of standard assigning operator
    inline lx_interval & operator = (const l_interval & a) throw();
    //! Implementation of standard assigning operator
    inline lx_interval & operator = (const l_real & a)     throw();
    //! Implementation of standard assigning operator
    inline lx_interval & operator = (const real & a)       throw();
    //! Implementation of standard assigning operator
    inline lx_interval & operator = (const interval & a)   throw();
    //! Implementation of standard assigning operator
    inline lx_interval & operator = (const lx_real & a)     throw();

    // l_interval & operator = (const lx_interval &a) throw(); declared
    // in l_interval.hpp and implemented in lx_interval.cpp; 
	 
	 // interval & operator = (const lx_interval &a) throw(); declared
    // in interval.hpp and implemented in lx_interval.cpp; 

    // ------------- Functions -----------------------------------------------

    //! Returns the actual precision of an interval of type lx_interval
    friend inline int StagPrec(const lx_interval&) throw(); 
    //! Returns the exponent of base 2 of an interval of type lx_interval
    friend inline real expo(const lx_interval&) throw();
    //! Returns the l_interval part of an interval of type lx_interval
    friend inline l_interval li_part(const lx_interval&) throw();
    //! The new interval a includs the old one and the old a.li is scaled downwards
    friend void scale_down(lx_interval &a);
    //! The new interval a includs the old one and the old a.li is scaled upwards
    friend void scale_up  (lx_interval &a);

    //! matches the precision of an interval to the actual stagprec value
    friend inline lx_interval adjust(const lx_interval &) throw();
    //! Returns the absolute value of an interval
    friend inline lx_interval abs(const lx_interval &) throw();
    //! Returns 1 if the interval is a point interval
    friend inline bool point_intv(const lx_interval &);
    //! Implementation of standard negation operation
    friend inline bool operator ! (const lx_interval &) throw(); 
    //! Returns 1 if an empty interval is given
    friend inline bool IsEmpty(const lx_interval &) throw(); 
    //! Multiplication of an interval with \f$ 2^n \f$
    friend inline void times2pown(lx_interval &, const real &) throw();
    //! Returns an inflated interval a, if abs(a) is too small
    friend inline lx_interval Blow(const lx_interval &) throw(); 
    //! Multiplication of an interval with \f$ 2^n, n<0, \f$ without integer overflow
    friend inline void times2pown_neg(lx_interval &, const real&) throw();
    //! Returns the relative diameter of an interval
    friend inline lx_real RelDiam( const lx_interval &);
    //! Returns the diameter of an interval
    friend inline lx_real diam(const lx_interval &) throw();
    //! Returns the rounded middle of an interval
    friend inline lx_real mid(const lx_interval &) throw();
    //! Returns the infimum of an interval
    friend inline lx_real Inf(const lx_interval &) throw();
    //! Returns the supremum of an interval
    friend inline lx_real Sup(const lx_interval &) throw();

    // ---------------- Monadic arithmetic operator ---------------------

    //! Implementation of standard algebraic negative sign operation
    friend inline lx_interval operator-(const lx_interval & a) throw();


    // ----------------------- Output -----------------------------------

    //! Implementation of standard output method
    friend std::ostream& operator << (std::ostream&, const lx_interval&) 
	throw();

    //! Implementation of standard output method
    friend std::string & operator << (std::string&, const lx_interval&) 
	throw();

}; // end of class lx_interval


// -------------------------------------------------------------------
//    Declaration of friend functions outside the class lx_interval 
// -------------------------------------------------------------------

inline int StagPrec(const lx_interval &a) throw();
inline real expo(const lx_interval &a) throw();
inline l_interval li_part(const lx_interval &a) throw();

       void scale_down(lx_interval &);
       void scale_up  (lx_interval &);
inline bool point_intv(const lx_interval &);
inline bool IsEmpty(const lx_interval &) throw();
inline void times2pown(lx_interval &, const real &) throw();
inline lx_interval Blow(const lx_interval &) throw();
inline void times2pown_neg(lx_interval &, const real&) throw();
inline lx_real RelDiam(const lx_interval &);
inline lx_real Inf(const lx_interval &) throw();
inline lx_real Sup(const lx_interval &) throw();
inline lx_interval abs(const lx_interval &) throw();
inline lx_interval adjust(const lx_interval &) throw();
inline lx_real diam(const lx_interval &) throw();
inline lx_real mid(const lx_interval &) throw();

// ------------------------ Input --------------------------------------

    //! Implementation of standard input method
    std::string & operator >> (std::string &s, lx_interval &a) throw();
    //! Implementation of standard input method
    void operator >> (const std::string &s, lx_interval &a) throw();
    //! Implementation of standard input method
    void operator >> (const char *s, lx_interval&) throw();

    //! Implementation of standard input method
    std::istream & operator >> (std::istream&, lx_interval&) throw();

// ------------------------ Output --------------------------------------

std::ostream& operator << (std::ostream& s,const lx_interval& a) throw();
// A value a of type lx_interval is written to the 
// output channel as decimal number in the form: 
// { exponent p to base 10, interval mantissa m } = 10^p * m;

std::string & operator << (std::string &s,const lx_interval& a) throw();
// The value of a variable a of type lx_interval is copied to a string s.
// s has the form:  {2**(ex), li} = 2^ex * li;


// -------------------------------------------------------------------
// ------- Function declarations outside the class lx_interval --------
// -------------------------------------------------------------------

void Bin2Dec(const lx_interval& a, real& p, l_interval& m);

//! b = expo2zero(a) returns \f$ a\subseteq b \f$ with \f$ \verb+b.ex+=0 \f$
lx_interval expo2zero(const lx_interval &) throw(OVERFLOW_ERROR);
//! Checks arguments for disjointness
inline int Disjoint(const lx_interval &, const lx_interval &);
//! Checks if first argument lies in the interior of second argument
inline int in (const lx_interval&, const lx_interval&);
//! Checks if first argument lies in the interior of second argument
inline int in (const l_interval&, const lx_interval&);
//! Checks if first argument is part of second argument
inline int in (const interval&, const lx_interval&);

//! Checks if first argument is part of second argument       
inline int in (const lx_real&, const lx_interval&);
//! Checks if first argument is part of second argument       
inline int in (const l_real&, const lx_interval&);
//! Checks if first argument is part of second argument       
inline int in (const real&, const lx_interval&);

//! Returns an epsilon inflation of the first argument
inline lx_interval Blow( const lx_interval&, const real& );
//! Computes the smallest absolute value \f$ \left< \left[x\right] \right> \f$
inline lx_real AbsMin (const lx_interval&); 
//! Computes the greatest absolute value \f$ \left|\left[x\right]\right| \f$  
inline lx_real AbsMax (const lx_interval&); 

// -----------------------------------------------------------------------
// ------------- set comparisons -----------------------------------------
// -----------------------------------------------------------------------

// ---- lx_interval--lx_interval

//! Implementation of standard less-than operation
inline bool operator <  (const lx_interval&, const lx_interval&) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const lx_interval&, const lx_interval&) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const lx_interval&, const lx_interval&) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const lx_interval&, const lx_interval&) throw();

// ---- lx_interval--l_interval

//! Implementation of standard less-than operation
inline bool operator <  (const lx_interval&, const l_interval&) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const lx_interval&, const l_interval&) throw();
//! Implementation of standard less-than operation
inline bool operator <  (const l_interval&, const lx_interval&) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const l_interval&, const lx_interval&) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const lx_interval&, const l_interval&) throw();
inline bool operator >= (const lx_interval&, const l_interval&) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const l_interval&, const lx_interval&) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const l_interval&, const lx_interval&) throw();

// ---- lx_interval--interval

//! Implementation of standard less-than operation
inline bool operator <  (const lx_interval&, const interval&) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const lx_interval&, const interval&) throw();
//! Implementation of standard less-than operation
inline bool operator <  (const interval&, const lx_interval&) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const interval&, const lx_interval&) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const lx_interval&, const interval&) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const lx_interval&, const interval&) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const interval&, const lx_interval&) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const interval&, const lx_interval&) throw();

// ---- lx_interval--real

//! Implementation of standard less-than operation
inline bool operator <  (const real &, const lx_interval &) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const real &, const lx_interval &) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const lx_interval &, const real &) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const lx_interval &, const real &) throw();

// ---- lx_interval--l_real

//! Implementation of standard less-than operation
inline bool operator <  (const l_real &, const lx_interval &) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const l_real &, const lx_interval &) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const lx_interval &, const l_real &) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const lx_interval &, const l_real &) throw();

// ---- lx_interval--lx_real

//! Implementation of standard less-than operation
inline bool operator <  (const lx_real &, const lx_interval &) throw();
//! Implementation of standard less-or-equal-than operation
inline bool operator <= (const lx_real &, const lx_interval &) throw();
//! Implementation of standard greater-than operation
inline bool operator >  (const lx_interval &, const lx_real &) throw();
//! Implementation of standard greater-or-equal-than operation
inline bool operator >= (const lx_interval &, const lx_real &) throw();


// -------------------------- comparisons --------------------------------

//! Implementation of standard negation operation
inline bool operator ! (const lx_interval &) throw();

//! Implementation of standard equality operation
inline bool operator == (const lx_interval &, const lx_interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const lx_interval &, const l_interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const l_interval &, const lx_interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const lx_interval &, const interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const interval &, const lx_interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const lx_interval &, const real &) throw();
//! Implementation of standard equality operation
inline bool operator == (const real &, const lx_interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const lx_interval &, const l_real &) throw();
//! Implementation of standard equality operation
inline bool operator == (const l_real &, const lx_interval &) throw();
//! Implementation of standard equality operation
inline bool operator == (const lx_interval &, const lx_real &) throw();
//! Implementation of standard equality operation
inline bool operator == (const lx_real &, const lx_interval &) throw();

//! Implementation of standard negated equality operation
inline bool operator != (const lx_interval &, const lx_interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const lx_interval &, const l_interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const l_interval &, const lx_interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const lx_interval &, const interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const interval &, const lx_interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const lx_interval &, const real &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const real &, const lx_interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const lx_interval &, const l_real &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const l_real &, const lx_interval &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const lx_interval &, const lx_real &) throw();
//! Implementation of standard negated equality operation
inline bool operator != (const lx_real &, const lx_interval &) throw();

//! Implementation of standard algebraic positive sign operation
inline lx_interval operator+(const lx_interval &) throw();
//! Implementation of standard algebraic negative sign operation
inline lx_interval operator-(const lx_interval &) throw();

//! Implementation of standard algebraic addition operation
lx_interval operator + (const lx_interval &, const lx_interval &) throw();

//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const lx_interval &, const l_interval &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const l_interval &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const lx_interval &, const l_real &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const l_real &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const lx_interval &, const lx_real &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const lx_real &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const lx_interval &, const real &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const real &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const lx_interval &, const interval &) 
                                                                 throw();
//! Implementation of standard algebraic addition operation
inline lx_interval operator + (const interval &, const lx_interval &)
                                                                 throw();

//! Implementation of standard algebraic addition and allocation operation
inline lx_interval & operator +=(lx_interval &, const lx_interval &) throw();
//! Implementation of standard algebraic addition and allocation operation
inline lx_interval & operator +=(lx_interval &, const l_interval &) throw();
//! Implementation of standard algebraic addition and allocation operation
inline lx_interval & operator +=(lx_interval &, const l_real     &) throw();
//! Implementation of standard algebraic addition and allocation operation
inline lx_interval & operator +=(lx_interval &, const lx_real     &) throw();
//! Implementation of standard algebraic addition and allocation operation
inline lx_interval & operator +=(lx_interval &, const real       &) throw();
//! Implementation of standard algebraic addition and allocation operation
inline lx_interval & operator +=(lx_interval &, const interval   &) throw();

//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_interval &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_interval &, const l_interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const l_interval &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_interval &, const l_real &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const l_real &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_interval &, const lx_real &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_real &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_interval &, const real &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const real &, const lx_interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const lx_interval &, const interval &) 
                                                                 throw();
//! Implementation of standard algebraic subtraction operation
inline lx_interval operator - (const interval &, const lx_interval &) 
                                                                 throw();

//! Implementation of standard algebraic subtraction and allocation operation
inline lx_interval & operator -=(lx_interval &, const lx_interval &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
inline lx_interval & operator -=(lx_interval &, const l_interval &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
inline lx_interval & operator -=(lx_interval &, const l_real     &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
inline lx_interval & operator -=(lx_interval &, const lx_real     &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
inline lx_interval & operator -=(lx_interval &, const real       &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
inline lx_interval & operator -=(lx_interval &, const interval   &) throw();

//! Implementation of standard algebraic multiplication operation
lx_interval operator * (const lx_interval &, const lx_interval &) 
    throw();

//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const lx_interval &, const l_interval &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const l_interval &, const lx_interval &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const lx_interval &, const l_real &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const l_real &, const lx_interval &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const lx_interval &, const lx_real &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const lx_real &, const lx_interval &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const lx_interval &, const real &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const real &, const lx_interval &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const lx_interval &, const interval &) 
    throw();
//! Implementation of standard algebraic multiplication operation
inline lx_interval operator * (const interval &, const lx_interval &) 
    throw();

//! Implementation of standard algebraic multiplication and allocation operation
inline lx_interval & operator *=(lx_interval &, const lx_interval &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
inline lx_interval & operator *=(lx_interval &, const l_interval &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
inline lx_interval & operator *=(lx_interval &, const l_real     &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
inline lx_interval & operator *=(lx_interval &, const lx_real     &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
inline lx_interval & operator *=(lx_interval &, const real       &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
inline lx_interval & operator *=(lx_interval &, const interval   &) throw();

//! Implementation of standard algebraic division operation
lx_interval operator / (const lx_interval &, const lx_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);

//! Implementation of standard algebraic division operation
inline lx_interval operator / (const lx_interval &, const l_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const l_interval &, const lx_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const lx_interval &, const l_real &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const l_real &, const lx_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const lx_interval &, const real &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const real &, const lx_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const lx_interval &, const interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const interval &, const lx_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const lx_interval &, const lx_real &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
inline lx_interval operator / (const lx_real &, const lx_interval &) 
    throw(ERROR_LINTERVAL_DIV_BY_ZERO);

//! Implementation of standard algebraic division and allocation operation
inline lx_interval & operator /=(lx_interval &, const lx_interval &) throw();
//! Implementation of standard algebraic division and allocation operation
inline lx_interval & operator /=(lx_interval &, const l_interval &) throw();
//! Implementation of standard algebraic division and allocation operation
inline lx_interval & operator /=(lx_interval &, const l_real     &) throw();
//! Implementation of standard algebraic division and allocation operation
inline lx_interval & operator /=(lx_interval &, const lx_real     &) throw();
//! Implementation of standard algebraic division and allocation operation
inline lx_interval & operator /=(lx_interval &, const real       &) throw();
//! Implementation of standard algebraic division and allocation operation
inline lx_interval & operator /=(lx_interval &, const interval   &) throw();

// ----------------------------- Convex hull -------------------------------

    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const lx_interval&, const lx_interval&) 
	throw();
    //! Returns the convex hull of the arguments 
    inline lx_interval operator | (const lx_interval&, const l_interval&) 
	throw();
    //! Returns the convex hull of the arguments 
    inline lx_interval operator | (const l_interval&, const lx_interval&) 
	throw();
    //! Returns the convex hull of the arguments 
    inline lx_interval operator | (const lx_interval&, const interval&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const interval&, const lx_interval&) 
	throw();
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_interval & operator |= (lx_interval&, const lx_interval&) 
	throw();
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_interval & operator |= (lx_interval&, const l_interval&) 
	throw();
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_interval & operator |= (lx_interval&, const interval&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const lx_real&, const lx_interval&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const real&, const lx_interval&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const lx_interval&, const lx_real&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const lx_interval&, const real&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const lx_interval&, const l_real&) 
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const l_real&, const lx_interval&) 
	throw();
    //! Allocates the convex hull of the arguments to the first argument 
    inline lx_interval & operator |= (lx_interval&, const real&)
	throw();
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_interval & operator |= (lx_interval&, const l_real&)
	throw();
    //! Allocates the convex hull of the arguments to the first argument
    inline lx_interval & operator |= (lx_interval&, const lx_real&)
	throw();
    //! Returns the convex hull of the arguments
    inline lx_interval operator | (const lx_real&, const lx_real&) 
	throw();

// --------------------------- Intersection -----------------------------

    //! Returns the intersection of the arguments
    inline lx_interval operator & (const lx_interval&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const lx_interval&, const l_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Allocates the intersection of the arguments to the first argument
    inline lx_interval & operator &= (lx_interval&, const l_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const l_interval&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Returns the intersection of the arguments 
    inline lx_interval operator & (const lx_interval&, const interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Allocates the intersection of the arguments to the first argument 
    inline lx_interval & operator &= (lx_interval &a, const interval &b) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const interval&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Allocates the intersection of the arguments to the first argument 
    inline lx_interval & operator &= (lx_interval&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const lx_interval&, const lx_real&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const lx_interval&, const l_real&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const lx_interval&, const real&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Returns the intersection of the arguments 
    inline lx_interval operator & (const lx_real&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const l_real&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Returns the intersection of the arguments
    inline lx_interval operator & (const real&, const lx_interval&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
    //! Allocates the intersection of the arguments to the first argument 
    inline lx_interval & operator &= (lx_interval&, const lx_real&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Allocates the intersection of the arguments to the first argument
    inline lx_interval & operator &= (lx_interval&, const l_real&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 
    //! Allocates the intersection of the arguments to the first argument
    inline lx_interval & operator &= (lx_interval&, const real&) 
	throw(ERROR_LINTERVAL_EMPTY_INTERVAL); 

// ------------------------- SetInf, SetSup -----------------------------

//! Returns the interval with the new given infimum value
inline lx_interval & SetInf(lx_interval&, const lx_real&) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
//! Returns the interval with the new given infimum value
inline lx_interval & SetInf(lx_interval&, const l_real&) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
//! Returns the interval with the new given infimum value
inline lx_interval & SetInf(lx_interval&, const real&) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL);

//! Returns the interval with the new given supremum value
inline lx_interval & SetSup(lx_interval&, const lx_real&) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
//! Returns the interval with the new given supremum value
inline lx_interval & SetSup(lx_interval&, const l_real&) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL);
//! Returns the interval with the new given supremum value
inline lx_interval & SetSup(lx_interval&, const real&) 
    throw(ERROR_LINTERVAL_EMPTY_INTERVAL);

// ------------------------- Help Functions: ----------------------------

// ----------------- Intervals for test calculations -------------------- 

//! Returns a point interval with maximum exponent 1020
l_interval point_max(void);
// returns a staggered point interval with maximum exponent 1020,
// whereby nearly all mantissa bits are set. 

//! Returns a point interval with exponent n
l_interval point_any(int n);
// returns a staggered point interval with exponent n,
// whereby nearly all mantissa bits are set.
// -1074 <= n <= +1020; 

//! Returns a wide interval with maximum exponent 1020
l_interval wide_max(void);
// returns a staggered interval a with maximum exponent 1020
// and diam(a)>0, whereby nearly all mantissa bits are set. 

//! Returns a wide interval with exponent n, \f$ -1074\le n \le +1020; \f$
l_interval wide_any(int n);
// returns a wide staggered interval a with exponent n,
// and diam(a)>0, whereby nearly all mantissa bits are set.
// -1074 <= n <= +1020;


// --------------------------------------------------------------------
// ----------------- lx_interval Constants ----------------------------
// --------------------------------------------------------------------

 //! Enclosure-Interval for \f$ \pi \f$
 lx_interval Pi_lx_interval() throw(); // pi
	 //! Enclosure-Interval for \f$ \pi^2 \f$
 lx_interval Pip2_lx_interval() throw(); // pi^2
	 //! Enclosure-Interval for \f$ 2\pi \f$
 lx_interval Pi2_lx_interval() throw(); // 2*pi
	 //! Enclosure-Interval for \f$ \frac{1}{2\pi} \f$
 lx_interval Pi2r_lx_interval() throw(); // 1/(2*pi)
	 //! Enclosure-Interval for \f$ \pi/4 \f$
 lx_interval Pid4_lx_interval() throw(); // pi/4
	 //! Enclosure-Interval for \f$ \pi/2 \f$
 lx_interval Pid2_lx_interval() throw(); // pi/2
    //! Enclosure-Interval for \f$ \ln(2) \f$
 lx_interval Ln2_lx_interval() throw(); // ln(2)
    //! Enclosure-Interval for \f$ \ln(10) \f$
 lx_interval Ln10_lx_interval() throw(); // ln(10)
    //! Enclosure-Interval for \f$ 1/\ln(10) \f$
 lx_interval Ln10r_lx_interval() throw(); // 1/ln(10)
    //! Enclosure-Interval for \f$ 1/\pi \f$
 lx_interval Pir_lx_interval() throw();  // 1/pi
    //! Enclosure-Interval for \f$ \sqrt{\pi} \f$
 lx_interval SqrtPi_lx_interval() throw();  // sqrt(pi)
    //! Enclosure-Interval for \f$ \sqrt{2\pi} \f$
 lx_interval Sqrt2Pi_lx_interval() throw(); // sqrt(2pi)
    //! Enclosure-Interval for \f$ \sqrt{2} \f$
 lx_interval Sqrt2_lx_interval() throw();  // sqrt(2)
	 //! Enclosure-Interval for \f$ \frac{1}{\sqrt{2}} \f$
 lx_interval Sqrt2r_lx_interval() throw();  // sqrt(2)
    //! Enclosure-Interval for \f$ \sqrt{3} \f$
 lx_interval Sqrt3_lx_interval() throw();  // sqrt(3)
	 //! Enclosure-Interval for \f$ \frac{\sqrt{3}}{2} \f$
 lx_interval Sqrt3d2_lx_interval() throw();  // sqrt(3)/2
	 //! Enclosure-Interval for \f$ \frac{1}{\sqrt{3}} \f$
 lx_interval Sqrt3r_lx_interval() throw();  // 1/sqrt(3)
    //! Enclosure-Interval for \f$ 1/\ln(2) \f$
 lx_interval Ln2r_lx_interval() throw();  // 1/ln(2)
    //! Enclosure-Interval for \f$ \pi/3 \f$
 lx_interval Pid3_lx_interval() throw();  // pi/3
    //! Enclosure-Interval for \f$ 1/\sqrt{\pi} \f$
 lx_interval SqrtPir_lx_interval() throw();  // 1/sqrt(pi)
    //! Enclosure-Interval for \f$ 1/\sqrt{2\pi} \f$
 lx_interval Sqrt2Pir_lx_interval() throw(); // 1/sqrt(2pi)
    //! Enclosure-Interval for \f$ \ln(\pi) \f$
 lx_interval LnPi_lx_interval() throw();  // ln(pi)
    //! Enclosure-Interval for \f$ \ln(2\pi) \f$
 lx_interval Ln2Pi_lx_interval() throw();  // ln(2pi)
    //! Enclosure-Interval for \f$ e=2.718... \f$
 lx_interval E_lx_interval() throw();  // e
	 //! Enclosure-Interval for \f$ e^2 \f$
 lx_interval Ep2_lx_interval() throw();  // e^2
	 //! Enclosure-Interval for \f$ e^{-2} \f$
 lx_interval Ep2r_lx_interval() throw();  // 1/e^2
	 //! Enclosure-Interval for \f$ \frac{1}{e} \f$
 lx_interval Er_lx_interval() throw();  // 1/e
    //! Enclosure-Interval for \f$ e^{\pi} \f$
 lx_interval EpPi_lx_interval() throw();  // e^pi
	 //! Enclosure-Interval for \f$ e^{\pi/2} \f$
 lx_interval EpPid2_lx_interval() throw();  // e^(pi/2)
	 //! Enclosure-Interval for \f$ e^{\pi/4} \f$
 lx_interval EpPid4_lx_interval() throw();  // e^(pi/4)
	 //! Enclosure-Interval for \f$ e^{2\pi} \f$
 lx_interval Ep2Pi_lx_interval() throw();  // e^(2*pi)
    //! Enclosure-Interval for \f$ \mbox{EulerGamma}=0.5772... \f$
 lx_interval EulerGamma_lx_interval() throw();
    //! Enclosure-Interval for \f$ \mbox{Catalan}=0.9159... \f$
 lx_interval Catalan_lx_interval() throw();
    //! Enclosure-Interval for \f$ \sqrt{5} \f$ 
 lx_interval sqrt5_lx_interval() throw();  // sqrt(5)
    //! Enclosure-Interval for \f$ \sqrt{7} \f$
 lx_interval sqrt7_lx_interval() throw();   // sqrt(7)
    //! Enclosure-Interval for \f$ 1-2^{-2097} \f$
 lx_interval One_m_lx_interval() throw();
    //! Enclosure-Interval for \f$ 1+2^{-2097} \f$
 lx_interval One_p_lx_interval() throw();

// -------------------------------------------------------------------------
// ---------------- lx_interval: elementary functions ----------------------
// -------------------------------------------------------------------------
 
//! Calculates \f$ \sqrt{[x]}  \f$
 lx_interval sqrt(const lx_interval&) throw();
//! Calculates \f$ [x]^2  \f$
 lx_interval sqr(const lx_interval&) throw();
//! Calculates \f$ \ln([x]) \f$
 lx_interval ln(const lx_interval &) throw();
//! Calculates \f$ \log2([x]) \f$
 lx_interval log2(const lx_interval &) throw();
//! Calculates \f$ \log10([x]) \f$
 lx_interval log10(const lx_interval &) throw();
//! Calculates \f$ \ln(1+[x]) \f$
 lx_interval lnp1(const lx_interval &) throw();
//! Calculates \f$ \exp([x]) \f$
 lx_interval exp(const lx_interval &) throw();
//! Calculates \f$ 2^{[x]} \f$
 lx_interval exp2(const lx_interval &) throw(); // 2^x
//! Calculates \f$ 10^{[x]} \f$
 lx_interval exp10(const lx_interval &) throw(); // 10^x
//! Calculates \f$ \exp([x])-1 \f$
 lx_interval expm1(const lx_interval &x) throw(); 
//! Calculates \f$ [x]^n \f$
 lx_interval power(const lx_interval &, const real &) throw();
//! Calculates \f$ [x]^{[y]} \f$
 lx_interval pow(const lx_interval &, const lx_interval &) throw();
//! Calculates \f$ (1+[x])^{[y]} \f$
 lx_interval xp1_pow_y(const lx_interval &, const lx_interval &) throw(); 
//! Calculates \f$ \sin([x]) \f$
 lx_interval sin(const lx_interval &)throw();
//! Calculates \f$ \sin(n\cdot\pi+[x]) \f$
 lx_interval sin_n(const lx_interval &x, const real& n) throw();
//! Calculates \f$ \cos([x]) \f$
 lx_interval cos(const lx_interval &) throw();
//! Calculates \f$ \cos((n+1/2)\cdot\pi+[x]) \f$
 lx_interval cos_n(const lx_interval &x, const real& n) throw();
//! Calculates \f$ \tan([x]) \f$
 lx_interval tan(const lx_interval &) throw();
//! Calculates \f$ \cot([x]) \f$
 lx_interval cot(const lx_interval &) throw();
//! Calculates \f$ \sqrt{1+[x]^2} \f$
 lx_interval sqrt1px2(const lx_interval &) throw();
//! Calculates \f$ \arctan([x]) \f$
 lx_interval atan(const lx_interval &) throw();
//! Calculates \f$ \sqrt{1-[x]^2} \f$
 lx_interval sqrt1mx2(const lx_interval &) throw();
//! Calculates \f$ \sqrt{[x]^2-1} \f$
 lx_interval sqrtx2m1(const lx_interval &) throw();
//! Calculates \f$ \arcsin([x]) \f$
 lx_interval asin(const lx_interval & ) throw();
//! Calculates \f$ \arccos([x]) \f$
 lx_interval acos(const lx_interval &) throw();
//! Calculates \f$ \mbox{arccot}([x]) \f$
 lx_interval acot(const lx_interval &) throw();
//! Calculates \f$ \sinh([x]) \f$
 lx_interval sinh(const lx_interval &) throw();
//! Calculates \f$ \cosh([x]) \f$
 lx_interval cosh(const lx_interval &) throw();
//! Calculates \f$ \tanh([x]) \f$
 lx_interval tanh(const lx_interval &) throw();
//! Calculates \f$ \coth([x]) \f$
 lx_interval coth(const lx_interval &) throw();
//! Calculates \f$ \sqrt{([x]+1)-1} \f$
 lx_interval sqrtp1m1(const lx_interval &) throw();
//! Calculates \f$ \mbox{arcsinh}([x]) \f$
 lx_interval asinh(const lx_interval &) throw();
//! Calculates \f$ \mbox{arccosh}([x]) \f$
 lx_interval acosh(const lx_interval &) throw();
//! Calculates \f$ \mbox{arccosh}(1+[x]) \f$
 lx_interval acoshp1(const lx_interval &) throw();
//! Calculates \f$ \mbox{arctanh}([x]) \f$
 lx_interval atanh(const lx_interval &) throw();
//! Calculates \f$ \mbox{arctanh}(1-[x]) \f$
 lx_interval atanh1m(const lx_interval &) throw();
//! Calculates \f$ \mbox{arctanh}(-1+[x]) \f$
 lx_interval atanhm1p(const lx_interval &) throw();
//! Calculates \f$ \mbox{arccoth}([x]) \f$
 lx_interval acoth(const lx_interval &) throw();
//! Calculates \f$ \mbox{arccoth}(+1+[x]) \f$
 lx_interval acothp1(const lx_interval &) throw();
//! Calculates \f$ \mbox{arctanh}(-1-[x]) \f$
 lx_interval acothm1m(const lx_interval &) throw();
//! Calculates \f$ \sqrt{[x]^2 + [y]^2} \f$
 lx_interval sqrtx2y2(const lx_interval &, const lx_interval &) throw();
//! Calculates \f$ \ln(\sqrt{[x]^2 + [y]^2}) \f$
 lx_interval ln_sqrtx2y2(const lx_interval &, const lx_interval &) throw();
//! Calculates \f$ \sqrt[n]{[x]} \f$
 lx_interval sqrt(const lx_interval &, int) throw();

} // end namespace cxsc

#include "lx_interval.inl"

#endif // _CXSC_LX_INTERVAL_HPP_INCLUDED
