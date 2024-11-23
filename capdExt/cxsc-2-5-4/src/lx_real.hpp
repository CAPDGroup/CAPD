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

/* CVS $Id: lx_real.hpp,v 1.10 2014/01/30 17:23:47 cxsc Exp $ */


/*
**  F. Blomquist, University of Wuppertal, 19.09.2007;
*/

/*
**  Implementation of the classes
**
**  lx_real      with all tools and elementary functions for real
**               point and interval aruments
**
*/

#ifndef _CXSC_LX_REAL_HPP_INCLUDED
#define _CXSC_LX_REAL_HPP_INCLUDED

#include <l_imath.hpp>
#include <sstream>
#include <cmath>
#include <iostream>

namespace cxsc {

	const real Max_Int_R =  9007199254740991.0; // = 2^(53) - 1
	const real Max_Int_N = -9007199254738891.0; // = -Max_Int_R+2100

   //! Returns 1 if x is an integer value and if \f$ |x|\le 2^{53}-1 \f$ 
	inline bool Is_Integer(const real& x);
    // returns 1 if x is an integer value and 
    //           if |x| <= 2^53 - 1 = 9007199254740991.0;
    // otherwise 0 is returnd

	class lx_real
	{

		private:
    // ------------- Datenelemente -------------------------------------------
			real ex;
			l_real lr;

		public:
    // ------------- Constructors --------------------------------------------

    //! Constructor of class lx_real
			lx_real(void)  throw() {}
    //! Constructor of class lx_real
			lx_real(const real& n, const l_real &a) throw() 
			{ 
				if ( !(Is_Integer(n)) ) 
					cxscthrow(REAL_NOT_ALLOWED("lx_real(const real&, const l_real&)"));
				else
				{
					ex = n; lr = a;
				}
			}

    //! Constructor of class lx_real
			lx_real(const real& n, const real &a) throw() 
			{ 
				if ( !(Is_Integer(n)) ) 
					cxscthrow(REAL_NOT_ALLOWED("lx_real(const real&, const real&)"));
				else
				{
					ex = n; lr = a;
				}
			}

    //! Constructor of class lx_real
			explicit lx_real(const l_real &a) throw() : ex(0), lr(a) { }
    //! Constructor of class lx_real
			explicit lx_real(const real &a)   throw() : ex(0), lr(a) { }
	 //! Constructor of class lx_real
			lx_real(const real&, const string&) throw();


    // ------------- Assignments ---------------------------------------------
	 
    //! Implementation of standard assigning operator
			inline lx_real & operator = (const lx_real &) throw();
	 //! Implementation of standard assigning operator
			inline lx_real & operator = (const l_real &)  throw();
	 //! Implementation of standard assigning operator
			inline lx_real & operator = (const real &)    throw();

    // ------------- Functions -----------------------------------------------

	 //! Returns the precision of a lx_real value
			friend inline int StagPrec(const lx_real&) throw(); 
	 //! Returns the exponent a.ex of a lx_real value a
			friend inline real expo(const lx_real&) throw();
	 //! Returns the sign of a lx_real value
			friend inline int sign(const lx_real&) throw();
	 //! Returns the l_real component a.lr of a lx_real value a
			friend inline l_real lr_part(const lx_real&) throw();
	 //! Returns the absolute value of a lx_real value a
			friend inline lx_real abs(const lx_real&) throw();
	 //! Scaling a.lr upwards and keeps the value of a unchanged
			friend void scale_up  (lx_real &) throw();
	 //! Scaling a.lr downwards and keeps the value of a unchanged
			friend void scale_down(lx_real &) throw();
	 //! matches the precision of an lx_real value to the actual stagprec value
			friend inline lx_real adjust(const lx_real &) throw();

	 //! Returns 1 if \f$ a=0 \f$ and 0 otherwise
			friend inline bool eq_zero (const lx_real &a) throw(); // a = 0;
	 //! Returns 1 if \f$ a>0 \f$ and 0 otherwise
			friend inline bool gr_zero (const lx_real &a) throw(); // a > 0;
	 //! Returns 1 if \f$ a\ge0 \f$ and 0 otherwise
			friend inline bool ge_zero (const lx_real &a) throw(); // a >=0;
	 //! Returns 1 if \f$ a<0 \f$ and 0 otherwise
			friend inline bool sm_zero (const lx_real &a) throw(); // a < 0;
	 //! Returns 1 if \f$ a\le0 \f$ and 0 otherwise
			friend inline bool se_zero (const lx_real &a) throw(); // a <=0;
	 
	 //! Multiplication with \f$ 2^r; \f$ the real value r must be an integer
			friend inline void times2pown(lx_real &, const real &) throw();
	 //! Multiplication with \f$ 2^r; \f$ the real value r must be a negative integer
			friend inline void times2pown_neg(lx_real& a, const real& n) throw();
	 //! Implementation of standard negation operation
			friend inline bool operator ! (const lx_real &) throw();
    //! Implementation of standard algebraic negative sign operation
			friend inline lx_real operator -(const lx_real &) throw();


// ----------------------- Output --------------------------------------------

//friend inline std::ostream& operator << (std::ostream& s, const lx_real& a) 
//    throw(); // A value a of type lx_real is written to the output channel.
// The above operator is declared and defined in
// lx_interval.hpp (outside the class lx_interval) , lx_interval.inl;
//! Implementation of standard output method
friend std::string & operator << (std::string &s, const lx_real& a) 
					throw();
// The value of a variable a of type lx_real is copied to a string s.
// s has the form:  {2**(ex)*li};  ex is the exponent to base 2.

}; // end of class lx_real

// -------------------------------------------------------------
// ------- friend functions declared in class lx_real: ---------
// -------------------------------------------------------------

inline int StagPrec(const lx_real&) throw();
inline real expo(const lx_real&) throw();
inline int sign(const lx_real&) throw();
inline l_real lr_part(const lx_real&) throw();
inline lx_real abs(const lx_real&) throw();
inline lx_real adjust(const lx_real &) throw();
inline void times2pown_neg(lx_real&, const real&) throw();

       void scale_up  (lx_real&)   throw();
		 void scale_down(lx_real &a) throw();

		 inline bool eq_zero (const lx_real &a) throw(); // a = 0;
		 inline bool gr_zero (const lx_real &a) throw(); // a > 0;
		 inline bool ge_zero (const lx_real &a) throw(); // a >=0;
		 inline bool sm_zero (const lx_real &a) throw(); // a < 0;
		 inline bool se_zero (const lx_real &a) throw(); // a <=0;

		 inline void times2pown(lx_real &, const real &) throw();
		 inline bool operator ! (const lx_real &) throw();

		 inline lx_real operator -(const lx_real &) throw();

// -------------------------- Input ------------------------------------

//! Implementation of standard input method
		 string & operator >> (string &s, lx_real &a) throw();
//! Implementation of standard input method
		 void operator >> (const string &s, lx_real &a) throw();
//! Implementation of standard input method
		 void operator >> (const char *s, lx_real&) throw();
//! Implementation of standard input method
		 std::istream & operator >> (std::istream &s, lx_real &a) throw();

// -------------------------- Output ------------------------------------

   std::string & operator << (std::string &s, const lx_real& a) 
					throw();
//inline std::ostream& operator << (std::ostream& s, const lx_real& a) 
//    throw(); // A value a of type lx_real is written to the output channel.
//! Implementation of standard output method
		 inline std::string & operator << (std::string &s, const lx_real& a) 
				 throw();
// The value of a variable a of type lx_real is copied to a string s.
// s has the form:  {2**(ex),lr}

// ---- function declarations outside the class lx_real ----

//! Returns \f$ a+b; \f$ a,b must be integers with \f$ |a|,|b|\le2^{53}. \f$
		 inline real add_real(const real &a, const real &b) throw();

//! Returns \f$ a-b; \f$ a,b must be integers with \f$ |a|,|b|\le2^{53}. \f$
		 inline real sub_real(const real &a, const real &b) throw();

//! Returns a rather small upper bound of x.
		 lx_real upper_bnd(const lx_real& x) throw();

//! Returns a rather great lower bound of x.
		 lx_real lower_bnd(const lx_real& x) throw();
//! Implementation of standard algebraic positive sign operation
		 inline lx_real operator +(const lx_real &) throw();

//! Implementation of standard algebraic addition operation
		 lx_real operator + (const lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic addition operation
		 inline lx_real operator + (const lx_real&, const l_real &) throw();
//! Implementation of standard algebraic addition operation
		 inline lx_real operator + (const l_real&, const lx_real &) throw();
//! Implementation of standard algebraic addition operation
		 inline lx_real operator + (const lx_real&, const real &) throw();
//! Implementation of standard algebraic addition operation
		 inline lx_real operator + (const real&, const lx_real &) throw();

//! Implementation of standard algebraic addition and allocation operation
		 inline lx_real & operator +=(lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic addition and allocation operation
		 inline lx_real & operator +=(lx_real &, const l_real &) throw();
//! Implementation of standard algebraic addition and allocation operation
		 inline lx_real & operator +=(lx_real &, const real &) throw();

//! Implementation of standard algebraic subtraction operation
		 inline lx_real operator - (const lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic subtraction operation
		 inline lx_real operator - (const lx_real &, const l_real &) throw();
//! Implementation of standard algebraic subtraction operation
		 inline lx_real operator - (const l_real &, const lx_real &) throw();
//! Implementation of standard algebraic subtraction operation
		 inline lx_real operator - (const lx_real &, const real &) throw();
//! Implementation of standard algebraic subtraction operation
		 inline lx_real operator - (const real &, const lx_real &) throw();

//! Implementation of standard algebraic subtraction and allocation operation
		 inline lx_real & operator -=(lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
		 inline lx_real & operator -=(lx_real &, const l_real &) throw();
//! Implementation of standard algebraic subtraction and allocation operation
		 inline lx_real & operator -=(lx_real &, const real &) throw();

//! Implementation of standard algebraic multiplication operation
		 lx_real operator * (const lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic multiplication operation
		 inline lx_real operator * (const lx_real&, const l_real &) throw();
//! Implementation of standard algebraic multiplication operation
		 inline lx_real operator * (const l_real&, const lx_real &) throw();
//! Implementation of standard algebraic multiplication operation
		 inline lx_real operator * (const lx_real&, const real &) throw();
//! Implementation of standard algebraic multiplication operation
		 inline lx_real operator * (const real&, const lx_real &) throw();

//! Implementation of standard algebraic multiplication and allocation operation
		 inline lx_real & operator *=(lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
		 inline lx_real & operator *=(lx_real &, const l_real &) throw();
//! Implementation of standard algebraic multiplication and allocation operation
		 inline lx_real & operator *=(lx_real &, const real &) throw();

//! Implementation of standard algebraic division operation
		 lx_real operator / (const lx_real &, const lx_real &) throw(DIV_BY_ZERO);
//! Implementation of standard algebraic division operation
		 inline lx_real operator / (const lx_real&, const l_real &) throw();
//! Implementation of standard algebraic division operation
		 inline lx_real operator / (const l_real&, const lx_real &) throw();
//! Implementation of standard algebraic division operation
		 inline lx_real operator / (const lx_real&, const real &) throw();
//! Implementation of standard algebraic division operation
		 inline lx_real operator / (const real&, const lx_real &) throw();

//! Implementation of standard algebraic division and allocation operation
		 inline lx_real & operator /=(lx_real &, const lx_real &) throw();
//! Implementation of standard algebraic division and allocation operation
		 inline lx_real & operator /=(lx_real &, const l_real &) throw();
//! Implementation of standard algebraic division and allocation operation
		 inline lx_real & operator /=(lx_real &, const real &) throw();

//! Implementation of standard equality operation
		 bool operator == (const lx_real &, const lx_real &) throw();
//! Implementation of standard equality operation
		 inline bool operator == (const lx_real &, const l_real &) throw();
//! Implementation of standard equality operation
		 inline bool operator == (const l_real &, const lx_real &) throw();
//! Implementation of standard equality operation
		 inline bool operator == (const lx_real &, const real &)   throw();
//! Implementation of standard equality operation
		 inline bool operator == (const real &,   const lx_real &) throw();

//! Implementation of standard negated equality operation
		 inline bool operator != (const lx_real &, const lx_real &) throw();
//! Implementation of standard negated equality operation
		 inline bool operator != (const lx_real &, const l_real &) throw();
//! Implementation of standard negated equality operation
		 inline bool operator != (const l_real &, const lx_real &) throw();
//! Implementation of standard negated equality operation
		 inline bool operator != (const lx_real &, const real &)   throw();
//! Implementation of standard negated equality operation
		 inline bool operator != (const real &,   const lx_real &) throw();

//! Implementation of standard greater-than operation
		 bool operator > (const lx_real &, const lx_real &) throw();

//! Implementation of standard less-or-equal-than operation
		 inline bool operator <= (const lx_real &, const lx_real &) throw();
//! Implementation of standard less-than operation
		 inline bool operator <  (const lx_real &, const lx_real &) throw();
//! Implementation of standard greater-or-equal-than operation
		 inline bool operator >= (const lx_real &, const lx_real &) throw();

//! Implementation of standard greater-than operation
		 inline bool operator >  (const real &, const lx_real &) throw();
//! Implementation of standard less-or-equal-than operation
		 inline bool operator <= (const real &, const lx_real &) throw();
//! Implementation of standard less-than operation
		 inline bool operator <  (const real &, const lx_real &) throw();
//! Implementation of standard greater-or-equal-than operation
		 inline bool operator >= (const real &, const lx_real &) throw();

//! Implementation of standard greater-than operation
		 inline bool operator >  (const lx_real &, const real &) throw();
//! Implementation of standard less-or-equal-than operation
		 inline bool operator <= (const lx_real &, const real &) throw();
//! Implementation of standard less-than operation
		 inline bool operator <  (const lx_real &, const real &) throw();
//! Implementation of standard greater-or-equal-than operation
		 inline bool operator >= (const lx_real &, const real &) throw();

//! Implementation of standard greater-than operation
		 inline bool operator >  (const l_real &, const lx_real &) throw();
//! Implementation of standard less-or-equal-than operation
		 inline bool operator <= (const l_real &, const lx_real &) throw();
//! Implementation of standard less-than operation
		 inline bool operator <  (const l_real &, const lx_real &) throw();
//! Implementation of standard greater-or-equal-than operation
		 inline bool operator >= (const l_real &, const lx_real &) throw();

//! Implementation of standard greater-than operation
		 inline bool operator >  (const lx_real &, const l_real &) throw();
//! Implementation of standard less-or-equal-than operation
		 inline bool operator <= (const lx_real &, const l_real &) throw();
//! Implementation of standard less-than operation
		 inline bool operator <  (const lx_real &, const l_real &) throw();
//! Implementation of standard greater-or-equal-than operation
		 inline bool operator >= (const lx_real &, const l_real &) throw();

//! Calculating the maximum of two lx_real values
		 inline lx_real max(const lx_real&, const lx_real&);
//! Calculating the minimum of two lx_real values
		 inline lx_real min(const lx_real&, const lx_real&);

//! Returns the truncated integer part of x.
		 inline real cutint(const real& x) throw();

// ------------------- lx_real Constants ------------------------------------

    //! lx_real approximation for \f$ \pi \f$
		 lx_real Pi_lx_real() throw(); // pi
	 //! lx_real approximation for \f$ \pi^2 \f$
		 lx_real Pip2_lx_real() throw(); // pi^2
	 //! lx_real approximation for \f$ 2\pi \f$
		 lx_real Pi2_lx_real() throw(); // 2*pi
	 //! lx_real approximation for \f$ \frac{\pi}{4} \f$
		 lx_real Pid4_lx_real() throw(); // Pi/4
	 //! lx_real approximation for \f$ \frac{\pi}{2} \f$
		 lx_real Pid2_lx_real() throw(); // Pi/2
    //! lx_real approximation for \f$ \ln(2) \f$
		 lx_real Ln2_lx_real() throw();
    //! lx_real approximation for \f$ \ln(10) \f$
		 lx_real Ln10_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{1}{\ln(10)} \f$
		 lx_real Ln10r_lx_real() throw();
    //! lx_real approximation for \f$ \frac{1}{\pi} \f$
		 lx_real Pir_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{1}{2\pi} \f$
		 lx_real Pi2r_lx_real() throw();  // 1/(2*pi)
    //! lx_real approximation for \f$ \sqrt{\pi} \f$
		 lx_real SqrtPi_lx_real() throw();
    //! lx_real approximation for \f$ \sqrt{2\pi} \f$
		 lx_real Sqrt2Pi_lx_real() throw();
    //! lx_real approximation for \f$ \sqrt{2} \f$
		 lx_real Sqrt2_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{1}{\sqrt{2}} \f$
		 lx_real Sqrt2r_lx_real() throw();
    //! lx_real approximation for \f$ \sqrt{3} \f$
		 lx_real Sqrt3_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{1}{\sqrt{3}} \f$
		 lx_real Sqrt3r_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{\sqrt{3}}{2} \f$
		 lx_real Sqrt3d2_lx_real() throw();
    //! lx_real approximation for \f$ 1/\ln(2) \f$
		 lx_real Ln2r_lx_real() throw();
    //! lx_real approximation for \f$ \pi/3 \f$
		 lx_real Pid3_lx_real() throw();
    //! lx_real approximation for \f$ 1/\sqrt{\pi} \f$
		 lx_real SqrtPir_lx_real() throw();
    //! lx_real approximation for \f$ 1/\sqrt{2\pi} \f$
		 lx_real Sqrt2Pir_lx_real() throw();
    //! lx_real approximation for \f$ \ln(\pi) \f$
		 lx_real LnPi_lx_real() throw();
    //! lx_real approximation for \f$ \ln(2\pi) \f$
		 lx_real Ln2Pi_lx_real() throw();
    //! lx_real approximation for \f$ e=2.718... \f$
		 lx_real E_lx_real() throw();
	 //! lx_real approximation for \f$ e^2 \f$
		 lx_real Ep2_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{1}{e^2} \f$
		 lx_real Ep2r_lx_real() throw();
	 //! lx_real approximation for \f$ \frac{1}{e} \f$
		 lx_real Er_lx_real() throw();
    //! lx_real approximation for \f$ e^{\pi} \f$
		 lx_real EpPi_lx_real() throw();
	 //! lx_real approximation for \f$ e^{\pi/2} \f$
		 lx_real EpPid2_lx_real() throw();
	 //! lx_real approximation for \f$ e^{\pi/4} \f$
		 lx_real EpPid4_lx_real() throw();
	 //! lx_real approximation for \f$ e^{2\pi} \f$
		 lx_real Ep2Pi_lx_real() throw();
    //! lx_real approximation for \f$ \mbox{EulerGamma}=0.5772... \f$
		 lx_real EulerGamma_lx_real() throw();
    //! lx_real approximation for \f$ \mbox{Catalan}=0.9159... \f$
		 lx_real Catalan_lx_real() throw();
    //! lx_real approximation for \f$ \sqrt{5} \f$
		 lx_real sqrt5_lx_real() throw();
    //! lx_real approximation for \f$ \sqrt{7} \f$
		 lx_real sqrt7_lx_real() throw();
    //! lx_real approximation for \f$ 1-2^{-2097} \f$
		 lx_real One_m_lx_real() throw();
    //! lx_real approximation for \f$ 1+2^{-2097} \f$
		 lx_real One_p_lx_real() throw();

// **********************************************************************
// **********************************************************************



// ------------------- Array of constants ----------------------

// const real ln_N[180];
// ln_N[0] = ln(2); ln_N[1] = ln(3); ... ln_N[179] = ln(181);

const real ln_N[180] =
{6243314768165359.0 / 9007199254740992.0,
   4947709893870347.0 / 4503599627370496.0,
	6243314768165359.0 / 4503599627370496.0,
   7248263982714163.0 / 4503599627370496.0,
	8069367277953026.0 / 4503599627370496.0,
   8763600222181975.0 / 4503599627370496.0,
	4682486076124019.0 / 2251799813685248.0,
   4947709893870347.0 / 2251799813685248.0,
	5184960683398422.0 / 2251799813685248.0,
   5399580128524108.0 / 2251799813685248.0,
   5595512331017853.0 / 2251799813685248.0,
   5775752485243985.0 / 2251799813685248.0,
   5942628803132327.0 / 2251799813685248.0,
   6097986938292255.0 / 2251799813685248.0,
   6243314768165359.0 / 2251799813685248.0,
   6379829280276346.0 / 2251799813685248.0,
   6508538585911686.0 / 2251799813685248.0,
   6630287144694572.0 / 2251799813685248.0,
   6745789375439761.0 / 2251799813685248.0,
   6855655058026161.0 / 2251799813685248.0,
   6960408820565448.0 / 2251799813685248.0,
   7060505291240432.0 / 2251799813685248.0,
   7156341023059193.0 / 2251799813685248.0,
   7248263982714163.0 / 2251799813685248.0,
   7336581177285325.0 / 2251799813685248.0,
   7421564840805520.0 / 2251799813685248.0,
   7503457495173667.0 / 2251799813685248.0,
   7582476122586655.0 / 2251799813685248.0,
   7658815630333595.0 / 2251799813685248.0,
   7732651747257178.0 / 2251799813685248.0,
   7804143460206699.0 / 2251799813685248.0,
   7873435075459281.0 / 2251799813685248.0,
   7940657972317686.0 / 2251799813685248.0,
   8005932102448069.0 / 2251799813685248.0,
   8069367277953026.0 / 2251799813685248.0,
   8131064282924989.0 / 2251799813685248.0,
   8191115836735912.0 / 2251799813685248.0,
   8249607432179158.0 / 2251799813685248.0,
   8306618067481101.0 / 2251799813685248.0,
   8362220887911502.0 / 2251799813685248.0,
   8416483750067501.0 / 2251799813685248.0,
   8469469719751697.0 / 2251799813685248.0,
   8521237512606787.0 / 2251799813685248.0,
   8571841885227428.0 / 2251799813685248.0,
   8621333983281772.0 / 2251799813685248.0,
   8669761652191414.0 / 2251799813685248.0,
   8717169715100533.0 / 2251799813685248.0,
   8763600222181975.0 / 2251799813685248.0,
   8809092674755503.0 / 2251799813685248.0,
   8853684227211519.0 / 2251799813685248.0,
   8897409869326665.0 / 2251799813685248.0,
   8940302591212715.0 / 2251799813685248.0,
   8982393532846860.0 / 2251799813685248.0,
   4511856059940595.0 / 1125899906842624.0,
	4532143093607504.0 / 1125899906842624.0,
   4552071045814873.0 / 1125899906842624.0,
	4571652407313997.0 / 1125899906842624.0,
 4590899028240761.0 / 1125899906842624.0,
 4609822161187467.0 / 1125899906842624.0,
 4628432500714509.0 / 1125899906842624.0,
 4646740219649259.0 / 1125899906842624.0,
 4664755002480667.0 / 1125899906842624.0,
 4682486076124019.0 / 1125899906842624.0,
 4699942238300533.0 / 1125899906842624.0,
 4717131883750310.0 / 1125899906842624.0,
 4734063028474157.0 / 1125899906842624.0,
 4750743332179513.0 / 1125899906842624.0,
 4767180119087803.0 / 1125899906842624.0,
 4783380397244705.0 / 1125899906842624.0,
 4799350876460745.0 / 1125899906842624.0,
 4815097984997183.0 / 1125899906842624.0,
 4830627885101030.0 / 1125899906842624.0,
 4845946487483165.0 / 1125899906842624.0,
 4861059464824668.0 / 1125899906842624.0,
 4875972264388626.0 / 1125899906842624.0,
 4890690119807548.0 / 1125899906842624.0,
 4905218062110249.0 / 1125899906842624.0,
 4919560930046323.0 / 1125899906842624.0,
 4933723379761220.0 / 1125899906842624.0,
 4947709893870347.0 / 1125899906842624.0,
 4961524789976421.0 / 1125899906842624.0,
 4975172228670594.0 / 1125899906842624.0,
 4988656221054420.0 / 1125899906842624.0,
 5001980635816714.0 / 1125899906842624.0,
 5015149205896518.0 / 1125899906842624.0,
 5028165534760914.0 / 1125899906842624.0,
 5041033102324064.0 / 1125899906842624.0,
 5053755270531830.0 / 1125899906842624.0,
 5066335288634384.0 / 1125899906842624.0,
 5078776298167486.0 / 1125899906842624.0,
 5091081337661556.0 / 1125899906842624.0,
 5103253347096176.0 / 1125899906842624.0,
 5115295172116377.0 / 1125899906842624.0,
 5127209568025827.0 / 1125899906842624.0,
 5138999203570936.0 / 1125899906842624.0,
 5150666664528888.0 / 1125899906842624.0,
 5162214457111658.0 / 1125899906842624.0,
 5173645011197227.0 / 1125899906842624.0,
 5184960683398422.0 / 1125899906842624.0,
 5196163759979057.0 / 1125899906842624.0,
 5207256459626429.0 / 1125899906842624.0,
 5218240936088556.0 / 1125899906842624.0,
 5229119280684002.0 / 1125899906842624.0,
 5239893524691621.0 / 1125899906842624.0,
 5250565641627027.0 / 1125899906842624.0,
 5261137549412187.0 / 1125899906842624.0,
 5271611112444100.0 / 1125899906842624.0,
 5281988143568134.0 / 1125899906842624.0,
 5292270405961265.0 / 1125899906842624.0,
 5302459614930081.0 / 1125899906842624.0,
 5312557439628173.0 / 1125899906842624.0,
 5322565504697180.0 / 1125899906842624.0,
 5332485391835543.0 / 1125899906842624.0,
 5342318641298757.0 / 1125899906842624.0,
 5352066753334667.0 / 1125899906842624.0,
 5361731189557166.0 / 1125899906842624.0,
 5371313374261431.0 / 1125899906842624.0,
 5380814695683667.0 / 1125899906842624.0,
 5390236507208137.0 / 1125899906842624.0,
 5399580128524108.0 / 1125899906842624.0,
 5408846846735179.0 / 1125899906842624.0,
 5418037917423337.0 / 1125899906842624.0,
 5427154565669929.0 / 1125899906842624.0,
 5436197987035623.0 / 1125899906842624.0,
 5445169348501337.0 / 1125899906842624.0,
 5454069789371970.0 / 1125899906842624.0,
 5462900422144689.0 / 1125899906842624.0,
 5471662333343435.0 / 1125899906842624.0,
 5480356584321203.0 / 1125899906842624.0,
 5488984212031586.0 / 1125899906842624.0,
 5497546229770980.0 / 1125899906842624.0,
 5506043627892780.0 / 1125899906842624.0,
 5514477374494827.0 / 1125899906842624.0,
 5522848416081301.0 / 1125899906842624.0,
 5531157678200183.0 / 1125899906842624.0,
 5539406066057373.0 / 1125899906842624.0,
 5547594465108473.0 / 1125899906842624.0,
 5555723741629202.0 / 1125899906842624.0,
 5563794743265374.0 / 1125899906842624.0,
 5571808299563294.0 / 1125899906842624.0,
 5579765222481415.0 / 1125899906842624.0,
 5587666306884046.0 / 1125899906842624.0,
 5595512331017853.0 / 1125899906842624.0,
 5603304056971868.0 / 1125899906842624.0,
 5611042231121700.0 / 1125899906842624.0,
 5618727584558574.0 / 1125899906842624.0,
 5626360833503834.0 / 1125899906842624.0,
 5633942679709485.0 / 1125899906842624.0,
 5641473810845338.0 / 1125899906842624.0,
 5648954900873299.0 / 1125899906842624.0,
 5656386610409296.0 / 1125899906842624.0,
 5663769587073346.0 / 1125899906842624.0,
 5671104465828218.0 / 1125899906842624.0,
 5678391869307130.0 / 1125899906842624.0,
 5685632408130919.0 / 1125899906842624.0,
 5692826681215068.0 / 1125899906842624.0,
 5699975276066993.0 / 1125899906842624.0,
 5707078769073944.0 / 1125899906842624.0,
 5714137725781890.0 / 1125899906842624.0,
 5721152701165710.0 / 1125899906842624.0,
 5728124239891016.0 / 1125899906842624.0,
 5735052876567931.0 / 1125899906842624.0,
 5741939135997091.0 / 1125899906842624.0,
 5748783533408181.0 / 1125899906842624.0,
 5755586574691264.0 / 1125899906842624.0,
 5762348756621151.0 / 1125899906842624.0,
 5769070567075090.0 / 1125899906842624.0,
 5775752485243985.0 / 1125899906842624.0,
 5782394981837384.0 / 1125899906842624.0,
 5788998519282460.0 / 1125899906842624.0,
 5795563551917188.0 / 1125899906842624.0,
 5802090526177927.0 / 1125899906842624.0,
 5808579880781584.0 / 1125899906842624.0,
 5815032046902576.0 / 1125899906842624.0,
 5821447448344733.0 / 1125899906842624.0,
 5827826501708347.0 / 1125899906842624.0,
 5834169616552500.0 / 1125899906842624.0,
 5840477195552856.0 / 1125899906842624.0,
 5846749634655054.0 / 1125899906842624.0,
 5852987323223851.0 / 1125899906842624.0};

// ------------------------------------------------------------------------
// --------------- lx_real elementary functions ---------------------------
// ------------------------------------------------------------------------

//! Calculates \f$ \sqrt{[x]}  \f$
 lx_real sqrt(const lx_real&) throw();
//! Calculates \f$ [x]^2  \f$
 lx_real sqr(const lx_real&) throw();
//! Calculates \f$ \ln([x]) \f$
 lx_real ln(const lx_real &) throw();
//! Calculates \f$ \log2([x]) \f$
 lx_real log2(const lx_real &) throw();
//! Calculates \f$ \log10([x]) \f$
 lx_real log10(const lx_real &) throw();
//! Calculates \f$ \ln(1+[x]) \f$
 lx_real lnp1(const lx_real &) throw();
//! Calculates \f$ \exp([x]) \f$
 lx_real exp(const lx_real &) throw();
//! Calculates \f$ 2^{[x]} \f$
 lx_real exp2(const lx_real &) throw(); // 2^x
//! Calculates \f$ 10^{[x]} \f$
 lx_real exp10(const lx_real &) throw(); // 10^x
//! Calculates \f$ \exp([x])-1 \f$
 lx_real expm1(const lx_real &x) throw(); 
//! Calculates \f$ [x]^n \f$
 lx_real power(const lx_real &, const real &) throw();
//! Calculates \f$ [x]^{[y]} \f$
 lx_real pow(const lx_real &, const lx_real &) throw();
//! Calculates \f$ (1+[x])^{[y]} \f$
 lx_real xp1_pow_y(const lx_real &, const lx_real &) throw(); 
//! Calculates \f$ \sin([x]) \f$
 lx_real sin(const lx_real &)throw();
//! Calculates \f$ \sin(n\cdot\pi+[x]) \f$
 lx_real sin_n(const lx_real &x, const real& n) throw();
//! Calculates \f$ \cos([x]) \f$
 lx_real cos(const lx_real &) throw();
//! Calculates \f$ \cos((n+1/2)\cdot\pi+[x]) \f$
 lx_real cos_n(const lx_real &x, const real& n) throw();
//! Calculates \f$ \tan([x]) \f$
 lx_real tan(const lx_real &) throw();
//! Calculates \f$ \cot([x]) \f$
 lx_real cot(const lx_real &) throw();
//! Calculates \f$ \sqrt{1+[x]^2} \f$
 lx_real sqrt1px2(const lx_real &) throw();
//! Calculates \f$ \arctan([x]) \f$
 lx_real atan(const lx_real &) throw();
//! Calculates \f$ \sqrt{1-[x]^2} \f$
 lx_real sqrt1mx2(const lx_real &) throw();
//! Calculates \f$ \sqrt{[x]^2-1} \f$
 lx_real sqrtx2m1(const lx_real &) throw();
//! Calculates \f$ \arcsin([x]) \f$
 lx_real asin(const lx_real & ) throw();
//! Calculates \f$ \arccos([x]) \f$
 lx_real acos(const lx_real &) throw();
//! Calculates \f$ \mbox{arccot}([x]) \f$
 lx_real acot(const lx_real &) throw();
//! Calculates \f$ \sinh([x]) \f$
 lx_real sinh(const lx_real &) throw();
//! Calculates \f$ \cosh([x]) \f$
 lx_real cosh(const lx_real &) throw();
//! Calculates \f$ \tanh([x]) \f$
 lx_real tanh(const lx_real &) throw();
//! Calculates \f$ \coth([x]) \f$
 lx_real coth(const lx_real &) throw();
//! Calculates \f$ \sqrt{([x]+1)-1} \f$
 lx_real sqrtp1m1(const lx_real &) throw();
//! Calculates \f$ \mbox{arcsinh}([x]) \f$
 lx_real asinh(const lx_real &) throw();
//! Calculates \f$ \mbox{arccosh}([x]) \f$
 lx_real acosh(const lx_real &) throw();
//! Calculates \f$ \mbox{arccosh}(1+[x]) \f$
 lx_real acoshp1(const lx_real &) throw();
//! Calculates \f$ \mbox{arctanh}([x]) \f$
 lx_real atanh(const lx_real &) throw();
//! Calculates \f$ \mbox{arctanh}(1-[x]) \f$
 lx_real atanh1m(const lx_real &) throw();
//! Calculates \f$ \mbox{arctanh}(-1+[x]) \f$
 lx_real atanhm1p(const lx_real &) throw();
//! Calculates \f$ \mbox{arccoth}([x]) \f$
 lx_real acoth(const lx_real &) throw();
//! Calculates \f$ \mbox{arccoth}(+1+[x]) \f$
 lx_real acothp1(const lx_real &) throw();
//! Calculates \f$ \mbox{arctanh}(-1-[x]) \f$
 lx_real acothm1m(const lx_real &) throw();
//! Calculates \f$ \sqrt{[x]^2 + [y]^2} \f$
 lx_real sqrtx2y2(const lx_real &, const lx_real &) throw();
//! Calculates \f$ \ln(\sqrt{[x]^2 + [y]^2}) \f$
 lx_real ln_sqrtx2y2(const lx_real &, const lx_real &) throw();
//! Calculates \f$ \sqrt[n]{[x]} \f$
 lx_real sqrt(const lx_real &, int) throw();

} // end namespace cxsc

#include "lx_real.inl"

#endif // _CXSC_LX_REAL_HPP_INCLUDED
