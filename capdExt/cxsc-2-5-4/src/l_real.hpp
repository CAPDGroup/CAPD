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

/* CVS $Id: l_real.hpp,v 1.43 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_L_REAL_HPP_INCLUDED
#define _CXSC_L_REAL_HPP_INCLUDED

#include <iostream>
#include <string>
#include <cstdlib>

#include "dot.hpp"
#include "real.hpp"
#include "interval.hpp" // Blomquist 03.10.02; for _interval(t) in times2pown
#include "except.hpp"

namespace cxsc {

class l_interval;
class interval;

#ifdef CXSC_USE_TLS_PREC

#ifdef _WIN32
extern __declspec(thread) int stagprec;
#else
extern __thread int stagprec;
#endif

#else

extern int stagprec;

#endif




//! The Multiple-Precision Data Type l_real
/*!
The multiple-precision data type l_real is a variant of the scalar type real, which provides support for longer numbers, thus increasing the accuracy of the data type.

The multiple-precision arithmetic is based on the staggered correction principle briefly described below.

\section basicideas Basic Ideas of the Staggered Correction Multiple-Precision Arithmetic

A real staggered correction mulitple-precision number \f$ X \f$ is a vector \f$ \vec X = (x_i)^{n}_{1} \f$ of C-XSC type real floating-point numbers. The
value of a multiple-precision number \f$ X \f$ is the number in \f$ R \f$ that is the exact sum of the components of the vector,

\f[
X = \sum \limits_{i=1}^n x_i , x_i \in R = R(b,l,e_{min},e{max})
\f]

\sa real
*/
class l_real
{
   private:
      // ---- Data elements ----
      int   prec;
      real  *data;

   public:
      // ---- Constructors ----
      //! Constructor of class l_real
      l_real(void) throw();
      //! Constructor of class l_real
      l_real(const l_real &) throw();
      

      //! Implementation of standard assigning operator
      l_real & operator =(const l_real &) throw();
      //! Implementation of standard assigning operator
      l_real & operator =(const real &) throw();
      //! Implementation of standard assigning operator
      l_real & operator =(const dotprecision &) throw();
      //! Implementation of standard assigning operator
      l_real & operator =(const lx_real &) throw(); // Blomquist, 12.11.2008;

      // ---- Destructors    ----
      ~l_real(void) throw();

      // ---- Type casting ----
      //! Constructor of class l_real
      explicit l_real(int)          throw();
      //! Constructor of class l_real
      explicit l_real(long)          throw();
      //! Constructor of class l_real
      explicit l_real(const real &) throw();
      //! Constructor of class l_real
      explicit l_real(const dotprecision &) throw();
      //! Constructor of class l_real
      explicit l_real(const double &) throw(); // Blomquist 10.09.02. {l_real.cpp}
      
      friend real::real(const l_real &) throw();
      friend dotprecision::dotprecision(const l_real &) throw();
      friend dotprecision & dotprecision::operator =(const l_real &) throw();

//      friend real _real(const l_real & a) throw() { return real(a); }
//      friend l_real _l_real(const real & a) throw()  { return l_real(a); }
//      friend dotprecision _dotprecision(const l_real & a) throw()
//                                               { return dotprecision(a); }
//      friend l_real _l_real(const dotprecision & a) throw() 
//                                               { return l_real(a); }
      friend l_interval _unchecked_l_interval(const l_real &, 
                                              const l_real &) throw();

      // The following are defined in the specific vector, matrix-files
#if(CXSC_INDEX_CHECK) 
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rvector &)       throw (ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rvector_slice &) throw (ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rmatrix &)       throw (ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rmatrix_slice &) throw (ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rvector &)
      */
      friend INLINE real _l_real(const l_rvector &)       throw (ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rvector_slice &)
      */
      friend INLINE real _l_real(const l_rvector_slice &) throw (ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rmatrix &)
      */
      friend INLINE real _l_real(const l_rmatrix &)       throw (ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rmatrix_slice &)
      */
      friend INLINE real _l_real(const l_rmatrix_slice &) throw (ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ);
#else
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rvector &)       throw ();
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rvector_slice &) throw ();
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rmatrix &)       throw ();
      //! Constructor of class l_real
      explicit INLINE l_real(const l_rmatrix_slice &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rvector &)
      */
      friend INLINE real _l_real(const l_rvector &)       throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rvector_slice &)
      */
      friend INLINE real _l_real(const l_rvector_slice &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rmatrix &)
      */
      friend INLINE real _l_real(const l_rmatrix &)       throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa l_real(const l_rmatrix_slice &)
      */
      friend INLINE real _l_real(const l_rmatrix_slice &) throw ();
#endif


      // ---- Output functions ----
      //! Implementation of standard output method
      friend std::ostream & operator <<(std::ostream &,const l_real &) throw();
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &,l_real &)       throw();
      //! Implementation of standard output method
      friend std::string & operator <<(std::string &,const l_real &)   throw();
      //! Implementation of standard input method
      friend std::string & operator >>(std::string &,l_real &)         throw();
      //! Implementation of standard input method
      friend void          operator >>(const std::string &,l_real &)   throw();
      //! Implementation of standard input method
      friend void          operator >>(const char *,l_real &)          throw();

      // ---- Standard functions ---- (arithmetic operators)
      //! Access to the single components used to store the long data type value
      real&             operator[](int) const throw();

      //! Implementation of standard algebraic negative sign operation
      friend     l_real operator -(const l_real& lr1) throw();
      //! Implementation of standard algebraic positive sign operation
      friend     l_real operator +(const l_real& lr1) throw();

      //! Implementation of standard algebraic addition operation
      friend     l_real operator +(const l_real &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend     l_real operator -(const l_real &,const l_real &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend     l_real operator *(const l_real &,const l_real &) throw();
      //! Implementation of standard algebraic division operation
      friend     l_real operator /(const l_real &,const l_real &) throw(DIV_BY_ZERO);
      //! Returns the convex hull of the arguments
      friend inline l_interval operator |(const l_real &,const l_real &) throw();

      //! Implementation of standard algebraic addition operation
      friend     l_real operator +(const l_real &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend     l_real operator +(const real &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend     l_real operator -(const l_real &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend     l_real operator -(const real &,const l_real &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend     l_real operator *(const l_real &,const real &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend     l_real operator *(const real &,const l_real &) throw();
      //! Implementation of standard algebraic division operation
      friend     l_real operator /(const l_real &,const real &) throw();
      //! Implementation of standard algebraic division operation
      friend     l_real operator /(const real &,const l_real &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_interval operator |(const real &,const l_real &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_interval operator |(const l_real &,const real &) throw();

      //! Implementation of standard algebraic addition operation
      friend  dotprecision operator +(const l_real &,const dotprecision &) throw();
      //! Implementation of standard algebraic addition operation
      friend  dotprecision operator +(const dotprecision &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend  dotprecision operator -(const l_real &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend  dotprecision operator -(const dotprecision &,const l_real &) throw();
      //! Returns the convex hull of the arguments
      friend inline idotprecision operator |(const dotprecision &,const l_real &) throw();
      //! Returns the convex hull of the arguments
      friend inline idotprecision operator |(const l_real &,const dotprecision &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend     l_real & operator +=(l_real &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend     l_real & operator -=(l_real &,const l_real &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend     l_real & operator *=(l_real &,const l_real &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend     l_real & operator /=(l_real &,const l_real &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend     l_real & operator +=(l_real &,const real &) throw();      
      //! Implementation of standard algebraic subtraction and allocation operation
      friend     l_real & operator -=(l_real &,const real &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend     l_real & operator *=(l_real &,const real &) throw();      
      //! Implementation of standard algebraic division and allocation operation
      friend     l_real & operator /=(l_real &,const real &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend     real   & operator +=(real &,const l_real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend     real   & operator -=(real &,const l_real &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend     real   & operator *=(real &,const l_real &) throw();
	//! Implementation of standard algebraic division and allocation operation
      friend     real   & operator /=(real &,const l_real &) throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend     inline dotprecision & operator +=(dotprecision &d,const l_real &lr) throw() { lr._akku_add(d); return d; }
      //! Implementation of standard algebraic subtraction and allocation operation
      friend     inline dotprecision & operator -=(dotprecision &d,const l_real &lr) throw() { lr._akku_sub(d); return d; }

      // ---- Compare operators ----
      //! Implementation of standard equality operation
      friend bool operator ==(const l_real &,const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const l_real &,const l_real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const l_real &,const l_real &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const l_real &,const l_real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const l_real &,const l_real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const l_real &,const l_real &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const real &,const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const real &,const l_real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const real &,const l_real &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const real &,const l_real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const real &,const l_real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const real &,const l_real &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const l_real &,const real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const l_real &,const real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const l_real &,const real &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const l_real &,const real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const l_real &,const real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const l_real &,const real &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const dotprecision &,const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const dotprecision &,const l_real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const dotprecision &,const l_real &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const dotprecision &,const l_real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const dotprecision &,const l_real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const dotprecision &,const l_real &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const l_real &,const dotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const l_real &,const dotprecision &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const l_real &,const dotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const l_real &,const dotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const l_real &,const dotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const l_real &,const dotprecision &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const interval &,const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const interval &,const l_real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const interval &,const l_real &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const interval &,const l_real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const interval &,const l_real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const interval &,const l_real &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const l_real &,const interval &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const l_real &,const interval &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const l_real &,const interval &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const l_real &,const interval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const l_real &,const interval &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const l_real &,const interval &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const idotprecision &,const l_real &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const idotprecision &,const l_real &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const idotprecision &,const l_real &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const idotprecision &,const l_real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const idotprecision &,const l_real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const idotprecision &,const l_real &) throw();

      //! Implementation of standard equality operation
      friend bool operator ==(const l_real &,const idotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend bool operator !=(const l_real &,const idotprecision &) throw();
      //! Implementation of standard less-than operation
      friend bool operator  <(const l_real &,const idotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend bool operator  >(const l_real &,const idotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend bool operator <=(const l_real &,const idotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend bool operator >=(const l_real &,const idotprecision &) throw();


      //! Implementation of standard negation operation
      friend bool operator!(const l_real& lr) throw();
      
      // ---- functions ----
	//! The accurate scalar product of the last two arguments added to the value of the first argument
      friend void accumulate(dotprecision&, const real&, const l_real&) throw();
	//! The accurate scalar product of the last two arguments added to the value of the first argument
      friend void accumulate(dotprecision&, const l_real&, const real&) throw();
	//! The accurate scalar product of the last two arguments added to the value of the first argument
      friend void accumulate(dotprecision&, const l_real&, const l_real&) throw();
      
	//! The accurate scalar product of the last two arguments added to the value of the first argument
      friend void accumulate(idotprecision&,const real&, const l_real&) throw();
	//! The accurate scalar product of the last two arguments added to the value of the first argument
      friend void accumulate(idotprecision&,const l_real&,const real&) throw();
	//! The accurate scalar product of the last two arguments added to the value of the first argument
      friend void accumulate(idotprecision&,const l_real&, const l_real&) throw();

      //! The absolute value of a l_real value
      friend l_real   abs  (const l_real&) throw();
      //! The sign of a l_real value
      friend int      sign (const l_real&) throw();
      //! Returns the precision of the long datatype value
      friend int      StagPrec(const l_real&) throw();
      //! Sets the precision of a specific long datatype value
      friend l_real   adjust(const l_real&) throw();

      //! Rounds the argument up to the next l_real value
      friend l_real   rnd_up(const dotprecision&);    // Blomquist, 20.11.2006;
      //! Rounds the argument down to the next l_real value
      friend l_real   rnd_down(const dotprecision&);  // Blomquist, 20.11.2006;

    //! Returns the minimum of the expo-function
    friend int expo_sm(const l_real&); // Blomquist, 25.03.2007;
    // Calculating expo(x[k]) of the smallest |x[k]| <> 0.

    //! Returns the maximum of the expo-function
    friend int expo_gr(const l_real&); // Blomquist, 25.03.2007;
    // Calculating expo(x[k]) of the greatest |x[k]|.

      // ---- Friends      -----
      
      //! Returns the infimum of an l_interval
      friend inline l_real Inf(const l_interval &) throw();
      //! Returns the supremum of an l_interval
      friend inline l_real Sup(const l_interval &) throw();
      //! Returns the rounded middle of the l_interval
      friend        l_real mid(const l_interval &) throw();
      //! Checks if the argument is zero
      friend inline bool zero_(const l_real &) throw(); // Blomquist,27.11.02
   private:
      void _clear(int) throw(); // filling a l_real number from element int p
                                //  up to the end with zero.
      void _akku_out(const dotprecision&) throw(); // The dotprecision value is rounded to the
                        // activated l_real number in its own precision.
      void _akku_out_up(const dotprecision&) throw(); // The dotprecision value is rounded up to 
                      // the activated l_real number in its own precision.
      void _akku_out_down(const dotprecision&) throw(); // The dotprecision value is rounded down 
                      // to the activated l_real number in its own precision.
      void _akku_add(dotprecision&) const throw(); // adding the activated
                                     // l_real number to the accumulator d.
      void _akku_sub(dotprecision&) const throw(); // subtracting the 
      // activated l_real number to the accumulator d of type dotprecision.
      inline real& elem(int i) const throw() {  return data[i-1];  }
};



      // ---- Compare operators ----
      bool operator ==(const l_real &,const l_real &) throw();
      bool operator !=(const l_real &,const l_real &) throw();
      bool operator  <(const l_real &,const l_real &) throw();
      bool operator  >(const l_real &,const l_real &) throw();
      bool operator <=(const l_real &,const l_real &) throw();
      bool operator >=(const l_real &,const l_real &) throw();

      bool operator ==(const real &,const l_real &) throw();
      bool operator !=(const real &,const l_real &) throw();
      bool operator  <(const real &,const l_real &) throw();
      bool operator  >(const real &,const l_real &) throw();
      bool operator <=(const real &,const l_real &) throw();
      bool operator >=(const real &,const l_real &) throw();

      bool operator ==(const l_real &,const real &) throw();
      bool operator !=(const l_real &,const real &) throw();
      bool operator  <(const l_real &,const real &) throw();
      bool operator  >(const l_real &,const real &) throw();
      bool operator <=(const l_real &,const real &) throw();
      bool operator >=(const l_real &,const real &) throw();

      bool operator ==(const dotprecision &,const l_real &) throw();
      bool operator !=(const dotprecision &,const l_real &) throw();
      bool operator  <(const dotprecision &,const l_real &) throw();
      bool operator  >(const dotprecision &,const l_real &) throw();
      bool operator <=(const dotprecision &,const l_real &) throw();
      bool operator >=(const dotprecision &,const l_real &) throw();

      bool operator ==(const l_real &,const dotprecision &) throw();
      bool operator !=(const l_real &,const dotprecision &) throw();
      bool operator  <(const l_real &,const dotprecision &) throw();
      bool operator  >(const l_real &,const dotprecision &) throw();
      bool operator <=(const l_real &,const dotprecision &) throw();
      bool operator >=(const l_real &,const dotprecision &) throw();

      bool operator ==(const interval &,const l_real &) throw();
      bool operator !=(const interval &,const l_real &) throw();
      bool operator  <(const interval &,const l_real &) throw();
      bool operator  >(const interval &,const l_real &) throw();
      bool operator <=(const interval &,const l_real &) throw();
      bool operator >=(const interval &,const l_real &) throw();

      bool operator ==(const l_real &,const interval &) throw();
      bool operator !=(const l_real &,const interval &) throw();
      bool operator  <(const l_real &,const interval &) throw();
      bool operator  >(const l_real &,const interval &) throw();
      bool operator <=(const l_real &,const interval &) throw();
      bool operator >=(const l_real &,const interval &) throw();

      bool operator ==(const idotprecision &,const l_real &) throw();
      bool operator !=(const idotprecision &,const l_real &) throw();
      bool operator  <(const idotprecision &,const l_real &) throw();
      bool operator  >(const idotprecision &,const l_real &) throw();
      bool operator <=(const idotprecision &,const l_real &) throw();
      bool operator >=(const idotprecision &,const l_real &) throw();

      bool operator ==(const l_real &,const idotprecision &) throw();
      bool operator !=(const l_real &,const idotprecision &) throw();
      bool operator  <(const l_real &,const idotprecision &) throw();
      bool operator  >(const l_real &,const idotprecision &) throw();
      bool operator <=(const l_real &,const idotprecision &) throw();
      bool operator >=(const l_real &,const idotprecision &) throw();


      bool operator!(const l_real& lr) throw();


inline real _real(const l_real & a) throw() { return real(a); }
inline l_real _l_real(const real & a) throw()  { return l_real(a); }
inline dotprecision _dotprecision(const l_real & a) throw()
                                               { return dotprecision(a); }
inline l_real _l_real(const dotprecision & a) throw() 
                                               { return l_real(a); }
l_interval _unchecked_l_interval(const l_real &, const l_real &) throw();

//! Rounds the argument up to the next l_real value
l_real   rnd_up(const dotprecision&);    // Blomquist, 20.11.2006;
//! Rounds the argument down to the next l_real value
l_real   rnd_down(const dotprecision&);  // Blomquist, 20.11.2006;

//! Fast calculation of \f$ 2^n \f$
    inline l_real l_pow2n(const int n) throw()
{   // Fast and exact calculation of 2^n; -1074 <= n <= 1023;
    // Blomquist 01.10.02.
    return l_real( comp(0.5,n+1) );
}

//! Fast multiplication of reference parameter lr with \f$ 2^n \f$
inline void times2pown(l_real& lr, const int n) throw() // Blomquist 03.10.02
{ // lr is multiplied with 2^n; if lr[i]*2^n are all normalized, the result
  // is exact. if one of the lr[i]*2^n are denormalized, the result is not
  // exact in general.
    int k = StagPrec(lr);
    for (int i=1; i<=k; i++)
    {
	times2pown(lr[i],n);
    }
}

//! Fast multiplication of reference parameter lr with \f$ 2^n \f$
/*! 
\param lr Value to multiply with \f$ 2^n \f$ and contains the result after the operation
\param z Enclosure-Interval for the error of the result
\param n The exponent of \f$ 2^n \f$

If we denote the old value of \f$ lr \f$ with \f$ y \f$, then times2pown(lr,z,n)
delivers with \f$ lr+z \f$ an inclusion of the exact value \f$ y*2^n \f$

\f$ z=0 \f$ signals that \f$ lr = y *2^n  \f$ was exactly calculated.

\f$ n \f$ is restricted to:  \f$ -1074 <= n <= +1023 \f$
*/
inline void times2pown(l_real& lr, interval& z, const int n) throw() 
{  // Blomquist 03.10.02;

    if ( n<-1074 || n>1023 ) 
    { std::cerr << "Error in:  " 
           << "times2pown(l_real& lr, interval& z, const int n): " << std::endl
           << " -1074 <= n <= +1023 not fulfilled" << std::endl; exit(0); } 
    int k = StagPrec(lr);
    z = 0;
    real mt,t;
    real F = comp(0.5,n+1);
    for (int i=1; i<=k; i++)
    {
	mt = mant(lr[i]);
	t = lr[i];
	times2pown(lr[i],n);
	if ( mt != mant(lr[i]) ) 
	{
	    lr[i] = 0;
	    z += _interval(t) * F;
        }
    }

}

//! Fast multiplication of reference parameter a with \f$ 2^n \f$
/*!
\param a Value to multiply with \f$ 2^n \f$ and contains the result after the operation
\param z Enclosure-Interval for the error of the result
\param n The exponent of \f$ 2^n \f$

Fast multiplication of a with \f$ 2^n \f$.

\f$ n \f$ is NOT restricted to: \f$ -1074 <= n <= +1023 \f$

If \f$ z = 0 \f$ then it holds \f$ a = y*2^n \f$ (exact multiplication!).
*/
inline void Times2pown(l_real& a, interval& z, int n) throw()
// If we denote the old value of a with y then with the
// new calculated values of a and z it holds:
//  -----  a+z is an inclusion of y*2^n;  ------
// If z==0 then it holds a = y*2^n; (exact multiplication!).
// Especially, if n>=0, the multiplication with 2^n is exact
// if no oveflow occurs.
{
    int fac,rest;
    interval z1;
    z=0;
    if (n>=0)
    {   // if (n>=0) a*2^n is exactly calculated
        // if no overflow occurs:
	fac = n/1023;  rest = n%1023;
	for (int k=1; k<=fac; k++)
	    times2pown(a,1023);
	times2pown(a,rest);
    }
    else // n < 0:
	if (n<-2100)
	{
	    if(a>0) z = interval(0,minreal);
	    else 
		if (a<0) z = interval(-minreal,0);
		else z=0;
	    a = 0;
	} else // -2100<=n<0
	{
	    fac = n/-1074;  rest = n%-1074;
	    for (int k=1; k<=fac; k++)
	    {
		times2pown(a,z1,-1074);
		times2pown(z,-1074);
		z += z1;
	    }
	    times2pown(a,z1,rest);
	    times2pown(z,rest);
	    z += z1;
	}
} // void Times2pown(...)

//! Fast multiplication of reference parameter a with \f$ 2^p \f$
/*!
\param a Value to multiply with \f$ 2^p \f$ and contains the result after the operation
\param p The exponent of \f$ 2^p \f$, p of type real must be an integer value!

Fast multiplication of a with \f$ 2^n \f$.

\f$ p \f$ is NOT restricted to: \f$ -1074 <= n <= +1023 \f$
*/
inline void Times2pown(l_real& a, const real& p) throw()
// The first parameter delivers an approximation of a * 2^p;
// For p in [-2100,+2100] p must be an integer value.
// This condition is NOT tested in this function!
// For p outside [-2100,+2100] an approximation of a * 2^p is
// calculated for any p of type real, unless an overflow occurs.
// If the function is activated with the second parameter of type int,
// then the first parameter delivers approximations of a * 2^p,
// unless an overflow occurs.
// Blomquist, 04.11.2008;
{
 const int c1 = -1000000, 
           c2 = 2100,
			  c3 = 1023,
	        c4 = -1074;
 int ex(expo_gr(a)),fac,rest,n;
 double dbl;

 if (ex > c1)
 {
	 if (p>=0)
		 if (p>c2)
			 times2pown(a,c2); // Produces an error
	 else // 0 <= p <= 2100
	 {  // By too great p-values overflow is possible!
		 dbl = _double(p);
		 n = (int) dbl;
		 fac = n/c3;
		 rest = n%c3;
		 for (int k=1; k<=fac; k++)
			 times2pown(a,c3);
		 times2pown(a,rest);
	 }
	 else // p<0; No overflow or underflow!
		 if (p<-c2) a = 0.0;
	 else // -2100 <= p < 0
	 {
		 dbl = _double(p);
		 n = (int) dbl;
		 fac = n/c4;
		 rest = n%c4;
		 for (int k=1; k<=fac; k++)
			 times2pown(a,c4);
		 times2pown(a,rest);
	 }
 }
} // Times2pown(...)


inline bool zero_(const l_real& lr) throw()
{  // returns only true if all lr.elem(i) == 0; Blomquist, 27.11.02; 
    int i=1, p=StagPrec(lr);
    bool tmp = true;
    do
    {
	if (sign(lr.elem(i))!=0) tmp = false;
	i++;
    }  while(tmp && i <= p );
    return tmp;
}

// real staggered constants:
//! Approximation of \f$ \ln 2 \f$
l_real Ln2_l_real()   throw();   // ln(2) 
//! Approximation of \f$ \ln 10 \f$
l_real Ln10_l_real()  throw();   // ln(10)
//! Approximation of \f$ \frac{1}{\ln 10} \f$
l_real Ln10r_l_real() throw();   // 1/ln(10)
//! Approximation of \f$ \frac{\pi}{4} \f$
l_real Pid4_l_real()  throw();   // Pi/4
//! Approximation of \f$ \sqrt{2} \f$
l_real Sqrt2_l_real() throw();   // sqrt(2)
//! Approximation of \f$ \sqrt{5} \f$
l_real Sqrt5_l_real() throw();   // sqrt(5)
//! Approximation of \f$ \sqrt{7} \f$
l_real Sqrt7_l_real() throw();   // sqrt(7)
//! Approximation of \f$ \frac{1}{\ln 2} \f$
l_real Ln2r_l_real() throw();     // 1/ln(2)
//! Approximation of \f$ \pi \f$
l_real Pi_l_real() throw();       // Pi
//! Approximation of \f$ \frac{\pi}{2} \f$
l_real Pid2_l_real() throw();     // Pi/2
//! Approximation of \f$ 2\pi \f$
l_real Pi2_l_real() throw();      // 2*Pi
//! Approximation of \f$ \frac{\pi}{3} \f$
l_real Pid3_l_real() throw();     // Pi/3
//! Approximation of \f$ \frac{1}{\pi} \f$
l_real Pir_l_real() throw();      // 1/Pi
//! Approximation of \f$ \frac{1}{2\pi} \f$
l_real Pi2r_l_real() throw();     // 1/(2*Pi)
//! Approximation of \f$ \sqrt{\pi} \f$
l_real SqrtPi_l_real() throw();   // sqrt(Pi)
//! Approximation of \f$ \sqrt{2\pi} \f$
l_real Sqrt2Pi_l_real() throw();  // sqrt(2*Pi)
//! Approximation of \f$ \frac{1}{\sqrt{\pi}} \f$
l_real SqrtPir_l_real() throw();  // 1/sqrt(Pi)
//! Approximation of \f$ \frac{1}{\sqrt{2\pi}} \f$
l_real Sqrt2Pir_l_real() throw(); // 1/sqrt(2*Pi)
//! Approximation of \f$ 2^\pi \f$
l_real Pip2_l_real() throw();     // Pi^2
//! Approximation of \f$ \frac{1}{\sqrt{2}} \f$
l_real Sqrt2r_l_real() throw();   // 1/sqrt(2)
//! Approximation of \f$ \sqrt{3} \f$
l_real Sqrt3_l_real() throw();    // sqrt(3)
//! Approximation of \f$ \frac{\sqrt{3}}{2} \f$
l_real Sqrt3d2_l_real() throw();  // sqrt(3)/2
//! Approximation of \f$ \frac{1}{\sqrt{3}} \f$
l_real Sqrt3r_l_real() throw();   // 1/sqrt(3)
//! Approximation of \f$ \ln \pi \f$
l_real LnPi_l_real() throw();     // ln(Pi)
//! Approximation of \f$ \ln 2\pi \f$
l_real Ln2Pi_l_real() throw();    // ln(2*Pi)
//! Approximation of \f$ e \f$
l_real E_l_real() throw();        // e = exp(1)
//! Approximation of \f$ \frac{1}{e} \f$
l_real Er_l_real() throw();       // 1/e
//! Approximation of \f$ e^2 \f$
l_real Ep2_l_real() throw();      // e^2
//! Approximation of \f$ \frac{1}{e^2} \f$
l_real Ep2r_l_real() throw();     // 1/e^2
//! Approximation of \f$ e^\pi \f$
l_real EpPi_l_real() throw();     // e^Pi
//! Approximation of \f$ e^{2\pi} \f$
l_real Ep2Pi_l_real() throw();    // e^(2*Pi)
//! Approximation of \f$ e^{\frac{\pi}{2}} \f$
l_real EpPid2_l_real() throw();   // e^(Pi/2)
//! Approximation of \f$ e^{\frac{\pi}{4}} \f$
l_real EpPid4_l_real() throw();   // e^(Pi/4)
//! Approximation of Euler Gamma
l_real EulerGa_l_real() throw();  // EulerGamma
//! Approximation of Catalan Numbers
l_real Catalan_l_real() throw();  // Catalan

} // namespace cxsc 

#endif // _CXSC_L_REAL_HPP_INCLUDED

