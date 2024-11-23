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

/* CVS $Id: cinterval.hpp,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CINTERVAL_HPP_INCLUDED
#define _CXSC_CINTERVAL_HPP_INCLUDED

#include <iostream>
#include <string>

// Headerfile for cinterval.

#include "except.hpp"
#include "real.hpp"
#include "complex.hpp"
#include "interval.hpp"

namespace cxsc {

class civector;
class cimatrix;
class civector_slice;
class cimatrix_slice;

//! The Scalar Type cinterval
/*!
The data type cinterval is used to store complex intervals. They are defined as rectangles of the form \f$ z = \left[ x \right] + i \left[ y \right] \in C \f$ (real part
\f$ \left[ x \right] \in R \f$ and imaginary part \f$ \left[ y \right] \in R \f$) with sides parallel to the
axes in the complex plane, as indicated in the figure below.

\image html "cinterval.png" "Complex Interval [3.0, 4.5] + [1.0, 2.0]i"
*/
class cinterval
{
   private:
      // ---- private data ----------------------------------------
      interval  re;
      interval  im;

   public:
      // ---- implicit constructors  ------------------------------
      //! Constructor of class cinterval
      inline cinterval(void)            {}
      //! Constructor of class cinterval
      inline cinterval(const interval & a,const interval &b); 
      //! Constructor of class cinterval
      inline cinterval(const complex & a,const complex & b) ; 

      // The following are defined in the specific vector, matrix-files
#if(CXSC_INDEX_CHECK) 
      //! Constructor of class cinterval
      INLINE cinterval(const civector &)      ;
      //! Constructor of class cinterval
      INLINE cinterval(const civector_slice &);
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix &)      ;
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix_slice &);
#else
      //! Constructor of class cinterval
      INLINE cinterval(const civector &)      ;
      //! Constructor of class cinterval
      INLINE cinterval(const civector_slice &);
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix &)      ;
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix_slice &);
#endif
      // ---- explicit constructors -------------------------------

      //! Constructor of class cinterval
      explicit inline cinterval(const real     & a) ;
      //! Constructor of class cinterval
      explicit inline cinterval(const interval & a) ;
      //! Constructor of class cinterval
      explicit inline cinterval(const complex  & a) ; 
      //! Constructor of class cinterval
      explicit        cinterval(const dotprecision &) ;
      //! Constructor of class cinterval
      explicit        cinterval(const cdotprecision &);
      //! Constructor of class cinterval
      explicit        cinterval(const idotprecision &);
      //! Constructor of class cinterval
      explicit        cinterval(const cidotprecision &);
      //! Constructor of class cinterval
      explicit        cinterval(const l_cinterval&);

      // ---- assignments -----------------------------------------

      //! Implementation of standard assigning operator
      inline cinterval & operator =(const real &);
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const interval &);
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const complex &);
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const cinterval &);
      
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const dotprecision &);
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const idotprecision &);
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const cdotprecision &);
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const cidotprecision &);
      //! Implementation of standard assigning operator
             cinterval & operator = (const l_cinterval&);      
      //! Implementation of standard assigning operator
             cinterval & operator = (const lx_cinterval&);

      // ---- compatiblility typecasts ----------------------------

      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const interval &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const dotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const idotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const cdotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const cidotprecision &);
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const interval &,const interval &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const real &,const interval &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const interval &,const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const complex &,const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const real &,const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const complex &,const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _unchecked_cinterval(const complex &,const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _unchecked_cinterval(const real &,const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _unchecked_cinterval(const complex &,const real &);

      // ---- Input/Output  ---------------------------------------

      //! Implementation of standard output method
      friend std::ostream & operator <<(std::ostream &,const cinterval &);
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &,cinterval &)      ;
      //! Implementation of standard output method
      friend std::string &  operator <<(std::string &,const cinterval &) ;
      //! Implementation of standard input method
      friend std::string &  operator >>(std::string &,cinterval &)       ;
      //! Implementation of standard input method
      friend void           operator >>(const char *,cinterval &)        ;
      //! Implementation of standard input method
      friend void           operator >>(const std::string &,cinterval &) ;

      // ---- Std.Operators ---------------------------------------

      //! Implementation of standard algebraic negative sign operation
      friend inline cinterval operator -(const cinterval &);
      //! Implementation of standard algebraic positive sign operation
      friend inline cinterval operator +(const cinterval &);
      //! Implementation of standard negation operation
      friend inline bool operator!  (const cinterval & a)                   ;

      // CI-CI

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const cinterval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const cinterval &);
      //! Implementation of standard algebraic multiplication operation
      friend        cinterval operator *(const cinterval &,const cinterval &);
      //! Implementation of standard algebraic division operation
      friend        cinterval operator /(const cinterval &,const cinterval &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const cinterval &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const cinterval &,const cinterval &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const cinterval &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const cinterval &);
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const cinterval &);
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const cinterval &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const cinterval &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const cinterval &);
      
      // CI-R
      
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const real &);
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const real &,const cinterval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const real &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const real &,const cinterval &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const cinterval &,const real &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const real &,const cinterval &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const cinterval &,const real &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const real &,const cinterval &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const real &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const real &,const cinterval &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const cinterval &,const real &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const real &,const cinterval &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const real &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const real &);
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const real &);
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const real &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const real &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const real &);
      
      // CI-I

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const interval &);
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const interval &,const cinterval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const interval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const interval &,const cinterval &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const cinterval &,const interval &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const interval &,const cinterval &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const cinterval &,const interval &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const interval &,const cinterval &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const interval &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const interval &,const cinterval &);
      //! Returns the intersection cinterval of the arguments
      friend inline cinterval operator &(const cinterval &,const interval &);
      //! Returns the intersection cinterval of the arguments
      friend inline cinterval operator &(const interval &,const cinterval &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const interval &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const interval &);
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const interval &);
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const interval &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const interval &);
      //! Allocates the intersection cinterval of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const interval &);

      // CI-C

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const complex &);
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const complex &,const cinterval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const complex &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const complex &,const cinterval &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const cinterval &,const complex &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const complex &,const cinterval &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const cinterval &,const complex &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const complex &,const cinterval &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const complex &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const complex &,const cinterval &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const cinterval &,const complex &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const complex &,const cinterval &);
      

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const complex &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const complex &);
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const complex &);
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const complex &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const complex &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const complex &);
      
      // C-R

      //! Returns the union cinterval of the arguments
      friend inline cinterval operator |(const complex &,const real &);
      //! Returns the union cinterval of the arguments
      friend inline cinterval operator |(const real &,const complex &);

      // C-I

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const complex &,const interval &);
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const interval &,const complex &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const complex &,const interval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const interval &,const complex &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const complex &,const interval &);
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const interval &,const complex &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const complex &,const interval &);
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const interval &,const complex &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const complex &,const interval &);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const interval &,const complex &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const complex &,const interval &);
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const interval &,const complex &);
      

      // C-C

      //! Returns the union cinterval of the arguments
      friend inline cinterval operator |(const complex &,const complex &);

      // ---- Comp.Operat.  ---------------------------------------
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const cinterval & b);
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const cinterval & b);
      
      // CI-R
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const real & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const real & a, const cinterval & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const real & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real & a, const cinterval & b)   ;

      // CI-I
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const interval & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const interval & a, const cinterval & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const interval & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const interval & a, const cinterval & b)   ;

      // CI-C
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const complex & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const cinterval & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const complex & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const cinterval & b)   ;

      // ---- Set Operators ----
      
      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const cinterval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const cinterval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const cinterval &);
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const cinterval &);

      // CI-R

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const real &,const cinterval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const real &,const cinterval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const real &,const cinterval &);
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const real &,const cinterval &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const real &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const real &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const real &);
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const real &);

      // CI-I

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const interval &,const cinterval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const interval &,const cinterval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const interval &,const cinterval &);
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const interval &,const cinterval &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const interval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const interval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const interval &);
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const interval &);

      // CI-C

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const complex &,const cinterval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const complex &,const cinterval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const complex &,const cinterval &);
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const complex &,const cinterval &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const complex &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const complex &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const complex &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const complex &);

      // ---- Others   -------------------------------------------
      //! Returns the infimum of a complex interval
      friend inline complex    Inf(const cinterval &);
      //! Returns the supremum of a complex interval
      friend inline complex    Sup(const cinterval &);
      
      //! Returns the complex interval with the new given infimum value
      friend inline cinterval & SetInf(cinterval &,const complex &);
      //! Returns the complex interval with the new given infimum value
      friend inline cinterval & SetInf(cinterval &,const real &);
      //! Returns the complex interval with the new given supremum value
      friend inline cinterval & SetSup(cinterval &,const complex &);
      //! Returns the complex interval with the new given supremum value
      friend inline cinterval & SetSup(cinterval &,const real &);
      //! Returns the complex interval with the unchecked new given infimum value
      friend inline cinterval & UncheckedSetInf(cinterval &,const complex &);
      //! Returns the complex interval with the unchecked new given infimum value
      friend inline cinterval & UncheckedSetInf(cinterval &,const real &)   ;
      //! Returns the complex interval with the unchecked new given supremum value
      friend inline cinterval & UncheckedSetSup(cinterval &,const complex &);
      //! Returns the cinterval with the unchecked new given supremum value
      friend inline cinterval & UncheckedSetSup(cinterval &,const real &)   ;
      
      //! Returns the real interval of the complex interval
      friend inline interval & Re(cinterval & a)      ;
      //! Returns the real interval of the complex interval
      friend inline interval   Re(const cinterval & a);
      //! Returns the imaginary interval of the complex interval
      friend inline interval & Im(cinterval & a)      ;
      //! Returns the imaginary interval of the complex interval
      friend inline interval   Im(const cinterval & a);
      
      //! Sets the real interval of the complex interval
      friend inline cinterval & SetRe(cinterval & a,const interval & b);
      //! Sets the imaginary interval of the complex interval
      friend inline cinterval & SetIm(cinterval & a,const interval & b); 
      //! Sets the real interval of the complex interval
      friend inline cinterval & SetRe(cinterval & a,const real     & b);
      //! Sets the imaginary interval of the complex interval
      friend inline cinterval & SetIm(cinterval & a,const real     & b);

      //! Returns the infimum of the real interval of the complex interval
      friend inline real InfRe(const cinterval &a);
      //! Returns the infimum of the imaginary interval of the complex interval
      friend inline real InfIm(const cinterval &a);
      //! Returns the supremum of the real interval of the complex interval
      friend inline real SupRe(const cinterval &a);
      //! Returns the supremum of the imaginary interval of the complex interval
      friend inline real SupIm(const cinterval &a);
      
      //! Returns the infimum of the real interval of the complex interval
      friend inline real & InfRe(cinterval &a);
      //! Returns the infimum of the imaginary interval of the complex interval
      friend inline real & InfIm(cinterval &a);
      //! Returns the supremum of the real interval of the complex interval
      friend inline real & SupRe(cinterval &a);
      //! Returns the supremum of the imaginary interval of the complex interval
      friend inline real & SupIm(cinterval &a);
      
      //! Returns the absolute value of a complex interval
      friend        interval  abs(const cinterval &);
      //! Returns the conjugated complex interval
      friend inline cinterval conj(const cinterval &);
      //! Returns the rounded middle of the complex interval
      friend inline   complex mid(const cinterval &);
      //! Returns the rounded diameter of the complex interval
      friend inline   complex diam(const cinterval &);
};

//! Checks if first argument is part of second argument
extern int       in   ( const cinterval&, const cinterval& );
//! Performs an epsilon inflation
extern cinterval Blow ( cinterval, const real& );

// Additional declaration of friend functions outside class cinterval
interval  abs(const cinterval &);

} // namespace cxsc 

#include "cinterval.inl"
#include "cimath.hpp"

#endif // _CXSC_CINTERVAL_HPP_INCLUDED
 
