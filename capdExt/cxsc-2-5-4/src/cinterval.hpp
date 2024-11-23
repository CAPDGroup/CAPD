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
      inline cinterval(void)  throw ()           {}
      //! Constructor of class cinterval
      inline cinterval(const interval & a,const interval &b) throw(); 
      //! Constructor of class cinterval
      inline cinterval(const complex & a,const complex & b)  throw(ERROR_CINTERVAL_EMPTY_INTERVAL); 

      // The following are defined in the specific vector, matrix-files
#if(CXSC_INDEX_CHECK) 
      //! Constructor of class cinterval
      INLINE cinterval(const civector &)       throw (ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class cinterval
      INLINE cinterval(const civector_slice &) throw (ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix &)       throw (ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix_slice &) throw (ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ);
#else
      //! Constructor of class cinterval
      INLINE cinterval(const civector &)       throw ();
      //! Constructor of class cinterval
      INLINE cinterval(const civector_slice &) throw ();
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix &)       throw ();
      //! Constructor of class cinterval
      INLINE cinterval(const cimatrix_slice &) throw ();
#endif
      // ---- explicit constructors -------------------------------

      //! Constructor of class cinterval
      explicit inline cinterval(const real     & a)  throw();
      //! Constructor of class cinterval
      explicit inline cinterval(const interval & a)  throw();
      //! Constructor of class cinterval
      explicit inline cinterval(const complex  & a)  throw(); 
      //! Constructor of class cinterval
      explicit        cinterval(const dotprecision &)  throw();
      //! Constructor of class cinterval
      explicit        cinterval(const cdotprecision &) throw();
      //! Constructor of class cinterval
      explicit        cinterval(const idotprecision &) throw();
      //! Constructor of class cinterval
      explicit        cinterval(const cidotprecision &) throw();
      //! Constructor of class cinterval
      explicit        cinterval(const l_cinterval&) throw();

      // ---- assignments -----------------------------------------

      //! Implementation of standard assigning operator
      inline cinterval & operator =(const real &) throw();
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const interval &) throw();
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const complex &) throw();
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const cinterval &) throw();
      
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const dotprecision &) throw();
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const idotprecision &) throw();
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const cdotprecision &) throw();
      //! Implementation of standard assigning operator
      inline cinterval & operator =(const cidotprecision &) throw();
      //! Implementation of standard assigning operator
             cinterval & operator = (const l_cinterval&) throw();      
      //! Implementation of standard assigning operator
             cinterval & operator = (const lx_cinterval&) throw();

      // ---- compatiblility typecasts ----------------------------

      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const real &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const interval &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const complex &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const dotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const idotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const cdotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const cidotprecision &) throw();
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const interval &,const interval &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const real &,const interval &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const interval &,const real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const complex &,const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const real &,const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _cinterval(const complex &,const real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _unchecked_cinterval(const complex &,const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _unchecked_cinterval(const real &,const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cinterval _unchecked_cinterval(const complex &,const real &) throw();

      // ---- Input/Output  ---------------------------------------

      //! Implementation of standard output method
      friend std::ostream & operator <<(std::ostream &,const cinterval &) throw();
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &,cinterval &)       throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Implementation of standard output method
      friend std::string &  operator <<(std::string &,const cinterval &)  throw();
      //! Implementation of standard input method
      friend std::string &  operator >>(std::string &,cinterval &)        throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Implementation of standard input method
      friend void           operator >>(const char *,cinterval &)         throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Implementation of standard input method
      friend void           operator >>(const std::string &,cinterval &)  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);

      // ---- Std.Operators ---------------------------------------

      //! Implementation of standard algebraic negative sign operation
      friend inline cinterval operator -(const cinterval &) throw ();
      //! Implementation of standard algebraic positive sign operation
      friend inline cinterval operator +(const cinterval &) throw ();
      //! Implementation of standard negation operation
      friend inline bool operator!  (const cinterval & a)                    throw();

      // CI-CI

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend        cinterval operator *(const cinterval &,const cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend        cinterval operator /(const cinterval &,const cinterval &) throw(DIV_BY_ZERO);
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const cinterval &,const cinterval &) throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const cinterval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const cinterval &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const cinterval &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const cinterval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const cinterval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const cinterval &) throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      
      // CI-R
      
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const real &,const cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const real &,const cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const cinterval &,const real &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const real &,const cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const cinterval &,const real &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const real &,const cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const real &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const real &,const cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const cinterval &,const real &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const real &,const cinterval &) throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const real &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const real &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const real &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const real &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const real &) throw();
      
      // CI-I

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const interval &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const interval &,const cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const interval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const interval &,const cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const cinterval &,const interval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const interval &,const cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const cinterval &,const interval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const interval &,const cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const interval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const interval &,const cinterval &) throw();
      //! Returns the intersection cinterval of the arguments
      friend inline cinterval operator &(const cinterval &,const interval &) throw();
      //! Returns the intersection cinterval of the arguments
      friend inline cinterval operator &(const interval &,const cinterval &) throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const interval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const interval &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const interval &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const interval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const interval &) throw();
      //! Allocates the intersection cinterval of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const interval &) throw();

      // CI-C

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const cinterval &,const complex &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const complex &,const cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const cinterval &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const complex &,const cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const cinterval &,const complex &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const complex &,const cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const cinterval &,const complex &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const complex &,const cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const cinterval &,const complex &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const complex &,const cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const cinterval &,const complex &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const complex &,const cinterval &) throw();
      

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cinterval & operator +=(cinterval &, const complex &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cinterval & operator -=(cinterval &, const complex &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline cinterval & operator *=(cinterval &, const complex &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline cinterval & operator /=(cinterval &, const complex &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cinterval & operator |=(cinterval &, const complex &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cinterval & operator &=(cinterval &, const complex &) throw();
      
      // C-R

      //! Returns the union cinterval of the arguments
      friend inline cinterval operator |(const complex &,const real &) throw();
      //! Returns the union cinterval of the arguments
      friend inline cinterval operator |(const real &,const complex &) throw();

      // C-I

      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const complex &,const interval &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cinterval operator +(const interval &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const complex &,const interval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cinterval operator -(const interval &,const complex &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const complex &,const interval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline cinterval operator *(const interval &,const complex &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const complex &,const interval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline cinterval operator /(const interval &,const complex &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const complex &,const interval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cinterval operator |(const interval &,const complex &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const complex &,const interval &) throw();
      //! Returns the intersection of the arguments
      friend inline cinterval operator &(const interval &,const complex &) throw();
      

      // C-C

      //! Returns the union cinterval of the arguments
      friend inline cinterval operator |(const complex &,const complex &) throw();

      // ---- Comp.Operat.  ---------------------------------------
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const cinterval & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const cinterval & b) throw();
      
      // CI-R
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const real & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real & a, const cinterval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const real & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real & a, const cinterval & b)    throw();

      // CI-I
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const interval & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const interval & a, const cinterval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const interval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const interval & a, const cinterval & b)    throw();

      // CI-C
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const complex & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const cinterval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const complex & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const cinterval & b)    throw();

      // ---- Set Operators ----
      
      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const cinterval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const cinterval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const cinterval &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const cinterval &) throw();

      // CI-R

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const real &,const cinterval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const real &,const cinterval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const real &,const cinterval &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const real &,const cinterval &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const real &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const real &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const real &) throw();

      // CI-I

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const interval &,const cinterval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const interval &,const cinterval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const interval &,const cinterval &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const interval &,const cinterval &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const interval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const interval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const interval &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const interval &) throw();

      // CI-C

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const complex &,const cinterval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const complex &,const cinterval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const complex &,const cinterval &) throw();
      //! Implementation of standard more-or-equal-than operation
      friend inline bool operator >=(const complex &,const cinterval &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const complex &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const complex &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const complex &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cinterval &,const complex &) throw();

      // ---- Others   -------------------------------------------
      //! Returns the infimum of a complex interval
      friend inline complex    Inf(const cinterval &) throw();
      //! Returns the supremum of a complex interval
      friend inline complex    Sup(const cinterval &) throw();
      
      //! Returns the complex interval with the new given infimum value
      friend inline cinterval & SetInf(cinterval &,const complex &) throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the complex interval with the new given infimum value
      friend inline cinterval & SetInf(cinterval &,const real &) throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the complex interval with the new given supremum value
      friend inline cinterval & SetSup(cinterval &,const complex &) throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the complex interval with the new given supremum value
      friend inline cinterval & SetSup(cinterval &,const real &) throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the complex interval with the unchecked new given infimum value
      friend inline cinterval & UncheckedSetInf(cinterval &,const complex &) throw();
      //! Returns the complex interval with the unchecked new given infimum value
      friend inline cinterval & UncheckedSetInf(cinterval &,const real &)    throw();
      //! Returns the complex interval with the unchecked new given supremum value
      friend inline cinterval & UncheckedSetSup(cinterval &,const complex &) throw();
      //! Returns the cinterval with the unchecked new given supremum value
      friend inline cinterval & UncheckedSetSup(cinterval &,const real &)    throw();
      
      //! Returns the real interval of the complex interval
      friend inline interval & Re(cinterval & a)       throw();
      //! Returns the real interval of the complex interval
      friend inline interval   Re(const cinterval & a) throw();
      //! Returns the imaginary interval of the complex interval
      friend inline interval & Im(cinterval & a)       throw();
      //! Returns the imaginary interval of the complex interval
      friend inline interval   Im(const cinterval & a) throw();
      
      //! Sets the real interval of the complex interval
      friend inline cinterval & SetRe(cinterval & a,const interval & b);
      //! Sets the imaginary interval of the complex interval
      friend inline cinterval & SetIm(cinterval & a,const interval & b); 
      //! Sets the real interval of the complex interval
      friend inline cinterval & SetRe(cinterval & a,const real     & b);
      //! Sets the imaginary interval of the complex interval
      friend inline cinterval & SetIm(cinterval & a,const real     & b);

      //! Returns the infimum of the real interval of the complex interval
      friend inline real InfRe(const cinterval &a) throw();
      //! Returns the infimum of the imaginary interval of the complex interval
      friend inline real InfIm(const cinterval &a) throw();
      //! Returns the supremum of the real interval of the complex interval
      friend inline real SupRe(const cinterval &a) throw();
      //! Returns the supremum of the imaginary interval of the complex interval
      friend inline real SupIm(const cinterval &a) throw();
      
      //! Returns the infimum of the real interval of the complex interval
      friend inline real & InfRe(cinterval &a) throw();
      //! Returns the infimum of the imaginary interval of the complex interval
      friend inline real & InfIm(cinterval &a) throw();
      //! Returns the supremum of the real interval of the complex interval
      friend inline real & SupRe(cinterval &a) throw();
      //! Returns the supremum of the imaginary interval of the complex interval
      friend inline real & SupIm(cinterval &a) throw();
      
      //! Returns the absolute value of a complex interval
      friend        interval  abs(const cinterval &) throw();
      //! Returns the conjugated complex interval
      friend inline cinterval conj(const cinterval &) throw();
      //! Returns the rounded middle of the complex interval
      friend inline   complex mid(const cinterval &) throw();
      //! Returns the rounded diameter of the complex interval
      friend inline   complex diam(const cinterval &) throw();
};

//! Checks if first argument is part of second argument
extern int       in   ( const cinterval&, const cinterval& );
//! Performs an epsilon inflation
extern cinterval Blow ( cinterval, const real& );

// Additional declaration of friend functions outside class cinterval
interval  abs(const cinterval &) throw();

} // namespace cxsc 

#include "cinterval.inl"
#include "cimath.hpp"

#endif // _CXSC_CINTERVAL_HPP_INCLUDED
 
