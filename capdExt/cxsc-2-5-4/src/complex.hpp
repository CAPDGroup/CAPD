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

/* CVS $Id: complex.hpp,v 1.32 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_COMPLEX_HPP_INCLUDED
#define _CXSC_COMPLEX_HPP_INCLUDED

#include <iostream>
#include <string>
#include <list>  // Blomquist, 02.12.2008;
#include "compiler.h"
#include "except.hpp"
#include "real.hpp"

namespace cxsc {

class cvector;
class cmatrix;
class cvector_slice;
class cmatrix_slice;


//!The Scalar Type complex
/*!
The data type complex is used to store complex numbers \f$ z = x + i y \in C \f$, where \f$ x \f$ denotes the real part, \f$ y \f$ denotes the imaginary
part of \f$ z \f$, and \f$ i \f$ denotes the imaginary unit \f$ \sqrt{-1} \f$.
*/
class complex
{  
   private:
      // ---- Datenelemente ---------------------------------------
      real  re;
      real  im;

   public:
      // ---- Constructors  ---------------------------------------
      //! Constructor of class complex
      complex(void)  throw ()           {}
      //! Constructor of class complex
      complex(const real & a,const real & b) throw () : re(a), im(b) { }
      
      //! Implementation of standard assigning operator
      inline complex & operator= (const real & r) throw();

      // ---- Type-Casts    ---------------------------------------

      //! Constructor of class complex
      explicit inline complex(const real &r) throw() : re(r),im(0.0)     { }


//     friend inline complex _complex(const real &a) throw ()             { return complex(a,0.0); }
//     friend inline complex _complex(const real &a,const real &b) throw(){ return complex(a,b); }


      // The following are defined in the specific vector, matrix-files
#if(CXSC_INDEX_CHECK) 
      //! Constructor of class complex
      explicit INLINE complex(const cvector &)       throw (ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class complex
      explicit INLINE complex(const cvector_slice &) throw (ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class complex
      explicit INLINE complex(const cmatrix &)       throw (ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class complex
      explicit INLINE complex(const cmatrix_slice &) throw (ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cvector &)
      */
      friend INLINE complex _complex(const cvector &)       throw (ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cvector_slice &)
      */
      friend INLINE complex _complex(const cvector_slice &) throw (ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cmatrix &)
      */
      friend INLINE complex _complex(const cmatrix &)       throw (ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cmatrix_slice &)
      */
      friend INLINE complex _complex(const cmatrix_slice &) throw (ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ);
#else
      //! Constructor of class complex
      explicit INLINE complex(const cvector &)       throw ();
      //! Constructor of class complex
      explicit INLINE complex(const cvector_slice &) throw ();
      //! Constructor of class complex
      explicit INLINE complex(const cmatrix &)       throw ();
      //! Constructor of class complex
      explicit INLINE complex(const cmatrix_slice &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cvector &)
      */
      friend INLINE complex _complex(const cvector &)       throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cvector_slice &)
      */
      friend INLINE complex _complex(const cvector_slice &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cmatrix &)
      */
      friend INLINE complex _complex(const cmatrix &)       throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cmatrix_slice &)
      */
      friend INLINE complex _complex(const cmatrix_slice &) throw ();
#endif
      //! Constructor of class complex
      explicit        complex(const cdotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::complex::complex(const cdotprecision &)
      */
      friend inline complex _complex(const cdotprecision &a) throw() { return complex(a); }
      
      //! Implementation of standard assigning operator
      complex & operator =(const cdotprecision &) throw();

      // ---- Input/Output  ---------------------------------------
      //! Implementation of standard output method
      friend std::ostream & operator <<(std::ostream &,const complex &) throw();
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &,complex &)       throw();
      //! Implementation of standard output method
      friend std::string &  operator <<(std::string &,const complex &)  throw();
      //! Implementation of standard input method
      friend std::string &  operator >>(std::string &,complex &)        throw();
      //! Implementation of standard input method
      friend void           operator >>(const char *,complex &)         throw();
      //! Implementation of standard input method
      friend void           operator >>(const std::string &,complex &)  throw();

      // ---- Std.Operators ---------------------------------------
      //! Implementation of standard algebraic negative sign operation
      friend inline complex operator -(const complex &) throw ();
      //! Implementation of standard algebraic positive sign operation
      friend inline complex operator +(const complex &) throw ();

      //! Implementation of standard algebraic addition operation
      friend inline complex operator +(const complex &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline complex operator -(const complex &,const complex &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend complex operator *(const complex &,const complex &) throw();
      //! Implementation of standard algebraic division operation
      friend complex operator /(const complex &,const complex &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline complex & operator +=(complex &, const complex &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline complex & operator -=(complex &, const complex &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline complex & operator *=(complex &, const complex &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline complex & operator /=(complex &, const complex &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline complex operator +(const complex &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline complex operator +(const real &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline complex operator -(const complex &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline complex operator -(const real &,const complex &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline complex operator *(const complex &,const real &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline complex operator *(const real &,const complex &) throw();
      //! Implementation of standard algebraic division operation
      friend inline complex operator /(const complex &,const real &) throw();
      //! Implementation of standard algebraic division operation
      friend inline complex operator /(const real &,const complex &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline complex & operator +=(complex &, const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline complex & operator -=(complex &, const real &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline complex & operator *=(complex &, const real &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline complex & operator /=(complex &, const real &) throw();

      // ---- Comp.Operat.  ---------------------------------------
//             inline       operator void *() const throw() { if(re) return (void *)1; if(im) return (void *)1; else return 0; }
      //! Implementation of standard negation operation
      friend inline bool operator!  (const complex & a)                    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const complex & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const complex & b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const real & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real & a, const complex & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const real & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real & a, const complex & b)    throw();
      
      //! Implementation of standard equality operation
      friend        bool operator== (const complex & a, const dotprecision & b)    throw();
      //! Implementation of standard equality operation
      friend        bool operator== (const dotprecision & a, const complex & b)    throw();
      //! Implementation of standard negated equality operation
      friend        bool operator!= (const complex & a, const dotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend        bool operator!= (const dotprecision & a, const complex & b)    throw();

      // ---- Others   -------------------------------------------
      //! Returns the real part of the complex value
      friend inline real & Re(complex & a); // { return a.re; }
      //! Returns the real part of the complex value
      friend inline real   Re(const complex & a); // { return a.re; }
      //! Returns the imaginary part of the complex value
      friend inline real & Im(complex & a); // { return a.im; }
      //! Returns the imaginary part of the complex value
      friend inline real   Im(const complex & a); // { return a.im; }
      
      //! Sets the real part of a complex value
      friend inline complex & SetRe(complex & a,const real & b); // { a.re=b; return a; }
      //! Sets the imaginary part of a complex value
      friend inline complex & SetIm(complex & a,const real & b); // { a.im=b; return a; } 

      //! Returns the absolute value of a complex value
      friend        real abs(complex) throw();
      //! Returns the absolute value of a complex value
      friend        real abs2(const complex &) throw();
      //! Returns the conjugated complex value
      friend inline complex conj(const complex &) throw();

// -------------- Directed rounding, Blomquist 07.11.02 --------------------

    //! Returns the nearest rounded result of the division operation
    friend complex divn (const complex &, const complex &);
    //! Returns the downward rounded result of the division operation
    friend complex divd (const complex &, const complex &);
    //! Returns the upward rounded result of the  operation
    friend complex divu (const complex &, const complex &);
    //! Returns the downward rounded result of the division operation
    friend inline complex divd(const complex &, const real &) throw();
    //! Returns the upward rounded result of the division operation
    friend inline complex divu(const complex &, const real &) throw();
    //! Returns the downward rounded result of the division operation
    friend inline complex divd(const real &, const complex &) throw();
    //! Returns the upward rounded result of the division operation
    friend inline complex divu(const real &, const complex &) throw();

    //! Returns the downward rounded result of the multiplication operation
    friend complex muld (const complex &, const complex &) throw();
    //! Returns the upward rounded result of the multiplication operation
    friend complex mulu (const complex &, const complex &) throw();
    //! Returns the downward rounded result of the multiplication operation
    friend inline complex muld(const complex &, const real &) throw();
    //! Returns the upward rounded result of the multiplication operation
    friend inline complex mulu(const complex &, const real &) throw();
    //! Returns the downward rounded result of the multiplication operation
    friend inline complex muld(const real &, const complex &) throw();
    //! Returns the upward rounded result of the multiplication operation
    friend inline complex mulu(const real &, const complex &) throw();


    //! Returns the downward rounded result of the addition operation
    friend inline complex addd(const complex &, const complex &) throw();
    //! Returns the upward rounded result of the addition operation
    friend inline complex addu(const complex &, const complex &) throw();
    //! Returns the downward rounded result of the addition operation
    friend inline complex addd(const complex &, const real &) throw();
    //! Returns the upward rounded result of the addition operation
    friend inline complex addu(const complex &, const real &) throw();
    //! Returns the downward rounded result of the addition operation
    friend inline complex addd(const real &, const complex &) throw();
    //! Returns the upward rounded result of the addition operation
    friend inline complex addu(const real &, const complex &) throw();


    //! Returns the downward rounded result of the subtraction operation
    friend inline complex subd(const complex &, const complex &) throw();
    //! Returns the upward rounded result of the subtraction operation
    friend inline complex subu(const complex &, const complex &) throw();
    //! Returns the downward rounded result of the subtraction operation
    friend inline complex subd(const complex &, const real &) throw();
    //! Returns the upward rounded result of the subtraction operation
    friend inline complex subu(const complex &, const real &) throw();
    //! Returns the downward rounded result of the subtraction operation
    friend inline complex subd(const real &, const complex &) throw();
    //! Returns the upward rounded result of the subtraction operation
    friend inline complex subu(const real &, const complex &) throw();
	 
    //! Assigning lx_complex to complex
    complex & operator = (const lx_complex&) throw(); // Blomquist, 12.11.2008;
    //! Assigning l_complex to complex
    complex & operator = (const l_complex&)  throw(); // Blomquist, 12.11.2008;
}; // end class complex


// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- friend functions of class real (not inline)                       ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

complex divn (const complex &, const complex &);
complex divd (const complex &, const complex &);
complex divu (const complex &, const complex &);
complex muld (const complex &, const complex &) throw();
complex mulu (const complex &, const complex &) throw();

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- global functions associated with class real                       ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::complex::complex(const real &r)
*/
inline complex _complex(const real &a) throw ()             { return complex(a,0.0); }
//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::complex::complex(const real & a,const real & b)
*/
inline complex _complex(const real &a,const real &b) throw(){ return complex(a,b); }

//! Returns the real part of a variable z of type complex
inline real & Re(complex & z) { return z.re; }
//! Returns the real part of a variable z of type complex
inline real   Re(const complex & z) { return z.re; }
//! Returns the imaginary part of a variable z of type complex
inline real & Im(complex & z) { return z.im; }
//! Returns the imaginary part of a variable z of type complex
inline real   Im(const complex & z) { return z.im; }

//! Sets a new real part of a variable z of type complex
inline complex & SetRe(complex & z,const real & b) { z.re=b; return z; }
//! Sets a new imaginary part of a variable z of type complex
inline complex & SetIm(complex & z,const real & b) { z.im=b; return z; } 

//! Calculates an approximation of \f$ z^2 \f$
inline complex sqr (const complex&) throw();
//! Calculates an approximation of \f$ \sqrt(z) \f$
complex sqrt(const complex&) throw();
//! Calculates an approximation of \f$ \sqrt(1+z)-1 \f$
complex sqrtp1m1(const complex&) throw();
//! Calculates an approximation of \f$ \sqrt(1+z^2) \f$
complex sqrt1px2(const complex&) throw();
//! Calculates an approximation of \f$ \sqrt(z^2-1) \f$
complex sqrtx2m1(const complex&) throw();
//! Calculates an approximation of \f$ \sqrt(1-z^2) \f$
complex sqrt1mx2(const complex&) throw();

//! Calculates an approximation of \f$ \exp(z) \f$
complex exp(const complex&) throw();
//! Calculates an approximation of \f$ \exp(z)-1 \f$
complex expm1(const complex&) throw();
//! Calculates an approximation of \f$ 2^z \f$
complex exp2(const complex&) throw();
//! Calculates an approximation of \f$ 10^z \f$
complex exp10(const complex&) throw();
//! Calculates an approximation of \f$ \sin(z) \f$
complex sin(const complex&) throw();
//! Calculates an approximation of \f$ \cos(z) \f$
complex cos(const complex&) throw();
//! Calculates an approximation of \f$ \tan(z) \f$
complex tan(const complex&) throw();
//! Calculates an approximation of \f$ \cot(z) \f$
complex cot(const complex&) throw();
//! Calculates an approximation of \f$ \arcsin(z) \f$
complex asin(const complex&) throw();
//! Calculates an approximation of \f$ \arccos(z) \f$
complex acos(const complex&) throw();
//! Calculates an approximation of \f$ \arctan(z) \f$
complex atan(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{arccot}(z) \f$
complex acot(const complex&) throw();
//! Calculates an approximation of \f$ \sinh(z) \f$
complex sinh(const complex&) throw();
//! Calculates an approximation of \f$ \cosh(z) \f$
complex cosh(const complex&) throw();
//! Calculates an approximation of \f$ \tanh(z) \f$
complex tanh(const complex&) throw();
//! Calculates an approximation of \f$ \coth(z) \f$
complex coth(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{arcsinh}(z) \f$
complex asinh(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{arccosh}(z) \f$
complex acosh(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{arctanh}(z) \f$
complex atanh(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{arccoth}(z) \f$
complex acoth(const complex&) throw();
//! Calculates an approximation of \f$ \sqrt{z}  \f$ and returns all possible solutions
std::list<complex>sqrt_all(const complex&);
//! Calculates an approximation of \f$ \sqrt[n]{z} \f$
complex sqrt(const complex&, int) throw();
//! Calculates an approximation of \f$ \mbox{arg}(z) \f$
real arg(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{arg}(z) \f$
real Arg(const complex&) throw();
//! Calculates an approximation of \f$ \sqrt[n]{z} \f$ and returns all possible solutions
std::list<complex>sqrt_all(const complex&, int);
//! Calculates an approximation of \f$ \ln(z) \f$
complex ln(const complex&) throw();
//! Calculates an approximation of \f$ \ln(1+z) \f$
complex lnp1(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{log2}(z) \f$
complex log2(const complex&) throw();
//! Calculates an approximation of \f$ \mbox{log10}(z) \f$
complex log10(const complex&) throw();
//! Calculates an approximation of \f$ z^n \f$
complex power(const complex&,int) throw();
//! Calculates an approximation of \f$ z^n \f$
complex power_fast(const complex&, int) throw();
//! Calculates an approximation of \f$ z^y \f$
complex pow(const complex&, const real&) throw();
//! Calculates an approximation of \f$ z_1^{z_2} \f$
complex pow(const complex&, const complex&) throw();

} // namespace cxsc 


#include "complex.inl"

#endif

