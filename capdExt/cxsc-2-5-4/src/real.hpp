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

/* CVS $Id: real.hpp,v 1.42 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_REAL_HPP_INCLUDED
#define _CXSC_REAL_HPP_INCLUDED

#include <iostream>
#include <string>

// namespace cxsc {

#include "compiler.h"
#include "RtsTyp.h"
#include "ioflags.hpp"

#include "except.hpp"
//! The namespace cxsc, providing all functionality of the class library C-XSC
/*!
In the namespace cxsc are all classes, data types, methods etc. defined, which
 are provided by the class library C-XSC.
*/
namespace cxsc {

class rvector;
class rmatrix;
class rvector_slice;
class rmatrix_slice;

#define  addu  addup
#define  addd  adddown
#define  subu  subup
#define  subd  subdown
#define  mulu  multup
#define  muld  multdown
#define  divu  divup
#define  divd  divdown

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- class real (declaration)                                          ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! The Scalar Type real
/*!
The arithmetic of C-XSC is based on the IEEE standard for binary floating-point
arithmetic. Data of the type real consist of all floating-point numbers and 
special values specified by the standard for floating-point numbers of double
mantissa length. Therefore, a number of the type real is a 64-bit value, and
the base \f$ b \f$ of the floating-point system 
\f$ R = R(b , l , e_min, e_max)\mbox{ is }2 \f$.

The C-XSC data type real differs only in some special aspects such as error 
handling from the C type double if used on a IEEE standard conforming 
arithmetic. If the C compiler on the host computer is not standard conforming,
the data type real uses its own IEEE software arithmetic. Hence, the 
introduction of a new data type is necessary for the portability of C-XSC.

\section structure Structure

A sketch of the real floating-point data format is given in the figure below.

\image html "ieee_real_float_format.png" "The real Floating-Point Format"

The most significant bit is the sign bit, denoted by \f$ s \f$ . If it is one,
the non-zero floating-point number is negative. Otherwise the floating-point 
number is positive.

The remaining 63 bits of the floating-point numbr are subdivided as follows:

\subsection mantissa Mantissa

The mantissa length \f$ l \f$ of a real value is 53 bits. Bit 53 is not
explicitly stored. Its value depends upon the normalization or denormalization
of the floating-point number represented. With this mantissa length \f$ l \f$,
decimal numbers with  a maximum of 15 fractional digits can be represented 
with high accuracy. Since a binary floating-point format is used, it is 
impossible to store all decimal numbers exactly. For example, the decimal
number 0.6 is not exactly representable as a binary floating-point number. 
However, we can read decimal constants with controlled rounding into a 
variable of the type real.

\subsection exponent Exponent

The exponent \f$ e \f$ of the floating-point number is stored in a field of 11
bit width using a biased notation. A value of 1023 is subtracted from the 
stored value to get the actual exponent value. This yields 
\f$ e_{min} = -1022 \f$ and \f$ e_{max} = 1023 \f$ .
Converted to the decimal system, this yields \f$ -308 \le e \le 308 \f$ .
*/
class real
{
   private:
      // ---- private data ----------------------------------------
      double w;

   public:
      // ---- implicit constructors -------------------------------
      //! Constructor of class real
      real(void)  throw ()                  { }
      //! Constructor of class real
      real(const float  &a) throw () : w(a) { }
      //! Constructor of class real
      real(const double &a) throw () : w(a) { }
      //! Constructor of class real
      real(const int     a) throw () : w(a) { }
      //! Constructor of class real
      real(const long    a) throw () : w(a) { }

      // ---- explicit constructors -------------------------------
      //! Constructor of class real
      explicit real(const l_real &) throw();

      // The following are defined in the specific vector, matrix-files
#if(CXSC_INDEX_CHECK) 
      //! Constructor of class real
      explicit INLINE real(const rvector &)       throw (ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class real
      explicit INLINE real(const rvector_slice &) throw (ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class real
      explicit INLINE real(const rmatrix &)       throw (ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Constructor of class real
      explicit INLINE real(const rmatrix_slice &) throw (ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ);
#else
      //! Constructor of class real
      explicit INLINE real(const rvector &)       throw ();
      //! Constructor of class real
      explicit INLINE real(const rvector_slice &) throw ();
      //! Constructor of class real
      explicit INLINE real(const rmatrix &)       throw ();
      //! Constructor of class real
      explicit INLINE real(const rmatrix_slice &) throw ();
#endif
      

      // ---- assignments -----------------------------------------

      // ---- compatibility typecasts -----------------------------

      //! Typecast to convert a real value into a double value
      friend inline double _double(const real &a) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline real   _real(const double &a) throw();

#if(CXSC_INDEX_CHECK)
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rvector &)
      */
      friend INLINE real _real(const rvector &)       throw (ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rvector_slice &)
      */
      friend INLINE real _real(const rvector_slice &) throw (ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rmatrix &)
      */
      friend INLINE real _real(const rmatrix &)       throw (ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rmatrix_slice &)
      */
      friend INLINE real _real(const rmatrix_slice &) throw (ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ);
#else
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rvector &)
      */
      friend INLINE real _real(const rvector &)       throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rvector_slice &)
      */
      friend INLINE real _real(const rvector_slice &) throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rmatrix &)
      */
      friend INLINE real _real(const rmatrix &)       throw ();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      /*!
      \deprecated use standard contructors for typecasting

      \sa cxsc::real::real(const rmatrix_slice &)
      */
      friend INLINE real _real(const rmatrix_slice &) throw ();
#endif
      
      // ---- Input/Output  ---------------------------------------

      //! Implementation of standard output method
      friend std::ostream & operator <<(std::ostream &,const real &) throw();
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &,real &)       throw();
      //! Implementation of standard output method
      friend std::string &  operator <<(std::string &,const real &)  throw();
      //! Implementation of standard input method
      friend std::string &  operator >>(std::string &,real &)        throw();
      //! Implementation of standard output method
      friend void           operator >>(const char *,real &)         throw();
      //! Implementation of standard input method
      friend void           operator >>(const std::string &,real &)  throw();

      // ---- Std.Operators ---------------------------------------
      // As the real-arithmetic should be as fast as double all
      // operators are inlined.                

      //! Implementation of standard algebraic negative sign operation
      friend inline real operator -(const real &) throw ();
      //! Implementation of standard algebraic positive sign operation
      friend inline real operator +(const real &) throw ();

      //! Implementation of standard algebraic addition operation
      friend inline real operator +(const real &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline real operator -(const real &,const real &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline real operator *(const real &,const real &) throw();
      //! Implementation of standard algebraic division operation
      friend inline real operator /(const real &,const real &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline real& operator +=(real &, const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline real& operator -=(real &, const real &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline real& operator *=(real &, const real &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline real& operator /=(real &, const real &) throw();

      // ---- Comp.Operat.  ---------------------------------------

      //! Implementation of standard negation operation
      friend inline bool operator!  (const real& a)                throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real& a, const real& b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real& a, const real& b) throw();
      //! Implementation of standard less-than operation
      friend inline bool operator<  (const real& a, const real& b) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator<= (const real& a, const real& b) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator>= (const real& a, const real& b) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator>  (const real& a, const real& b) throw();

      //! Implementation of standard equality operation
      friend inline bool operator== (const real& a, const int & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real& a, const int & b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const int & a, const real& b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const int & a, const real& b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real& a, const long & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real& a, const long & b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const long & a, const real& b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const long & a, const real& b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real& a, const float & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real& a, const float & b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const float & a, const real& b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const float & a, const real& b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real& a, const double & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real& a, const double & b) throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const double & a, const real& b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const double & a, const real& b) throw();
   
      // ---- Rounding Operators ---------------------------------

      //! Addition of two real values and rounding the result upwards
      friend inline real addup(const real &, const real &);
      //! Addition of two real values and rounding the result downwards
      friend inline real adddown(const real &, const real &);
      //! Subtraction of two real values and rounding the result upwards
      friend inline real subup(const real &, const real &);
      //! Subtraction of two real values and rounding the result downwards
      friend inline real subdown(const real &, const real &);
      //! Multiplication of two real values and rounding the result upwards
      friend inline real multup(const real &, const real &);
      //! Multiplication of two real values and rounding the result downwards
      friend inline real multdown(const real &, const real &);
      //! Division of two real values and rounding the result upwards
      friend inline real divup(const real &, const real &);
      //! Division of two real values and rounding the result downwards
      friend inline real divdown(const real &, const real &);

      // ---- Others   -------------------------------------------

      //! The absolute value of a real value
      friend inline real abs(const real &a) throw();
      //! The sign of a real value
      friend inline int  sign(const real &) throw();
      
      //! The predecessor of a real value
      friend inline real pred(const real &) throw();
      //! The successor of a real value
      friend inline real succ(const real &) throw();
      //! The exponent of a real value
      friend inline a_intg expo(const real &) throw();
      //! Composes an IEEE floating point value out of an given mantissa and exponent
      friend inline real comp(const real &,a_intg) throw();
      //! The mantissa of a real value
      friend inline real mant(const real &) throw();
		
      real & operator = (const lx_real&) throw();  // Blomquist, 12.11.2008;
      real & operator = (const l_real&)  throw();  // Blomquist, 12.11.2008;
}; // end of class real

inline real comp(const real &,a_intg) throw();

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- friend functions of class real (not inline)                       ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

std::ostream & operator <<(std::ostream &,const real &) throw();
std::istream & operator >>(std::istream &,real &)       throw();
std::string &  operator <<(std::string &,const real &)  throw();
std::string &  operator >>(std::string &,real &)        throw();
void           operator >>(const char *,real &)         throw();
void           operator >>(const std::string &,real &)  throw();

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- global functions associated with class real                       ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! Returns the greater value of two real values
inline real max(const real & a, const real & b);
//! Returns the smaller value of two real values
inline real min(const real & a, const real & b);
//! Returns the greater value of two real values (for Compatibility with former r_util.hpp)
inline real Max(const real & a, const real & b); 

//----------------------------------------------------------------------
// MakeHexReal - erstellt aus den binaer angegebenen Einzelteilen einer
//               IEEE 64-bit Gleitkommazahl eine solche

//! Produces an IEEE 64-bit floating-point number from given binary coded parts of an IEEE 64-bit floating-point number
const real& MakeHexReal ( 
   int sign, unsigned int expo, a_btyp manthigh, a_btyp mantlow);

//! Returns if the given real value represents the value infinity
inline bool IsInfinity(const real &a);
//! Returns if the given real value represents the value of a quiet NaN 
inline bool IsQuietNaN(const real &a);
//! Returns if the given real value represents the value of a signaling NaN
inline bool IsSignalingNaN(const real &a);

//-------------------------------------------------------------------------
// times2pown - Fast multiplication of reference parameter r with 2^n -----
//! Fast multiplication of reference parameter r with \f$ 2^n \f$
inline void times2pown(real& r,const int n); // Blomquist 1.10.02. {real.inl}

//! Returns the value of \f$ 2^n \f$
inline real pow2n(const int n) throw(); // returns 2^n; 

//!Returns a real number in hexadecimal format as string
string realToHex(const real& a);

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- special constants associated with class real                      ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

// ----  special constants -----------------------------------------
//! Smallest normalized representable floating-point number
extern const real MinReal;
//! Smallest positive denormalized representable floating-point number
extern const real minreal;
//! Greatest representable floating-point number
extern const real MaxReal;
//! Representation of positive infinity in floating-point format
extern const real Infinity;
//! Not defined result in floating-point format
extern const real SignalingNaN;
//! Representation of Not-a-Number in floating-point format
extern const real QuietNaN; 
//! Machine epsilon
extern const real Epsilon;
extern const real Factor;

//! Constant for \f$ \pi \f$ rounded to the nearest machine number
extern const real Pi_real;        // Pi
//! Constant for \f$ 2 \pi \f$ rounded to the nearest machine number
extern const real Pi2_real;       // 2*Pi
//! Constant for \f$ 3 \pi \f$ rounded to the nearest machine number
extern const real Pi3_real;       // 3*Pi
//! Constant for \f$ \frac{\pi}{2} \f$ rounded to the nearest machine number
extern const real Pid2_real;      // Pi/2
//! Constant for \f$ \frac{\pi}{3} \f$ rounded to the nearest machine number
extern const real Pid3_real;      // Pi/3
//! Constant for \f$ \frac{\pi}{4} \f$ rounded to the nearest machine number
extern const real Pid4_real;      // Pi/4
//! Constant for \f$ \frac{1}{\pi} \f$ rounded to the nearest machine number
extern const real Pir_real;       // 1/Pi
//! Constant for \f$ \frac{1}{2 \cdot \pi} \f$ rounded to the nearest machine number
extern const real Pi2r_real;      // 1/(2*Pi)
//! Constant for \f$ \pi^2 \f$ rounded to the nearest machine number
extern const real Pip2_real;      // Pi^2
//! Constant for \f$ \sqrt{\pi} \f$ rounded to the nearest machine number
extern const real SqrtPi_real;    // sqrt(Pi)
//! Constant for \f$ \sqrt{2 \pi} \f$ rounded to the nearest machine number
extern const real Sqrt2Pi_real;   // sqrt(2Pi)
//! Constant for \f$ \frac{1}{\sqrt{\pi}} \f$ rounded to the nearest machine number
extern const real SqrtPir_real;   // 1/sqrt(Pi)
//! Constant for \f$ \frac{1}{\sqrt{2 \pi}} \f$ rounded to the nearest machine number
extern const real Sqrt2Pir_real;  // 1/sqrt(2Pi)
//! Constant for \f$ \sqrt{2} \f$ rounded to the nearest machine number
extern const real Sqrt2_real;     // sqrt(2)
//! Constant for \f$ \sqrt{5} \f$ rounded to the nearest machine number
extern const real Sqrt5_real;     // sqrt(5)
//! Constant for \f$ \sqrt{7} \f$ rounded to the nearest machine number
extern const real Sqrt7_real;     // sqrt(7)
//! Constant for \f$ \frac{1}{\sqrt{2}} \f$ rounded to the nearest machine number
extern const real Sqrt2r_real;    // 1/sqrt(2)
//! Constant for \f$ \sqrt{3} \f$ rounded to the nearest machine number
extern const real Sqrt3_real;     // sqrt(3)
//! Constant for \f$ \frac{\sqrt{3}}{2} \f$ rounded to the nearest machine number
extern const real Sqrt3d2_real;   // sqrt(3)/2
//! Constant for \f$ \frac{1}{\sqrt{3}} \f$ rounded to the nearest machine number
extern const real Sqrt3r_real;    // 1/sqrt(3)
//! Constant for \f$ \ln 2 \f$ rounded to the nearest machine number
extern const real Ln2_real;       // ln(2)
//! Constant for \f$ \frac{1}{\ln 2} \f$ rounded to the nearest machine number
extern const real Ln2r_real;      // 1/ln(2)
//! Constant for \f$ \ln 10 \f$ rounded to the nearest machine number
extern const real Ln10_real;      // ln(10)
//! Constant for \f$ \frac{1}{\ln 10} \f$ rounded to the nearest machine number
extern const real Ln10r_real;     // 1/ln(10)
//! Constant for \f$ \ln \pi \f$ rounded to the nearest machine number
extern const real LnPi_real;      // ln(Pi)
//! Constant for \f$ \ln (2 \pi) \f$ rounded to the nearest machine number
extern const real Ln2Pi_real;     // ln(2Pi)
//! Constant for \f$ \mbox{e} \f$ rounded to the nearest machine number
extern const real E_real;         // e
//! Constant for \f$ \frac{1}{\mbox{e}} \f$ rounded to the nearest machine number
extern const real Er_real;        // 1/e
//! Constant for \f$ \mbox{e}^2 \f$ rounded to the nearest machine number
extern const real Ep2_real;       // e^2
//! Constant for \f$ \frac{1}{\mbox{e}^2} \f$ rounded to the nearest machine number
extern const real Ep2r_real;      // 1/e^2
//! Constant for \f$ \mbox{e}^\pi \f$ rounded to the nearest machine number
extern const real EpPi_real;      // e^(Pi)
//! Constant for \f$ \mbox{€}^{2 \pi} \f$ rounded to the nearest machine number
extern const real Ep2Pi_real;     // e^(2Pi)
//! Constant for \f$ \mbox{e}^{\frac{\pi}{2}} \f$ rounded to the nearest machine number
extern const real EpPid2_real;    // e^(Pi/2)
//! Constant for \f$ \mbox{e}^{\frac{\pi}{4}} \f$ rounded to the nearest machine number
extern const real EpPid4_real;    // e^(Pi/4)

} // namespace cxsc 

#include "real.inl"
#include "rmath.hpp"
// }
#endif // _CXSC_REAL_HPP_INCLUDED

