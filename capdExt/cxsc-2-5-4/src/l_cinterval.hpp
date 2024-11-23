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

/* CVS $Id: l_cinterval.hpp,v 1.20 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_L_CINTERVAL_HPP_INCLUDED
#define _CXSC_L_CINTERVAL_HPP_INCLUDED

#include <iostream>
#include <string>

// Headerfile for l_cinterval.

#include "except.hpp"
#include "l_real.hpp"
#include "l_complex.hpp"
#include "l_interval.hpp"
#include "cinterval.hpp"

namespace cxsc {

// class civector;
// class cimatrix;
// class civector_slice;
// class cimatrix_slice;

//! The Multiple-Precision Data Type l_cinterval
/*!
The multiple-precision data type l_cinterval is a variant of the scalar type cinterval, which provides support for longer numbers, thus increasing the accuracy of the data type.

\sa cinterval
*/
class l_cinterval
{
   private:
      // ---- private data ----------------------------------------
      l_interval  re;
      l_interval  im;

   public:
      // ---- implicit constructors  ------------------------------
      //! Constructor of class l_cinterval
      inline l_cinterval(void)  throw ()           {}
      //! Constructor of class l_cinterval
      inline l_cinterval(const interval & a,const interval &b) throw();
      //! Constructor of class l_cinterval
      inline l_cinterval(const l_interval & a,const l_interval &b) throw(); 
      //! Constructor of class l_cinterval
      inline l_cinterval(const complex & a, const complex & b) 
                                     throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Constructor of class l_cinterval
      inline l_cinterval(const l_complex & a, const l_complex & b) 
                                     throw(ERROR_CINTERVAL_EMPTY_INTERVAL); 

      friend cinterval::cinterval(const l_cinterval &) throw();

      // ---- explicit constructors -------------------------------

      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const real       & a)  throw();
      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const l_real     & a)  throw();
      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const interval   & a)  throw();
      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const l_interval & a)  throw();
      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const complex    & a)  throw();
      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const l_complex  & a)  throw();
      //! Constructor of class l_cinterval
      explicit inline l_cinterval(const cinterval  & a)  throw();  

      //! Constructor of class l_cinterval
      explicit        l_cinterval(const dotprecision &)   throw();
      //! Constructor of class l_cinterval
      explicit        l_cinterval(const cdotprecision &)  throw();
      //! Constructor of class l_cinterval
      explicit        l_cinterval(const idotprecision &)  throw();
      //! Constructor of class l_cinterval
      explicit        l_cinterval(const cidotprecision &) throw();

      // ---- assignments -----------------------------------------

      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const real        &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const l_real      &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const interval    &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const l_interval  &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const complex     &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const l_complex   &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const cinterval   &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const l_cinterval &) throw();
      //! Implementation of standard assigning operator
             l_cinterval & operator =(const lx_cinterval &) throw();
      
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const dotprecision   &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const idotprecision  &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const cdotprecision  &) throw();
      //! Implementation of standard assigning operator
      inline l_cinterval & operator =(const cidotprecision &) throw();

    friend cinterval & cinterval::operator = (const l_cinterval &) throw();

      // ---- Input/Output  ---------------------------------------

      //! Implementation of standard output method
      friend std::ostream & operator << (std::ostream &,const l_cinterval &) 
                                                                    throw();
      //! Implementation of standard output method
      friend std::string &  operator << (std::string &, const l_cinterval &)  
                                                                    throw();
      //! Implementation of standard input method
      friend std::string &  operator >> (std::string &, l_cinterval &)
                                                      throw(EMPTY_INTERVAL);
      //! Implementation of standard input method
      friend std::istream & operator >>(std::istream &, l_cinterval &) 
                                                      throw(EMPTY_INTERVAL);
      //! Implementation of standard input method
      friend void operator >> (const char *, l_cinterval &) 
                                                      throw(EMPTY_INTERVAL);
      //! Implementation of standard input method
      friend void operator >> (const std::string &, l_cinterval &)  
                                                      throw(EMPTY_INTERVAL);


// ----------------------------- Std.Operators ------------------------------

      //! Implementation of standard algebraic negative sign operation
      friend inline l_cinterval operator -(const l_cinterval &) throw ();
      //! Implementation of standard algebraic positive sign operation
      friend inline l_cinterval operator +(const l_cinterval &) throw ();
      //! Implementation of standard negation operation
      friend inline bool operator! (const l_cinterval & a)      throw ();

// LCI <--> LCI

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &, 
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend        l_cinterval operator *(const l_cinterval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend        l_cinterval operator /(const l_cinterval &,
                                           const l_cinterval &) 
                                                     throw(DIV_BY_ZERO);
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &,
                                           const l_cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &,
                                           const l_cinterval &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, 
                                              const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, 
                                              const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, 
                                              const l_cinterval &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, 
                                              const l_cinterval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, 
                                              const l_cinterval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, 
                                              const l_cinterval &) 
                                     throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      
      // LCI <--> R
      
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const real &, const l_cinterval &) 
                                                                     throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const real &, const l_cinterval &) 
                                                                     throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const real &, const l_cinterval &) 
                                                                     throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const real &, const l_cinterval &) 
                                                                     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &, const real &) 
                                                                     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const real &, const l_cinterval &) 
                                                                     throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &, const real &) 
                                                                     throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const real &, const l_cinterval &) 
                                                                     throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, const real &) 
                                                                     throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, const real &) 
                                                                     throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, const real &) 
                                                                     throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, const real &) 
                                                                     throw();

      // LCI <--> LR
      
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &, const l_real &)
                                                                       throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_real &, const l_cinterval &)
                                                                       throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &, const l_real &)
                                                                       throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_real &, const l_cinterval &)
                                                                       throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &, const l_real &)
                                                                       throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_real &, const l_cinterval &)
                                                                       throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &, const l_real &)
                                                                       throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_real &, const l_cinterval &)
                                                                       throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &, const l_real &)
                                                                       throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_real &, const l_cinterval &)
                                                                       throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &, const l_real &)
                                                                       throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_real &, const l_cinterval &)
                                                                       throw();
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, const l_real &) 
                                                                     throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, const l_real &) 
                                                                     throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, const l_real &) 
                                                                     throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, const l_real &) 
                                                                     throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, const l_real &) 
                                                                     throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, const l_real &) 
                                                                     throw();

      // LCI <--> C

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &,
                                           const complex &)     throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const complex &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &,
                                           const complex &)     throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const complex &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &,
                                           const complex &)     throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const complex &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &,
                                           const complex &)     throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const complex &,
                                           const l_cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &,
                                           const complex &)     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const complex &,
                                           const l_cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &,
                                           const complex &)     throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const complex &,
                                           const l_cinterval &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, const complex &) 
                                                                      throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, const complex &) 
                                                                      throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, const complex &)
                                                                      throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, const complex &) 
                                                                      throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, const complex &) 
                                                                      throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, const complex &) 
                                                                      throw();

      // LCI <--> LC

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &,
                                           const l_complex &)   throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_complex &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &,
                                           const l_complex &)   throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_complex &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &,
                                           const l_complex &)   throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_complex &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &,
                                           const l_complex &)   throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_complex &,
                                           const l_cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &,
                                           const l_complex &)   throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &,
                                           const l_cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &,
                                           const l_complex &)   throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_complex &,
                                           const l_cinterval &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, const l_complex &)
                                                                      throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, const l_complex &)
                                                                      throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, const l_complex &)
                                                                      throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, const l_complex &)
                                                                      throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, const l_complex &)
                                                                      throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, const l_complex &)
                                                                      throw();
      
      // LCI <--> I

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &,
                                           const interval &)    throw(); 
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const interval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &,
                                           const interval &)    throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const interval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &,
                                           const interval &)    throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const interval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &,
                                           const interval &)    throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const interval &,
                                           const l_cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &,
                                           const interval &)    throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const interval &,
                                           const l_cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &,
                                           const interval &)    throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const interval &,
                                           const l_cinterval &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, 
                                              const interval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, 
                                              const interval &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, 
                                              const interval &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, 
                                              const interval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, 
                                              const interval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, 
                                              const interval &) throw();

      // LCI <--> LI

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &,
                                           const l_interval &)  throw(); 
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_interval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &,
                                           const l_interval &)  throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_interval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &,
                                           const l_interval &)  throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_interval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &,
                                           const l_interval &)  throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_interval &,
                                           const l_cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &,
                                           const l_interval &)  throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_interval &,
                                           const l_cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &,
                                           const l_interval &)  throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_interval &,
                                           const l_cinterval &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, 
                                              const l_interval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, 
                                              const l_interval &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, 
                                              const l_interval &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, 
                                              const l_interval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, 
                                              const l_interval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, 
                                              const l_interval &) throw();

    // LCI <--> CI

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_cinterval &,
                                           const cinterval &)   throw(); 
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const cinterval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_cinterval &,
                                           const cinterval &)   throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const cinterval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_cinterval &,
                                           const cinterval &)   throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const cinterval &,
                                           const l_cinterval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_cinterval &,
                                           const cinterval &)   throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const cinterval &,
                                           const l_cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_cinterval &,
                                           const cinterval &)   throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const cinterval &,
                                           const l_cinterval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_cinterval &,
                                           const cinterval &)   throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const cinterval &,
                                           const l_cinterval &) throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline l_cinterval & operator +=(l_cinterval &, 
                                              const cinterval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline l_cinterval & operator -=(l_cinterval &, 
                                              const cinterval &) throw();
      //! Implementation of standard algebraic multiplication and allocation operation
      friend inline l_cinterval & operator *=(l_cinterval &, 
                                              const cinterval &) throw();
      //! Implementation of standard algebraic division and allocation operation
      friend inline l_cinterval & operator /=(l_cinterval &, 
                                              const cinterval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline l_cinterval & operator |=(l_cinterval &, 
                                              const cinterval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline l_cinterval & operator &=(l_cinterval &, 
                                              const cinterval &) throw();



// Convex Hull: LC <--> R; C <--> LR; LC <--> LR; CI <--> LR; CI <--> LC;

      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, const real &) 
                                                                   throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const real &, const l_complex &) 
                                                                   throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const complex &, const l_real &) 
                                                                   throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_real &, const complex &) 
                                                                   throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, const l_real &) 
                                                                     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_real &, const l_complex &) 
                                                                     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const cinterval &, const l_real &) 
                                                                     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_real &, const cinterval &) 
                                                                     throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const cinterval &, 
                                           const l_complex &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, 
                                           const cinterval &) throw();

      // LC <--> I

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_complex &, const interval &)
                                                                       throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const interval &, const l_complex &)
                                                                       throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_complex &, const interval &)
                                                                       throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const interval &, const l_complex &)
                                                                       throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_complex &, const interval &)
                                                                       throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const interval &, const l_complex &)
                                                                       throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_complex &, const interval &)
                                                                       throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const interval &, const l_complex &)
                                                                       throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, const interval &)
                                                                       throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const interval &, const l_complex &)
                                                                       throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_complex &, const interval &)
                                                                       throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const interval &, const l_complex &)
                                                                       throw();
      
      // C <--> LI

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const complex &, const l_interval &)
                                                                       throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_interval &, const complex &)
                                                                       throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const complex &, const l_interval &)
                                                                       throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_interval &, const complex &)
                                                                       throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const complex &, const l_interval &)
                                                                       throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_interval &, const complex &)
                                                                       throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const complex &, const l_interval &)
                                                                       throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_interval &, const complex &)
                                                                       throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const complex &, const l_interval &)
                                                                       throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_interval &, const complex &)
                                                                       throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const complex &, const l_interval &)
                                                                       throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_interval &, const complex &)
                                                                       throw();
      
      // LC <--> LI

      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_complex &, 
                                           const l_interval &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline l_cinterval operator +(const l_interval &, 
                                           const l_complex &)  throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_complex &, 
                                           const l_interval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline l_cinterval operator -(const l_interval &, 
                                           const l_complex &)  throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_complex &, 
                                           const l_interval &) throw();
      //! Implementation of standard algebraic multiplication operation
      friend inline l_cinterval operator *(const l_interval &, 
                                           const l_complex &)  throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_complex &, 
                                           const l_interval &) throw();
      //! Implementation of standard algebraic division operation
      friend inline l_cinterval operator /(const l_interval &, 
                                           const l_complex &)  throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, 
                                           const l_interval &) throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_interval &, 
                                           const l_complex &)  throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_complex &, 
                                           const l_interval &) throw();
      //! Returns the intersection of the arguments
      friend inline l_cinterval operator &(const l_interval &, 
                                           const l_complex &)  throw();
      

      // LC <--> C; LC <--> LC;  

      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, const complex &) 
                                                                      throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const complex &, const l_complex &) 
                                                                      throw();
      //! Returns the convex hull of the arguments
      friend inline l_cinterval operator |(const l_complex &, 
                                           const l_complex &) throw();

// --------------------- Comp.Operat.  ---------------------------------------
      // LCI <--> LCI      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, 
                                     const l_cinterval & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, 
                                     const l_cinterval & b) throw();
   
      // LCI <--> R
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, const real & b)
                                                                   throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real & a, const l_cinterval & b)
                                                                   throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, const real & b)
                                                                   throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real & a, const l_cinterval & b)
                                                                   throw();

      // LCI <--> LR
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, const l_real & b)
                                                                     throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_real & a, const l_cinterval & b)
                                                                     throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, const l_real & b)
                                                                     throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_real & a, const l_cinterval & b)
                                                                     throw();

      // LCI <--> I
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, 
                                     const interval & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const interval & a, 
                                     const l_cinterval & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, 
                                     const interval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const interval & a, 
                                     const l_cinterval & b) throw();

      // LCI <--> LI
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, 
                                     const l_interval & b)  throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_interval & a, 
                                     const l_cinterval & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, 
                                     const l_interval & b)  throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_interval & a, 
                                     const l_cinterval & b) throw();

      // LCI <--> C
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, const complex & b)
                                                                      throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const l_cinterval & b)
                                                                      throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, const complex & b)
                                                                      throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const l_cinterval & b)
                                                                      throw();

      // LCI <--> LC
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, 
                                     const l_complex & b)   throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_complex & a, 
                                     const l_cinterval & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, 
                                     const l_complex & b)   throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_complex & a, 
                                     const l_cinterval & b) throw();

      // LCI <--> CI
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const l_cinterval & a, 
                                     const cinterval & b)   throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, 
                                     const l_cinterval & b) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const l_cinterval & a, 
                                     const cinterval & b)   throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, 
                                     const l_cinterval & b) throw();

      // ---- Set Operators ----
      
      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const l_cinterval &)
                                                                      throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const l_cinterval &)
                                                                      throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const l_cinterval &)
                                                                      throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const l_cinterval &)
                                                                      throw();
      // LCI <--> R

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const real &, const l_cinterval &)
                                                               throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const real &, const l_cinterval &)
                                                               throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const real &, const l_cinterval &)
                                                               throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const real &, const l_cinterval &) 
                                                               throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const real &)
                                                               throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const real &)
                                                               throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const real &)
                                                               throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const real &)
                                                               throw();

      // LCI <--> LR

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_real &, const l_cinterval &)
                                                                 throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_real &, const l_cinterval &)
                                                                 throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_real &, const l_cinterval &)
                                                                 throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_real &, const l_cinterval &) 
                                                                 throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const l_real &)
                                                                 throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const l_real &)
                                                                 throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const l_real &)
                                                                 throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const l_real &)
                                                                 throw();

      // LCI <--> I

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const interval &, const l_cinterval &)
                                                                   throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const interval &, const l_cinterval &)
                                                                   throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const interval &, const l_cinterval &)
                                                                   throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const interval &, const l_cinterval &)
                                                                   throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const interval &)
                                                                   throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const interval &)
                                                                   throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const interval &)
                                                                   throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const interval &)
                                                                   throw();

      // LCI <--> LI

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_interval &, const l_cinterval &)
                                                                     throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_interval &, const l_cinterval &)
                                                                     throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_interval &, const l_cinterval &)
                                                                     throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_interval &, const l_cinterval &)
                                                                     throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const l_interval &)
                                                                     throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const l_interval &)
                                                                     throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const l_interval &)
                                                                     throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const l_interval &)
                                                                     throw();

      // LCI <--> C

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const complex &, const l_cinterval &)
                                                                  throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const complex &, const l_cinterval &)
                                                                  throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const complex &, const l_cinterval &)
                                                                  throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const complex &, const l_cinterval &)
                                                                  throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const complex &)
                                                                  throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const complex &)
                                                                  throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const complex &)
                                                                  throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const complex &)
                                                                  throw();

      // LCI <--> LC

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_complex &, const l_cinterval &)
                                                                    throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_complex &, const l_cinterval &)
                                                                    throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_complex &, const l_cinterval &)
                                                                    throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_complex &, const l_cinterval &)
                                                                    throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const l_complex &)
                                                                    throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const l_complex &)
                                                                    throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const l_complex &)
                                                                    throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const l_complex &)
                                                                    throw();

      // LCI <--> CI

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &, const l_cinterval &)
                                                                    throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &, const l_cinterval &)
                                                                    throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &, const l_cinterval &)
                                                                    throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cinterval &, const l_cinterval &)
                                                                    throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const l_cinterval &, const cinterval &)
                                                                    throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const l_cinterval &, const cinterval &)
                                                                    throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const l_cinterval &, const cinterval &)
                                                                    throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const l_cinterval &, const cinterval &)
                                                                    throw();

      // ---- Others   -------------------------------------------

      //! Returns the infimum of an interval
      friend inline l_complex    Inf(const l_cinterval &) throw();
      //! Returns the supremum of an interval
      friend inline l_complex    Sup(const l_cinterval &) throw();
      
      //! Returns the interval with the new given infimum value
      friend inline l_cinterval & SetInf(l_cinterval &, const complex &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the interval with the new given supremum value
      friend inline l_cinterval & SetSup(l_cinterval &, const complex &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the interval with the new given infimum value
      friend inline l_cinterval & SetInf(l_cinterval &, const l_complex &) 
                                    throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the interval with the new given supremum value
      friend inline l_cinterval & SetSup(l_cinterval &, const l_complex &) 
                                    throw(ERROR_CINTERVAL_EMPTY_INTERVAL);

      //! Returns the interval with the new given infimum value
      friend inline l_cinterval & SetInf(l_cinterval &, const real &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);

      //! Returns the interval with the new given supremum value
      friend inline l_cinterval & SetSup(l_cinterval &, const real &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);
      //! Returns the interval with the new given infimum value
      friend inline l_cinterval & SetInf(l_cinterval &, const l_real &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);

      //! Returns the interval with the new given supremum value
      friend inline l_cinterval & SetSup(l_cinterval &, const l_real &) 
                                  throw(ERROR_CINTERVAL_EMPTY_INTERVAL);

      //! Returns the interval with the unchecked new given infimum value
      friend inline l_cinterval & UncheckedSetInf(l_cinterval &, 
                                                  const complex &) throw();
      //! Returns the interval with the unchecked new given infimum value
      friend inline l_cinterval & UncheckedSetInf(l_cinterval &, 
                                                  const real &)    throw();
      //! Returns the interval with the unchecked new given supremum value
      friend inline l_cinterval & UncheckedSetSup(l_cinterval &, 
                                                  const complex &) throw();
      //! Returns the interval with the unchecked new given supremum value
      friend inline l_cinterval & UncheckedSetSup(l_cinterval &, 
                                                  const real &)    throw();

      //! Returns the interval with the unchecked new given infimum value
      friend inline l_cinterval & UncheckedSetInf(l_cinterval &, 
                                                  const l_complex &) throw();
      //! Returns the interval with the unchecked new given infimum value
      friend inline l_cinterval & UncheckedSetInf(l_cinterval &, 
                                                  const l_real &)    throw();
      //! Returns the interval with the unchecked new given supremum value
      friend inline l_cinterval & UncheckedSetSup(l_cinterval &, 
                                                  const l_complex &) throw();
      //! Returns the interval with the unchecked new given supremum value
      friend inline l_cinterval & UncheckedSetSup(l_cinterval &, 
                                                  const l_real &)    throw();
     
      //! Returns the real part of the complex interval
      friend l_interval & Re(l_cinterval & a)       throw() { return a.re; }
      //! Returns the real part of the complex interval
      friend l_interval   Re(const l_cinterval & a) throw() { return a.re; }
      //! Returns the imaginary part of the complex interval
      friend l_interval & Im(l_cinterval & a)       throw() { return a.im; }
      //! Returns the imaginary part of the complex interval
      friend l_interval   Im(const l_cinterval & a) throw() { return a.im; }
      
      //! Sets the real part of a complex interval
      friend l_cinterval & SetRe(l_cinterval & a, const interval & b) 
                                                       { a.re=b; return a; }
      //! Sets the imaginary part of a complex interval
      friend l_cinterval & SetIm(l_cinterval & a, const interval & b) 
                                                       { a.im=b; return a; } 
      //! Sets the real part of a complex interval
      friend l_cinterval & SetRe(l_cinterval & a, const real     & b) 
                                                       { a.re=b; return a; }
      //! Sets the imaginary part of a complex interval
      friend l_cinterval & SetIm(l_cinterval & a, const real     & b) 
                                                       { a.im=b; return a; }
 
      //! Sets the real part of a complex interval
      friend l_cinterval & SetRe(l_cinterval & a, const l_interval & b) 
                                                       { a.re=b; return a; }
      //! Sets the imaginary part of a complex interval
      friend l_cinterval & SetIm(l_cinterval & a, const l_interval & b) 
                                                       { a.im=b; return a; } 
      //! Sets the real part of a complex interval
      friend l_cinterval & SetRe(l_cinterval & a, const l_real     & b) 
                                                       { a.re=b; return a; }
      //! Sets the imaginary part of a complex interval
      friend l_cinterval & SetIm(l_cinterval & a, const l_real     & b) 
                                                       { a.im=b; return a; } 

      //! Returns the infimum of the real interval of the complex interval
      friend inline l_real InfRe(const l_cinterval &a) throw() 
                                                     { return Inf(a.re); }
      //! Returns the infimum of the imaginary interval of the complex interval
      friend inline l_real InfIm(const l_cinterval &a) throw() 
                                                     { return Inf(a.im); }
      //! Returns the supremum of the real interval of the complex interval
      friend inline l_real SupRe(const l_cinterval &a) throw() 
                                                     { return Sup(a.re); }
      //! Returns the supremum of the imaginary interval of the complex interval
      friend inline l_real SupIm(const l_cinterval &a) throw() 
                                                     { return Sup(a.im); }
    //! Multiplication of interval with \f$ 2^n \f$
    friend inline void times2pown(l_cinterval& x, 
                              const int& n) throw(); 
    //! Multiplication of interval with \f$ 2^n \f$
    friend inline void Times2pown(l_cinterval& x, 
                              const int& n) throw(); 

     
    //! Returns the absolute value of the complex interval
    friend        l_interval  abs (const l_cinterval &) throw();
    //! Returns the conjugated complex interval
    friend inline l_cinterval conj(const l_cinterval &) throw();
    //! Returns the rounded middle of the complex interval
    friend inline l_complex   mid (const l_cinterval &) throw();
    //! Returns the rounded diameter of the complex interval
    friend inline l_complex   diam(const l_cinterval &) throw();

    //! Sets the precision of a specific long datatype value
    friend inline l_cinterval adjust(const l_cinterval &) throw();

};

} // namespace cxsc 

#include "l_cinterval.inl"
#include "l_cimath.hpp"

#endif // _CXSC_L_CINTERVAL_HPP_INCLUDED
