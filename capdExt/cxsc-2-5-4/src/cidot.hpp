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

/* CVS $Id: cidot.hpp,v 1.33 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_CIDOT_HPP_INCLUDED
#define _CXSC_CIDOT_HPP_INCLUDED

#include <iostream>
#include <string>
#include "dot.hpp"
#include "idot.hpp"
#include "cdot.hpp"
#include "cinterval.hpp"

namespace cxsc {

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- class cidotprecision (declaration)                                ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! The Data Type cidotprecision
/*!
The data types dotprecision, idotprecision, cdotprecision and cidotprecision
are based on the scalar data types real, interval, complex, and cinterval, 
respectively. Variables of these data types permit the exact representation 
of products of two arbitrary numbers of the corresponding scalar base type 
and the exact summation of an arbitrary number of such products in a 
dotprecision accumulator, i.e. in a fixed-point format of suitable size.

Since C-XSC Version 2.3.0, the precision for dot products computed with the dotprecision types can be choosen by the user. The default is precision k=0, which means maximum precision (the behaviour of older C-XSC Versions). For k=1, pure floating point operations are used and an error bound is computed using manipulation of the rounding mode of the processor. For k>=2, the so called DotK algorithm is used, simulating higher precision computations and also computing an error bound. When calling the rnd() function, the error bound will be added to the result interval. The resulting intervals will be wider for lower k, but computations will be significantly faster than with maximum precision (k=0).

\sa cxsc::dotprecision
*/
class cidotprecision
{
   private:
      // ---- Datenelemente ---------------------------------------
      dotprecision reinf,resup,iminf,imsup;
      int k;

   public:
      // ---- Constructors  ---------------------------------------
      //! Constructor of class cidotprecision
      cidotprecision() : k(0) {}
      //! Constructor of class cidotprecision
               inline cidotprecision(const cidotprecision &) throw();
      
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const real &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const dotprecision &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const interval &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const idotprecision &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const complex &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const cdotprecision &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const cinterval &) throw();
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const idotprecision &, const idotprecision &) throw();

      //! Get currently set precision for computation of dot products
      inline int get_k() const { return k; }
      //! Set precision for computation of dot products
      inline void set_k(unsigned int i) { k=i; reinf.set_k(i); resup.set_k(i); iminf.set_k(i); imsup.set_k(i);}
      //! Get currently set precision for computation of dot products
      inline int get_dotprec() const { return k; }
      //! Set precision for computation of dot products
      inline void set_dotprec(unsigned int i) { k=i; reinf.set_k(i); resup.set_k(i); iminf.set_k(i); imsup.set_k(i);}
                  
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const real & a)         throw() { reinf=resup=a; iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const complex & a)      throw() { reinf=resup=Re(a); iminf=imsup=Im(a); return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const interval & a)     throw() { reinf=Inf(a),resup=Sup(a),iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const cinterval & a)    throw() { reinf=Inf(Re(a)),resup=Sup(Re(a)),iminf=Inf(Im(a)),imsup=Sup(Im(a)); return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const dotprecision & a) throw() { reinf=resup=a; iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const cdotprecision & a)throw() { reinf=resup=Re(a),iminf=imsup=Im(a); return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const idotprecision & a)throw() { reinf=Inf(a),resup=Sup(a),iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const cidotprecision& a)throw() { reinf=a.reinf,resup=a.resup,iminf=a.iminf,imsup=a.imsup; return *this; }

      // ---- Destruktor    ----
      // ~cidotprecision() {} unnoetig

      // ---- Typwandlungen ----
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const complex &,const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const complex &,const real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const real &,const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const interval &,const interval &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const interval &,const real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const real &,const interval &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const interval &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cinterval &) throw();
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const idotprecision &,const idotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cdotprecision &,const cdotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const idotprecision &,const dotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cdotprecision &,const dotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const dotprecision &,const idotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const dotprecision &,const cdotprecision&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cdotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const idotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const dotprecision &) throw();
      
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const complex &, const complex &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const complex &, const real &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const real &, const complex &) throw();
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const cdotprecision &, const cdotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const cdotprecision &, const dotprecision &) throw();
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const dotprecision &, const cdotprecision &) throw();

      // ---- Ausgabefunkt. ---------------------------------------
      //! Implementation of standard input method
      friend std::istream& operator >> (std::istream& s, cidotprecision& a)       throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Implementation of standard output method
      friend std::ostream& operator << (std::ostream& s, const cidotprecision& a) throw();
      //! Implementation of standard input method
      friend std::string&  operator >> (std::string&  s, cidotprecision& a)       throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Implementation of standard output method
      friend std::string&  operator << (std::string&  s, const cidotprecision& a) throw();
      //! Implementation of standard input method
      friend void          operator >> (const std::string &s,cidotprecision& a)   throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Implementation of standard input method
      friend void          operator >> (const char *s       ,cidotprecision& a)   throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);

      // ---- Standardfunkt ---- (arithmetische Operatoren)
      //! Implementation of standard algebraic negative sign operation
      friend inline cidotprecision operator -(cidotprecision) throw();
      //! Implementation of standard algebraic positive sign operation
      friend inline cidotprecision operator +(const cidotprecision &) throw();

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);

      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cdotprecision &,const cdotprecision &)  throw();

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const cidotprecision &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const cidotprecision &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);

      // CID-R
      
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const real &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const real &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const real &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const real &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const real &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const real &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const real &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const real &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const real &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const real &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const real &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const real &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      // CID-C

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const complex &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const complex &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const complex &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const complex &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const complex &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const complex &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const complex &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const complex &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const complex &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const complex &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const complex &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const complex &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      // CID-I

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const interval &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const interval &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const interval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const interval &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const interval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const interval &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const interval &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const interval &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const interval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const interval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const interval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const interval &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);

      // CID-CI

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const cinterval &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cinterval &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const cinterval &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cinterval &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const cinterval &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cinterval &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const cinterval &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cinterval &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const cinterval &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const cinterval &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const cinterval &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const cinterval &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      // CID-D
      
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const dotprecision &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const dotprecision &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const dotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const dotprecision &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const dotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const dotprecision &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const dotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const dotprecision &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const dotprecision &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const dotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      // CID-CD

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cdotprecision &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cdotprecision &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const cdotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cdotprecision &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const cdotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cdotprecision &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const cdotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const cdotprecision &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const cdotprecision &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const cdotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      // CID-ID

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const idotprecision &) throw();
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const idotprecision &,const cidotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const idotprecision &) throw();
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const idotprecision &,const cidotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const idotprecision &) throw();
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const idotprecision &,const cidotprecision &) throw();
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const idotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const idotprecision &,const cidotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const idotprecision &) throw();
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const idotprecision &) throw();
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const idotprecision &) throw();
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const idotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);

      // ---- Vergleichsop. ----
      //! Implementation of standard negation operation
      friend inline bool operator !(const cidotprecision &) throw();
//             inline      operator void *() const throw();

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cidotprecision &,const cidotprecision &) throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cidotprecision &,const cidotprecision &) throw();

      // CID-R
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const real & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const real & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const real & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real & a, const cidotprecision & b)    throw();

      // CID-C
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const complex & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const complex & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const cidotprecision & b)    throw();

      // CID-I
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const interval & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const interval & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const interval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const interval & a, const cidotprecision & b)    throw();

      // CID-CI
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const cinterval & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const cinterval & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const cidotprecision & b)    throw();
      
      // CID-D
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const dotprecision & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const dotprecision & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const dotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const dotprecision & a, const cidotprecision & b)    throw();

      // CID-CD
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const cdotprecision & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const cdotprecision & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const cdotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cdotprecision & a, const cidotprecision & b)    throw();

      // CID-ID
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const idotprecision & b)    throw();
      //! Implementation of standard equality operation
      friend inline bool operator== (const idotprecision & a, const cidotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const idotprecision & b)    throw();
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const idotprecision & a, const cidotprecision & b)    throw();

      // ---- Set Operators ----
      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const cidotprecision &) throw();

      // CID-R

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const real &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const real &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const real &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const real &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const real &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const real &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const real &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const real &) throw();

      // CID-C

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const complex &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const complex &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const complex &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const complex &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const complex &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const complex &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const complex &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const complex &) throw();

      // CID-I

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const interval &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const interval &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const interval &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const interval &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const interval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const interval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const interval &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const interval &) throw();

      // CID-CI

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const cidotprecision &) throw();
      friend inline bool operator >=(const cinterval &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const cinterval &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const cinterval &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const cinterval &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const cinterval &) throw();

      // CID-D

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const dotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const dotprecision &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const dotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const dotprecision &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const dotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const dotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const dotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const dotprecision &) throw();

      // CID-CD

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cdotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cdotprecision &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cdotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cdotprecision &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const cdotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const cdotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const cdotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const cdotprecision &) throw();

      // CID-ID

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const idotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const idotprecision &,const cidotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const idotprecision &,const cidotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const idotprecision &,const cidotprecision &) throw();

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const idotprecision &) throw();
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const idotprecision &) throw();
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const idotprecision &) throw();
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const idotprecision &) throw();


      // ---- Funktionen    ----
      
      //! Returns the infimum of a complex dotprecison interval
      friend inline cdotprecision   Inf(const cidotprecision&) throw();
      //! Returns the supremum of a complex dotprecison interval
      friend inline cdotprecision   Sup(const cidotprecision&) throw();
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const cdotprecision&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const cdotprecision&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const dotprecision&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const dotprecision&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const complex&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const complex&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const real&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const real&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const cdotprecision&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const cdotprecision&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const dotprecision&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const dotprecision&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const complex&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const complex&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const real&) throw();
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const real&) throw();
   
      //! Returns the real part of the complex dotprecision interval
      friend inline idotprecision   Re(const cidotprecision &) throw();
      //! Returns the imaginary part of the complex dotprecision interval
      friend inline idotprecision   Im(const cidotprecision &) throw();
      
      //! Returns the infimum of the real part of the complex dotprecision interval
      friend inline const dotprecision & InfRe(const cidotprecision &) throw();
      //! Returns the infimum of the imaginary part of the complex dotprecision interval
      friend inline const dotprecision & InfIm(const cidotprecision &) throw();
      //! Returns the supremum of the real part of the complex dotprecision interval
      friend inline const dotprecision & SupRe(const cidotprecision &) throw();
      //! Returns the supremum of the imaginary part of the complex dotprecision interval
      friend inline const dotprecision & SupIm(const cidotprecision &) throw(); 
      
      //! Returns the infimum of the real part of the complex dotprecision interval
      friend inline       dotprecision & InfRe(cidotprecision &) throw();
      //! Returns the infimum of the imaginary part of the complex dotprecision interval
      friend inline       dotprecision & InfIm(cidotprecision &) throw();
      //! Returns the supremum of the real part of the complex dotprecision interval
      friend inline       dotprecision & SupRe(cidotprecision &) throw();
      //! Returns the supremum of the imaginary part of the complex dotprecision interval
      friend inline       dotprecision & SupIm(cidotprecision &) throw();
      
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const idotprecision&) throw();
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const idotprecision&) throw();
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const dotprecision&) throw();
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const dotprecision&) throw();
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const interval&) throw();
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const interval&) throw();
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const real&) throw();
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const real&) throw();

      
      friend inline void rnd(const cidotprecision &,cinterval &) throw();
      friend inline cinterval rnd(const cidotprecision &) throw();
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend        void accumulate  (cidotprecision&, const cinterval&, const cinterval&) throw();

      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const cinterval&, const interval&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const cinterval&, const complex&) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const cinterval&, const real&) throw();
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const cinterval &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const cinterval &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const cinterval&) throw();
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const interval &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const complex &) throw();

      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const interval &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const real &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const interval &) throw();

      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const complex &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const complex &) throw();
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const real &) throw();
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const real &) throw();
};

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- friend functions of class cdotprecision (not inline)              ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

std::istream& operator >> (std::istream& s, cidotprecision& a)       throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
std::ostream& operator << (std::ostream& s, const cidotprecision& a) throw();
std::string&  operator >> (std::string&  s, cidotprecision& a)       throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
std::string&  operator << (std::string&  s, const cidotprecision& a) throw();
void          operator >> (const std::string &s,cidotprecision& a)   throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
void          operator >> (const char *s       ,cidotprecision& a)   throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);

void accumulate  (cidotprecision&, const cinterval&, const cinterval&) throw();

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ----  global CIDotprecision Akku's                                     ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//#define MAXCIDOTAKKU     (MAXDOTAKKU / 2)
//extern cidotprecision cidotakku[MAXCIDOTAKKU];
  
//----------------------------------------------------------------------
} // namespace cxsc 

#include "cidot.inl"

#endif // _CXSC_CIDOT_HPP_INCLUDED

