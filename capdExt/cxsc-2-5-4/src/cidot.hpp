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
               inline cidotprecision(const cidotprecision &);
      
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const real &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const dotprecision &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const interval &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const idotprecision &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const complex &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const cdotprecision &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const cinterval &);
      //! Constructor of class cidotprecision
      explicit inline cidotprecision(const idotprecision &, const idotprecision &);

      //! Get currently set precision for computation of dot products
      inline int get_k() const { return k; }
      //! Set precision for computation of dot products
      inline void set_k(unsigned int i) { k=i; reinf.set_k(i); resup.set_k(i); iminf.set_k(i); imsup.set_k(i);}
      //! Get currently set precision for computation of dot products
      inline int get_dotprec() const { return k; }
      //! Set precision for computation of dot products
      inline void set_dotprec(unsigned int i) { k=i; reinf.set_k(i); resup.set_k(i); iminf.set_k(i); imsup.set_k(i);}
                  
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const real & a)         { reinf=resup=a; iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const complex & a)      { reinf=resup=Re(a); iminf=imsup=Im(a); return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const interval & a)     { reinf=Inf(a),resup=Sup(a),iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const cinterval & a)    { reinf=Inf(Re(a)),resup=Sup(Re(a)),iminf=Inf(Im(a)),imsup=Sup(Im(a)); return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const dotprecision & a) { reinf=resup=a; iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const cdotprecision & a) { reinf=resup=Re(a),iminf=imsup=Im(a); return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const idotprecision & a) { reinf=Inf(a),resup=Sup(a),iminf=imsup=0.0; return *this; }
      //! Implementation of standard assigning operator
      inline cidotprecision & operator= (const cidotprecision& a) { reinf=a.reinf,resup=a.resup,iminf=a.iminf,imsup=a.imsup; return *this; }

      // ---- Destruktor    ----
      // ~cidotprecision() {} unnoetig

      // ---- Typwandlungen ----
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const complex &,const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const complex &,const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const real &,const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const interval &,const interval &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const interval &,const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const real &,const interval &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const interval &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cinterval &);
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const idotprecision &,const idotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cdotprecision &,const cdotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const idotprecision &,const dotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cdotprecision &,const dotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const dotprecision &,const idotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const dotprecision &,const cdotprecision&);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const cdotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const idotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _cidotprecision(const dotprecision &);
      
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const complex &, const complex &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const complex &, const real &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const real &, const complex &);
      
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const cdotprecision &, const cdotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const cdotprecision &, const dotprecision &);
      //! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
      friend inline cidotprecision _unchecked_cidotprecision(const dotprecision &, const cdotprecision &);

      // ---- Ausgabefunkt. ---------------------------------------
      //! Implementation of standard input method
      friend std::istream& operator >> (std::istream& s, cidotprecision& a)      ;
      //! Implementation of standard output method
      friend std::ostream& operator << (std::ostream& s, const cidotprecision& a);
      //! Implementation of standard input method
      friend std::string&  operator >> (std::string&  s, cidotprecision& a)      ;
      //! Implementation of standard output method
      friend std::string&  operator << (std::string&  s, const cidotprecision& a);
      //! Implementation of standard input method
      friend void          operator >> (const std::string &s,cidotprecision& a)  ;
      //! Implementation of standard input method
      friend void          operator >> (const char *s       ,cidotprecision& a)  ;

      // ---- Standardfunkt ---- (arithmetische Operatoren)
      //! Implementation of standard algebraic negative sign operation
      friend inline cidotprecision operator -(cidotprecision);
      //! Implementation of standard algebraic positive sign operation
      friend inline cidotprecision operator +(const cidotprecision &);

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const cidotprecision &);

      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cdotprecision &,const cdotprecision &) ;

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const cidotprecision &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const cidotprecision &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const cidotprecision &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const cidotprecision &);

      // CID-R
      
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const real &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const real &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const real &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const real &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const real &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const real &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const real &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const real &,const cidotprecision &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const real &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const real &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const real &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const real &);
      
      // CID-C

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const complex &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const complex &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const complex &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const complex &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const complex &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const complex &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const complex &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const complex &,const cidotprecision &);
      

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const complex &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const complex &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const complex &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const complex &);
      
      // CID-I

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const interval &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const interval &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const interval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const interval &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const interval &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const interval &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const interval &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const interval &,const cidotprecision &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const interval &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const interval &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const interval &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const interval &);

      // CID-CI

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const cinterval &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cinterval &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const cinterval &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cinterval &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const cinterval &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cinterval &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const cinterval &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cinterval &,const cidotprecision &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const cinterval &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const cinterval &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const cinterval &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const cinterval &);
      
      // CID-D
      
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const dotprecision &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const dotprecision &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const dotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const dotprecision &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const dotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const dotprecision &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const dotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const dotprecision &,const cidotprecision &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const dotprecision &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const dotprecision &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const dotprecision &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const dotprecision &);
      
      // CID-CD

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const cdotprecision &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cdotprecision &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const cdotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cdotprecision &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const cdotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cdotprecision &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const cdotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cdotprecision &,const cidotprecision &);
      

      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const cdotprecision &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const cdotprecision &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const cdotprecision &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const cdotprecision &);
      
      // CID-ID

      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const cidotprecision &,const idotprecision &);
      //! Implementation of standard algebraic addition operation
      friend inline cidotprecision operator +(const idotprecision &,const cidotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const cidotprecision &,const idotprecision &);
      //! Implementation of standard algebraic subtraction operation
      friend inline cidotprecision operator -(const idotprecision &,const cidotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const cidotprecision &,const idotprecision &);
      //! Returns the convex hull of the arguments
      friend inline cidotprecision operator |(const idotprecision &,const cidotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const cidotprecision &,const idotprecision &);
      //! Returns the intersection of the arguments
      friend inline cidotprecision operator &(const idotprecision &,const cidotprecision &);
      
      //! Implementation of standard algebraic addition and allocation operation
      friend inline cidotprecision & operator +=(cidotprecision &, const idotprecision &);
      //! Implementation of standard algebraic subtraction and allocation operation
      friend inline cidotprecision & operator -=(cidotprecision &, const idotprecision &);
      //! Allocates the convex hull of the arguments to the first argument
      friend inline cidotprecision & operator |=(cidotprecision &, const idotprecision &);
      //! Allocates the intersection of the arguments to the first argument
      friend inline cidotprecision & operator &=(cidotprecision &, const idotprecision &);

      // ---- Vergleichsop. ----
      //! Implementation of standard negation operation
      friend inline bool operator !(const cidotprecision &);
//             inline      operator void *() const;

      //! Implementation of standard equality operation
      friend inline bool operator ==(const cidotprecision &,const cidotprecision &);
      //! Implementation of standard negated equality operation
      friend inline bool operator !=(const cidotprecision &,const cidotprecision &);

      // CID-R
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const real & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const real & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const real & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const real & a, const cidotprecision & b)   ;

      // CID-C
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const complex & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const complex & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const complex & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const complex & a, const cidotprecision & b)   ;

      // CID-I
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const interval & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const interval & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const interval & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const interval & a, const cidotprecision & b)   ;

      // CID-CI
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const cinterval & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const cinterval & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const cinterval & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cinterval & a, const cidotprecision & b)   ;
      
      // CID-D
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const dotprecision & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const dotprecision & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const dotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const dotprecision & a, const cidotprecision & b)   ;

      // CID-CD
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const cdotprecision & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const cdotprecision & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const cdotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cdotprecision & a, const cidotprecision & b)   ;

      // CID-ID
      
      //! Implementation of standard equality operation
      friend inline bool operator== (const cidotprecision & a, const idotprecision & b)   ;
      //! Implementation of standard equality operation
      friend inline bool operator== (const idotprecision & a, const cidotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const cidotprecision & a, const idotprecision & b)   ;
      //! Implementation of standard negated equality operation
      friend inline bool operator!= (const idotprecision & a, const cidotprecision & b)   ;

      // ---- Set Operators ----
      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const cidotprecision &);

      // CID-R

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const real &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const real &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const real &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const real &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const real &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const real &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const real &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const real &);

      // CID-C

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const complex &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const complex &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const complex &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const complex &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const complex &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const complex &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const complex &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const complex &);

      // CID-I

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const interval &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const interval &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const interval &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const interval &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const interval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const interval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const interval &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const interval &);

      // CID-CI

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cinterval &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cinterval &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cinterval &,const cidotprecision &);
      friend inline bool operator >=(const cinterval &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const cinterval &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const cinterval &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const cinterval &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const cinterval &);

      // CID-D

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const dotprecision &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const dotprecision &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const dotprecision &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const dotprecision &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const dotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const dotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const dotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const dotprecision &);

      // CID-CD

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cdotprecision &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cdotprecision &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cdotprecision &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cdotprecision &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const cdotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const cdotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const cdotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const cdotprecision &);

      // CID-ID

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const idotprecision &,const cidotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const idotprecision &,const cidotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const idotprecision &,const cidotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const idotprecision &,const cidotprecision &);

      //! Implementation of standard less-than operation
      friend inline bool operator  <(const cidotprecision &,const idotprecision &);
      //! Implementation of standard greater-than operation
      friend inline bool operator  >(const cidotprecision &,const idotprecision &);
      //! Implementation of standard less-or-equal-than operation
      friend inline bool operator <=(const cidotprecision &,const idotprecision &);
      //! Implementation of standard greater-or-equal-than operation
      friend inline bool operator >=(const cidotprecision &,const idotprecision &);


      // ---- Funktionen    ----
      
      //! Returns the infimum of a complex dotprecison interval
      friend inline cdotprecision   Inf(const cidotprecision&);
      //! Returns the supremum of a complex dotprecison interval
      friend inline cdotprecision   Sup(const cidotprecision&);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const cdotprecision&);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const cdotprecision&);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const dotprecision&);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const dotprecision&);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const complex&);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const complex&);
      //! Returns the complex dotprecison interval with the new given infimum value
      friend inline cidotprecision& SetInf(cidotprecision&, const real&);
      //! Returns the complex dotprecison interval with the new given supremum value
      friend inline cidotprecision& SetSup(cidotprecision&, const real&);
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const cdotprecision&);
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const cdotprecision&);
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const dotprecision&);
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const dotprecision&);
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const complex&);
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const complex&);
      //! Returns the complex dotprecison interval with the unchecked new given infimum value
      friend inline cidotprecision& UncheckedSetInf(cidotprecision&, const real&);
      //! Returns the complex dotprecison interval with the unchecked new given supremum value
      friend inline cidotprecision& UncheckedSetSup(cidotprecision&, const real&);
   
      //! Returns the real part of the complex dotprecision interval
      friend inline idotprecision   Re(const cidotprecision &);
      //! Returns the imaginary part of the complex dotprecision interval
      friend inline idotprecision   Im(const cidotprecision &);
      
      //! Returns the infimum of the real part of the complex dotprecision interval
      friend inline const dotprecision & InfRe(const cidotprecision &);
      //! Returns the infimum of the imaginary part of the complex dotprecision interval
      friend inline const dotprecision & InfIm(const cidotprecision &);
      //! Returns the supremum of the real part of the complex dotprecision interval
      friend inline const dotprecision & SupRe(const cidotprecision &);
      //! Returns the supremum of the imaginary part of the complex dotprecision interval
      friend inline const dotprecision & SupIm(const cidotprecision &); 
      
      //! Returns the infimum of the real part of the complex dotprecision interval
      friend inline       dotprecision & InfRe(cidotprecision &);
      //! Returns the infimum of the imaginary part of the complex dotprecision interval
      friend inline       dotprecision & InfIm(cidotprecision &);
      //! Returns the supremum of the real part of the complex dotprecision interval
      friend inline       dotprecision & SupRe(cidotprecision &);
      //! Returns the supremum of the imaginary part of the complex dotprecision interval
      friend inline       dotprecision & SupIm(cidotprecision &);
      
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const idotprecision&);
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const idotprecision&);
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const dotprecision&);
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const dotprecision&);
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const interval&);
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const interval&);
      //! Sets the real part of the complex dotprecision interval
      friend inline cidotprecision& SetRe(cidotprecision&, const real&);
      //! Sets the imaginary part of the complex dotprecision interval
      friend inline cidotprecision& SetIm(cidotprecision&, const real&);

      
      friend inline void rnd(const cidotprecision &,cinterval &);
      friend inline cinterval rnd(const cidotprecision &);
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend        void accumulate  (cidotprecision&, const cinterval&, const cinterval&);

      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const cinterval&, const interval&);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const cinterval&, const complex&);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const cinterval&, const real&);
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const cinterval &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const cinterval &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const cinterval&);
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const interval &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const complex &);

      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const interval &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const interval &,const real &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const interval &);

      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const complex &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const complex &);
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const complex &,const real &);
      
      //! The accurate scalar product of the last two arguments added to the value of the first argument
      friend inline void accumulate  (cidotprecision&, const real &,const real &);
};

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- friend functions of class cdotprecision (not inline)              ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

std::istream& operator >> (std::istream& s, cidotprecision& a)      ;
std::ostream& operator << (std::ostream& s, const cidotprecision& a);
std::string&  operator >> (std::string&  s, cidotprecision& a)      ;
std::string&  operator << (std::string&  s, const cidotprecision& a);
void          operator >> (const std::string &s,cidotprecision& a)  ;
void          operator >> (const char *s       ,cidotprecision& a)  ;

void accumulate  (cidotprecision&, const cinterval&, const cinterval&);

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

