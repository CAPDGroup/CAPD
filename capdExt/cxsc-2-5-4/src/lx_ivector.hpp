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

/* CVS $Id: lx_ivector.hpp,v 1.9 2014/01/30 17:23:47 cxsc Exp $ */

#ifndef _CXSC_LX_IVECTOR_HPP_INCLUDED
#define _CXSC_LX_IVECTOR_HPP_INCLUDED 

#include <except.hpp>
#include "lx_interval.hpp"
#include <iostream>

namespace cxsc {
	
//! The Multiple-Precision Data Type lx_ivector
/*!
	The vectors of C-XSC are one dimensional arrays of the corresponding scalar base type. 

	\sa lx_ivector
 */

class lx_ivector
{	
	private:
		lx_interval *dat;
		int l,u,size;
		
	public:
	//------ Konstruktoren ----------------------------------------------------
		
	//! Constructor of class lx_ivector
	inline lx_ivector () throw();
		
	/*!
	\param i Dimension of vector
	Creation of a variable of type lx_ivector with length \f$ n = i \f$ and index bounds \f$ lb = 1 \f$, and \f$ ub = i \f$. The values of the elements are undefined.
	*/
	explicit inline lx_ivector(int i) throw();
		
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(int i1, int i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_WRONG_BOUNDARIES,ERROR_IVECTOR_NO_MORE_MEMORY);
#else
		throw();
#endif
		
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(const lx_interval &) throw();
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(const l_interval &)  throw();
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(const interval &)    throw();
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(const lx_real &)     throw();
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(const l_real &)      throw();
	//! Constructor of class lx_ivector
	explicit inline lx_ivector(const real &)        throw();
	
	//! Constructor of class lx_ivector
	inline lx_ivector(const lx_ivector &) throw();
		
	
	//! Implementation of standard assigning operator
	inline lx_ivector & operator = (const lx_ivector &) throw();
	//! Implementation of standard assigning operator
	inline lx_ivector & operator =(const lx_interval &) throw();
	//! Implementation of standard assigning operator
	inline lx_ivector & operator =(const l_interval &) throw();
	//! Implementation of standard assigning operator
	inline lx_ivector & operator =(const interval &) throw();
	//! Implementation of standard assigning operator
	inline lx_ivector & operator =(const lx_real &) throw();
	//! Implementation of standard assigning operator
	inline lx_ivector & operator =(const l_real &) throw();
	//! Implementation of standard assigning operator
	inline lx_ivector & operator =(const real &) throw();
	
	//--------- Destruktor ----------------------------------------------------
	inline ~lx_ivector() { delete [] dat; }
	
	
	//! Operator for accessing the single elements of the vector
	inline lx_interval & operator [](const int &i)
#if(CXSC_INDEX_CHECK)
			throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	//! Operator for accessing the single elements of the vector
	inline const lx_interval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	
//------ Standardfunktionen -----------------------------------------------
	
//! Returns the lower bound of the vector
friend inline int Lb(const lx_ivector &a) throw() { return a.l; }
//! Returns the upper bound of the vector
friend inline int Ub(const lx_ivector &a) throw() { return a.u; }
//! Returns the dimension of the vector
friend inline int VecLen(const lx_ivector &a) throw() { return a.size; }
//! Sets the lower bound of the vector
friend inline lx_ivector& SetLb(lx_ivector &a, int l) throw() 
{ a.l=l; a.u=l+a.size-1; return a; }
//! Sets the upper bound of the vector
friend inline lx_ivector & SetUb(lx_ivector &a, int u) throw()
{ a.u=u; a.l=u-a.size+1; return a; }

//! Resizes the vector
friend inline void Resize(lx_ivector &rv, int lb, int ub)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<lx_ivector>);
#else
		throw();
#endif

//! Resizes the vector
friend inline void Resize(lx_ivector &rv, int len)
#if(CXSC_INDEX_CHECK)
				throw(ERROR__WRONG_BOUNDARIES<lx_ivector>);
#else
		throw();
#endif		
	
}; // End of class lx_ivector

//! Doubles the vector size
inline void DoubleSize(lx_ivector&) throw();

inline void Resize(lx_ivector &rv, int lb, int ub)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<lx_ivector>);
#else
		throw();
#endif
		
		inline void Resize(lx_ivector &rv, int len)
#if(CXSC_INDEX_CHECK)
				throw(ERROR__WRONG_BOUNDARIES<lx_ivector>);
#else
		throw();
#endif		

} // End namespace cxsc
 

#include "lx_ivector.inl"
					
#endif
