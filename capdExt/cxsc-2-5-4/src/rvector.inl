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

/* CVS $Id: rvector.inl,v 1.26 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_RVECTOR_INL_INCLUDED
#define _CXSC_RVECTOR_INL_INCLUDED

#include "rvector.hpp"
#include "intvector.hpp"

namespace cxsc {

	/*!
	Creation of a variable of type rvector with length \f$ n = 1 \f$ and index bounds \f$ lb = ub = 1 \f$. The value of the element is undefined.
	*/
	INLINE rvector::rvector () throw():dat(NULL),l(1),u(0),size(0)
	{
	}

	/*!
	\param i Dimension of vector

	Creation of a variable of type rvector with length \f$ n = i \f$ and index bounds \f$ lb = 1 \f$, and \f$ ub = i \f$. The values of the elements are undefined.
	*/
	INLINE rvector::rvector(const int &i) throw():l(1),u(i),size(i)
	{
		dat=new real[i];
	}

#ifdef OLD_CXSC  
	INLINE rvector::rvector(const class index &i) throw():l(1),u(i._int()),size(i._int())
	{
		dat=new real[i._int()];
	}
#endif
	/*!
	\param i1 Starting dimension of vector
	\param i2 Ending dimension of vector

	Creation of a variable of type rvector with length \f$ n = i2 - i1 + 1 \f$ and index bounds \f$ lb = i1 \f$, and \f$ ub = i2 \f$. The values of the elements are undefined.
	*/
	INLINE rvector::rvector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_WRONG_BOUNDARIES,ERROR_RVECTOR_NO_MORE_MEMORY):l(i1),u(i2),size(i2-i1+1)
#else
	throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_RVECTOR_WRONG_BOUNDARIES("rvector::rvector(const int &i1,const int &i2)"));
#endif
		dat=new real[size];
	}

	INLINE rvector::rvector(const rvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new real[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE rvector::rvector(const rvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new real[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE rvector::rvector(const real &r) throw():l(1),u(1),size(1)
	{
		dat=new real[1];
		*dat=r;
	}

        INLINE rvector::rvector(const intvector& v) : l(Lb(v)),u(Ub(v)),size(VecLen(v)) {
		dat=new real[size];
		for(int i=0;i<size;i++)
		  dat[i] = v[i+l];
        }

	INLINE real & rvector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC("real & rvector::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}
	
	INLINE real & rvector::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC("real & rvector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}

	INLINE real & rvector_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC("real & rvector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE real & rvector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC("real & rvector_slice::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}
	

	/*!
	\param i The maximum dimension of the wanted part of the vector
	\return The wanted part of the vector

	Example:

	You have the vector \f$ a = \left( \begin {array} {c} x_1 \\ x_2 \\ x_3 \\ x_4 \\ x_5 \end {array} \right) 
	\f$  and then use the operation \f$ b = a(3) \f$ which results in \f$ b = \left( \begin {array} {c} x_1 \\ x_2 \\ x_3 \end {array} \right) \f$
	*/
	INLINE rvector_slice rvector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG("rvector_slice rvector::operator ()(const int &i)"));
#endif
		return rvector_slice(*this,1,i);
	}

	/*!
	\param i1 The starting dimension of the wanted part of the vector
	\param i2 The ending dimension of the wanted part of the vector
	\return The wanted part of the vector

	Example:

	You have the vector \f$ a = \left( \begin {array} {c} x_1 \\ x_2 \\ x_3 \\ x_4 \\ x_5 \end {array} \right) 
	\f$  and then use the operation \f$ b = a(2,4) \f$ which results in \f$ b = \left( \begin {array} {c} x_2 \\ x_3 \\ x_4 \end {array} \right) \f$
	*/
	INLINE rvector_slice rvector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG("rvector_slice rvector::operator ()(const int &i1,const int &i2)"));
#endif
		return rvector_slice(*this,i1,i2);
	}
	
	INLINE rvector_slice rvector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG("rvector_slice rvector_slice::operator ()(const int &i)"));
#endif
		return rvector_slice(*this,1,i);
	}
	
   INLINE rvector_slice rvector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG("rvector_slice rvector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return rvector_slice(*this,i1,i2);
	}

	INLINE real::real(const rvector &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size>1) cxscthrow(ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ("real::real(const rvector &rv)"));
		else if(rv.size<1) cxscthrow(ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ("real::real(const rvector &rv)"));
#endif
		*this=rv.dat[0];
	}
	
	INLINE real::real(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size>1) cxscthrow(ERROR_RVECTOR_TYPE_CAST_OF_THICK_OBJ("real::real(const rvector_slice &sl)"));
		else if(sl.size<1) cxscthrow(ERROR_RVECTOR_USE_OF_UNINITIALIZED_OBJ("real::real(const rvector_slice &sl)"));
#endif
		*this=sl.dat[sl.start-sl.l];
	}

	INLINE rvector &rvector::operator =(const rvector &rv) throw() { return _vvassign<rvector,rvector,real>(*this,rv); }
	INLINE rvector &rvector::operator =(const real &r) throw() { return _vsassign<rvector,real>(*this,r); }
	INLINE rvector::operator void*() throw() { return _vvoid(*this); }

	INLINE rvector_slice & rvector_slice::operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsassign<rvector_slice,rvector_slice>(*this,sl); }
	INLINE rvector_slice & rvector_slice::operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvassign<rvector_slice,rvector>(*this,rv); }
	INLINE rvector_slice & rvector_slice::operator =(const real &r) throw() { return _vssassign<rvector_slice,real>(*this,r); }
	INLINE rvector_slice::operator void*() throw() { return _vsvoid(*this); }

//======================== Vector Functions =============================
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::rvector::rvector(const real &)
	*/
	INLINE rvector _rvector(const real &r) throw() { return rvector(r); }

	INLINE void Resize(rvector &rv) throw() { _vresize(rv); } 
	INLINE void Resize(rvector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<rvector>)
#else
	throw()
#endif
	{ _vresize<class rvector,class real>(rv,len); }
	INLINE void Resize(rvector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<rvector>)
#else
	throw()
#endif
	{ _vresize<class rvector,class real>(rv,lb,ub); }
	
	INLINE rvector abs(const rvector &rv) throw() { return _vabs<rvector,rvector>(rv); }
	INLINE rvector abs(const rvector_slice &sl) throw() { return _vsabs<rvector_slice,rvector>(sl); }
	INLINE bool operator !(const rvector &rv) throw() { return _vnot(rv); }
	INLINE bool operator !(const rvector_slice &sl) throw() { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

	INLINE rvector operator *(const rvector &rv, const real &s) throw() { return _vsmult<rvector,real,rvector>(rv,s); }
	INLINE rvector operator *(const rvector_slice &sl, const real &s) throw() { return _vssmult<rvector_slice,real,rvector>(sl,s); }
	INLINE rvector operator *(const real &s, const rvector &rv) throw() { return _vsmult<rvector,real,rvector>(rv,s); }
	INLINE rvector operator *(const real &s, const rvector_slice &sl) throw() { return _vssmult<rvector_slice,real,rvector>(sl,s); }
	INLINE rvector &operator *=(rvector &rv,const real &r) throw() { return _vsmultassign(rv,r); }
	INLINE rvector_slice &rvector_slice::operator *=(const real &r) throw() { return _vssmultassign(*this,r); }

	INLINE rvector operator /(const rvector &rv, const real &s) throw() { return _vsdiv<rvector,real,rvector>(rv,s); }
	INLINE rvector operator /(const rvector_slice &sl, const real &s) throw() { return _vssdiv<rvector_slice,real,rvector>(sl,s); }
	INLINE rvector &operator /=(rvector &rv,const real &r) throw() { return _vsdivassign(rv,r); }
	INLINE rvector_slice &rvector_slice::operator /=(const real &r) throw() { return _vssdivassign(*this,r); }

//======================= Vector / Vector ===============================

	INLINE rvector &rvector::operator =(const rvector_slice &sl) throw() { return _vvsassign<rvector,rvector_slice,real>(*this,sl); }


		
	INLINE real operator *(const rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvmult<rvector,rvector,real>(rv1,rv2); }
	INLINE real operator *(const rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvmult<rvector_slice,rvector,real>(sl,rv); }
	INLINE real operator *(const rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvmult<rvector_slice,rvector,real>(sl,rv); }
	INLINE real operator *(const rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsmult<rvector_slice,rvector_slice,real>(sl1,sl2); }
	
	INLINE const rvector &operator +(const rvector &rv) throw() { return rv; }
	INLINE rvector operator +(const rvector_slice &sl) throw() { return sl; }
	INLINE rvector operator +(const rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,rvector,rvector>(rv1,rv2); }
	INLINE rvector operator +(const rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,rvector_slice,rvector>(rv,sl); }
	INLINE rvector operator +(const rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,rvector_slice,rvector>(rv,sl); }
	INLINE rvector operator +(const rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,rvector_slice,rvector>(sl1,sl2); }
	INLINE rvector & operator +=(rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE rvector &operator +=(rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE rvector_slice &rvector_slice::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE rvector_slice &rvector_slice::operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE rvector operator -(const rvector &rv) throw() { return _vminus(rv); }
	INLINE rvector operator -(const rvector_slice &sl) throw() { return _vsminus<rvector_slice,rvector>(sl); }
	INLINE rvector operator -(const rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvminus<rvector,rvector,rvector>(rv1,rv2); }
	INLINE rvector operator -(const rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsminus<rvector,rvector_slice,rvector>(rv,sl); }
	INLINE rvector operator -(const rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvminus<rvector_slice,rvector,rvector>(sl,rv); }
	INLINE rvector operator -(const rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsminus<rvector_slice,rvector_slice,rvector>(sl1,sl2); }
	INLINE rvector & operator -=(rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE rvector &operator -=(rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE rvector_slice &rvector_slice::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE rvector_slice &rvector_slice::operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE bool operator ==(const rvector &rv1, const rvector &rv2) throw() { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const rvector_slice &sl1, const rvector_slice &sl2) throw() { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const rvector_slice &sl, const rvector &rv) throw() { return _vsveq(sl,rv); }
	INLINE bool operator ==(const rvector &rv, const rvector_slice &sl) throw() { return _vsveq(sl,rv); }
	INLINE bool operator !=(const rvector &rv1, const rvector &rv2) throw() { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const rvector_slice &sl1, const rvector_slice &sl2) throw() { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const rvector_slice &sl, const rvector &rv) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const rvector &rv, const rvector_slice &sl) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator <(const rvector &rv1, const rvector &rv2) throw() { return _vvless(rv1,rv2); }
	INLINE bool operator <(const rvector_slice &sl1, const rvector_slice &sl2) throw() { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const rvector_slice &sl, const rvector &rv) throw() { return _vsvless(sl,rv); }
	INLINE bool operator < (const rvector &rv, const rvector_slice &sl) throw() { return _vvsless(rv,sl); }
	INLINE bool operator <=(const rvector &rv1, const rvector &rv2) throw() { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const rvector_slice &sl1, const rvector_slice &sl2) throw() { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const rvector_slice &sl, const rvector &rv) throw() { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const rvector &rv, const rvector_slice &sl) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >(const rvector &rv1, const rvector &rv2) throw() { return _vvless(rv2,rv1); }
	INLINE bool operator >(const rvector_slice &sl1, const rvector_slice &sl2) throw() { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const rvector_slice &sl, const rvector &rv) throw() { return _vvsless(rv,sl); }
	INLINE bool operator >(const rvector &rv, const rvector_slice &sl) throw() { return _vsvless(sl,rv); }
	INLINE bool operator >=(const rvector &rv1, const rvector &rv2) throw() { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const rvector_slice &sl1, const rvector_slice &sl2) throw() { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const rvector_slice &sl, const rvector &rv) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const rvector &rv, const rvector_slice &sl) throw() { return _vsvleq(sl,rv); }

	INLINE std::ostream &operator <<(std::ostream &s, const rvector &rv) throw() { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const rvector_slice &sl) throw() { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, rvector &rv) throw() { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, rvector_slice &rv) throw() { return _vsin(s,rv); }

        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE rvector rvector::operator()(const intvector& p) {
          rvector x(*this);
          for(int i=0 ; i<VecLen(x) ; i++)
              x[i+Lb(x)] = (*this)[p[i+Lb(p)]+Lb(*this)];
          return x;
        }




} // namespace cxsc

#endif

