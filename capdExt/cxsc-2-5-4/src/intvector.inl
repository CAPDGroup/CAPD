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

/* CVS $Id: intvector.inl,v 1.19 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_INTVECTOR_INL_INCLUDED
#define _CXSC_INTVECTOR_INL_INCLUDED

#include "intvector.hpp"

namespace cxsc {

	INLINE intvector::intvector () throw():dat(NULL),l(1),u(0),size(0)
	{
	}

	INLINE intvector::intvector(const int &i) throw():l(1),u(i),size(i)
	{
		dat=new int[i];
	}

#ifdef OLD_CXSC  
	INLINE intvector::intvector(const class index &i) throw():l(1),u(i._int()),size(i._int())
	{
		dat=new int[i._int()];
	}
#endif

	INLINE intvector::intvector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_WRONG_BOUNDARIES,ERROR_INTVECTOR_NO_MORE_MEMORY):l(i1),u(i2),size(i2-i1+1)
#else
	throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_INTVECTOR_WRONG_BOUNDARIES("intvector::intvector(const int &i1,const int &i2)"));
#endif
		dat=new int[size];
	}

	INLINE intvector::intvector(const intvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new int[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE intvector::intvector(const intvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new int[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE int & intvector_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC("int & intvector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE const int & intvector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC("int & intvector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE const int & intvector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC("int & intvector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE int & intvector::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC("int & intvector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	/*!
	\param i The maximum dimension of the wanted part of the vector
	\return The wanted part of the vector

	\sa rvector::operator ()(const int &i)
	*/
	INLINE intvector_slice intvector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG("intvector_slice intvector::operator ()(const int &i)"));
#endif
		return intvector_slice(*this,1,i);
	}
	
	/*!
	\param i1 The starting dimension of the wanted part of the vector
	\param i2 The ending dimension of the wanted part of the vector
	\return The wanted part of the vector

	\sa rvector::operator ()(const int &i1,const int &i2)
	*/
	INLINE intvector_slice intvector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG("intvector_slice intvector::operator ()(const int &i1,const int &i2)"));
#endif
		return intvector_slice(*this,i1,i2);
	}
	
	INLINE intvector_slice intvector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG("intvector_slice intvector_slice::operator ()(const int &i)"));
#endif
		return intvector_slice(*this,1,i);
	}
	
   INLINE intvector_slice intvector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG("intvector_slice intvector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return intvector_slice(*this,i1,i2);
	}

	INLINE intvector &intvector::operator =(const intvector &rv) throw() { return _vvassign<intvector,intvector,int>(*this,rv); }
	INLINE intvector &intvector::operator =(const int &r) throw() { return _vsassign<intvector,int>(*this,r); }
	INLINE intvector::operator void*() throw() { return _vvoid(*this); }

	INLINE intvector_slice & intvector_slice::operator =(const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvsassign<intvector_slice,intvector_slice>(*this,sl); }
	INLINE intvector_slice & intvector_slice::operator =(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvassign<intvector_slice,intvector>(*this,rv); }
	INLINE intvector_slice & intvector_slice::operator =(const int &r) throw() { return _vssassign<intvector_slice,int>(*this,r); }
	INLINE intvector_slice::operator void*() throw() { return _vsvoid(*this); }

//======================== Vector Functions =============================
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::intvector::intvector(const int &)
	*/
	INLINE intvector _intvector(const int &r) throw() { return intvector(r); }

	INLINE void Resize(intvector &rv) throw() { _vresize(rv); } 
	INLINE void Resize(intvector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<intvector>)
#else
	throw()
#endif
	{ _vresize<intvector,int>(rv,len); }
	INLINE void Resize(intvector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<intvector>)
#else
	throw()
#endif
	{ _vresize<intvector,int>(rv,lb,ub); }
	
	INLINE intvector abs(const intvector &rv) throw() { return _vabs<intvector,intvector>(rv); }
	INLINE intvector abs(const intvector_slice &sl) throw() { return _vsabs<intvector_slice,intvector>(sl); }
	INLINE bool operator !(const intvector &rv) throw() { return _vnot(rv); }
	INLINE bool operator !(const intvector_slice &sl) throw() { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

	INLINE intvector operator *(const intvector &rv, const int &s) throw() { return _vsmult<intvector,int,intvector>(rv,s); }
	INLINE intvector operator *(const intvector_slice &sl, const int &s) throw() { return _vssmult<intvector_slice,int,intvector>(sl,s); }
	INLINE intvector operator *(const int &s, const intvector &rv) throw() { return _vsmult<intvector,int,intvector>(rv,s); }
	INLINE intvector operator *(const int &s, const intvector_slice &sl) throw() { return _vssmult<intvector_slice,int,intvector>(sl,s); }
	INLINE intvector &operator *=(intvector &rv,const int &r) throw() { return _vsmultassign(rv,r); }
	INLINE intvector_slice &intvector_slice::operator *=(const int &r) throw() { return _vssmultassign(*this,r); }

	INLINE intvector operator /(const intvector &rv, const int &s) throw() { return _vsdiv<intvector,int,intvector>(rv,s); }
	INLINE intvector operator /(const intvector_slice &sl, const int &s) throw() { return _vssdiv<intvector_slice,int,intvector>(sl,s); }
	INLINE intvector &operator /=(intvector &rv,const int &r) throw() { return _vsdivassign(rv,r); }
	INLINE intvector_slice &intvector_slice::operator /=(const int &r) throw() { return _vssdivassign(*this,r); }

//======================= Vector / Vector ===============================

	INLINE intvector &intvector::operator =(const intvector_slice &sl) throw() { return _vvsassign<intvector,intvector_slice,int>(*this,sl); }


	INLINE void accumulate(dotprecision &dp, const intvector & rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv1,rv2); }
//	INLINE void accumulate(dotprecision &dp, const intvector & rv1, const intmatrix_subv &rv2) throw(OP_WITH_WRONG_DIM);
//	INLINE void accumulate(dotprecision &dp, const intmatrix_subv & rv1, const intvector &rv2) throw(OP_WITH_WRONG_DIM);
	INLINE void accumulate(dotprecision &dp,const intvector_slice &sl,const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp,const intvector &rv,const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp, const intvector_slice & sl1, const intvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE const intvector &operator +(const intvector &rv) throw() { return rv; }
	INLINE intvector operator +(const intvector_slice &sl) throw() { return sl; }
	INLINE intvector operator +(const intvector &rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvplus<intvector,intvector,intvector>(rv1,rv2); }
	INLINE intvector operator +(const intvector &rv, const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvsplus<intvector,intvector_slice,intvector>(rv,sl); }
	INLINE intvector operator +(const intvector_slice &sl, const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvsplus<intvector,intvector_slice,intvector>(rv,sl); }
	INLINE intvector operator +(const intvector_slice &sl1, const intvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvsplus<intvector_slice,intvector_slice,intvector>(sl1,sl2); }
	INLINE intvector & operator +=(intvector &rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE intvector &operator +=(intvector &rv, const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE intvector_slice &intvector_slice::operator +=(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE intvector_slice &intvector_slice::operator +=(const intvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE intvector operator -(const intvector &rv) throw() { return _vminus(rv); }
	INLINE intvector operator -(const intvector_slice &sl) throw() { return _vsminus<intvector_slice,intvector>(sl); }
	INLINE intvector operator -(const intvector &rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvminus<intvector,intvector,intvector>(rv1,rv2); }
	INLINE intvector operator -(const intvector &rv, const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvsminus<intvector,intvector_slice,intvector>(rv,sl); }
	INLINE intvector operator -(const intvector_slice &sl, const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvminus<intvector_slice,intvector,intvector>(sl,rv); }
	INLINE intvector operator -(const intvector_slice &sl1, const intvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvsminus<intvector_slice,intvector_slice,intvector>(sl1,sl2); }
	INLINE intvector & operator -=(intvector &rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE intvector &operator -=(intvector &rv, const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE intvector_slice &intvector_slice::operator -=(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE intvector_slice &intvector_slice::operator -=(const intvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE bool operator ==(const intvector &rv1, const intvector &rv2) throw() { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const intvector_slice &sl1, const intvector_slice &sl2) throw() { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const intvector_slice &sl, const intvector &rv) throw() { return _vsveq(sl,rv); }
	INLINE bool operator ==(const intvector &rv, const intvector_slice &sl) throw() { return _vsveq(sl,rv); }
	INLINE bool operator !=(const intvector &rv1, const intvector &rv2) throw() { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const intvector_slice &sl1, const intvector_slice &sl2) throw() { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const intvector_slice &sl, const intvector &rv) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const intvector &rv, const intvector_slice &sl) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator <(const intvector &rv1, const intvector &rv2) throw() { return _vvless(rv1,rv2); }
	INLINE bool operator <(const intvector_slice &sl1, const intvector_slice &sl2) throw() { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const intvector_slice &sl, const intvector &rv) throw() { return _vsvless(sl,rv); }
	INLINE bool operator < (const intvector &rv, const intvector_slice &sl) throw() { return _vvsless(rv,sl); }
	INLINE bool operator <=(const intvector &rv1, const intvector &rv2) throw() { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const intvector_slice &sl1, const intvector_slice &sl2) throw() { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const intvector_slice &sl, const intvector &rv) throw() { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const intvector &rv, const intvector_slice &sl) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >(const intvector &rv1, const intvector &rv2) throw() { return _vvless(rv2,rv1); }
	INLINE bool operator >(const intvector_slice &sl1, const intvector_slice &sl2) throw() { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const intvector_slice &sl, const intvector &rv) throw() { return _vvsless(rv,sl); }
	INLINE bool operator >(const intvector &rv, const intvector_slice &sl) throw() { return _vsvless(sl,rv); }
	INLINE bool operator >=(const intvector &rv1, const intvector &rv2) throw() { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const intvector_slice &sl1, const intvector_slice &sl2) throw() { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const intvector_slice &sl, const intvector &rv) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const intvector &rv, const intvector_slice &sl) throw() { return _vsvleq(sl,rv); }

	INLINE std::ostream &operator <<(std::ostream &s, const intvector &rv) throw() { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const intvector_slice &sl) throw() { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, intvector &rv) throw() { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, intvector_slice &rv) throw() { return _vsin(s,rv); }

        INLINE intvector perminv(const intvector& x) {
          intvector p(0,VecLen(x));
          for(int i=0 ; i<VecLen(x) ; i++)
            p[x[i+Lb(x)]] = i;
          return p;
        }


} // namespace cxsc

#endif

