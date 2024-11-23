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

/* CVS $Id: ivector.inl,v 1.27 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_IVECTOR_INL_INCLUDED
#define _CXSC_IVECTOR_INL_INCLUDED

namespace cxsc {

	INLINE ivector::ivector () throw():dat(NULL),l(1),u(0),size(0)
	{
	}

	INLINE ivector::ivector(const int &i) throw():l(1),u(i),size(i)
	{
		dat=new interval[i];
	}

#ifdef OLD_CXSC
	INLINE ivector::ivector(const class index &i) throw():l(1),u(i._int()),size(i._int())
	{
		dat=new interval[i._int()];
	}
#endif

	INLINE ivector::ivector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_WRONG_BOUNDARIES,ERROR_IVECTOR_NO_MORE_MEMORY):l(i1),u(i2),size(i2-i1+1)
#else
	throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_IVECTOR_WRONG_BOUNDARIES("ivector::ivector(const int &i1,const int &i2)"));
#endif
		dat=new interval[size];
	}

	INLINE ivector::ivector(const ivector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new interval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE ivector::ivector(const ivector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new interval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE ivector::ivector(const interval &r) throw():l(1),u(1),size(1)
	{
		dat=new interval[1];
		*dat=r;
	}
	
	INLINE ivector::ivector(const rvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new interval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE ivector::ivector(const rvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new interval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE ivector::ivector(const real &r) throw():l(1),u(1),size(1)
	{
		dat=new interval[1];
		*dat=r;
	}
	
	INLINE interval & ivector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC("interval & ivector::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}

	INLINE interval & ivector::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC("interval & ivector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}

	INLINE interval & ivector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC("interval & ivector_slice::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}

	INLINE interval & ivector_slice::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC("interval & ivector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}

	
	INLINE ivector_slice ivector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG("ivector_slice ivector::operator ()(const int &i)"));
#endif
		return ivector_slice(*this,1,i);
	}
	
   INLINE ivector_slice ivector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG("ivector_slice ivector::operator ()(const int &i1,const int &i2)"));
#endif
		return ivector_slice(*this,i1,i2);
	}
	
	INLINE ivector_slice ivector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG("ivector_slice ivector_slice::operator ()(const int &i)"));
#endif
		return ivector_slice(*this,1,i);
	}
	
   INLINE ivector_slice ivector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG("ivector_slice ivector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return ivector_slice(*this,i1,i2);
	}
	
	INLINE interval::interval(const ivector &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size>1) cxscthrow(ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ("interval::interval(const ivector &rv)"));
		else if(rv.size<1) cxscthrow(ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ("interval::interval(const ivector &rv)"));
#endif
		*this=rv.dat[0];
	}
	
	INLINE interval::interval(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size>1) cxscthrow(ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ("interval::interval(const ivector_slice &sl)"));
		else if(sl.size<1) cxscthrow(ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ("interval::interval(const ivector_slice &sl)"));
#endif
		*this=sl.dat[sl.start-sl.l];
	}

	/*!
	\deprecated use standard contructors for typecasting

	\sa ???
	*/
	INLINE ivector _ivector(const interval &r) throw() { return ivector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa ???
	*/
	INLINE ivector _ivector(const real &r) throw() { return ivector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa ???
	*/
	INLINE ivector _ivector(const rvector_slice &rs) throw() { return ivector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa ???
	*/
	INLINE ivector _ivector(const rvector &rs) throw() { return ivector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa ???
	*/
	INLINE ivector _ivector(const rmatrix_subv &rs) throw() { return ivector(rs); }
	INLINE ivector &ivector::operator =(const ivector &rv) throw() { return _vvassign<ivector,ivector,interval>(*this,rv); }
	INLINE ivector &ivector::operator =(const interval &r) throw() { return _vsassign<ivector,interval>(*this,r); }
	INLINE ivector &ivector::operator =(const rvector &rv) throw() { return _vvassign<ivector,rvector,interval>(*this,rv); }
	INLINE ivector &ivector::operator =(const real &r) throw() { return _vsassign<ivector,real>(*this,r); }
	INLINE ivector::operator void*() throw() { return _vvoid(*this); }
	INLINE ivector_slice & ivector_slice::operator =(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE ivector_slice & ivector_slice::operator =(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE ivector_slice & ivector_slice::operator =(const interval &r) throw() { return _vssassign<ivector_slice,interval>(*this,r); }
	INLINE ivector_slice & ivector_slice::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,ivector(m)); }
	INLINE ivector_slice & ivector_slice::operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE ivector_slice & ivector_slice::operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE ivector_slice & ivector_slice::operator =(const real &r) throw() { return _vssassign(*this,r); }
	INLINE ivector_slice::operator void*() throw() { return _vsvoid(*this); }

//=======================================================================
//======================== Vector Functions =============================


	INLINE ivector &SetInf(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetinf(iv,rv); }
	INLINE ivector_slice &SetInf(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetinf(iv,rv); }
	INLINE ivector &SetInf(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetinf(iv,rv); }
	INLINE ivector_slice &SetInf(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetinf(iv,rv); }
	INLINE ivector &UncheckedSetInf(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvusetinf(iv,rv); }
	INLINE ivector_slice &UncheckedSetInf(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvusetinf(iv,rv); }
	INLINE ivector &UncheckedSetInf(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsusetinf(iv,rv); }
	INLINE ivector_slice &UncheckedSetInf(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsusetinf(iv,rv); }

	INLINE ivector &SetSup(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetsup(iv,rv); }
	INLINE ivector_slice &SetSup(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetsup(iv,rv); }
	INLINE ivector &SetSup(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetsup(iv,rv); }
	INLINE ivector_slice &SetSup(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetsup(iv,rv); }
	INLINE ivector &UncheckedSetSup(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvusetsup(iv,rv); }
	INLINE ivector_slice &UncheckedSetSup(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvusetsup(iv,rv); }
	INLINE ivector &UncheckedSetSup(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsusetsup(iv,rv); }
	INLINE ivector_slice &UncheckedSetSup(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsusetsup(iv,rv); }

	INLINE ivector &SetSup(ivector &iv,const real &r) throw() { return _vssetsup(iv,r); }
	INLINE ivector &SetInf(ivector &iv,const real &r) throw() { return _vssetinf(iv,r); }
	INLINE ivector &UncheckedSetSup(ivector &iv,const real &r) throw() { return _vsusetsup(iv,r); }
	INLINE ivector &SetUncheckedInf(ivector &iv,const real &r) throw() { return _vsusetinf(iv,r); }
	INLINE ivector_slice &SetSup(ivector_slice &iv,const real &r) throw() { return _vsssetsup(iv,r); }
	INLINE ivector_slice &SetInf(ivector_slice &iv,const real &r) throw() { return _vsssetinf(iv,r); }
	INLINE ivector_slice &UncheckedSetSup(ivector_slice &iv,const real &r) throw() { return _vssusetsup(iv,r); }
	INLINE ivector_slice &SetUncheckedInf(ivector_slice &iv,const real &r) throw() { return _vssusetinf(iv,r); }

	INLINE void Resize(ivector &rv) throw() { _vresize(rv); } 
	INLINE void Resize(ivector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<ivector>)
#else
	throw()
#endif
	{ _vresize<class ivector,class interval>(rv,len); }
	INLINE void Resize(ivector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<ivector>)
#else
	throw()
#endif
	{ _vresize<class ivector,class interval>(rv,lb,ub); }
	
	INLINE ivector abs(const ivector &rv) throw() { return _vabs<ivector,ivector>(rv); }
	INLINE ivector abs(const ivector_slice &sl) throw() { return _vsabs<ivector_slice,ivector>(sl); }
	INLINE rvector absmin(const ivector &rv) throw() {
          rvector x(Lb(rv),Ub(rv));
          for(int i=Lb(rv) ; i<=Ub(rv) ; i++)
            x[i] = AbsMin(rv[i]);
          return x;
        }
	INLINE rvector absmin(const ivector_slice &sl) throw() { 
          rvector x(Lb(sl),Ub(sl));
          for(int i=Lb(sl) ; i<=Ub(sl) ; i++)
            x[i] = AbsMin(sl[i]);
          return x;
        }
	INLINE rvector absmax(const ivector &rv) throw() {
          rvector x(Lb(rv),Ub(rv));
          for(int i=Lb(rv) ; i<=Ub(rv) ; i++)
            x[i] = AbsMax(rv[i]);
          return x;
        }
	INLINE rvector absmax(const ivector_slice &sl) throw() { 
          rvector x(Lb(sl),Ub(sl));
          for(int i=Lb(sl) ; i<=Ub(sl) ; i++)
            x[i] = AbsMax(sl[i]);
          return x;
        }
	INLINE rvector diam(const ivector &v) throw() { return _vdiam<ivector,rvector>(v); }
	INLINE rvector diam(const ivector_slice &v) throw() { return _vsdiam<ivector_slice,rvector>(v); }
	INLINE rvector mid(const ivector &v) throw() { return _vmid<ivector,rvector>(v); }
	INLINE rvector mid(const ivector_slice &v) throw() { return _vsmid<ivector_slice,rvector>(v); }
	INLINE rvector Inf(const ivector &v) throw() { return _vinf<ivector,rvector>(v); }
	INLINE rvector Inf(const ivector_slice &v) throw() { return _vsinf<ivector_slice,rvector>(v); }
	INLINE rvector Sup(const ivector &v) throw() { return _vsup<ivector,rvector>(v); }
	INLINE rvector Sup(const ivector_slice &v) throw() { return _vssup<ivector_slice,rvector>(v); }
	INLINE bool operator !(const ivector &rv) throw() { return _vnot(rv); }
	INLINE bool operator !(const ivector_slice &sl) throw() { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

//----------------------------- Interval ---------------------------

	INLINE ivector operator *(const ivector &rv, const interval &s) throw() { return _vsmult<ivector,interval,ivector>(rv,s); }
	INLINE ivector operator *(const ivector_slice &sl, const interval &s) throw() { return _vssmult<ivector_slice,interval,ivector>(sl,s); }
	INLINE ivector operator *(const interval &s, const ivector &rv) throw() { return _vsmult<ivector,interval,ivector>(rv,s); }
	INLINE ivector operator *(const interval &s, const ivector_slice &sl) throw() { return _vssmult<ivector_slice,interval,ivector>(sl,s); }
	INLINE ivector &operator *=(ivector &rv,const interval &r) throw() { return _vsmultassign(rv,r); }
	INLINE ivector_slice &ivector_slice::operator *=(const interval &r) throw() { return _vssmultassign(*this,r); }

	INLINE ivector operator /(const ivector &rv, const interval &s) throw() { return _vsdiv<ivector,interval,ivector>(rv,s); }
	INLINE ivector operator /(const ivector_slice &sl, const interval &s) throw() { return _vssdiv<ivector_slice,interval,ivector>(sl,s); }
	INLINE ivector &operator /=(ivector &rv,const interval &r) throw() { return _vsdivassign(rv,r); }
	INLINE ivector_slice &ivector_slice::operator /=(const interval &r) throw() { return _vssdivassign(*this,r); }

//---------------------------- Real --------------------------------------

	INLINE ivector operator *(const ivector &rv, const real &s) throw() { return _vsmult<ivector,real,ivector>(rv,s); }
	INLINE ivector operator *(const ivector_slice &sl, const real &s) throw() { return _vssmult<ivector_slice,real,ivector>(sl,s); }
	INLINE ivector operator *(const real &s, const ivector &rv) throw() { return _vsmult<ivector,real,ivector>(rv,s); }
	INLINE ivector operator *(const real &s, const ivector_slice &sl) throw() { return _vssmult<ivector_slice,real,ivector>(sl,s); }
	INLINE ivector &operator *=(ivector &rv,const real &r) throw() { return _vsmultassign(rv,r); }
	INLINE ivector_slice &ivector_slice::operator *=(const real &r) throw() { return _vssmultassign(*this,r); }

	INLINE ivector operator /(const ivector &rv, const real &s) throw() { return _vsdiv<ivector,real,ivector>(rv,s); }
	INLINE ivector operator /(const ivector_slice &sl, const real &s) throw() { return _vssdiv<ivector_slice,real,ivector>(sl,s); }
	INLINE ivector &operator /=(ivector &rv,const real &r) throw() { return _vsdivassign(rv,r); }
	INLINE ivector_slice &ivector_slice::operator /=(const real &r) throw() { return _vssdivassign(*this,r); }

	INLINE ivector operator *(const rvector &rv, const interval &s) throw() { return _vsmult<rvector,interval,ivector>(rv,s); }
	INLINE ivector operator *(const rvector_slice &sl, const interval &s) throw() { return _vssmult<rvector_slice,interval,ivector>(sl,s); }
	INLINE ivector operator *(const interval &s, const rvector &rv) throw() { return _vsmult<rvector,interval,ivector>(rv,s); }
	INLINE ivector operator *(const interval &s, const rvector_slice &sl) throw() { return _vssmult<rvector_slice,interval,ivector>(sl,s); }

	INLINE ivector operator /(const rvector &rv, const interval &s) throw() { return _vsdiv<rvector,interval,ivector>(rv,s); }
	INLINE ivector operator /(const rvector_slice &sl, const interval &s) throw() { return _vssdiv<rvector_slice,interval,ivector>(sl,s); }

//======================= Vector / Vector ===============================


	INLINE std::ostream &operator <<(std::ostream &s, const ivector &rv) throw() { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const ivector_slice &sl) throw() { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, ivector &rv) throw() { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, ivector_slice &rv) throw() { return _vsin(s,rv); }
	
//----------------------- Interval / Interval ---------------------------
	INLINE ivector & ivector::operator =(const ivector_slice &sl) throw() { return _vvsassign<ivector,ivector_slice,interval>(*this,sl); }


	INLINE interval operator *(const ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvimult<ivector,ivector,interval>(rv1,rv2); }
	INLINE interval operator *(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvimult<ivector_slice,ivector,interval>(sl,rv); }
	INLINE interval operator *(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvimult<ivector_slice,ivector,interval>(sl,rv); }
	INLINE interval operator *(const ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsimult<ivector_slice,ivector_slice,interval>(sl1,sl2); }
	
	INLINE const ivector &operator +(const ivector &rv) throw() { return rv; }
	INLINE ivector operator +(const ivector_slice &sl) throw() { return sl; }
	INLINE ivector operator +(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvplus<ivector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator +(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator +(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator +(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsplus<ivector_slice,ivector_slice,ivector>(sl1,sl2); }
	INLINE ivector & operator +=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE ivector &operator +=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator +=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE ivector operator -(const ivector &rv) throw() { return _vminus(rv); }
	INLINE ivector operator -(const ivector_slice &sl) throw() { return _vsminus<ivector_slice,ivector>(sl); }
	INLINE ivector operator -(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvminus<ivector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator -(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsminus<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator -(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvminus<ivector_slice,ivector,ivector>(sl,rv); }
	INLINE ivector operator -(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsminus<ivector_slice,ivector_slice,ivector>(sl1,sl2); }
	INLINE ivector & operator -=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE ivector &operator -=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator -=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE ivector operator |(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvconv<ivector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator |(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsconv<ivector_slice,ivector_slice,ivector>(sl1,sl2); }
	INLINE ivector & operator |=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE ivector &operator |=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator |=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE ivector operator &(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsect<ivector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator &(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator &(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator &(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvssect<ivector_slice,ivector_slice,ivector>(sl1,sl2); }
	INLINE ivector & operator &=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE ivector &operator &=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator &=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvssectassign(*this,sl2); }

	INLINE bool operator ==(const ivector &rv1, const ivector &rv2) throw() { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const ivector_slice &sl1, const ivector_slice &sl2) throw() { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const ivector_slice &sl, const ivector &rv) throw() { return _vsveq(sl,rv); }
	INLINE bool operator ==(const ivector &rv, const ivector_slice &sl) throw() { return _vsveq(sl,rv); }
	INLINE bool operator !=(const ivector &rv1, const ivector &rv2) throw() { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const ivector_slice &sl1, const ivector_slice &sl2) throw() { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const ivector_slice &sl, const ivector &rv) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const ivector &rv, const ivector_slice &sl) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator <(const ivector &rv1, const ivector &rv2) throw() { return _vvless(rv1,rv2); }
	INLINE bool operator <(const ivector_slice &sl1, const ivector_slice &sl2) throw() { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const ivector_slice &sl, const ivector &rv) throw() { return _vsvless(sl,rv); }
	INLINE bool operator < (const ivector &rv, const ivector_slice &sl) throw() { return _vvsless(rv,sl); }
	INLINE bool operator <=(const ivector &rv1, const ivector &rv2) throw() { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const ivector_slice &sl1, const ivector_slice &sl2) throw() { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const ivector_slice &sl, const ivector &rv) throw() { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const ivector &rv, const ivector_slice &sl) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >(const ivector &rv1, const ivector &rv2) throw() { return _vvless(rv2,rv1); }
	INLINE bool operator >(const ivector_slice &sl1, const ivector_slice &sl2) throw() { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const ivector_slice &sl, const ivector &rv) throw() { return _vvsless(rv,sl); }
	INLINE bool operator >(const ivector &rv, const ivector_slice &sl) throw() { return _vsvless(sl,rv); }
	INLINE bool operator >=(const ivector &rv1, const ivector &rv2) throw() { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const ivector_slice &sl1, const ivector_slice &sl2) throw() { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const ivector_slice &sl, const ivector &rv) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const ivector &rv, const ivector_slice &sl) throw() { return _vsvleq(sl,rv); }

//-------------------------------- Interval / Real --------------------------------

	INLINE ivector & ivector::operator =(const rvector_slice &sl) throw() { return _vvsassign<ivector,rvector_slice,interval>(*this,sl); }


	INLINE interval operator *(const rvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvimult<rvector,ivector,interval>(rv1,rv2); }
	INLINE interval operator *(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvimult<rvector_slice,ivector,interval>(sl,rv); }
	INLINE interval operator *(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvimult<ivector_slice,rvector,interval>(sl,rv); }
	INLINE interval operator *(const rvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsimult<rvector_slice,ivector_slice,interval>(sl1,sl2); }
	
	INLINE interval operator *(const ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvimult<rvector,ivector,interval>(rv2,rv1); }
	INLINE interval operator *(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvimult<ivector_slice,rvector,interval>(sl,rv); }
	INLINE interval operator *(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvimult<rvector_slice,ivector,interval>(sl,rv); }
	INLINE interval operator *(const ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsimult<rvector_slice,ivector_slice,interval>(sl2,sl1); }
	
	INLINE ivector operator +(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator +(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator +(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator +(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,ivector_slice,ivector>(sl1,sl2); }

	INLINE ivector operator +(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,ivector,ivector>(rv2,rv1); }
	INLINE ivector operator +(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator +(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator +(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,ivector_slice,ivector>(sl2,sl1); }

	INLINE ivector & operator +=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE ivector &operator +=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE ivector operator -(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvminus<rvector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator -(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsminus<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator -(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvminus<rvector_slice,ivector,ivector>(sl,rv); }
	INLINE ivector operator -(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsminus<rvector_slice,ivector_slice,ivector>(sl1,sl2); }

	INLINE ivector operator -(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvminus<ivector,rvector,ivector>(rv1,rv2); }
	INLINE ivector operator -(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsminus<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator -(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvminus<ivector_slice,rvector,ivector>(sl,rv); }
	INLINE ivector operator -(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsminus<ivector_slice,rvector_slice,ivector>(sl1,sl2); }

	INLINE ivector & operator -=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE ivector &operator -=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE ivector operator |(const rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,rvector,ivector>(rv1,rv2); }
	INLINE ivector operator |(const rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,rvector_slice,ivector>(sl1,sl2); }
	INLINE ivector operator |(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator |(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,ivector_slice,ivector>(sl1,sl2); }

	INLINE ivector operator |(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,ivector,ivector>(rv2,rv1); }
	INLINE ivector operator |(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator |(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,ivector_slice,ivector>(sl2,sl1); }

	INLINE ivector & operator |=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE ivector &operator |=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE ivector operator &(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsect<rvector,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator &(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssect<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator &(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator &(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvssect<rvector_slice,ivector_slice,ivector>(sl1,sl2); }

	INLINE ivector operator &(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsect<rvector,ivector,ivector>(rv2,rv1); }
	INLINE ivector operator &(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,rvector_slice,ivector>(rv,sl); }
	INLINE ivector operator &(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssect<rvector,ivector_slice,ivector>(rv,sl); }
	INLINE ivector operator &(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvssect<rvector_slice,ivector_slice,ivector>(sl2,sl1); }

	INLINE ivector & operator &=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE ivector &operator &=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE ivector_slice &ivector_slice::operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE ivector_slice &ivector_slice::operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>)
#else
	throw()
#endif
	{ return _vsvssectassign(*this,sl2); }


        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE ivector ivector::operator()(const intvector& p) {
          ivector x(*this);
          for(int i=0 ; i<VecLen(x) ; i++)
              x[i+Lb(x)] = (*this)[p[i+Lb(p)]+Lb(*this)];
          return x;
        }

} // namespace cxsc

#endif // _CXSC_IVECTOR_INL_INCLUDED
