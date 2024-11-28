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

/* CVS $Id: l_ivector.inl,v 1.20 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_LIVECTOR_INL_INCLUDED
#define _CXSC_LIVECTOR_INL_INCLUDED

namespace cxsc {

	INLINE l_ivector::l_ivector ():dat(NULL),l(1),u(0),size(0)
	{
	}

	INLINE l_ivector::l_ivector(const int &i):l(1),u(i),size(i)
	{
		dat=new l_interval[i];
	}

#ifdef OLD_CXSC
	INLINE l_ivector::l_ivector(const class index &i):l(1),u(i._int()),size(i._int())
	{
		dat=new l_interval[i._int()];
	}
#endif

	INLINE l_ivector::l_ivector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		:l(i1),u(i2),size(i2-i1+1)
#else
	:l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_LIVECTOR_WRONG_BOUNDARIES("l_ivector::l_ivector(const int &i1,const int &i2)"));
#endif
		dat=new l_interval[size];
	}

	INLINE l_ivector::l_ivector(const l_ivector_slice &rs):l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new l_interval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE l_ivector::l_ivector(const l_rvector_slice &rs):l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new l_interval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE l_ivector::l_ivector(const ivector_slice &rs):l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new l_interval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE l_ivector::l_ivector(const rvector_slice &rs):l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new l_interval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE l_ivector::l_ivector(const l_ivector &v):l(v.l),u(v.u),size(v.size)
	{
		dat=new l_interval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE l_ivector::l_ivector(const l_interval &r):l(1),u(1),size(1)
	{
		dat=new l_interval[1];
		*dat=r;
	}
	
	INLINE l_ivector::l_ivector(const l_rvector &v):l(v.l),u(v.u),size(v.size)
	{
		dat=new l_interval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE l_ivector::l_ivector(const ivector &v):l(v.l),u(v.u),size(v.size)
	{
		dat=new l_interval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE l_ivector::l_ivector(const rvector &v):l(v.l),u(v.u),size(v.size)
	{
		dat=new l_interval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE l_ivector::l_ivector(const real &r):l(1),u(1),size(1)
	{
		dat=new l_interval[1];
		*dat=r;
	}
	
	INLINE l_ivector::l_ivector(const interval &r):l(1),u(1),size(1)
	{
		dat=new l_interval[1];
		*dat=r;
	}
	
	INLINE l_ivector::l_ivector(const l_real &r):l(1),u(1),size(1)
	{
		dat=new l_interval[1];
		*dat=r;
	}
	

	INLINE l_interval & l_ivector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC("l_interval & l_ivector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE l_interval & l_ivector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC("l_interval & l_ivector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE l_ivector_slice l_ivector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG("l_ivector_slice l_ivector::operator ()(const int &i)"));
#endif
		return l_ivector_slice(*this,1,i);
	}
	
   INLINE l_ivector_slice l_ivector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG("l_ivector_slice l_ivector::operator ()(const int &i1,const int &i2)"));
#endif
		return l_ivector_slice(*this,i1,i2);
	}
	
	INLINE l_ivector_slice l_ivector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG("l_ivector_slice l_ivector_slice::operator ()(const int &i)"));
#endif
		return l_ivector_slice(*this,1,i);
	}
	
   INLINE l_ivector_slice l_ivector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG("l_ivector_slice l_ivector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return l_ivector_slice(*this,i1,i2);
	}
	
	INLINE l_interval::l_interval(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size>1) cxscthrow(ERROR_LIVECTOR_TYPE_CAST_OF_THICK_OBJ("l_interval::l_interval(const l_ivector &rv)"));
		else if(rv.size<1) cxscthrow(ERROR_LIVECTOR_USE_OF_UNINITIALIZED_OBJ("l_interval::l_interval(const l_ivector &rv)"));
#endif
		*this=rv.dat[0];
	}
	
	INLINE l_interval::l_interval(const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
		
#else
	
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size>1) cxscthrow(ERROR_LIVECTOR_TYPE_CAST_OF_THICK_OBJ("l_interval::l_interval(const l_ivector_slice &sl)"));
		else if(sl.size<1) cxscthrow(ERROR_LIVECTOR_USE_OF_UNINITIALIZED_OBJ("l_interval::l_interval(const l_ivector_slice &sl)"));
#endif
		*this=sl.dat[sl.start-sl.l];
	}
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const l_interval &)
	*/
	INLINE l_ivector _l_ivector(const l_interval &r) { return l_ivector(r); }

	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const real &)
	*/
	INLINE l_ivector _l_ivector(const real &r) { return l_ivector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const rvector_slice &rs)
	*/
	INLINE l_ivector _l_ivector(const rvector_slice &rs) { return l_ivector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const rvector &v)
	*/
	INLINE l_ivector _l_ivector(const rvector &rs) { return l_ivector(rs); }
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const l_real &)
	*/
	INLINE l_ivector _l_ivector(const l_real &r) { return l_ivector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const l_rvector_slice &rs)
	*/
	INLINE l_ivector _l_ivector(const l_rvector_slice &rs) { return l_ivector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const l_rvector &v)
	*/
	INLINE l_ivector _l_ivector(const l_rvector &rs) { return l_ivector(rs); }
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const interval &)
	*/
	INLINE l_ivector _l_ivector(const interval &r) { return l_ivector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const ivector_slice &rs)
	*/
	INLINE l_ivector _l_ivector(const ivector_slice &rs) { return l_ivector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const ivector &v)
	*/
	INLINE l_ivector _l_ivector(const ivector &rs) { return l_ivector(rs); }
	
	INLINE l_ivector &l_ivector::operator =(const l_ivector &rv) { return _vvassign<l_ivector,l_ivector,l_interval>(*this,rv); }
	INLINE l_ivector &l_ivector::operator =(const l_interval &r) { return _vsassign<l_ivector,l_interval>(*this,r); }
	INLINE l_ivector::operator void*() { return _vvoid(*this); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvassign(*this,rv); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const l_interval &r) { return _vssassign<l_ivector_slice,l_interval>(*this,r); }
	INLINE l_ivector_slice::operator void*() { return _vsvoid(*this); }
	
//=======================================================================
//======================== Vector Functions =============================


	INLINE l_ivector &SetInf(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsetinf(iv,rv); }
	INLINE l_ivector_slice &SetInf(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsetinf(iv,rv); }
	INLINE l_ivector &SetInf(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssetinf(iv,rv); }
	INLINE l_ivector_slice &SetInf(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssetinf(iv,rv); }
	INLINE l_ivector &UncheckedSetInf(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvusetinf(iv,rv); }
	INLINE l_ivector_slice &UncheckedSetInf(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvusetinf(iv,rv); }
	INLINE l_ivector &UncheckedSetInf(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsusetinf(iv,rv); }
	INLINE l_ivector_slice &UncheckedSetInf(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsusetinf(iv,rv); }

	INLINE l_ivector &SetSup(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsetsup(iv,rv); }
	INLINE l_ivector_slice &SetSup(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsetsup(iv,rv); }
	INLINE l_ivector &SetSup(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssetsup(iv,rv); }
	INLINE l_ivector_slice &SetSup(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssetsup(iv,rv); }
	INLINE l_ivector &UncheckedSetSup(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvusetsup(iv,rv); }
	INLINE l_ivector_slice &UncheckedSetSup(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvusetsup(iv,rv); }
	INLINE l_ivector &UncheckedSetSup(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsusetsup(iv,rv); }
	INLINE l_ivector_slice &UncheckedSetSup(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsusetsup(iv,rv); }

	INLINE l_ivector &SetSup(l_ivector &iv,const l_real &r) { return _vssetsup(iv,r); }
	INLINE l_ivector &SetInf(l_ivector &iv,const l_real &r) { return _vssetinf(iv,r); }
	INLINE l_ivector &UncheckedSetSup(l_ivector &iv,const l_real &r) { return _vsusetsup(iv,r); }
	INLINE l_ivector &SetUncheckedInf(l_ivector &iv,const l_real &r) { return _vsusetinf(iv,r); }

	INLINE l_ivector_slice &SetSup(l_ivector_slice &iv,const l_real &r) { return _vsssetsup(iv,r); }
	INLINE l_ivector_slice &SetInf(l_ivector_slice &iv,const l_real &r) { return _vsssetinf(iv,r); }
	INLINE l_ivector_slice &UncheckedSetSup(l_ivector_slice &iv,const l_real &r) { return _vssusetsup(iv,r); }
	INLINE l_ivector_slice &SetUncheckedInf(l_ivector_slice &iv,const l_real &r) { return _vssusetinf(iv,r); }

	INLINE void Resize(l_ivector &rv) { _vresize(rv); } 
	INLINE void Resize(l_ivector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vresize<class l_ivector,class l_interval>(rv,len); }
	INLINE void Resize(l_ivector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vresize<class l_ivector,class l_interval>(rv,lb,ub); }
	
	INLINE l_ivector abs(const l_ivector &rv) { return _vabs<l_ivector,l_ivector>(rv); }
	INLINE l_ivector abs(const l_ivector_slice &sl) { return _vsabs<l_ivector_slice,l_ivector>(sl); }
	INLINE l_rvector diam(const l_ivector &v) { return _vdiam<l_ivector,l_rvector>(v); }
	INLINE l_rvector diam(const l_ivector_slice &v) { return _vsdiam<l_ivector_slice,l_rvector>(v); }
	INLINE l_rvector mid(const l_ivector &v) { return _vmid<l_ivector,l_rvector>(v); }
	INLINE l_rvector mid(const l_ivector_slice &v) { return _vsmid<l_ivector_slice,l_rvector>(v); }
	INLINE l_rvector Inf(const l_ivector &v) { return _vinf<l_ivector,l_rvector>(v); }
	INLINE l_rvector Inf(const l_ivector_slice &v) { return _vsinf<l_ivector_slice,l_rvector>(v); }
	INLINE l_rvector Sup(const l_ivector &v) { return _vsup<l_ivector,l_rvector>(v); }
	INLINE l_rvector Sup(const l_ivector_slice &v) { return _vssup<l_ivector_slice,l_rvector>(v); }
	INLINE bool operator !(const l_ivector &rv) { return _vnot(rv); }
	INLINE bool operator !(const l_ivector_slice &sl) { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

//----------------------------- l_interval ---------------------------

	INLINE l_ivector operator *(const l_ivector &rv, const l_interval &s) { return _vsmult<l_ivector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_ivector_slice &sl, const l_interval &s) { return _vssmult<l_ivector_slice,l_interval,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const l_interval &s, const l_ivector &rv) { return _vsmult<l_ivector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_interval &s, const l_ivector_slice &sl) { return _vssmult<l_ivector_slice,l_interval,l_ivector>(sl,s); }
	INLINE l_ivector &operator *=(l_ivector &rv,const l_interval &r) { return _vsmultassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const l_interval &r) { return _vssmultassign(*this,r); }

	INLINE l_ivector operator /(const l_ivector &rv, const l_interval &s) { return _vsdiv<l_ivector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const l_ivector_slice &sl, const l_interval &s) { return _vssdiv<l_ivector_slice,l_interval,l_ivector>(sl,s); }
	INLINE l_ivector &operator /=(l_ivector &rv,const l_interval &r) { return _vsdivassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator /=(const l_interval &r) { return _vssdivassign(*this,r); }

//---------------------------- Real --------------------------------------

	INLINE l_ivector operator *(const l_ivector &rv, const real &s) { return _vsmult<l_ivector,real,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_ivector_slice &sl, const real &s) { return _vssmult<l_ivector_slice,real,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const real &s, const l_ivector &rv) { return _vsmult<l_ivector,real,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const real &s, const l_ivector_slice &sl) { return _vssmult<l_ivector_slice,real,l_ivector>(sl,s); }
	INLINE l_ivector &operator *=(l_ivector &rv,const real &r) { return _vsmultassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const real &r) { return _vssmultassign(*this,r); }

	INLINE l_ivector operator /(const l_ivector &rv, const real &s) { return _vsdiv<l_ivector,real,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const l_ivector_slice &sl, const real &s) { return _vssdiv<l_ivector_slice,real,l_ivector>(sl,s); }
	INLINE l_ivector &operator /=(l_ivector &rv,const real &r) { return _vsdivassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator /=(const real &r) { return _vssdivassign(*this,r); }

	INLINE l_ivector operator *(const rvector &rv, const l_interval &s) { return _vsmult<rvector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const rvector_slice &sl, const l_interval &s) { return _vssmult<rvector_slice,l_interval,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const l_interval &s, const rvector &rv) { return _vsmult<rvector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_interval &s, const rvector_slice &sl) { return _vssmult<rvector_slice,l_interval,l_ivector>(sl,s); }

	INLINE l_ivector operator /(const rvector &rv, const l_interval &s) { return _vsdiv<rvector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const rvector_slice &sl, const l_interval &s) { return _vssdiv<rvector_slice,l_interval,l_ivector>(sl,s); }

//---------------------------- l_real --------------------------------------

	INLINE l_ivector operator *(const l_ivector &rv, const l_real &s) { return _vsmult<l_ivector,l_real,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_ivector_slice &sl, const l_real &s) { return _vssmult<l_ivector_slice,l_real,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const l_real &s, const l_ivector &rv) { return _vsmult<l_ivector,l_real,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_real &s, const l_ivector_slice &sl) { return _vssmult<l_ivector_slice,l_real,l_ivector>(sl,s); }
	INLINE l_ivector &operator *=(l_ivector &rv,const l_real &r) { return _vsmultassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const l_real &r) { return _vssmultassign(*this,r); }

	INLINE l_ivector operator /(const l_ivector &rv, const l_real &s) { return _vsdiv<l_ivector,l_real,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const l_ivector_slice &sl, const l_real &s) { return _vssdiv<l_ivector_slice,l_real,l_ivector>(sl,s); }
	INLINE l_ivector &operator /=(l_ivector &rv,const l_real &r) { return _vsdivassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator /=(const l_real &r) { return _vssdivassign(*this,r); }

	INLINE l_ivector operator *(const l_rvector &rv, const l_interval &s) { return _vsmult<l_rvector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_rvector_slice &sl, const l_interval &s) { return _vssmult<l_rvector_slice,l_interval,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const l_interval &s, const l_rvector &rv) { return _vsmult<l_rvector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_interval &s, const l_rvector_slice &sl) { return _vssmult<l_rvector_slice,l_interval,l_ivector>(sl,s); }

	INLINE l_ivector operator /(const l_rvector &rv, const l_interval &s) { return _vsdiv<l_rvector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const l_rvector_slice &sl, const l_interval &s) { return _vssdiv<l_rvector_slice,l_interval,l_ivector>(sl,s); }

//---------------------------- interval --------------------------------------

	INLINE l_ivector operator *(const l_ivector &rv, const interval &s) { return _vsmult<l_ivector,interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_ivector_slice &sl, const interval &s) { return _vssmult<l_ivector_slice,interval,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const interval &s, const l_ivector &rv) { return _vsmult<l_ivector,interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const interval &s, const l_ivector_slice &sl) { return _vssmult<l_ivector_slice,interval,l_ivector>(sl,s); }
	INLINE l_ivector &operator *=(l_ivector &rv,const interval &r) { return _vsmultassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const interval &r) { return _vssmultassign(*this,r); }

	INLINE l_ivector operator /(const l_ivector &rv, const interval &s) { return _vsdiv<l_ivector,interval,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const l_ivector_slice &sl, const interval &s) { return _vssdiv<l_ivector_slice,interval,l_ivector>(sl,s); }
	INLINE l_ivector &operator /=(l_ivector &rv,const interval &r) { return _vsdivassign(rv,r); }
	INLINE l_ivector_slice &l_ivector_slice::operator /=(const interval &r) { return _vssdivassign(*this,r); }

	INLINE l_ivector operator *(const ivector &rv, const l_interval &s) { return _vsmult<ivector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const ivector_slice &sl, const l_interval &s) { return _vssmult<ivector_slice,l_interval,l_ivector>(sl,s); }
	INLINE l_ivector operator *(const l_interval &s, const ivector &rv) { return _vsmult<ivector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_interval &s, const ivector_slice &sl) { return _vssmult<ivector_slice,l_interval,l_ivector>(sl,s); }

	INLINE l_ivector operator /(const ivector &rv, const l_interval &s) { return _vsdiv<ivector,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator /(const ivector_slice &sl, const l_interval &s) { return _vssdiv<ivector_slice,l_interval,l_ivector>(sl,s); }

//======================= Vector / Vector ===============================


	INLINE std::ostream &operator <<(std::ostream &s, const l_ivector &rv) { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const l_ivector_slice &sl) { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, l_ivector &rv) { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, l_ivector_slice &rv) { return _vsin(s,rv); }
	
//----------------------- l_interval / l_interval ---------------------------
	INLINE l_ivector & l_ivector::operator =(const l_ivector_slice &sl) { return _vvsassign<l_ivector,l_ivector_slice,l_interval>(*this,sl); }

	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE l_interval operator *(const l_ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<l_ivector,l_ivector,l_interval>(rv1,rv2); }
	INLINE l_interval operator *(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<l_ivector_slice,l_ivector_slice,l_interval>(sl1,sl2); }
	
	INLINE const l_ivector &operator +(const l_ivector &rv) { return rv; }
	INLINE l_ivector operator +(const l_ivector_slice &sl) { return sl; }
	INLINE l_ivector operator +(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<l_ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<l_ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }
	INLINE l_ivector & operator +=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE l_ivector &operator +=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE l_ivector operator -(const l_ivector &rv) { return _vminus(rv); }
	INLINE l_ivector operator -(const l_ivector_slice &sl) { return _vsminus<l_ivector_slice,l_ivector>(sl); }
	INLINE l_ivector operator -(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<l_ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<l_ivector_slice,l_ivector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<l_ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }
	INLINE l_ivector & operator -=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE l_ivector &operator -=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE l_ivector operator |(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<l_ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<l_ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }
	INLINE l_ivector & operator |=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE l_ivector &operator |=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE l_ivector operator &(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<l_ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<l_ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }
	INLINE l_ivector & operator &=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE l_ivector &operator &=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssectassign(*this,sl2); }

	INLINE bool operator ==(const l_ivector &rv1, const l_ivector &rv2) { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const l_ivector_slice &sl1, const l_ivector_slice &sl2) { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const l_ivector_slice &sl, const l_ivector &rv) { return _vsveq(sl,rv); }
	INLINE bool operator ==(const l_ivector &rv, const l_ivector_slice &sl) { return _vsveq(sl,rv); }
	INLINE bool operator !=(const l_ivector &rv1, const l_ivector &rv2) { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const l_ivector_slice &sl1, const l_ivector_slice &sl2) { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const l_ivector_slice &sl, const l_ivector &rv) { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const l_ivector &rv, const l_ivector_slice &sl) { return _vsvneq(sl,rv); }
	INLINE bool operator <(const l_ivector &rv1, const l_ivector &rv2) { return _vvless(rv1,rv2); }
	INLINE bool operator <(const l_ivector_slice &sl1, const l_ivector_slice &sl2) { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const l_ivector_slice &sl, const l_ivector &rv) { return _vsvless(sl,rv); }
	INLINE bool operator < (const l_ivector &rv, const l_ivector_slice &sl) { return _vvsless(rv,sl); }
	INLINE bool operator <=(const l_ivector &rv1, const l_ivector &rv2) { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const l_ivector_slice &sl1, const l_ivector_slice &sl2) { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const l_ivector_slice &sl, const l_ivector &rv) { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const l_ivector &rv, const l_ivector_slice &sl) { return _vvsleq(rv,sl); }
	INLINE bool operator >(const l_ivector &rv1, const l_ivector &rv2) { return _vvless(rv2,rv1); }
	INLINE bool operator >(const l_ivector_slice &sl1, const l_ivector_slice &sl2) { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const l_ivector_slice &sl, const l_ivector &rv) { return _vvsless(rv,sl); }
	INLINE bool operator >(const l_ivector &rv, const l_ivector_slice &sl) { return _vsvless(sl,rv); }
	INLINE bool operator >=(const l_ivector &rv1, const l_ivector &rv2) { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const l_ivector_slice &sl1, const l_ivector_slice &sl2) { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const l_ivector_slice &sl, const l_ivector &rv) { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const l_ivector &rv, const l_ivector_slice &sl) { return _vsvleq(sl,rv); }

//-------------------------------- l_interval / Real --------------------------------

	INLINE l_ivector &l_ivector::operator =(const rvector &rv) { return _vvassign<l_ivector,rvector,l_interval>(*this,rv); }
	INLINE l_ivector &l_ivector::operator =(const real &r) { return _vsassign<l_ivector,real>(*this,r); }
	INLINE l_ivector & l_ivector::operator =(const rvector_slice &sl) { return _vvsassign<l_ivector,rvector_slice,l_interval>(*this,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const real &r) { return _vssassign<l_ivector_slice,real>(*this,r); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsassign(*this,sl); }

	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp,const l_ivector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp,const l_ivector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl2,sl1); }
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE l_interval operator *(const rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<rvector,l_ivector,l_interval>(rv1,rv2); }
	INLINE l_interval operator *(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<rvector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,rvector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<rvector_slice,l_ivector_slice,l_interval>(sl1,sl2); }
	
	INLINE l_interval operator *(const l_ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<rvector,l_ivector,l_interval>(rv2,rv1); }
	INLINE l_interval operator *(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,rvector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<rvector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<rvector_slice,l_ivector_slice,l_interval>(sl2,sl1); }
	
	INLINE l_ivector operator +(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator +(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<rvector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator +(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<rvector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator +=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE l_ivector &operator +=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE l_ivector operator -(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<rvector_slice,l_ivector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator -(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<l_ivector,rvector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<l_ivector_slice,rvector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<l_ivector_slice,rvector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector & operator -=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE l_ivector &operator -=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE l_ivector operator |(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator |(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<rvector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<rvector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator |=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE l_ivector &operator |=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE l_ivector operator &(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator &(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<rvector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator &(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<rvector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator &=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE l_ivector &operator &=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssectassign(*this,sl2); }

//-------------------------------- l_interval / l_real --------------------------------

	INLINE l_ivector &l_ivector::operator =(const l_rvector &rv) { return _vvassign<l_ivector,l_rvector,l_interval>(*this,rv); }
	INLINE l_ivector &l_ivector::operator =(const l_real &r) { return _vsassign<l_ivector,l_real>(*this,r); }
	INLINE l_ivector & l_ivector::operator =(const l_rvector_slice &sl) { return _vvsassign<l_ivector,l_rvector_slice,l_interval>(*this,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const l_real &r) { return _vssassign<l_ivector_slice,l_real>(*this,r); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsassign(*this,sl); }

	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp,const l_ivector_slice &sl,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp,const l_ivector &rv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl2,sl1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE l_interval operator *(const l_rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<l_rvector,l_ivector,l_interval>(rv1,rv2); }
	INLINE l_interval operator *(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_rvector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,l_rvector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<l_rvector_slice,l_ivector_slice,l_interval>(sl1,sl2); }
	
	INLINE l_interval operator *(const l_ivector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<l_rvector,l_ivector,l_interval>(rv2,rv1); }
	INLINE l_interval operator *(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,l_rvector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_rvector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<l_rvector_slice,l_ivector_slice,l_interval>(sl2,sl1); }
	
	INLINE l_ivector operator +(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<l_rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<l_rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator +(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<l_rvector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator +(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<l_rvector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator +=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE l_ivector &operator +=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE l_ivector operator -(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<l_rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<l_rvector_slice,l_ivector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<l_rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator -(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<l_ivector,l_rvector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<l_ivector_slice,l_rvector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<l_ivector_slice,l_rvector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector & operator -=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE l_ivector &operator -=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE l_ivector operator |(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<l_rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<l_rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator |(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<l_rvector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<l_rvector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator |=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE l_ivector &operator |=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE l_ivector operator &(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<l_rvector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<l_rvector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator &(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<l_rvector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator &(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_rvector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<l_rvector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator &=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE l_ivector &operator &=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssectassign(*this,sl2); }

//-------------------------------- l_interval / interval --------------------------------

	INLINE l_ivector &l_ivector::operator =(const ivector &rv) { return _vvassign<l_ivector,ivector,l_interval>(*this,rv); }
	INLINE l_ivector &l_ivector::operator =(const interval &r) { return _vsassign<l_ivector,interval>(*this,r); }
	INLINE l_ivector & l_ivector::operator =(const ivector_slice &sl) { return _vvsassign<l_ivector,ivector_slice,l_interval>(*this,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const interval &r) { return _vssassign<l_ivector_slice,interval>(*this,r); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsassign(*this,sl); }

	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp,const l_ivector_slice &sl,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp,const l_ivector &rv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl2,sl1); }
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE l_interval operator *(const ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<ivector,l_ivector,l_interval>(rv1,rv2); }
	INLINE l_interval operator *(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<ivector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<ivector_slice,l_ivector_slice,l_interval>(sl1,sl2); }
	
	INLINE l_interval operator *(const l_ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvlimult<ivector,l_ivector,l_interval>(rv2,rv1); }
	INLINE l_interval operator *(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<l_ivector_slice,ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvlimult<ivector_slice,l_ivector,l_interval>(sl,rv); }
	INLINE l_interval operator *(const l_ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvslimult<ivector_slice,l_ivector_slice,l_interval>(sl2,sl1); }
	
	INLINE l_ivector operator +(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator +(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<ivector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator +(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<ivector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator +=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE l_ivector &operator +=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator +=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE l_ivector operator -(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<ivector_slice,l_ivector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator -(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<l_ivector,ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<l_ivector_slice,ivector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<l_ivector_slice,ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector & operator -=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE l_ivector &operator -=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator -=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE l_ivector operator |(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator |(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<ivector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<ivector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator |=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE l_ivector &operator |=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator |=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE l_ivector operator &(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<ivector,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<ivector_slice,l_ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator &(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<ivector,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator &(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_ivector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<ivector,l_ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<ivector_slice,l_ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector & operator &=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE l_ivector &operator &=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE l_ivector_slice &l_ivector_slice::operator &=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssectassign(*this,sl2); }

//------------- real x l_real ------------------------
	INLINE l_ivector operator |(const rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<rvector,l_rvector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const l_rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<rvector,l_rvector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const rvector_slice &sl,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<rvector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const rvector &rv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<rvector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<rvector_slice,l_rvector_slice,l_ivector>(sl2,sl1); }
	INLINE l_ivector operator |(const rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<rvector_slice,l_rvector_slice,l_ivector>(sl1,sl2); }

//------------- l_real x l_real ------------------------
	INLINE l_ivector operator |(const l_rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<l_rvector,l_rvector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<l_rvector_slice,l_rvector_slice,l_ivector>(sl1,sl2); }

//-------------------------------- interval / l_real --------------------------------

	
	INLINE l_ivector operator +(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<l_rvector,ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<l_rvector_slice,ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator +(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvplus<l_rvector,ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator +(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsplus<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator +(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsplus<l_rvector_slice,ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector operator -(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<l_rvector,ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<l_rvector_slice,ivector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<l_rvector_slice,ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator -(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvminus<ivector,l_rvector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsminus<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator -(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvminus<ivector_slice,l_rvector,l_ivector>(sl,rv); }
	INLINE l_ivector operator -(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsminus<ivector_slice,l_rvector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator |(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<l_rvector,ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<l_rvector_slice,ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator |(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvconv<l_rvector,ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsconv<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator |(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvsconv<l_rvector_slice,ivector_slice,l_ivector>(sl2,sl1); }

	INLINE l_ivector operator &(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<l_rvector,ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<l_rvector_slice,ivector_slice,l_ivector>(sl1,sl2); }

	INLINE l_ivector operator &(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvsect<l_rvector,ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator &(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<ivector,l_rvector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vvssect<l_rvector,ivector_slice,l_ivector>(rv,sl); }
	INLINE l_ivector operator &(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvssect<l_rvector_slice,ivector_slice,l_ivector>(sl2,sl1); }

} // namespace cxsc

#endif // _CXSC_LIVECTOR_INL_INCLUDED
