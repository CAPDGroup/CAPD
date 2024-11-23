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

/* CVS $Id: l_rvector.inl,v 1.19 2014/01/30 17:23:47 cxsc Exp $ */

#ifndef _CXSC_LRVECTOR_INL_INCLUDED
#define _CXSC_LRVECTOR_INL_INCLUDED

namespace cxsc {

	INLINE l_rvector::l_rvector () throw():dat(NULL),l(1),u(0),size(0)
	{
	}

	INLINE l_rvector::l_rvector(const int &i) throw():l(1),u(i),size(i)
	{
		dat=new l_real[i];
	}

#ifdef OLD_CXSC
	INLINE l_rvector::l_rvector(const class index &i) throw():l(1),u(i._int()),size(i._int())
	{
		dat=new l_real[i._int()];
	}
#endif

	INLINE l_rvector::l_rvector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_WRONG_BOUNDARIES,ERROR_LRVECTOR_NO_MORE_MEMORY):l(i1),u(i2),size(i2-i1+1)
#else
	throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_LRVECTOR_WRONG_BOUNDARIES("l_rvector::l_rvector(const int &i1,const int &i2)"));
#endif
		dat=new l_real[size];
	}

	INLINE l_rvector::l_rvector(const l_rvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new l_real[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE l_rvector::l_rvector(const l_rvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new l_real[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE l_rvector::l_rvector(const l_real &r) throw():l(1),u(1),size(1)
	{
		dat=new l_real[1];
		*dat=r;
	}
	
	INLINE l_rvector::l_rvector(const rvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new l_real[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE l_rvector::l_rvector(const rvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new l_real[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE l_rvector::l_rvector(const real &r) throw():l(1),u(1),size(1)
	{
		dat=new l_real[1];
		*dat=r;
	}
	
	INLINE l_real & l_rvector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC("l_real & l_rvector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE l_real & l_rvector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC("l_real & l_rvector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
	
	INLINE l_rvector_slice l_rvector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG("l_rvector_slice l_rvector::operator ()(const int &i)"));
#endif
		return l_rvector_slice(*this,1,i);
	}
	
   INLINE l_rvector_slice l_rvector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG("l_rvector_slice l_rvector::operator ()(const int &i1,const int &i2)"));
#endif
		return l_rvector_slice(*this,i1,i2);
	}
	
	INLINE l_rvector_slice l_rvector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG("l_rvector_slice l_rvector_slice::operator ()(const int &i)"));
#endif
		return l_rvector_slice(*this,1,i);
	}
	
   INLINE l_rvector_slice l_rvector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG("l_rvector_slice l_rvector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return l_rvector_slice(*this,i1,i2);
	}
	
	INLINE l_real::l_real(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size>1) cxscthrow(ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ("l_real::l_real(const l_rvector &rv)"));
		else if(rv.size<1) cxscthrow(ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ("l_real::l_real(const l_rvector &rv)"));
#endif
		*this=rv.dat[0];
	}
	
	INLINE l_real::l_real(const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size>1) cxscthrow(ERROR_LRVECTOR_TYPE_CAST_OF_THICK_OBJ("l_real::l_real(const l_rvector_slice &sl)"));
		else if(sl.size<1) cxscthrow(ERROR_LRVECTOR_USE_OF_UNINITIALIZED_OBJ("l_real::l_real(const l_rvector_slice &sl)"));
#endif
		*this=sl.dat[sl.start-sl.l];
	}
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const l_real &)
	*/
	INLINE l_rvector _l_rvector(const l_real &r) throw() { return l_rvector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const real &)
	*/
	INLINE l_rvector _l_rvector(const real &r) throw() { return l_rvector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const rvector_slice &rs)
	*/
	INLINE l_rvector _l_rvector(const rvector_slice &rs) throw() { return l_rvector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const rvector &v)
	*/
	INLINE l_rvector _l_rvector(const rvector &rs) throw() { return l_rvector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const rmatrix_subv &)
	*/
	INLINE l_rvector &l_rvector::operator =(const l_rvector &rv) throw() { return _vvassign<l_rvector,l_rvector,l_real>(*this,rv); }
	INLINE l_rvector &l_rvector::operator =(const l_real &r) throw() { return _vsassign<l_rvector,l_real>(*this,r); }
	INLINE l_rvector &l_rvector::operator =(const rvector &rv) throw() { return _vvassign<l_rvector,rvector,l_real>(*this,rv); }
	INLINE l_rvector &l_rvector::operator =(const real &r) throw() { return _vsassign<l_rvector,real>(*this,r); }
	INLINE l_rvector::operator void*() throw() { return _vvoid(*this); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const l_real &r) throw() { return _vssassign<l_rvector_slice,l_real>(*this,r); }

	INLINE l_rvector_slice & l_rvector_slice::operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const real &r) throw() { return _vssassign(*this,r); }
	INLINE l_rvector_slice::operator void*() throw() { return _vsvoid(*this); }

//=======================================================================
//======================== Vector Functions =============================

	INLINE void Resize(l_rvector &rv) throw() { _vresize(rv); } 
	INLINE void Resize(l_rvector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<l_rvector>)
#else
	throw()
#endif
	{ _vresize<class l_rvector,class l_real>(rv,len); }
	INLINE void Resize(l_rvector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<l_rvector>)
#else
	throw()
#endif
	{ _vresize<class l_rvector,class l_real>(rv,lb,ub); }
	
	INLINE l_rvector abs(const l_rvector &rv) throw() { return _vabs<l_rvector,l_rvector>(rv); }
	INLINE l_rvector abs(const l_rvector_slice &sl) throw() { return _vsabs<l_rvector_slice,l_rvector>(sl); }
	INLINE bool operator !(const l_rvector &rv) throw() { return _vnot(rv); }
	INLINE bool operator !(const l_rvector_slice &sl) throw() { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

//----------------------------- l_real ---------------------------

	INLINE l_rvector operator *(const l_rvector &rv, const l_real &s) throw() { return _vsmult<l_rvector,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const l_rvector_slice &sl, const l_real &s) throw() { return _vssmult<l_rvector_slice,l_real,l_rvector>(sl,s); }
	INLINE l_rvector operator *(const l_real &s, const l_rvector &rv) throw() { return _vsmult<l_rvector,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const l_real &s, const l_rvector_slice &sl) throw() { return _vssmult<l_rvector_slice,l_real,l_rvector>(sl,s); }
	INLINE l_rvector &operator *=(l_rvector &rv,const l_real &r) throw() { return _vsmultassign(rv,r); }
	INLINE l_rvector_slice &l_rvector_slice::operator *=(const l_real &r) throw() { return _vssmultassign(*this,r); }

	INLINE l_rvector operator /(const l_rvector &rv, const l_real &s) throw() { return _vsdiv<l_rvector,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator /(const l_rvector_slice &sl, const l_real &s) throw() { return _vssdiv<l_rvector_slice,l_real,l_rvector>(sl,s); }
	INLINE l_rvector &operator /=(l_rvector &rv,const l_real &r) throw() { return _vsdivassign(rv,r); }
	INLINE l_rvector_slice &l_rvector_slice::operator /=(const l_real &r) throw() { return _vssdivassign(*this,r); }

//---------------------------- Real --------------------------------------

	INLINE l_rvector operator *(const l_rvector &rv, const real &s) throw() { return _vsmult<l_rvector,real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const l_rvector_slice &sl, const real &s) throw() { return _vssmult<l_rvector_slice,real,l_rvector>(sl,s); }
	INLINE l_rvector operator *(const real &s, const l_rvector &rv) throw() { return _vsmult<l_rvector,real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const real &s, const l_rvector_slice &sl) throw() { return _vssmult<l_rvector_slice,real,l_rvector>(sl,s); }
	INLINE l_rvector &operator *=(l_rvector &rv,const real &r) throw() { return _vsmultassign(rv,r); }
	INLINE l_rvector_slice &l_rvector_slice::operator *=(const real &r) throw() { return _vssmultassign(*this,r); }

	INLINE l_rvector operator /(const l_rvector &rv, const real &s) throw() { return _vsdiv<l_rvector,real,l_rvector>(rv,s); }
	INLINE l_rvector operator /(const l_rvector_slice &sl, const real &s) throw() { return _vssdiv<l_rvector_slice,real,l_rvector>(sl,s); }
	INLINE l_rvector &operator /=(l_rvector &rv,const real &r) throw() { return _vsdivassign(rv,r); }
	INLINE l_rvector_slice &l_rvector_slice::operator /=(const real &r) throw() { return _vssdivassign(*this,r); }

	INLINE l_rvector operator *(const rvector &rv, const l_real &s) throw() { return _vsmult<rvector,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const rvector_slice &sl, const l_real &s) throw() { return _vssmult<rvector_slice,l_real,l_rvector>(sl,s); }
	INLINE l_rvector operator *(const l_real &s, const rvector &rv) throw() { return _vsmult<rvector,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const l_real &s, const rvector_slice &sl) throw() { return _vssmult<rvector_slice,l_real,l_rvector>(sl,s); }

	INLINE l_rvector operator /(const rvector &rv, const l_real &s) throw() { return _vsdiv<rvector,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator /(const rvector_slice &sl, const l_real &s) throw() { return _vssdiv<rvector_slice,l_real,l_rvector>(sl,s); }

//======================= Vector / Vector ===============================


	INLINE std::ostream &operator <<(std::ostream &s, const l_rvector &rv) throw() { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const l_rvector_slice &sl) throw() { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, l_rvector &rv) throw() { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, l_rvector_slice &rv) throw() { return _vsin(s,rv); }
	
//----------------------- l_real / l_real ---------------------------
	INLINE l_rvector & l_rvector::operator =(const l_rvector_slice &sl) throw() { return _vvsassign<l_rvector,l_rvector_slice,l_real>(*this,sl); }

	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rvector_slice & sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp, const l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rvector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl1,sl2); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl1,sl2); }


	INLINE l_real operator *(const l_rvector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvlmult<l_rvector,l_rvector,l_real>(rv1,rv2); }
	INLINE l_real operator *(const l_rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvlmult<l_rvector_slice,l_rvector,l_real>(sl,rv); }
	INLINE l_real operator *(const l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvlmult<l_rvector_slice,l_rvector,l_real>(sl,rv); }
	INLINE l_real operator *(const l_rvector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvslmult<l_rvector_slice,l_rvector_slice,l_real>(sl1,sl2); }
	
	INLINE const l_rvector &operator +(const l_rvector &rv) throw() { return rv; }
	INLINE l_rvector operator +(const l_rvector_slice &sl) throw() { return sl; }
	INLINE l_rvector operator +(const l_rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvplus<l_rvector,l_rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator +(const l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplus<l_rvector,l_rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator +(const l_rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplus<l_rvector,l_rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator +(const l_rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsplus<l_rvector_slice,l_rvector_slice,l_rvector>(sl1,sl2); }
	INLINE l_rvector & operator +=(l_rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE l_rvector &operator +=(l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE l_rvector_slice &l_rvector_slice::operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE l_rvector_slice &l_rvector_slice::operator +=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE l_rvector operator -(const l_rvector &rv) throw() { return _vminus(rv); }
	INLINE l_rvector operator -(const l_rvector_slice &sl) throw() { return _vsminus<l_rvector_slice,l_rvector>(sl); }
	INLINE l_rvector operator -(const l_rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvminus<l_rvector,l_rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator -(const l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsminus<l_rvector,l_rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator -(const l_rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvminus<l_rvector_slice,l_rvector,l_rvector>(sl,rv); }
	INLINE l_rvector operator -(const l_rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsminus<l_rvector_slice,l_rvector_slice,l_rvector>(sl1,sl2); }
	INLINE l_rvector & operator -=(l_rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE l_rvector &operator -=(l_rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE l_rvector_slice &l_rvector_slice::operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE l_rvector_slice &l_rvector_slice::operator -=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE bool operator ==(const l_rvector &rv1, const l_rvector &rv2) throw() { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const l_rvector_slice &sl1, const l_rvector_slice &sl2) throw() { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const l_rvector_slice &sl, const l_rvector &rv) throw() { return _vsveq(sl,rv); }
	INLINE bool operator ==(const l_rvector &rv, const l_rvector_slice &sl) throw() { return _vsveq(sl,rv); }
	INLINE bool operator !=(const l_rvector &rv1, const l_rvector &rv2) throw() { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const l_rvector_slice &sl1, const l_rvector_slice &sl2) throw() { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const l_rvector_slice &sl, const l_rvector &rv) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const l_rvector &rv, const l_rvector_slice &sl) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator <(const l_rvector &rv1, const l_rvector &rv2) throw() { return _vvless(rv1,rv2); }
	INLINE bool operator <(const l_rvector_slice &sl1, const l_rvector_slice &sl2) throw() { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const l_rvector_slice &sl, const l_rvector &rv) throw() { return _vsvless(sl,rv); }
	INLINE bool operator < (const l_rvector &rv, const l_rvector_slice &sl) throw() { return _vvsless(rv,sl); }
	INLINE bool operator <=(const l_rvector &rv1, const l_rvector &rv2) throw() { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const l_rvector_slice &sl1, const l_rvector_slice &sl2) throw() { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const l_rvector_slice &sl, const l_rvector &rv) throw() { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const l_rvector &rv, const l_rvector_slice &sl) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >(const l_rvector &rv1, const l_rvector &rv2) throw() { return _vvless(rv2,rv1); }
	INLINE bool operator >(const l_rvector_slice &sl1, const l_rvector_slice &sl2) throw() { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const l_rvector_slice &sl, const l_rvector &rv) throw() { return _vvsless(rv,sl); }
	INLINE bool operator >(const l_rvector &rv, const l_rvector_slice &sl) throw() { return _vsvless(sl,rv); }
	INLINE bool operator >=(const l_rvector &rv1, const l_rvector &rv2) throw() { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const l_rvector_slice &sl1, const l_rvector_slice &sl2) throw() { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const l_rvector_slice &sl, const l_rvector &rv) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const l_rvector &rv, const l_rvector_slice &sl) throw() { return _vsvleq(sl,rv); }

//-------------------------------- l_real / Real --------------------------------

	INLINE l_rvector & l_rvector::operator =(const rvector_slice &sl) throw() { return _vvsassign<l_rvector,rvector_slice,l_real>(*this,sl); }

	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv2,rv1); }
	INLINE void accumulate(dotprecision &dp, const rvector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const rvector_slice & sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp,const l_rvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp, const rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp, const rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp,const l_rvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(dotprecision &dp, const l_rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl2,sl1); }
	INLINE void accumulate(dotprecision &dp, const rvector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp,const l_rvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp,const l_rvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvaccu(dp,sl,rv); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	;
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl2,sl1); }
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsaccu(dp,sl1,sl2); }

	INLINE l_real operator *(const rvector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvlmult<rvector,l_rvector,l_real>(rv1,rv2); }
	INLINE l_real operator *(const rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvlmult<rvector_slice,l_rvector,l_real>(sl,rv); }
	INLINE l_real operator *(const rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvlmult<l_rvector_slice,rvector,l_real>(sl,rv); }
	INLINE l_real operator *(const rvector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvslmult<rvector_slice,l_rvector_slice,l_real>(sl1,sl2); }
	
	INLINE l_real operator *(const l_rvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvlmult<rvector,l_rvector,l_real>(rv2,rv1); }
	INLINE l_real operator *(const l_rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvlmult<l_rvector_slice,rvector,l_real>(sl,rv); }
	INLINE l_real operator *(const l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvlmult<rvector_slice,l_rvector,l_real>(sl,rv); }
	INLINE l_real operator *(const l_rvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvslmult<rvector_slice,l_rvector_slice,l_real>(sl2,sl1); }
	
	INLINE l_rvector operator +(const rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,l_rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator +(const rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,l_rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator +(const rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplus<l_rvector,rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator +(const rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,l_rvector_slice,l_rvector>(sl1,sl2); }

	INLINE l_rvector operator +(const l_rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,l_rvector,l_rvector>(rv2,rv1); }
	INLINE l_rvector operator +(const l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplus<l_rvector,rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator +(const l_rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,l_rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator +(const l_rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,l_rvector_slice,l_rvector>(sl2,sl1); }

	INLINE l_rvector & operator +=(l_rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE l_rvector &operator +=(l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE l_rvector_slice &l_rvector_slice::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE l_rvector_slice &l_rvector_slice::operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE l_rvector operator -(const rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvminus<rvector,l_rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator -(const rvector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsminus<rvector,l_rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator -(const rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvminus<rvector_slice,l_rvector,l_rvector>(sl,rv); }
	INLINE l_rvector operator -(const rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsminus<rvector_slice,l_rvector_slice,l_rvector>(sl1,sl2); }

	INLINE l_rvector operator -(const l_rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvminus<l_rvector,rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator -(const l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsminus<l_rvector,rvector_slice,l_rvector>(rv,sl); }
	INLINE l_rvector operator -(const l_rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvminus<l_rvector_slice,rvector,l_rvector>(sl,rv); }
	INLINE l_rvector operator -(const l_rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsminus<l_rvector_slice,rvector_slice,l_rvector>(sl1,sl2); }

	INLINE l_rvector & operator -=(l_rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE l_rvector &operator -=(l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE l_rvector_slice &l_rvector_slice::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE l_rvector_slice &l_rvector_slice::operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

} // namespace cxsc

#endif

