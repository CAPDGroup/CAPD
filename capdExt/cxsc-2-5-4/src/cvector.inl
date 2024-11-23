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

/* CVS $Id: cvector.inl,v 1.28 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CVECTOR_INL_INCLUDED
#define _CXSC_CVECTOR_INL_INCLUDED

namespace cxsc {

	INLINE cvector::cvector () throw():dat(NULL),l(1),u(0),size(0)
	{
	}

	INLINE cvector::cvector(const int &i) throw():l(1),u(i),size(i)
	{
		dat=new complex[i];
	}

#ifdef OLD_CXSC
	INLINE cvector::cvector(const class index &i) throw():l(1),u(i._int()),size(i._int())
	{
		dat=new complex[i._int()];
	}
#endif

	INLINE cvector::cvector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_WRONG_BOUNDARIES,ERROR_CVECTOR_NO_MORE_MEMORY):l(i1),u(i2),size(i2-i1+1)
#else
	throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_CVECTOR_WRONG_BOUNDARIES("cvector::cvector(const int &i1,const int &i2)"));
#endif
		dat=new complex[size];
	}

	INLINE cvector::cvector(const cvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new complex[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE cvector::cvector(const cvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new complex[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE cvector::cvector(const complex &r) throw():l(1),u(1),size(1)
	{
		dat=new complex[1];
		*dat=r;
	}
	
	INLINE cvector::cvector(const rvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new complex[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE cvector::cvector(const rvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new complex[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE cvector::cvector(const real &r) throw():l(1),u(1),size(1)
	{
		dat=new complex[1];
		*dat=r;
	}
	
	INLINE complex & cvector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC("complex & cvector_slice::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}

	INLINE complex & cvector_slice::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC("complex & cvector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}

	INLINE complex & cvector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC("complex & cvector::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}

	INLINE complex & cvector::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC("complex & cvector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}


	/*!
	\param i The maximum dimension of the wanted part of the vector
	\return The wanted part of the vector

	\sa rvector::operator ()(const int &i)
	*/
	INLINE cvector_slice cvector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG("cvector_slice cvector::operator ()(const int &i)"));
#endif
		return cvector_slice(*this,1,i);
	}

	/*!
	\param i1 The starting dimension of the wanted part of the vector
	\param i2 The ending dimension of the wanted part of the vector
	\return The wanted part of the vector

	\sa rvector::operator ()(const int &i1,const int &i2)
	*/
   INLINE cvector_slice cvector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG("cvector_slice cvector::operator ()(const int &i1,const int &i2)"));
#endif
		return cvector_slice(*this,i1,i2);
	}
	
	INLINE cvector_slice cvector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG("cvector_slice cvector_slice::operator ()(const int &i)"));
#endif
		return cvector_slice(*this,1,i);
	}
	
   INLINE cvector_slice cvector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG("cvector_slice cvector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return cvector_slice(*this,i1,i2);
	}
	
	INLINE complex::complex(const cvector &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size>1) cxscthrow(ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ("complex::complex(const cvector &rv)"));
		else if(rv.size<1) cxscthrow(ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ("complex::complex(const cvector &rv)"));
#endif
		*this=rv.dat[0];
	}
	
	INLINE complex::complex(const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size>1) cxscthrow(ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ("complex::complex(const cvector_slice &sl)"));
		else if(sl.size<1) cxscthrow(ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ("complex::complex(const cvector_slice &sl)"));
#endif
		*this=sl.dat[sl.start-sl.l];
	}
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cvector::cvector(const complex& r)
	*/
	INLINE cvector _cvector(const complex &r) throw() { return cvector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cvector::cvector(const real &)
	*/
	INLINE cvector _cvector(const real &r) throw() { return cvector(r); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cvector::cvector(const rvector_slice &rs)
	*/
	INLINE cvector _cvector(const rvector_slice &rs) throw() { return cvector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cvector::cvector(const rvector &v)
	*/
	INLINE cvector _cvector(const rvector &rs) throw() { return cvector(rs); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cvector::cvector(const rmatrix &)
	*/
	INLINE cvector _cvector(const rmatrix_subv &rs) throw() { return cvector(rs); }
	INLINE cvector &cvector::operator =(const cvector &rv) throw() { return _vvassign<cvector,cvector,complex>(*this,rv); }
	INLINE cvector &cvector::operator =(const complex &r) throw() { return _vsassign<cvector,complex>(*this,r); }
	INLINE cvector &cvector::operator =(const rvector &rv) throw() { return _vvassign<cvector,rvector,complex>(*this,rv); }
	INLINE cvector &cvector::operator =(const real &r) throw() { return _vsassign<cvector,real>(*this,r); }
	INLINE cvector::operator void*() throw() { return _vvoid(*this); }
	INLINE cvector_slice & cvector_slice::operator =(const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE cvector_slice & cvector_slice::operator =(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE cvector_slice & cvector_slice::operator =(const complex &r) throw() { return _vssassign<cvector_slice,complex>(*this,r); }
	INLINE cvector_slice & cvector_slice::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,cvector(m)); }
	INLINE cvector_slice & cvector_slice::operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE cvector_slice & cvector_slice::operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE cvector_slice & cvector_slice::operator =(const real &r) throw() { return _vssassign(*this,r); }
	INLINE cvector_slice::operator void*() throw() { return _vsvoid(*this); }

//=======================================================================
//======================== Vector Functions =============================

	INLINE cvector &SetRe(cvector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return iv;} // S.W. temp return _vvsetre(iv,rv); }
	INLINE cvector_slice &SetRe(cvector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetre(iv,rv); }
	INLINE cvector &SetRe(cvector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetre(iv,rv); }
	INLINE cvector_slice &SetRe(cvector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetre(iv,rv); }

	INLINE cvector &SetIm(cvector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetim(iv,rv); }
	INLINE cvector_slice &SetIm(cvector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetim(iv,rv); }
	INLINE cvector &SetIm(cvector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetim(iv,rv); }
	INLINE cvector_slice &SetIm(cvector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetim(iv,rv); }

	INLINE cvector &SetRe(cvector &iv,const real &r) throw() { return _vssetre(iv,r); }
	INLINE cvector &SetIm(cvector &iv,const real &r) throw() { return _vssetim(iv,r); }
	INLINE cvector_slice &SetRe(cvector_slice &iv,const real &r) throw() { return _vsssetre(iv,r); }
	INLINE cvector_slice &SetIm(cvector_slice &iv,const real &r) throw() { return _vsssetim(iv,r); }

	INLINE void Resize(cvector &rv) throw() { _vresize(rv); } 
	INLINE void Resize(cvector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<cvector>)
#else
	throw()
#endif
	{ _vresize<class cvector,class complex>(rv,len); }
	INLINE void Resize(cvector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<cvector>)
#else
	throw()
#endif
	{ _vresize<class cvector,class complex>(rv,lb,ub); }
	
	INLINE cvector conj(const cvector &rv) throw() { return _vconj<cvector>(rv); }
	INLINE cvector conj(const cvector_slice &sl) throw() { return _vsconj<cvector_slice,cvector>(sl); }
	
	INLINE rvector abs(const cvector &rv) throw() { return _vabs<cvector,rvector>(rv); }
	INLINE rvector abs(const cvector_slice &sl) throw() { return _vsabs<cvector_slice,rvector>(sl); }
	INLINE rvector Im(const cvector &v) throw() { return _vim<cvector,rvector>(v); }
	INLINE rvector Im(const cvector_slice &v) throw() { return _vsim<cvector_slice,rvector>(v); }
	INLINE rvector Re(const cvector &v) throw() { return _vre<cvector,rvector>(v); }
	INLINE rvector Re(const cvector_slice &v) throw() { return _vsre<cvector_slice,rvector>(v); }
	INLINE bool operator !(const cvector &rv) throw() { return _vnot(rv); }
	INLINE bool operator !(const cvector_slice &sl) throw() { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

//----------------------------- complex ---------------------------

	INLINE cvector operator *(const cvector &rv, const complex &s) throw() { return _vsmult<cvector,complex,cvector>(rv,s); }
	INLINE cvector operator *(const cvector_slice &sl, const complex &s) throw() { return _vssmult<cvector_slice,complex,cvector>(sl,s); }
	INLINE cvector operator *(const complex &s, const cvector &rv) throw() { return _vsmult<cvector,complex,cvector>(rv,s); }
	INLINE cvector operator *(const complex &s, const cvector_slice &sl) throw() { return _vssmult<cvector_slice,complex,cvector>(sl,s); }
	INLINE cvector &operator *=(cvector &rv,const complex &r) throw() { return _vsmultassign(rv,r); }
	INLINE cvector_slice &cvector_slice::operator *=(const complex &r) throw() { return _vssmultassign(*this,r); }

	INLINE cvector operator /(const cvector &rv, const complex &s) throw() { return _vsdiv<cvector,complex,cvector>(rv,s); }
	INLINE cvector operator /(const cvector_slice &sl, const complex &s) throw() { return _vssdiv<cvector_slice,complex,cvector>(sl,s); }
	INLINE cvector &operator /=(cvector &rv,const complex &r) throw() { return _vsdivassign(rv,r); }
	INLINE cvector_slice &cvector_slice::operator /=(const complex &r) throw() { return _vssdivassign(*this,r); }

//---------------------------- Real --------------------------------------

	INLINE cvector operator *(const cvector &rv, const real &s) throw() { return _vsmult<cvector,real,cvector>(rv,s); }
	INLINE cvector operator *(const cvector_slice &sl, const real &s) throw() { return _vssmult<cvector_slice,real,cvector>(sl,s); }
	INLINE cvector operator *(const real &s, const cvector &rv) throw() { return _vsmult<cvector,real,cvector>(rv,s); }
	INLINE cvector operator *(const real &s, const cvector_slice &sl) throw() { return _vssmult<cvector_slice,real,cvector>(sl,s); }
	INLINE cvector &operator *=(cvector &rv,const real &r) throw() { return _vsmultassign(rv,r); }
	INLINE cvector_slice &cvector_slice::operator *=(const real &r) throw() { return _vssmultassign(*this,r); }

	INLINE cvector operator /(const cvector &rv, const real &s) throw() { return _vsdiv<cvector,real,cvector>(rv,s); }
	INLINE cvector operator /(const cvector_slice &sl, const real &s) throw() { return _vssdiv<cvector_slice,real,cvector>(sl,s); }
	INLINE cvector &operator /=(cvector &rv,const real &r) throw() { return _vsdivassign(rv,r); }
	INLINE cvector_slice &cvector_slice::operator /=(const real &r) throw() { return _vssdivassign(*this,r); }

	INLINE cvector operator *(const rvector &rv, const complex &s) throw() { return _vsmult<rvector,complex,cvector>(rv,s); }
	INLINE cvector operator *(const rvector_slice &sl, const complex &s) throw() { return _vssmult<rvector_slice,complex,cvector>(sl,s); }
	INLINE cvector operator *(const complex &s, const rvector &rv) throw() { return _vsmult<rvector,complex,cvector>(rv,s); }
	INLINE cvector operator *(const complex &s, const rvector_slice &sl) throw() { return _vssmult<rvector_slice,complex,cvector>(sl,s); }

	INLINE cvector operator /(const rvector &rv, const complex &s) throw() { return _vsdiv<rvector,complex,cvector>(rv,s); }
	INLINE cvector operator /(const rvector_slice &sl, const complex &s) throw() { return _vssdiv<rvector_slice,complex,cvector>(sl,s); }

//======================= Vector / Vector ===============================


	INLINE std::ostream &operator <<(std::ostream &s, const cvector &rv) throw() { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const cvector_slice &sl) throw() { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, cvector &rv) throw() { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, cvector_slice &rv) throw() { return _vsin(s,rv); }
	
//----------------------- complex / complex ---------------------------
	INLINE cvector & cvector::operator =(const cvector_slice &sl) throw() { return _vvsassign<cvector,cvector_slice,complex>(*this,sl); }



	INLINE complex operator *(const cvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvcmult<cvector,cvector,complex>(rv1,rv2); }
	INLINE complex operator *(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvcmult<cvector_slice,cvector,complex>(sl,rv); }
	INLINE complex operator *(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvcmult<cvector_slice,cvector,complex>(sl,rv); }
	INLINE complex operator *(const cvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvscmult<cvector_slice,cvector_slice,complex>(sl1,sl2); }
	
	INLINE const cvector &operator +(const cvector &rv) throw() { return rv; }
	INLINE cvector operator +(const cvector_slice &sl) throw() { return sl; }
	INLINE cvector operator +(const cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvplus<cvector,cvector,cvector>(rv1,rv2); }
	INLINE cvector operator +(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,cvector_slice,cvector>(rv,sl); }
	INLINE cvector operator +(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,cvector_slice,cvector>(rv,sl); }
	INLINE cvector operator +(const cvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsplus<cvector_slice,cvector_slice,cvector>(sl1,sl2); }
	INLINE cvector & operator +=(cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE cvector &operator +=(cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE cvector_slice &cvector_slice::operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE cvector_slice &cvector_slice::operator +=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE cvector operator -(const cvector &rv) throw() { return _vminus(rv); }
	INLINE cvector operator -(const cvector_slice &sl) throw() { return _vsminus<cvector_slice,cvector>(sl); }
	INLINE cvector operator -(const cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvminus<cvector,cvector,cvector>(rv1,rv2); }
	INLINE cvector operator -(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsminus<cvector,cvector_slice,cvector>(rv,sl); }
	INLINE cvector operator -(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvminus<cvector_slice,cvector,cvector>(sl,rv); }
	INLINE cvector operator -(const cvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsminus<cvector_slice,cvector_slice,cvector>(sl1,sl2); }
	INLINE cvector & operator -=(cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE cvector &operator -=(cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE cvector_slice &cvector_slice::operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE cvector_slice &cvector_slice::operator -=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE bool operator ==(const cvector &rv1, const cvector &rv2) throw() { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const cvector_slice &sl1, const cvector_slice &sl2) throw() { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const cvector_slice &sl, const cvector &rv) throw() { return _vsveq(sl,rv); }
	INLINE bool operator ==(const cvector &rv, const cvector_slice &sl) throw() { return _vsveq(sl,rv); }
	INLINE bool operator !=(const cvector &rv1, const cvector &rv2) throw() { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const cvector_slice &sl1, const cvector_slice &sl2) throw() { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const cvector_slice &sl, const cvector &rv) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const cvector &rv, const cvector_slice &sl) throw() { return _vsvneq(sl,rv); }
/*	INLINE bool operator <(const cvector &rv1, const cvector &rv2) throw() { return _vvless(rv1,rv2); }
	INLINE bool operator <(const cvector_slice &sl1, const cvector_slice &sl2) throw() { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const cvector_slice &sl, const cvector &rv) throw() { return _vsvless(sl,rv); }
	INLINE bool operator < (const cvector &rv, const cvector_slice &sl) throw() { return _vvsless(rv,sl); }
	INLINE bool operator <=(const cvector &rv1, const cvector &rv2) throw() { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const cvector_slice &sl1, const cvector_slice &sl2) throw() { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const cvector_slice &sl, const cvector &rv) throw() { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const cvector &rv, const cvector_slice &sl) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >(const cvector &rv1, const cvector &rv2) throw() { return _vvless(rv2,rv1); }
	INLINE bool operator >(const cvector_slice &sl1, const cvector_slice &sl2) throw() { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const cvector_slice &sl, const cvector &rv) throw() { return _vvsless(rv,sl); }
	INLINE bool operator >(const cvector &rv, const cvector_slice &sl) throw() { return _vsvless(sl,rv); }
	INLINE bool operator >=(const cvector &rv1, const cvector &rv2) throw() { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const cvector_slice &sl1, const cvector_slice &sl2) throw() { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const cvector_slice &sl, const cvector &rv) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const cvector &rv, const cvector_slice &sl) throw() { return _vsvleq(sl,rv); }
*/
//-------------------------------- complex / Real --------------------------------

	INLINE cvector & cvector::operator =(const rvector_slice &sl) throw() { return _vvsassign<cvector,rvector_slice,complex>(*this,sl); }


	INLINE complex operator *(const rvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvcmult<rvector,cvector,complex>(rv1,rv2); }
	INLINE complex operator *(const rvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvcmult<rvector_slice,cvector,complex>(sl,rv); }
	INLINE complex operator *(const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvcmult<cvector_slice,rvector,complex>(sl,rv); }
	INLINE complex operator *(const rvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvscmult<rvector_slice,cvector_slice,complex>(sl1,sl2); }
	
	INLINE complex operator *(const cvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvcmult<rvector,cvector,complex>(rv2,rv1); }
	INLINE complex operator *(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvcmult<cvector_slice,rvector,complex>(sl,rv); }
	INLINE complex operator *(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvcmult<rvector_slice,cvector,complex>(sl,rv); }
	INLINE complex operator *(const cvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvscmult<rvector_slice,cvector_slice,complex>(sl2,sl1); }
	
	INLINE cvector operator +(const rvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,cvector,cvector>(rv1,rv2); }
	INLINE cvector operator +(const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,cvector_slice,cvector>(rv,sl); }
	INLINE cvector operator +(const rvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,rvector_slice,cvector>(rv,sl); }
	INLINE cvector operator +(const rvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,cvector_slice,cvector>(sl1,sl2); }

	INLINE cvector operator +(const cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,cvector,cvector>(rv2,rv1); }
	INLINE cvector operator +(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,rvector_slice,cvector>(rv,sl); }
	INLINE cvector operator +(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,cvector_slice,cvector>(rv,sl); }
	INLINE cvector operator +(const cvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,cvector_slice,cvector>(sl2,sl1); }

	INLINE cvector & operator +=(cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE cvector &operator +=(cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE cvector_slice &cvector_slice::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE cvector_slice &cvector_slice::operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE cvector operator -(const rvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvminus<rvector,cvector,cvector>(rv1,rv2); }
	INLINE cvector operator -(const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsminus<rvector,cvector_slice,cvector>(rv,sl); }
	INLINE cvector operator -(const rvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvminus<rvector_slice,cvector,cvector>(sl,rv); }
	INLINE cvector operator -(const rvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsminus<rvector_slice,cvector_slice,cvector>(sl1,sl2); }

	INLINE cvector operator -(const cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvminus<cvector,rvector,cvector>(rv1,rv2); }
	INLINE cvector operator -(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsminus<cvector,rvector_slice,cvector>(rv,sl); }
	INLINE cvector operator -(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvminus<cvector_slice,rvector,cvector>(sl,rv); }
	INLINE cvector operator -(const cvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsminus<cvector_slice,rvector_slice,cvector>(sl1,sl2); }

	INLINE cvector & operator -=(cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE cvector &operator -=(cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE cvector_slice &cvector_slice::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE cvector_slice &cvector_slice::operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE cvector cvector::operator()(const intvector& p) {
          cvector x(*this);
          for(int i=0 ; i<VecLen(x) ; i++)
              x[i+Lb(x)] = (*this)[p[i+Lb(p)]+Lb(*this)];
          return x;
        }

} // namespace cxsc

#endif // _CXSC_CVECTOR_INL_INCLUDED

