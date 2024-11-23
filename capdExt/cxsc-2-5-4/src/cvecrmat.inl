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

/* CVS $Id: cvecrmat.inl,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

// Here are definitions for cvector x rmatrix-Functions
#ifndef _CXSC_CVECRMAT_INL_INCLUDED
#define _CXSC_CVECRMAT_INL_INCLUDED

namespace cxsc {

	INLINE cvector::cvector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<cvector,rmatrix,complex>(*this,sl); }
	INLINE cvector::cvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<cvector,rmatrix_slice,complex>(*this,sl); }
	INLINE cvector::cvector(const rmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new complex[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	INLINE cvector _cvector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return cvector(sl); }
	INLINE cvector _cvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return cvector(sl); }

// 	INLINE void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cdotprecision,cvector,rmatrix_subv>(dp,rv2,rv1); }
// 	INLINE void accumulate(cdotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cdotprecision,cvector,rmatrix_subv>(dp,rv1,rv2); }
// 	INLINE void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cidotprecision,cvector,rmatrix_subv>(dp,rv2,rv1); }
// 	INLINE void accumulate(cidotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cidotprecision,cvector,rmatrix_subv>(dp,rv1,rv2); }
// 	
// 	INLINE void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cdotprecision,cvector,rmatrix_subv>(dp,cvector(rv2),rv1); }
// 	INLINE void accumulate(cdotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cdotprecision,cvector,rmatrix_subv>(dp,cvector(rv1),rv2); }
// 	INLINE void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cidotprecision,cvector,rmatrix_subv>(dp,cvector(rv2),rv1); }
// 	INLINE void accumulate(cidotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2)
// #if(CXSC_INDEX_CHECK)
// 	throw(OP_WITH_WRONG_DIM)
// #else
// 	throw()
// #endif
// 	{ _vmvaccu<cidotprecision,cvector,rmatrix_subv>(dp,cvector(rv1),rv2); }

	INLINE cvector &cvector::operator =(const rmatrix_subv &mv) throw() { return _vmvassign<cvector,rmatrix_subv,complex>(*this,mv); }
	INLINE cvector_slice &cvector_slice::operator =(const rmatrix_subv &mv) throw() { return _vsvassign(*this,rvector(mv)); }
	INLINE cvector &cvector::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<cvector,rmatrix,complex>(*this,m); }
	INLINE cvector &cvector::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<cvector,rmatrix,complex>(*this,rmatrix(m)); }
	INLINE cvector_slice &cvector_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,rvector(m)); }
	INLINE cvector_slice & cvector_slice::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,cvector(rmatrix(m))); }

	INLINE cvector operator *(const rmatrix &m,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcmult<rmatrix,cvector,cvector>(m,v); }
	INLINE cvector operator *(const rmatrix_slice &ms,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcmult<rmatrix_slice,cvector,cvector>(ms,v); }
	INLINE cvector operator *(const cvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmult<cvector,rmatrix,cvector>(v,m); }
	INLINE cvector operator *(const cvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscmult<cvector,rmatrix_slice,cvector>(v,ms); }
	INLINE cvector &operator *=(cvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmultassign<cvector,rmatrix,complex>(v,m); }
	INLINE cvector &operator *=(cvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscmultassign<cvector,rmatrix_slice,complex>(v,ms); }

	INLINE cvector operator *(const cvector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmult<cvector,rmatrix,cvector>(cvector(v),m); }
	INLINE cvector_slice &cvector_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmcmultassign<cvector_slice,rmatrix,complex>(*this,m); }
	
} // namespace cxsc

#endif

