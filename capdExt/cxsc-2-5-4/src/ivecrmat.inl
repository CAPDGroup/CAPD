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

/* CVS $Id: ivecrmat.inl,v 1.26 2014/01/30 17:23:45 cxsc Exp $ */

// Here are definitions for ivector x rmatrix-Functions
#ifndef _CXSC_IVECRMAT_INL_INCLUDED
#define _CXSC_IVECRMAT_INL_INCLUDED

namespace cxsc {

	INLINE ivector::ivector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<ivector,rmatrix,interval>(*this,sl); }
	INLINE ivector::ivector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<ivector,rmatrix_slice,interval>(*this,sl); }
	INLINE ivector::ivector(const rmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new interval[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::ivector::ivector(const rmatrix &)
	*/
	INLINE ivector _ivector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return ivector(sl); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::ivector::ivector(const rmatrix_slice &sl)
	*/
	INLINE ivector _ivector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return ivector(sl); }


	INLINE void SetInf(ivector &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvsetinf(iv,rv); }
	INLINE void SetSup(ivector &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvsetsup(iv,rv); }
	INLINE void SetInf(ivector_slice &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsetinf(iv,rvector(rv)); }
	INLINE void SetSup(ivector_slice &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsetsup(iv,rvector(rv)); }

	INLINE void UncheckedSetInf(ivector &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvusetinf(iv,rv); }
	INLINE void UncheckedSetSup(ivector &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvusetsup(iv,rv); }
	INLINE void UncheckedSetInf(ivector_slice &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvusetinf(iv,rvector(rv)); }
	INLINE void UncheckedSetSup(ivector_slice &iv,const rmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvusetsup(iv,rvector(rv)); }

	INLINE ivector &ivector::operator =(const rmatrix_subv &mv) throw() { return _vmvassign<ivector,rmatrix_subv,interval>(*this,mv); }
	INLINE ivector_slice &ivector_slice::operator =(const rmatrix_subv &mv) throw() { return _vsvassign(*this,rvector(mv)); }
	INLINE ivector &ivector::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<ivector,rmatrix,interval>(*this,m); }
	INLINE ivector &ivector::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<ivector,rmatrix,interval>(*this,rmatrix(m)); }
	INLINE ivector_slice &ivector_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,rvector(m)); }
	INLINE ivector_slice & ivector_slice::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,ivector(rmatrix(m))); }

	INLINE ivector operator *(const rmatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvimult<rmatrix,ivector,ivector>(m,v); }
	INLINE ivector operator *(const rmatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvimult<rmatrix_slice,ivector,ivector>(ms,v); }
	INLINE ivector operator *(const ivector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimult<ivector,rmatrix,ivector>(v,m); }
	INLINE ivector operator *(const ivector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsimult<ivector,rmatrix_slice,ivector>(v,ms); }
	INLINE ivector &operator *=(ivector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimultassign<ivector,rmatrix,interval>(v,m); }
	INLINE ivector &operator *=(ivector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsimultassign<ivector,rmatrix_slice,interval>(v,ms); }

	INLINE ivector operator *(const ivector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimult<ivector,rmatrix,ivector>(ivector(v),m); }
	INLINE ivector_slice &ivector_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmimultassign<ivector_slice,rmatrix,interval>(*this,m); }

} // namespace cxsc

#endif

