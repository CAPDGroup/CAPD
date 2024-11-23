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

/* CVS $Id: livecimat.inl,v 1.25 2014/01/30 17:23:47 cxsc Exp $ */

// Here are definitions for l_ivector x imatrix-Functions
#ifndef _CXSC_LIVECIMAT_INL_INCLUDED
#define _CXSC_LIVECIMAT_INL_INCLUDED

namespace cxsc {

	INLINE l_ivector::l_ivector(const imatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<l_ivector,imatrix,l_interval>(*this,sl); }
	INLINE l_ivector::l_ivector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<l_ivector,imatrix_slice,l_interval>(*this,sl); }
	INLINE l_ivector::l_ivector(const imatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new l_interval[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const imatrix &)
	*/
	INLINE l_ivector _l_ivector(const imatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return l_ivector(sl); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const imatrix_slice &sl)
	*/
	INLINE l_ivector _l_ivector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return l_ivector(sl); }

	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_ivector,imatrix_subv>(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_ivector,imatrix_subv>(dp,rv1,rv2); }
	
	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_ivector,imatrix_subv>(dp,l_ivector(rv2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_ivector,imatrix_subv>(dp,l_ivector(rv1),rv2); }

	INLINE l_ivector &l_ivector::operator =(const imatrix_subv &mv) throw() { return _vmvassign<l_ivector,imatrix_subv,l_interval>(*this,mv); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const imatrix_subv &mv) throw() { return _vsvassign(*this,ivector(mv)); }
	INLINE l_ivector &l_ivector::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_ivector,imatrix,l_interval>(*this,m); }
	INLINE l_ivector &l_ivector::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_ivector,imatrix,l_interval>(*this,imatrix(m)); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,ivector(m)); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,l_ivector(imatrix(m))); }

	INLINE l_ivector operator *(const imatrix &m,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlimult<imatrix,l_ivector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const imatrix_slice &ms,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlimult<imatrix_slice,l_ivector,l_ivector>(ms,v); }
	INLINE l_ivector operator *(const l_ivector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,imatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const l_ivector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<l_ivector,imatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector &operator *=(l_ivector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimultassign<l_ivector,imatrix,l_interval>(v,m); }
	INLINE l_ivector &operator *=(l_ivector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimultassign<l_ivector,imatrix_slice,l_interval>(v,ms); }

	INLINE l_ivector operator *(const l_ivector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,imatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmlimultassign<l_ivector_slice,imatrix,l_interval>(*this,m); }
	
	INLINE l_ivector operator *(const l_rvector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_rvector,imatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const l_rvector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<l_rvector,imatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector operator *(const l_rvector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,imatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector operator *(const imatrix &m,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlimult<imatrix,l_rvector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const imatrix_slice &ms,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlimult<imatrix_slice,l_rvector,l_ivector>(ms,v); }

} // namespace cxsc

#endif

