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

/* CVS $Id: lrvecrmat.inl,v 1.25 2014/01/30 17:23:47 cxsc Exp $ */

// Here are definitions for l_rvector x rmatrix-Functions
#ifndef _CXSC_LRVECRMAT_INL_INCLUDED
#define _CXSC_LRVECRMAT_INL_INCLUDED

namespace cxsc {

	INLINE l_rvector::l_rvector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<l_rvector,rmatrix,l_real>(*this,sl); }
	INLINE l_rvector::l_rvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<l_rvector,rmatrix_slice,l_real>(*this,sl); }
	INLINE l_rvector::l_rvector(const rmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new l_real[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const l_rmatrix &)
	*/
	INLINE l_rvector _l_rvector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return l_rvector(sl); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rvector::l_rvector(const l_rmatrix_slice &sl)
	*/
	INLINE l_rvector _l_rvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return l_rvector(sl); }

	INLINE void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<dotprecision,l_rvector,rmatrix_subv>(dp,rv2,rv1); }
	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<dotprecision,l_rvector,rmatrix_subv>(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,rmatrix_subv>(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,rmatrix_subv>(dp,rv1,rv2); }
	
	INLINE void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const l_rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<dotprecision,l_rvector,rmatrix_subv>(dp,l_rvector(rv2),rv1); }
	INLINE void accumulate(dotprecision &dp, const l_rvector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<dotprecision,l_rvector,rmatrix_subv>(dp,l_rvector(rv1),rv2); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,rmatrix_subv>(dp,l_rvector(rv2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu<idotprecision,l_rvector,rmatrix_subv>(dp,l_rvector(rv1),rv2); }

	INLINE l_rvector &l_rvector::operator =(const rmatrix_subv &mv) throw() { return _vmvassign<l_rvector,rmatrix_subv,l_real>(*this,mv); }
	INLINE l_rvector_slice &l_rvector_slice::operator =(const rmatrix_subv &mv) throw() { return _vsvassign(*this,rvector(mv)); }
	INLINE l_rvector &l_rvector::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_rvector,rmatrix,l_real>(*this,m); }
	INLINE l_rvector &l_rvector::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_rvector,rmatrix,l_real>(*this,rmatrix(m)); }
	INLINE l_rvector_slice &l_rvector_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,rvector(m)); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,l_rvector(rmatrix(m))); }

	INLINE l_rvector operator *(const rmatrix &m,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlmult<rmatrix,l_rvector,l_rvector>(m,v); }
	INLINE l_rvector operator *(const rmatrix_slice &ms,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlmult<rmatrix_slice,l_rvector,l_rvector>(ms,v); }
	INLINE l_rvector operator *(const l_rvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmult<l_rvector,rmatrix,l_rvector>(v,m); }
	INLINE l_rvector operator *(const l_rvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslmult<l_rvector,rmatrix_slice,l_rvector>(v,ms); }
	INLINE l_rvector &operator *=(l_rvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmultassign<l_rvector,rmatrix,l_real>(v,m); }
	INLINE l_rvector &operator *=(l_rvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslmultassign<l_rvector,rmatrix_slice,l_real>(v,ms); }

	INLINE l_rvector operator *(const l_rvector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmult<l_rvector,rmatrix,l_rvector>(l_rvector(v),m); }
	INLINE l_rvector_slice &l_rvector_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ *this=*this*m; return *this; }
	
} // namespace cxsc

#endif

