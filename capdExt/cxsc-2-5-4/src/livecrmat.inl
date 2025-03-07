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

/* CVS $Id: livecrmat.inl,v 1.25 2014/01/30 17:23:47 cxsc Exp $ */

// Here are definitions for l_ivector x rmatrix-Functions
#ifndef _CXSC_LIVECRMAT_INL_INCLUDED
#define _CXSC_LIVECRMAT_INL_INCLUDED

namespace cxsc {

	INLINE l_ivector::l_ivector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vmconstr<l_ivector,rmatrix,l_interval>(*this,sl); }
	INLINE l_ivector::l_ivector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vmsconstr<l_ivector,rmatrix_slice,l_interval>(*this,sl); }
	INLINE l_ivector::l_ivector(const rmatrix_subv &v):l(v.lb),u(v.ub),size(v.size)
	{
		dat=new l_interval[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const rmatrix &)
	*/
	INLINE l_ivector _l_ivector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return l_ivector(sl); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_ivector::l_ivector(const rmatrix_slice &sl)
	*/
	INLINE l_ivector _l_ivector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return l_ivector(sl); }

	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vmvaccu<idotprecision,l_ivector,rmatrix_subv>(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vmvaccu<idotprecision,l_ivector,rmatrix_subv>(dp,rv1,rv2); }
	
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vmvaccu<idotprecision,l_ivector,rmatrix_subv>(dp,l_ivector(rv2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ _vmvaccu<idotprecision,l_ivector,rmatrix_subv>(dp,l_ivector(rv1),rv2); }

	INLINE l_ivector &l_ivector::operator =(const rmatrix_subv &mv) { return _vmvassign<l_ivector,rmatrix_subv,l_interval>(*this,mv); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const rmatrix_subv &mv) { return _vsvassign(*this,rvector(mv)); }
	INLINE l_ivector &l_ivector::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmassign<l_ivector,rmatrix,l_interval>(*this,m); }
	INLINE l_ivector &l_ivector::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmassign<l_ivector,rmatrix,l_interval>(*this,rmatrix(m)); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvassign(*this,rvector(m)); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsvassign(*this,l_ivector(rmatrix(m))); }

	INLINE l_ivector operator *(const rmatrix &m,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _mvlimult<rmatrix,l_ivector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const rmatrix_slice &ms,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _msvlimult<rmatrix_slice,l_ivector,l_ivector>(ms,v); }
	INLINE l_ivector operator *(const l_ivector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmlimult<l_ivector,rmatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const l_ivector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmslimult<l_ivector,rmatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector &operator *=(l_ivector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmlimultassign<l_ivector,rmatrix,l_interval>(v,m); }
	INLINE l_ivector &operator *=(l_ivector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmslimultassign<l_ivector,rmatrix_slice,l_interval>(v,ms); }

	INLINE l_ivector operator *(const l_ivector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vmlimult<l_ivector,rmatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	
#else
	
#endif
	{ return _vsmlimultassign<l_ivector_slice,rmatrix,l_interval>(*this,m); }

} // namespace cxsc

#endif

