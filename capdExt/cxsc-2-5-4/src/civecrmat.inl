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

/* CVS $Id: civecrmat.inl,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

// Here are definitions for civector x rmatrix-Functions
#ifndef _CXSC_CIVECRMAT_INL_INCLUDED
#define _CXSC_CIVECRMAT_INL_INCLUDED

namespace cxsc {

	INLINE civector::civector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<civector,rmatrix,cinterval>(*this,sl); }
	INLINE civector::civector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<civector,rmatrix_slice,cinterval>(*this,sl); }
	INLINE civector::civector(const rmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	INLINE civector _civector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return civector(sl); }
	INLINE civector _civector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return civector(sl); }

	INLINE civector &civector::operator =(const rmatrix_subv &mv) throw() { return _vmvassign<civector,rmatrix_subv,cinterval>(*this,mv); }
	INLINE civector_slice &civector_slice::operator =(const rmatrix_subv &mv) throw() { return _vsvassign(*this,rvector(mv)); }
	INLINE civector &civector::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,rmatrix,cinterval>(*this,m); }
	INLINE civector &civector::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,rmatrix,cinterval>(*this,rmatrix(m)); }
	INLINE civector_slice &civector_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,rvector(m)); }
	INLINE civector_slice & civector_slice::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,civector(rmatrix(m))); }

	INLINE civector operator *(const rmatrix &m,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<rmatrix,civector,civector>(m,v); }
	INLINE civector operator *(const rmatrix_slice &ms,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<rmatrix_slice,civector,civector>(ms,v); }
	INLINE civector operator *(const civector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,rmatrix,civector>(v,m); }
	INLINE civector operator *(const civector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<civector,rmatrix_slice,civector>(v,ms); }
	INLINE civector &operator *=(civector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimultassign<civector,rmatrix,cinterval>(v,m); }
	INLINE civector &operator *=(civector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimultassign<civector,rmatrix_slice,cinterval>(v,ms); }

	INLINE civector operator *(const civector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,rmatrix,civector>(civector(v),m); }
	INLINE civector_slice &civector_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmcimultassign<civector_slice,rmatrix,cinterval>(*this,m); }

} // namespace cxsc

#endif

