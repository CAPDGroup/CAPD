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

/* CVS $Id: civecimat.inl,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

// Here are definitions for civector x imatrix-Functions
#ifndef _CXSC_CIVECIMAT_INL_INCLUDED
#define _CXSC_CIVECIMAT_INL_INCLUDED

namespace cxsc {

	INLINE civector::civector(const imatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<civector,imatrix,cinterval>(*this,sl); }
	INLINE civector::civector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<civector,imatrix_slice,cinterval>(*this,sl); }
	INLINE civector::civector(const imatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	INLINE civector _civector(const imatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return civector(sl); }
	INLINE civector _civector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return civector(sl); }


	INLINE void SetIm(civector &iv,const imatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvsetim(iv,rv); }
	INLINE void SetRe(civector &iv,const imatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvsetre(iv,rv); }
	INLINE void SetIm(civector_slice &iv,const imatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsetim(iv,ivector(rv)); }
	INLINE void SetRe(civector_slice &iv,const imatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsetre(iv,ivector(rv)); }

	INLINE civector &civector::operator =(const imatrix_subv &mv) throw() { return _vmvassign<civector,imatrix_subv,cinterval>(*this,mv); }
	INLINE civector_slice &civector_slice::operator =(const imatrix_subv &mv) throw() { return _vsvassign(*this,ivector(mv)); }
	INLINE civector &civector::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,imatrix,cinterval>(*this,m); }
	INLINE civector &civector::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,imatrix,cinterval>(*this,imatrix(m)); }
	INLINE civector_slice &civector_slice::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,ivector(m)); }
	INLINE civector_slice & civector_slice::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,civector(imatrix(m))); }

	INLINE civector operator *(const imatrix &m,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<imatrix,civector,civector>(m,v); }
	INLINE civector operator *(const imatrix_slice &ms,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<imatrix_slice,civector,civector>(ms,v); }
	INLINE civector operator *(const civector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,imatrix,civector>(v,m); }
	INLINE civector operator *(const civector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<civector,imatrix_slice,civector>(v,ms); }
	INLINE civector &operator *=(civector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimultassign<civector,imatrix,cinterval>(v,m); }
	INLINE civector &operator *=(civector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimultassign<civector,imatrix_slice,cinterval>(v,ms); }

	INLINE civector operator *(const civector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,imatrix,civector>(civector(v),m); }
	INLINE civector_slice &civector_slice::operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmcimultassign<civector_slice,imatrix,cinterval>(*this,m); }
	
	INLINE civector operator *(const cvector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<cvector,imatrix,civector>(v,m); }
	INLINE civector operator *(const cvector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<cvector,imatrix_slice,civector>(v,ms); }
	INLINE civector operator *(const cvector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,imatrix,civector>(civector(v),m); }
	INLINE civector operator *(const imatrix &m,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<imatrix,cvector,civector>(m,v); }
	INLINE civector operator *(const imatrix_slice &ms,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<imatrix_slice,cvector,civector>(ms,v); }

} // namespace cxsc

#endif

