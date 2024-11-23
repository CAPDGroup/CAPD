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

/* CVS $Id: civeccmat.inl,v 1.24 2014/01/30 17:23:44 cxsc Exp $ */

// Here are definitions for civector x cmatrix-Functions
#ifndef _CXSC_CIVECCMAT_INL_INCLUDED
#define _CXSC_CIVECCMAT_INL_INCLUDED

namespace cxsc {

	INLINE civector::civector(const cmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<civector,cmatrix,cinterval>(*this,sl); }
	INLINE civector::civector(const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<civector,cmatrix_slice,cinterval>(*this,sl); }
	INLINE civector::civector(const cmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
			dat[i]=v.dat[j];
	}
	INLINE civector _civector(const cmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return civector(sl); }
	INLINE civector _civector(const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return civector(sl); }

	INLINE void SetInf(civector &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvsetinf(iv,rv); }
	INLINE void SetSup(civector &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvsetsup(iv,rv); }
	INLINE void SetInf(civector_slice &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsetinf(iv,cvector(rv)); }
	INLINE void SetSup(civector_slice &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvsetsup(iv,cvector(rv)); }

	INLINE void UncheckedSetInf(civector &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvusetinf(iv,rv); }
	INLINE void UncheckedSetSup(civector &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvusetsup(iv,rv); }
	INLINE void UncheckedSetInf(civector_slice &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvusetinf(iv,cvector(rv)); }
	INLINE void UncheckedSetSup(civector_slice &iv,const cmatrix_subv &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vsvusetsup(iv,cvector(rv)); }

	INLINE civector &civector::operator =(const cmatrix_subv &mv) throw() { return _vmvassign<civector,cmatrix_subv,cinterval>(*this,mv); }
	INLINE civector_slice &civector_slice::operator =(const cmatrix_subv &mv) throw() { return _vsvassign(*this,cvector(mv)); }
	INLINE civector &civector::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,cmatrix,cinterval>(*this,m); }
	INLINE civector &civector::operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,cmatrix,cinterval>(*this,cmatrix(m)); }
	INLINE civector_slice &civector_slice::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,cvector(m)); }
	INLINE civector_slice & civector_slice::operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,civector(cmatrix(m))); }

	INLINE civector operator *(const cmatrix &m,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<cmatrix,civector,civector>(m,v); }
	INLINE civector operator *(const cmatrix_slice &ms,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<cmatrix_slice,civector,civector>(ms,v); }
	INLINE civector operator *(const civector &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cmatrix,civector>(v,m); }
	INLINE civector operator *(const civector &v,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<civector,cmatrix_slice,civector>(v,ms); }
	INLINE civector &operator *=(civector &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimultassign<civector,cmatrix,cinterval>(v,m); }
	INLINE civector &operator *=(civector &v,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimultassign<civector,cmatrix_slice,cinterval>(v,ms); }

	INLINE civector operator *(const civector_slice &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cmatrix,civector>(civector(v),m); }
	INLINE civector_slice &civector_slice::operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmcimultassign<civector_slice,cmatrix,cinterval>(*this,m); }
	
	INLINE civector operator *(const ivector &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<ivector,cmatrix,civector>(v,m); }
	INLINE civector operator *(const ivector &v,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<ivector,cmatrix_slice,civector>(v,ms); }
	INLINE civector operator *(const ivector_slice &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cmatrix,civector>(civector(v),m); }
	INLINE civector operator *(const cmatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<cmatrix,ivector,civector>(m,v); }
	INLINE civector operator *(const cmatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<cmatrix_slice,ivector,civector>(ms,v); }

} // namespace cxsc

#endif

