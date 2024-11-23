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

/* CVS $Id: iveccvec.inl,v 1.25 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_IVECCVEC_INL_INCLUDED
#define _CXSC_IVECCVEC_INL_INCLUDED

#include "cinterval.hpp"

namespace cxsc {


	INLINE cinterval operator *(const cvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<cvector,ivector,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<cvector_slice,ivector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<ivector_slice,cvector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const cvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<cvector_slice,ivector_slice,cinterval>(sl1,sl2); }
	
	INLINE cinterval operator *(const ivector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<cvector,ivector,cinterval>(rv2,rv1); }
	INLINE cinterval operator *(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<ivector_slice,cvector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<cvector_slice,ivector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const ivector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<cvector_slice,ivector_slice,cinterval>(sl2,sl1); }
	
} // namespace cxsc

#endif

