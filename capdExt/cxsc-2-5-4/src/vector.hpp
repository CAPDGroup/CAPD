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

/* CVS $Id: vector.hpp,v 1.24 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_VECTOR_HPP_INCLUDED
#define _CXSC_VECTOR_HPP_INCLUDED

#include "except.hpp"
#include "matrix.hpp" // there are the definitions matrix x vector

namespace cxsc {


	template <class V>
	TINLINE void _vresize(V &rv) throw();

	template <class V,class S>
	TINLINE void _vresize(V &rv, const int &len)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>);
#else
	throw();
#endif

	template <class V,class S>
	TINLINE void _vresize(V &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>);
#else
	throw();
#endif

	template <class V1,class V2,class S>
	TINLINE V1 &_vvassign(V1 &rv1,const V2 &rv2) throw();

	template <class V,class S>
	TINLINE V & _vsassign(V &rv,const S &r) throw();

	template <class VS1,class VS2>
	TINLINE VS1 & _vsvsassign(VS1 &sl1,const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
	
	template <class V,class VS,class S>
	TINLINE V & _vvsassign(V &rv,const VS &sl) throw();

	template <class VS,class V>
	TINLINE VS & _vsvassign(VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS,class S>
	TINLINE VS & _vssassign(VS &sl,const S &r) throw();

	template <class DP,class V1,class V2>
	TINLINE void _vvaccu(DP &dp, const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	template <class DP,class VS,class V>
	TINLINE void _vsvaccu(DP &dp, const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	template <class V,class S,class E>
	TINLINE E _vsdiv(const V &rv, const S &s) throw();

	template <class V,class S>
	TINLINE V &_vsdivassign(V &rv,const S &r) throw();

	template <class VS,class S,class E>
	TINLINE E _vssdiv(const VS &sl, const S &s) throw();

	template <class V,class S,class E>
	TINLINE E _vsmult(const V &rv, const S &s) throw();

	template <class VS,class S,class E>
	TINLINE E _vssmult(const VS &sl, const S &s) throw();

	template <class V1,class V2,class E>
	TINLINE E _vvlmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class VS,class V,class E>
	TINLINE E _vsvlmult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvslmult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvlimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class VS,class V,class E>
	TINLINE E _vsvlimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvslimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class VS,class V,class E>
	TINLINE E _vsvmult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsmult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class VS,class V,class E>
	TINLINE E _vsvimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvcmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class VS,class V,class E>
	TINLINE E _vsvcmult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvscmult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvcimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class VS,class V,class E>
	TINLINE E _vsvcimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V,class S>
	TINLINE V &_vsmultassign(V &rv,const S &r) throw();

	template <class VS,class S>
	TINLINE VS &_vssmultassign(VS &rv,const S &r) throw();

	template <class VS,class S>
	TINLINE VS &_vssdivassign(VS &rv,const S &r) throw();

	template <class V1,class V2,class E>
	TINLINE E _vvplus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS,class E>
	TINLINE E _vvsplus(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsplus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsminus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvplusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvsplusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvplusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsplusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsminusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvminusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvsminusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvminusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class V>
	TINLINE V _vminus(const V &rv) throw();

	template <class VS,class V>
	TINLINE V _vsminus(const VS &sl) throw();

	template <class V1,class V2,class E>
	TINLINE E _vvminus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class V,class VS,class E>
	TINLINE E _vvsminus(const V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class VS,class V,class E>
	TINLINE E _vsvminus(const VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvconv(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class V,class VS,class E>
	TINLINE E _vvsconv(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsconv(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvconvassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvsconvassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvconvassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsconvassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2,class E>
	TINLINE E _vvsect(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS,class E>
	TINLINE E _vvssect(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class VS1,class VS2,class E>
	TINLINE E _vsvssect(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvsectassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvssectassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvsectassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssectassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE bool _vveq(const V1 &rv1, const V2 &rv2) throw();

	template <class VS,class V>
	TINLINE bool _vsveq(const VS &sl, const V &rv) throw();

	template <class V1,class V2>
	TINLINE bool _vvneq(const V1 &rv1, const V2 &rv2) throw();

	template <class VS,class V>
	TINLINE bool _vsvneq(const VS &sl, const V &rv) throw();

	template <class V1,class V2>
	TINLINE bool _vvless(const V1 &rv1, const V2 &rv2) throw();

	template <class VS,class V>
	TINLINE bool _vsvless(const VS &sl, const V &rv) throw();

	template <class V1,class V2>
	TINLINE bool _vvleq(const V1 &rv1, const V2 &rv2) throw();

	template <class VS,class V>
	TINLINE bool _vsvleq(const VS &sl, const V &rv) throw();

	template <class V,class VS>
	TINLINE bool _vvsless(const V &rv, const VS &sl) throw();

	template <class V,class VS>
	TINLINE bool _vvsleq(const V &rv, const VS &sl) throw();

	template <class V>
	TINLINE bool _vnot(const V &rv) throw();

	template <class V>
	TINLINE void *_vvoid(const V &rv) throw();

	template <class V>
	TINLINE V _vconj(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsconj(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vabs(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsabs(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vdiam(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsdiam(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vmid(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsmid(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vinf(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsinf(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vsup(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vssup(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vre(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsre(const VS &sl) throw();

	template <class V,class E>
	TINLINE E _vim(const V &rv) throw();

	template <class VS,class E>
	TINLINE E _vsim(const VS &sl) throw();

	template <class V,class S>
	TINLINE V &_vsusetsup(V &v, const S &s) throw();

	template <class V,class S>
	TINLINE V &_vsusetinf(V &v, const S &s) throw();

	template <class V,class S>
	TINLINE V &_vssetinf(V &v, const S &s) throw();

	template <class V,class S>
	TINLINE V &_vssetsup(V &v, const S &s) throw();

	template <class V,class S>
	TINLINE V &_vssetre(V &v, const S &s) throw();

	template <class V,class S>
	TINLINE V &_vssetim(V &v, const S &s) throw();

	template <class VS,class S>
	TINLINE VS &_vssusetsup(VS &vs, const S &s) throw();

	template <class VS,class S>
	TINLINE VS &_vssusetinf(VS &vs, const S &s) throw();

	template <class VS,class S>
	TINLINE VS &_vsssetinf(VS &vs, const S &s) throw();

	template <class VS,class S>
	TINLINE VS &_vsssetsup(VS &vs, const S &s) throw();

	template <class VS,class S>
	TINLINE VS &_vsssetre(VS &vs, const S &s) throw();

	template <class VS,class S>
	TINLINE VS &_vsssetim(VS &vs, const S &s) throw();

	template <class V1,class V2>
	TINLINE V1 &_vvsetinf(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvssetinf(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvsetinf(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetinf(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvsetsup(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvssetsup(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvsetsup(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetsup(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvusetinf(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvsusetinf(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvusetinf(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsusetinf(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvusetsup(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvsusetsup(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvusetsup(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsusetsup(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvsetim(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvssetim(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvsetim(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetim(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class V1,class V2>
	TINLINE V1 &_vvsetre(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif

	template <class V,class VS>
	TINLINE V &_vvssetre(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class VS,class V>
	TINLINE VS &_vsvsetre(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetre(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	template <class DP,class VS1,class VS2>
	TINLINE void _vsvsaccu(DP &dp, const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	template <class VS1,class VS2>
	TINLINE bool _vsvseq(const VS1 &sl1, const VS2 &sl2) throw();

	template <class VS1,class VS2>
	TINLINE bool _vsvsneq(const VS1 &sl1, const VS2 &sl2) throw();

	template <class VS1,class VS2>
	TINLINE bool _vsvsless(const VS1 &sl1, const VS2 &sl2) throw();

	template <class VS1,class VS2>
	TINLINE bool _vsvsleq(const VS1 &sl1, const VS2 &sl2) throw();

	template <class VS>
	TINLINE bool _vsnot(const VS &sl) throw();

	template <class VS>
	TINLINE void *_vsvoid(const VS &sl) throw();

	template <class V>
	std::ostream &_vout(std::ostream &s, const V &rv) throw();

	template <class V>
	std::istream &_vin(std::istream &s, V &rv) throw();

	template <class V>
	std::ostream &_vsout(std::ostream &s, const V &rv) throw();

	template <class V>
	std::istream &_vsin(std::istream &s, V &rv) throw();

} // namespace cxsc 

#endif

