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

/* CVS $Id: civector.hpp,v 1.39 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CIVECTOR_HPP_INCLUDED
#define _CXSC_CIVECTOR_HPP_INCLUDED

#include "xscclass.hpp"
#include "except.hpp"
#include "cidot.hpp"
#include "cinterval.hpp" // used for declaration of Inf, Sup,...
//#include "cxscmatr.hpp"
#include "rvector.hpp"
#include "ivector.hpp"
#include "cvector.hpp"
#include "vector.hpp"


#include <iostream>

//#include "matrix.hpp" // hat hier eigentlich nichts zu suchen, sonst aber Internal Compiler Error #9

namespace cxsc {

class civector_slice;
class scivector;
class scivector_slice;

//! The Data Type civector
/*!
The vectors of C-XSC are one dimensional arrays of the corresponding scalar base type. 

\sa rvector
*/
class civector
{
	friend class civector_slice;
	friend class cimatrix;
	friend class cimatrix_subv;
	private:
	cinterval *dat;
	int l,u,size;

	public:
//#if(CXSC_INDEX_CHECK)
#ifdef _CXSC_FRIEND_TPL
	//------------ Templates --------------------------------------------------
	// cinterval
template <class V,class MS,class S> friend  void _vmsconstr(V &v,const MS &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>);
#else
	throw();
#endif
template <class V,class M,class S> friend  void _vmconstr(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
 template <class V> friend 	 void _vresize(V &rv) throw();
 template <class V,class S> friend 	 void _vresize(V &rv, const int &len)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>);
#else
	throw();
#endif
 template <class V,class S> friend 	 void _vresize(V &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>);
#else
	throw();
#endif
 template <class V1,class V2,class S> friend 	 V1 &_vvassign(V1 &rv1,const V2 &rv2) throw();
 template <class V,class S> friend 	 V & _vsassign(V &rv,const S &r) throw();
 template <class V,class VS,class S> friend 	 V & _vvsassign(V &rv,const VS &sl) throw();
 template <class VS,class V> friend 	 VS & _vsvassign(VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
template <class V,class M,class S> friend  V &_vmassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
template <class M,class V,class S> friend  M &_mvassign(M &m,const V &v) throw();
 template <class V1,class V2> friend 	 V1 &_vvsetinf(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvsetsup(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssetinf(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssetsup(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvsetinf(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvsetsup(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvusetinf(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvusetsup(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsusetinf(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsusetsup(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvusetinf(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvusetsup(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V,class S> friend 	 V &_vssetinf(V &v, const S &s) throw();
 template <class V,class S> friend 	 V &_vssetsup(V &v, const S &s) throw();
 template <class V,class S> friend 	 V &_vsusetinf(V &v, const S &s) throw();
 template <class V,class S> friend 	 V &_vsusetsup(V &v, const S &s) throw();
 template <class V,class S> friend 	 V &_vssetim(V &v, const S &s) throw();
 template <class V,class S> friend 	 V &_vssetre(V &v, const S &s) throw();
 template <class V> friend 	 V _vconj(const V &rv) throw();
 template <class VS,class E> friend 	 E _vsconj(const VS &sl) throw();
 template <class V,class E> friend 	 E _vabs(const V &rv) throw();
 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();
template <class MV,class V> friend  V _mvabs(const MV &mv) throw();
 template <class V,class E> friend 	 E _vdiam(const V &rv) throw();
 template <class V,class E> friend 	 E _vmid(const V &rv) throw();
 template <class V,class E> friend 	 E _vinf(const V &rv) throw();
 template <class V,class E> friend 	 E _vsup(const V &rv) throw();
 template <class V,class E> friend 	 E _vim(const V &rv) throw();
 template <class V,class E> friend 	 E _vre(const V &rv) throw();
	friend INLINE ivector Re(const civector &v) throw();
	friend INLINE ivector Im(const civector &v) throw();
	friend INLINE cvector Inf(const civector &v) throw();
	friend INLINE cvector Sup(const civector &v) throw();
	friend INLINE rvector SupRe(const civector &v) throw();
	friend INLINE rvector SupIm(const civector &v) throw();
	friend INLINE rvector InfRe(const civector &v) throw();
	friend INLINE rvector InfIm(const civector &v) throw();
 template <class V1,class V2> friend 	 V1 &_vvsetim(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvsetre(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssetim(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssetre(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvsetim(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvsetre(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

//-------- vector-vector -----------------------
 template <class DP,class V1,class V2> friend 	 void _vvaccu(DP &dp, const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
 template <class DP,class VS,class V> friend 	 void _vsvaccu(DP &dp, const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

 template <class V1,class V2,class E> friend 	 E _vvcimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class VS,class V,class E> friend 	 E _vsvcimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
	
 template <class V,class S> friend 	 V &_vsmultassign(V &rv,const S &r) throw();
 template <class V1,class V2,class E> friend 	 E _vvplus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS,class E> friend 	 E _vvsplus(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvsplus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvsminus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvplusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsplusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvplusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsplusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsminusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvminusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsminusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvminusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class V> friend 	 V _vminus(const V &rv) throw();
 template <class VS,class V> friend 	 V _vsminus(const VS &sl) throw();
 template <class V1,class V2,class E> friend 	 E _vvminus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V,class VS,class E> friend 	 E _vvsminus(const V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class VS,class V,class E> friend 	 E _vsvminus(const VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvconv(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V,class VS,class E> friend 	 E _vvsconv(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvsconv(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvconvassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsconvassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvconvassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsconvassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvsect(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS,class E> friend 	 E _vvssect(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvssect(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvsectassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssectassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvsectassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssectassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class MV1,class MV2,class E> friend 	 E _mvmvsect(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MV,class V,class E> friend 	 E _mvvsect(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvsectassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvsectassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class MV1,class MV2,class E> friend 	 E _mvmvconv(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MV,class V,class E> friend 	 E _mvvconv(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvconvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvconvassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V,class MV,class S> friend 	 S _vmvcimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
	//--------- vector-scalar -----------------
 template <class V,class S,class E> friend 	 E _vsdiv(const V &rv, const S &s) throw();
 template <class V,class S> friend 	 V &_vsdivassign(V &rv,const S &r) throw();
 template <class VS,class S,class E> friend 	 E _vssdiv(const VS &sl, const S &s) throw();
 template <class V,class S,class E> friend 	 E _vsmult(const V &rv, const S &s) throw();
 template <class VS,class S,class E> friend 	 E _vssmult(const VS &sl, const S &s) throw();
 template <class MV,class S,class E> friend 	 E _mvsmult(const MV &rv, const S &s) throw();
 template <class MV1,class MV2,class E> friend 	 E _mvmvplus(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MV,class V,class E> friend 	 E _mvvplus(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MV,class V,class E> friend 	 E _mvvminus(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V,class MV,class E> friend 	 E _vmvminus(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MV1,class MV2,class E> friend 	 E _mvmvminus(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvplusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvminusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
 template <class MV,class S,class E> friend 	 E _mvsdiv(const MV &rv, const S &s) throw();
template <class MV,class V> friend  MV &_mvvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

 template <class V1,class V2> friend 	 bool _vveq(const V1 &rv1, const V2 &rv2) throw();
 template <class VS,class V> friend 	 bool _vsveq(const VS &sl, const V &rv) throw();
 template <class V1,class V2> friend 	 bool _vvneq(const V1 &rv1, const V2 &rv2) throw();
 template <class VS,class V> friend 	 bool _vsvneq(const VS &sl, const V &rv) throw();
 template <class V1,class V2> friend 	 bool _vvless(const V1 &rv1, const V2 &rv2) throw();
 template <class VS,class V> friend 	 bool _vsvless(const VS &sl, const V &rv) throw();
 template <class V1,class V2> friend 	 bool _vvleq(const V1 &rv1, const V2 &rv2) throw();
 template <class VS,class V> friend 	 bool _vsvleq(const VS &sl, const V &rv) throw();
 template <class V,class VS> friend 	 bool _vvsless(const V &rv, const VS &sl) throw();
 template <class V,class VS> friend 	 bool _vvsleq(const V &rv, const VS &sl) throw();
 template <class V> friend 	 bool _vnot(const V &rv) throw();
 template <class V> friend 	 void *_vvoid(const V &rv) throw();
 template <class VS1,class VS2> friend 	 bool _vsvseq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsneq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsless(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsleq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS> friend 	 bool _vsnot(const VS &sl) throw();
 template <class VS> friend 	 void *_vsvoid(const VS &sl) throw();
 template <class V> friend 	std::ostream &_vout(std::ostream &s, const V &rv) throw();
 template <class V> friend 	std::istream &_vin(std::istream &s, V &rv) throw();

	//------------- vector-matrix ---------------
template <class DP,class V,class SV> friend 	 void _vmvaccu(DP &dp, const V & rv1, const SV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
template <class V,class MV2,class S> friend  V &_vmvassign(V &v,const MV2 &rv) throw();
 template <class M,class V,class E> friend 	 E _mvcimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmcimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class S> friend 	 V &_vmimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class S> friend 	 V &_vmcimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class V,class E> friend 	 E _msvcimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmscimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class S> friend 	 V &_vmscimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	// Real
 template <class DP,class VS1,class VS2> friend 	 void _vsvsaccu(DP &dp, const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//--Real -------- vector-scalar ------------
	//--Real--------- Vector-vector---------
	//-- Real -------- Vector-matrix ----------
	// complex
	//--complex -------- vector-scalar ------------
	//--complex--------- Vector-vector---------
	//-- complex -------- Vector-matrix ----------
	// interval
	//--interval -------- vector-scalar ------------
	//--interval--------- Vector-vector---------
	//-- interval -------- Vector-matrix ----------
	// cvector x ivector ----------------
	// vector - scalar -------
	// vector - vector -------------
	// vector - matrix ------------

/* template<class T1,class T2,class T3>   friend T3 _mvscimult(const T1 &m,const T2 &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif  
template<class T1,class T2,class T3> friend T3 _vsmcimult(const T1 &v,const T2 &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif */
	
	//  real x complex --------------
	// vector - vector --------------

#endif
	

	//------ Konstruktoren ----------------------------------------------------
	//! Constructor of class civector
	INLINE civector () throw();
	//! Constructor of class civector
	explicit INLINE civector(const int &i) throw();
#ifdef OLD_CXSC
	//! Constructor of class civector
	explicit INLINE civector(const class index &i) throw(); // for backwards compatibility
#endif
	//! Constructor of class civector
	explicit INLINE civector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_WRONG_BOUNDARIES,ERROR_CIVECTOR_NO_MORE_MEMORY);
#else
	throw();
#endif
	//! Constructor of class civector
	INLINE civector(const cimatrix_subv &) throw();
	//! Constructor of class civector
	explicit INLINE civector(const cinterval &) throw();
	//! Constructor of class civector
//	explicit INLINE civector(const cimatrix &)
	explicit        civector(const cimatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const cimatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	INLINE civector(const civector_slice &rs) throw();
	//! Constructor of class civector
	INLINE civector(const civector &v) throw();
	//! Constructor of class civector
	INLINE civector(const scivector_slice &rs);
	//! Constructor of class civector
	INLINE civector(const scivector &v);
	// Real
	//! Constructor of class civector
	explicit INLINE civector(const real &) throw();
	//! Constructor of class civector
	explicit INLINE civector(const rvector_slice &rs) throw();
	//! Constructor of class civector
	explicit INLINE civector(const rvector &v) throw();
	//! Constructor of class civector
	explicit INLINE civector(const srvector_slice &rs);
	//! Constructor of class civector
	explicit INLINE civector(const srvector &v);
	//! Constructor of class civector
	explicit INLINE civector(const rmatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const rmatrix_subv &) throw();
	
	// complex
	//! Constructor of class civector
	explicit INLINE civector(const complex &) throw();
	//! Constructor of class civector
	explicit INLINE civector(const cvector_slice &rs) throw();
	//! Constructor of class civector
	explicit INLINE civector(const cvector &v) throw();
	//! Constructor of class civector
	explicit INLINE civector(const scvector_slice &rs);
	//! Constructor of class civector
	explicit INLINE civector(const scvector &v);
	//! Constructor of class civector
	explicit INLINE civector(const cmatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const cmatrix_subv &) throw();
	
	// interval
	//! Constructor of class civector
	explicit INLINE civector(const interval &) throw();
	//! Constructor of class civector
	explicit INLINE civector(const ivector_slice &rs) throw();
	//! Constructor of class civector
	explicit INLINE civector(const ivector &v) throw();
	//! Constructor of class civector
	explicit INLINE civector(const sivector_slice &rs);
	//! Constructor of class civector
	explicit INLINE civector(const sivector &v);
	//! Constructor of class civector
	explicit INLINE civector(const imatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class civector
	explicit INLINE civector(const imatrix_subv &) throw();
	
	// cinterval
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const civector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const civector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const scivector &rv) ;
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const scivector_slice &sl) ;
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cinterval &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cimatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cimatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const rvector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const rvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const srvector &rv);
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const srvector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const rmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const rmatrix_subv &) throw();

	// complex
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cvector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const scvector &rv);
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const scvector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const cmatrix_subv &) throw();

	// interval
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const ivector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const ivector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const sivector &rv) ;
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const sivector_slice &sl) ;
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const imatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector &operator =(const imatrix_subv &) throw();


        civector& operator+=(const srvector&);
        civector& operator+=(const scvector&);
        civector& operator+=(const sivector&);
        civector& operator+=(const scivector&);
        civector& operator-=(const srvector&);
        civector& operator-=(const scvector&);
        civector& operator-=(const sivector&);
        civector& operator-=(const scivector&);
        civector& operator|=(const srvector&);
        civector& operator|=(const scvector&);
        civector& operator|=(const sivector&);
        civector& operator|=(const scivector&);
        civector& operator&=(const sivector&);
        civector& operator&=(const scivector&);
        civector& operator+=(const srvector_slice&);
        civector& operator+=(const scvector_slice&);
        civector& operator+=(const sivector_slice&);
        civector& operator+=(const scivector_slice&);
        civector& operator-=(const srvector_slice&);
        civector& operator-=(const scvector_slice&);
        civector& operator-=(const sivector_slice&);
        civector& operator-=(const scivector_slice&);
        civector& operator|=(const srvector_slice&);
        civector& operator|=(const scvector_slice&);
        civector& operator|=(const sivector_slice&);
        civector& operator|=(const scivector_slice&);
        civector& operator&=(const sivector_slice&);
        civector& operator&=(const scivector_slice&);

        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE civector operator()(const intvector& p);
        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE civector operator()(const intmatrix& P);

	//--------- Destruktor ----------------------------------------------------
	INLINE ~civector() { delete [] dat; }

	//------ Standardfunktionen -----------------------------------------------
	
	friend INLINE cinterval::cinterval(const civector &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const civector &rv) throw() { return rv.l; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const civector &rv) throw() { return rv.u; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const civector &rv) throw() { return rv.size; }
	//! Sets the lower bound of the vector
	friend INLINE civector & SetLb(civector &rv, const int &l) throw() { rv.l=l; rv.u=l+rv.size-1; return rv;}
	//! Sets the upper bound of the vector
	friend INLINE civector & SetUb(civector &rv, const int &u) throw() { rv.u=u; rv.l=u-rv.size+1; return rv;}
	//! Operator for accessing the single elements of the vector (read-only)
	INLINE cinterval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector
	INLINE cinterval & operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the whole vector
	INLINE civector & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE civector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	civector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	INLINE operator void*() throw();
//#else
//#endif
};

//! The Data Type civector_slice
/*!
This data type represents a partial civector.

\sa civector
*/
class civector_slice
{
	friend class civector;
	friend class cimatrix;
	private:
	cinterval *dat;
	int l,u,size;
	int start,end;

	public:
//#if(CXSC_INDEX_CHECK)	
#ifdef _CXSC_FRIEND_TPL
//------------------------- Templates -------------------------------------------
// cinterval / cinterval

 template <class VS1,class VS2> friend 	 VS1 & _vsvsassign(VS1 &sl1,const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V,class VS,class S> friend 	 V & _vvsassign(V &rv,const VS &sl) throw();
 template <class VS,class V> friend 	 VS & _vsvassign(VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class S> friend 	 VS & _vssassign(VS &sl,const S &r) throw();
 template <class VS,class V> friend 	 VS &_vsvsetinf(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvsetsup(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssetinf(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssetsup(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvusetinf(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvusetsup(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsusetinf(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsusetsup(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

 template <class VS,class E> friend 	 E _vsconj(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsdiam(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsmid(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsinf(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vssup(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsim(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsre(const VS &sl) throw();
	friend INLINE ivector Re(const civector_slice &v) throw();
	friend INLINE ivector Im(const civector_slice &v) throw();
	friend INLINE cvector Inf(const civector_slice &v) throw();
	friend INLINE cvector Sup(const civector_slice &v) throw();
	friend INLINE rvector SupRe(const civector_slice &v) throw();
	friend INLINE rvector SupIm(const civector_slice &v) throw();
	friend INLINE rvector InfRe(const civector_slice &v) throw();
	friend INLINE rvector InfIm(const civector_slice &v) throw();
 template <class VS,class V> friend 	 VS &_vsvsetim(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvsetre(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssetim(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssetre(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
  /*	friend TINLINE civector_slice &_vsmvsetim(civector_slice &,const imatrix_subv &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	friend TINLINE civector_slice &_vsmvsetre(civector_slice &,const imatrix_subv &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
	#endif */
 template <class VS,class S> friend 	 VS &_vsssetinf(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vsssetsup(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vssusetinf(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vssusetsup(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vsssetim(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vsssetre(VS &vs, const S &s) throw();

 template <class DP,class VS,class V> friend 	 void _vsvaccu(DP &dp, const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
 template <class DP,class VS1,class VS2> friend 	 void _vsvsaccu(DP &dp, const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
 template <class VS,class S,class E> friend 	 E _vssdiv(const VS &sl, const S &s) throw();
 template <class VS,class S,class E> friend 	 E _vssmult(const VS &sl, const S &s) throw();

 template <class VS,class V,class E> friend 	 E _vsvcimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class S> friend 	 VS &_vssmultassign(VS &rv,const S &r) throw();
 template <class VS,class S> friend 	 VS &_vssdivassign(VS &rv,const S &r) throw();
 template <class V,class VS,class E> friend 	 E _vvsplus(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvsplus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvsminus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsplusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvplusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsplusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsminusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsminusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvminusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 V _vsminus(const VS &sl) throw();
 template <class V,class VS,class E> friend 	 E _vvsminus(const V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class VS,class V,class E> friend 	 E _vsvminus(const VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V,class VS,class E> friend 	 E _vvssect(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvssect(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssectassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvsectassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssectassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V,class VS,class E> friend 	 E _vvsconv(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class VS1,class VS2,class E> friend 	 E _vsvsconv(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvsconvassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvconvassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsconvassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS,class M,class S> friend 	 VS &_vsmcimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

 template <class VS,class V> friend 	 bool _vsveq(const VS &sl, const V &rv) throw();
 template <class VS,class V> friend 	 bool _vsvneq(const VS &sl, const V &rv) throw();
 template <class VS,class V> friend 	 bool _vsvless(const VS &sl, const V &rv) throw();
 template <class VS,class V> friend 	 bool _vsvleq(const VS &sl, const V &rv) throw();
 template <class V,class VS> friend 	 bool _vvsless(const V &rv, const VS &sl) throw();
 template <class V,class VS> friend 	 bool _vvsleq(const V &rv, const VS &sl) throw();
 template <class VS1,class VS2,class E> friend 	 E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 bool _vsvseq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsneq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsless(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsleq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS> friend 	 bool _vsnot(const VS &sl) throw();
 template <class VS> friend 	 void *_vsvoid(const VS &sl) throw();
 template <class V> friend 	std::ostream &_vsout(std::ostream &s, const V &rv) throw();
 template <class V> friend 	std::istream &_vsin(std::istream &s, V &rv) throw();
	
	// cinterval / Real
 template <class V,class MS,class E> friend 	 E _vmscimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
	// cinterval / complex
 template <class DP,class V1,class V2> friend 	 void _vvaccu(DP &dp, const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	// cinterval / interval
	// complex
 template <class V1,class V2,class S> friend 	 V1 &_vvassign(V1 &rv1,const V2 &rv2) throw();
 template <class V,class S> friend 	 V & _vsassign(V &rv,const S &r) throw();

template <class V,class M,class S> friend  V &_vmassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
template <class M,class V,class S> friend  M &_mvassign(M &m,const V &v) throw();
template <class V,class MV2,class S> friend  V &_vmvassign(V &v,const MV2 &rv) throw();

 template <class V1,class V2,class E> friend 	 E _vvconv(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	//--complex -------- vector-scalar ------------
 template <class MV,class S,class E> friend 	 E _mvsmult(const MV &rv, const S &s) throw();
 template <class V,class S,class E> friend 	 E _vsmult(const V &rv, const S &s) throw();
 template <class V,class S,class E> friend 	 E _vsdiv(const V &rv, const S &s) throw();
 template <class V,class S> friend 	 V &_vsdivassign(V &rv,const S &r) throw();
 template <class V,class S> friend 	 V &_vsmultassign(V &rv,const S &r) throw();
	//--complex--------- Vector-vector---------
 template <class V1,class V2,class E> friend 	 E _vvcimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvplus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvplusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvminusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvminus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvconvassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvsect(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class V1,class V2> friend 	 V1 &_vvsectassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
	
	//-- complex -------- Vector-matrix ----------
template <class V,class MS,class S> friend  void _vmsconstr(V &v,const MS &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>);
#else
	throw();
#endif
template <class V,class M,class S> friend  void _vmconstr(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
 template <class M,class V,class E> friend 	 E _mvcimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class V,class E> friend 	 E _msvcimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmcimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class MS,class S> friend 	 V &_vmscimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class M,class S> friend 	 V &_vmcimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	// interval
  /*	friend TINLINE civector &_vsmassign<civector_slice,imatrix,cinterval>(civector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif */

	//--interval -------- vector-scalar ------------
	//--interval--------- Vector-vector---------
	//-- interval -------- Vector-matrix ----------
/*  friend TINLINE civector _mvscimult<imatrix,civector_slice,civector>(const imatrix &m,const civector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif */
  /*   friend TINLINE civector _msvscimult<imatrix_slice,civector_slice,civector>(const imatrix_slice &ms,const civector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
  /*   friend TINLINE civector _vsmcimult<civector_slice,imatrix,civector>(const civector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif */
  /*   friend TINLINE civector _vsmscimult<civector_slice,imatrix_slice,civector>(const civector_slice &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif */
  /*  friend TINLINE civector &_vsmscimultassign<civector_slice,imatrix_slice,cinterval>(civector_slice &v,const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif */


#endif
	
	//--------------------- Konstruktoren -----------------------------------
	//! Constructor of class civector_slice
	explicit INLINE civector_slice(civector &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	//! Constructor of class civector_slice
	explicit INLINE civector_slice(civector_slice &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	public: 
	//! Constructor of class civector_slice
	INLINE civector_slice(const civector_slice &a) throw():dat(a.dat),l(a.l),u(a.u),size(a.size),start(a.start),end(a.end) { }
	public:
	// cinterval
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const scivector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const scivector &sl);

	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cinterval &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>,ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>,ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice &operator =(const cimatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const srvector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const srvector &sl);

	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice &operator =(const rmatrix_subv &mv) throw();

	// complex
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const scvector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const scvector &sl);

	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice &operator =(const cmatrix_subv &mv) throw();

	// interval
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const sivector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const sivector &sl);

	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice & operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE civector_slice &operator =(const imatrix_subv &mv) throw();

	//--------------------- Standardfunktionen ------------------------------

	friend INLINE cinterval::cinterval(const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const civector_slice &sl) throw() { return sl.start; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const civector_slice &sl) throw() { return sl.end; }
	//! Returns the size of the vector
	friend INLINE int VecLen(const civector_slice &sl) throw() { return sl.size; }

	//! Operator for accessing the single elements of the vector (read-only)
	INLINE cinterval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector
	INLINE cinterval & operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the whole vector
	INLINE civector_slice & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE civector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	civector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();
	
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const cinterval &r) throw();
	//! Implementation of division and allocation operation
	INLINE civector_slice &operator /=(const cinterval &r) throw();
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const real &r) throw();
	//! Implementation of division and allocation operation
	INLINE civector_slice &operator /=(const real &r) throw();
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const complex &r) throw();
	//! Implementation of division and allocation operation
	INLINE civector_slice &operator /=(const complex &r) throw();
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const interval &r) throw();
	//! Implementation of division and allocation operation
	INLINE civector_slice &operator /=(const interval &r) throw();
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector_slice &operator +=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector_slice &operator -=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector_slice &operator |=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector_slice &operator &=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector_slice &operator *=(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
//#else
//#endif

        civector_slice& operator+=(const srvector&);
        civector_slice& operator+=(const scvector&);
        civector_slice& operator+=(const sivector&);
        civector_slice& operator+=(const scivector&);
        civector_slice& operator-=(const srvector&);
        civector_slice& operator-=(const scvector&);
        civector_slice& operator-=(const sivector&);
        civector_slice& operator-=(const scivector&);
        civector_slice& operator|=(const srvector&);
        civector_slice& operator|=(const scvector&);
        civector_slice& operator|=(const sivector&);
        civector_slice& operator|=(const scivector&);
        civector_slice& operator&=(const sivector&);
        civector_slice& operator&=(const scivector&);
        civector_slice& operator+=(const srvector_slice&);
        civector_slice& operator+=(const scvector_slice&);
        civector_slice& operator+=(const sivector_slice&);
        civector_slice& operator+=(const scivector_slice&);
        civector_slice& operator-=(const srvector_slice&);
        civector_slice& operator-=(const scvector_slice&);
        civector_slice& operator-=(const sivector_slice&);
        civector_slice& operator-=(const scivector_slice&);
        civector_slice& operator|=(const srvector_slice&);
        civector_slice& operator|=(const scvector_slice&);
        civector_slice& operator|=(const sivector_slice&);
        civector_slice& operator|=(const scivector_slice&);
        civector_slice& operator&=(const sivector_slice&);
        civector_slice& operator&=(const scivector_slice&);

};

//=======================================================================
//======================== Vector Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE civector _civector(const cinterval &r) throw();
//	INLINE civector _civector(const cimatrix &m) throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE civector _civector(const cimatrix_slice &sl) throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE civector _civector(const real &r) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE civector _civector(const rvector_slice &rs) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE civector _civector(const rvector &rs) throw();
//	INLINE civector _civector(const rmatrix &m) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE civector _civector(const rmatrix_slice &sl) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE civector _civector(const rmatrix_subv &rs) throw();

	//! Returns the vector with the new given infimum vector
	INLINE civector &SetInf(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE civector_slice &SetInf(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE civector &SetInf(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE civector_slice &SetInf(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE civector &UncheckedSetInf(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE civector_slice &UncheckedSetInf(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE civector &UncheckedSetInf(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE civector_slice &UncheckedSetInf(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given supremum vector
	INLINE civector &SetSup(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE civector_slice &SetSup(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE civector &SetSup(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE civector_slice &SetSup(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE civector &UncheckedSetSup(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE civector_slice &UncheckedSetSup(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE civector &UncheckedSetSup(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE civector_slice &UncheckedSetSup(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given real part vector
	INLINE civector &SetRe(civector &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given real part vector
	INLINE civector_slice &SetRe(civector_slice &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given real part vector
	INLINE civector &SetRe(civector &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given real part vector
	INLINE civector_slice &SetRe(civector_slice &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given imaginary part vector
	INLINE civector &SetIm(civector &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given imaginary part vector
	INLINE civector_slice &SetIm(civector_slice &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given imaginary part vector
	INLINE civector &SetIm(civector &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given imaginary part vector
	INLINE civector_slice &SetIm(civector_slice &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given supremum value
	INLINE civector &SetSup(civector &iv,const complex &r) throw();
	//! Returns the vector with the new given infimum value
	INLINE civector &SetInf(civector &iv,const complex &r) throw();
	//! Returns the vector with the new unchecked given supremum value
	INLINE civector &UncheckedSetSup(civector &iv,const complex &r) throw();
	//! Returns the vector with the new unchecked given infimum value
	INLINE civector &SetUncheckedInf(civector &iv,const complex &r) throw();
	//! Sets componentwise the real parts of the civector
	INLINE civector &SetRe(civector &iv,const interval &r) throw();
	//! Sets componentwise the imaginary parts of the civector
	INLINE civector &SetIm(civector &iv,const interval &r) throw();

	//! Returns the vector with the new given supremum value
	INLINE civector_slice &SetSup(civector_slice &iv,const complex &r) throw();
	//! Returns the vector with the new given infimum value
	INLINE civector_slice &SetInf(civector_slice &iv,const complex &r) throw();
	//! Returns the vector with the new unchecked given supremum value
	INLINE civector_slice &UncheckedSetSup(civector_slice &iv,const complex &r) throw();
	//! Returns the vector with the new unchecked given infimum value
	INLINE civector_slice &SetUncheckedInf(civector_slice &iv,const complex &r) throw();
	//! Sets componentwise the real parts of the civector
	INLINE civector_slice &SetRe(civector_slice &iv,const interval &r) throw();
	//! Sets componentwise the imaginary parts of the civector
	INLINE civector_slice &SetIm(civector_slice &iv,const interval &r) throw();

	//! Resizes the vector
	INLINE void Resize(civector &rv) throw();
	//! Resizes the vector
	INLINE void Resize(civector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<civector>);
#else
	throw();
#endif
	//! Resizes the vector
	INLINE void Resize(civector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<civector>);
#else
	throw();
#endif
	
	//! Returns the conjugated civector
	INLINE civector conj(const civector &rv) throw();
	//! Returns the conjugated civector
	INLINE civector conj(const civector_slice &sl) throw();
	
	//! Returns the absolute value of the vector
	INLINE ivector abs(const civector &rv) throw();
	//! Returns the absolute value of the vector
	INLINE ivector abs(const civector_slice &sl) throw();
	//! Returns the diameter of the vector
	INLINE cvector diam(const civector &v) throw();
	//! Returns the diameter of the vector
	INLINE cvector diam(const civector_slice &v) throw();
	//! Returns the middle of the vector
	INLINE cvector mid(const civector &v) throw();
	//! Returns the middle of the vector
	INLINE cvector mid(const civector_slice &v) throw();
	//! Returns the infimum of the vector
	INLINE cvector Inf(const civector &v) throw();
	//! Returns the infimum of the vector
	INLINE cvector Inf(const civector_slice &v) throw();
	//! Returns the supremum of the vector
	INLINE cvector Sup(const civector &v) throw();
	//! Returns the supremum of the vector
	INLINE cvector Sup(const civector_slice &v) throw();
	//! Returns the supremum of real part of the vector
	INLINE rvector SupRe(const civector &v) throw();
	//! Returns the supremum of imaginary part of the vector
	INLINE rvector SupIm(const civector &v) throw();
	//! Returns the infimum of real part of the vector
	INLINE rvector InfRe(const civector &v) throw();
	//! Returns the infimum of imaginary part of the vector
	INLINE rvector InfIm(const civector &v) throw();
	//! Returns the supremum of real part of the vector
	INLINE rvector SupRe(const civector_slice &v) throw();
	//! Returns the supremum of imaginary part of the vector
	INLINE rvector SupIm(const civector_slice &v) throw();
	//! Returns the infimum of real part of the vector
	INLINE rvector InfRe(const civector_slice &v) throw();
	//! Returns the infimum of imaginary part of the vector
	INLINE rvector InfIm(const civector_slice &v) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const civector &rv) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const civector_slice &sl) throw();

//======================= Vector / Scalar ===============================

//----------------------------- cinterval ---------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const civector &rv, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector_slice &sl, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const civector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const civector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE civector &operator *=(civector &rv,const cinterval &r) throw();

	//! Implementation of division operation
	INLINE civector operator /(const civector &rv, const cinterval &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const civector_slice &sl, const cinterval &s) throw();
	//! Implementation of division and allocation operation
	INLINE civector &operator /=(civector &rv,const cinterval &r) throw();

//---------------------------- Real --------------------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const civector &rv, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector_slice &sl, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const real &s, const civector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const real &s, const civector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE civector &operator *=(civector &rv,const real &r) throw();

	//! Implementation of division operation
	INLINE civector operator /(const civector &rv, const real &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const civector_slice &sl, const real &s) throw();
	//! Implementation of division and allocation operation
	INLINE civector &operator /=(civector &rv,const real &r) throw();

	//! Implementation of multiplication operation
	INLINE civector operator *(const rvector &rv, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const rvector_slice &sl, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const rvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const rvector_slice &sl) throw();

	//! Implementation of division operation
	INLINE civector operator /(const rvector &rv, const cinterval &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const rvector_slice &sl, const cinterval &s) throw();

//---------------------------- Complex --------------------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const civector &rv, const complex &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector_slice &sl, const complex &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const complex &s, const civector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const complex &s, const civector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE civector &operator *=(civector &rv,const complex &r) throw();

	//! Implementation of division operation
	INLINE civector operator /(const civector &rv, const complex &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const civector_slice &sl, const complex &s) throw();
	//! Implementation of division and allocation operation
	INLINE civector &operator /=(civector &rv,const complex &r) throw();

	//! Implementation of multiplication operation
	INLINE civector operator *(const cvector &rv, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cvector_slice &sl, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const cvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const cvector_slice &sl) throw();

	//! Implementation of division operation
	INLINE civector operator /(const cvector &rv, const cinterval &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const cvector_slice &sl, const cinterval &s) throw();

//---------------------------- interval --------------------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const civector &rv, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector_slice &sl, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const interval &s, const civector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const interval &s, const civector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE civector &operator *=(civector &rv,const interval &r) throw();

	//! Implementation of division operation
	INLINE civector operator /(const civector &rv, const interval &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const civector_slice &sl, const interval &s) throw();
	//! Implementation of division and allocation operation
	INLINE civector &operator /=(civector &rv,const interval &r) throw();

	//! Implementation of multiplication operation
	INLINE civector operator *(const ivector &rv, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const ivector_slice &sl, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const ivector_slice &sl) throw();

	//! Implementation of division operation
	INLINE civector operator /(const ivector &rv, const cinterval &s) throw();
	//! Implementation of division operation
	INLINE civector operator /(const ivector_slice &sl, const cinterval &s) throw();

//======================= Vector / Vector ===============================


	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s, const civector &rv) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &o, const civector_slice &sl) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, civector &rv) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, civector_slice &rv) throw();
	
//----------------------- cinterval / cinterval ---------------------------

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(cidotprecision &dp, const cvector &);

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(cidotprecision &dp, const rvector &);

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(cidotprecision &dp, const civector &);

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(cidotprecision &dp, const ivector &);


	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const cvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of positive sign operation
	INLINE const civector &operator +(const civector &rv) throw();
	//! Implementation of positive sign operation
	INLINE civector operator +(const civector_slice &sl) throw();

	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector & operator +=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector &operator +=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of negative sign operation
	INLINE civector operator -(const civector &rv) throw();
	//! Implementation of negative sign operation
	INLINE civector operator -(const civector_slice &sl) throw();
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector & operator -=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector &operator -=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector & operator |=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector &operator |=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector & operator &=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector &operator &=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of standard equality operation
	INLINE bool operator ==(const civector &rv1, const civector &rv2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const civector_slice &sl1, const civector_slice &sl2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const civector_slice &sl, const civector &rv) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const civector &rv, const civector_slice &sl) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const civector &rv1, const civector &rv2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const civector_slice &sl1, const civector_slice &sl2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const civector_slice &sl, const civector &rv) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const civector &rv, const civector_slice &sl) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const civector &rv1, const civector &rv2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const civector_slice &sl1, const civector_slice &sl2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator < (const civector_slice &sl, const civector &rv) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator < (const civector &rv, const civector_slice &sl) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const civector &rv1, const civector &rv2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const civector_slice &sl1, const civector_slice &sl2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const civector_slice &sl, const civector &rv) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const civector &rv, const civector_slice &sl) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const civector &rv1, const civector &rv2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const civector_slice &sl1, const civector_slice &sl2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const civector_slice &sl, const civector &rv) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const civector &rv, const civector_slice &sl) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const civector &rv1, const civector &rv2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const civector_slice &sl1, const civector_slice &sl2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const civector_slice &sl, const civector &rv) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const civector &rv, const civector_slice &sl) throw();

//-------------------------------- cinterval / Real --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const civector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const civector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE cinterval operator *(const rvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const rvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE civector operator +(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE civector & operator +=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector &operator +=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE civector operator -(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE civector operator -(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE civector & operator -=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector &operator -=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector & operator |=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector &operator |=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE civector & operator &=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector &operator &=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
//-------------------------------- cinterval / complex --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const civector_slice &sl,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const civector &rv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif




	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE civector operator +(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE civector & operator +=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector &operator +=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE civector operator -(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE civector & operator -=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector &operator -=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector & operator |=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector &operator |=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE civector & operator &=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector &operator &=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

//-------------------------------- cinterval / interval --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const civector_slice &sl,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const civector &rv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const civector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const civector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const civector_slice & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const rmatrix_subv &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rmatrix_subv & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE cinterval operator *(const ivector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const ivector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE civector operator +(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE civector & operator +=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE civector &operator +=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector operator -(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE civector & operator -=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE civector &operator -=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector & operator |=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE civector &operator |=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE civector & operator &=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE civector &operator &=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

//------------- real x complex ------------------------
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector_slice &sl,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector &rv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const rvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

//------------- complex x complex ------------------------
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

//-------------------------------- interval / complex --------------------------------

// multiplication in iveccvec.hpp
	
	//! Implementation of addition operation
	INLINE civector operator +(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE civector operator +(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif


	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE civector operator -(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif


	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE civector operator |(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE civector operator &(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>);
#else
	throw();
#endif

        //! Checks if v1 lies in the interior of v2
	INLINE bool in(const civector& v1, const civector& v2);

} // namespace cxsc 

#ifdef _CXSC_INCL_INL
#include "vector.inl"
#include "civector.inl"
#endif

#ifdef _CXSC_RMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "civecrmat.inl"
# else
#  include "civecrmat.hpp"
# endif
#endif

#ifdef _CXSC_CMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "civeccmat.inl"
# else
#  include "civeccmat.hpp"
# endif
#endif

#ifdef _CXSC_IMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "civecimat.inl"
# else
#  include "civecimat.hpp"
# endif
#endif

#ifdef CXSC_USE_BLAS
#define _CXSC_BLAS_CIVECTOR
#include "cxsc_blas.inl"
#endif


#endif
