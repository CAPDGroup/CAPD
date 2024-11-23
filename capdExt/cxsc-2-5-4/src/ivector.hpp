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

/* CVS $Id: ivector.hpp,v 1.36 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_IVECTOR_HPP_INCLUDED
#define _CXSC_IVECTOR_HPP_INCLUDED

#include "xscclass.hpp"
#include "except.hpp"
#include "idot.hpp"
#include "cidot.hpp"
#include "interval.hpp" // used for declaration of Inf, Sup,...
#include "rvector.hpp"
#include "vector.hpp"


#include <iostream>

namespace cxsc {

class srvector;
class srvector_slice;
class sivector;
class sivector_slice;
class ivector_slice;

//! The Data Type ivector
/*!
The vectors of C-XSC are one dimensional arrays of the corresponding scalar base type. 

\sa rvector
*/
class ivector
{
	friend class ivector_slice;
	friend class imatrix;
	friend class imatrix_subv;
	friend class civector;
	friend class cimatrix;
	friend class l_ivector;
	friend class l_imatrix;
	private:
	interval *dat;
	int l,u,size;

	public:
//#if(CXSC_INDEX_CHECK)
#ifdef _CXSC_FRIEND_TPL
	//------------ Templates --------------------------------------------------
	// Interval
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
 template <class V,class E> friend 	 E _vabs(const V &rv) throw();
 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();
template <class MV,class V> friend  V _mvabs(const MV &mv) throw();
 template <class V,class E> friend 	 E _vdiam(const V &rv) throw();
 template <class V,class E> friend 	 E _vmid(const V &rv) throw();
 template <class V,class E> friend 	 E _vinf(const V &rv) throw();
 template <class V,class E> friend 	 E _vsup(const V &rv) throw();

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
 template <class V1,class V2,class E> friend 	 E _vvimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class VS,class V,class E> friend 	 E _vsvimult(const VS & sl, const V &rv)
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
 template <class V,class MV,class S> friend 	 S _vmvimult(const V &rv1, const MV &rv2)
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
 template <class M,class V,class E> friend 	 E _mvimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmimult(const V &v,const M &m)
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
 template <class MS,class V,class E> friend 	 E _msvimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmsimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class S> friend 	 V &_vmsimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	// Real
  	//--Real -------- vector-scalar ------------
	//--Real--------- Vector-vector---------
	//-- Real -------- Vector-matrix ----------
	// complex
	//--complex -------- vector-scalar ------------
	//--complex--------- Vector-vector---------

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
	
	//-- complex -------- Vector-matrix ----------
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
 template <class V,class MS,class E> friend 	 E _vmscimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
	
	// civector ---------------------
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
 template <class V,class E> friend 	 E _vim(const V &rv) throw();
 template <class V,class E> friend 	 E _vre(const V &rv) throw();

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
template <class MV,class V> friend  MV &_mvvsetre(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvsetim(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  V _mvim(const MV &mv) throw();
template <class MV,class V> friend  V _mvre(const MV &mv) throw();

	// vector-scalar
	// vector-vector
 template <class VS1,class VS2,class E> friend 	 E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	// vector-matrix
	// l_real
	//--l_real -------- vector-scalar ------------
	//--l_real--------- Vector-vector---------

 template <class V1,class V2,class E> friend 	 E _vvlimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class VS,class V,class E> friend 	 E _vsvlimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
	
	//-- l_real -------- Vector-matrix ----------
 template <class M,class V,class E> friend 	 E _mvlimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class V,class E> friend 	 E _msvlimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmlimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmslimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
	
	// l_ivector ---------------------
	// vector-scalar
	// vector-vector
 template <class VS1,class VS2,class E> friend 	 E _vsvslimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
	// vector-matrix

#endif

	//------ Konstruktoren ----------------------------------------------------
	//! Constructor of class ivector
	INLINE ivector () throw();
	//! Constructor of class ivector
	explicit INLINE ivector(const int &i) throw();
#ifdef OLD_CXSC
	//! Constructor of class ivector
	explicit INLINE ivector(const class index &i) throw(); // for backwards compatibility
#endif
	//! Constructor of class ivector
	explicit INLINE ivector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_WRONG_BOUNDARIES,ERROR_IVECTOR_NO_MORE_MEMORY);
#else
	throw();
#endif
	//! Constructor of class ivector
	INLINE ivector(const imatrix_subv &) throw();
	//! Constructor of class ivector
	explicit INLINE ivector(const interval &) throw();
	//! Constructor of class ivector
//	explicit INLINE ivector(const imatrix &)
	explicit        ivector(const imatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class ivector
	explicit INLINE ivector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class ivector
	INLINE ivector(const ivector_slice &rs) throw();
	//! Constructor of class ivector
	INLINE ivector(const ivector &v) throw();
	//! Constructor of class ivector
	INLINE ivector(const sivector &rs);
	//! Constructor of class ivector
	INLINE ivector(const sivector_slice &rs);

	// Real
	//! Constructor of class ivector
	explicit INLINE ivector(const real &) throw();
	//! Constructor of class ivector
	explicit INLINE ivector(const rvector_slice &rs) throw();
	//! Constructor of class ivector
	explicit INLINE ivector(const srvector &rs);
	//! Constructor of class ivector
	explicit INLINE ivector(const srvector_slice &rs);

	//! Constructor of class ivector
	explicit INLINE ivector(const rvector &v) throw();
	//! Constructor of class ivector
	explicit INLINE ivector(const rmatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class ivector
	explicit INLINE ivector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class ivector
//	explicit INLINE ivector(const rmatrix_subv &) throw();
        explicit        ivector(const rmatrix_subv &) throw();	
	
	// Interval
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const ivector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const ivector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const sivector &rv);
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const sivector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const srvector &rv);
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const srvector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const imatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const imatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const rvector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const rvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const rmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector &operator =(const rmatrix_subv &) throw();

	//! Implementation of assignment and addition operator	
        INLINE ivector &operator+=(const srvector&);
	//! Implementation of assignment and addition operator	
        INLINE ivector &operator+=(const sivector&);
	//! Implementation of assignment and addition operator	
        INLINE ivector &operator+=(const srvector_slice&);
	//! Implementation of assignment and addition operator	
        INLINE ivector &operator+=(const sivector_slice&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector &operator-=(const srvector&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector &operator-=(const sivector&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector &operator-=(const srvector_slice&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector &operator-=(const sivector_slice&);
	//! Implementation of assignment and hull operator	
        INLINE ivector &operator|=(const srvector&);
	//! Implementation of assignment and hull operator	
        INLINE ivector &operator|=(const sivector&);
	//! Implementation of assignment and hull operator	
        INLINE ivector &operator|=(const srvector_slice&);
	//! Implementation of assignment and hull operator	
        INLINE ivector &operator|=(const sivector_slice&);
	//! Implementation of assignment and intersection operator	
        INLINE ivector &operator&=(const sivector&);
	//! Implementation of assignment and intersection operator	
        INLINE ivector &operator&=(const sivector_slice&);

        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE ivector operator()(const intvector& p);
        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE ivector operator()(const intmatrix& P);

	//--------- Destruktor ----------------------------------------------------
	INLINE ~ivector() { delete [] dat; }

	//------ Standardfunktionen -----------------------------------------------
	
	friend INLINE interval::interval(const ivector &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const ivector &rv) throw() { return rv.l; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const ivector &rv) throw() { return rv.u; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const ivector &rv) throw() { return rv.size; }
	//! Sets the lower bound of the vector
	friend INLINE ivector &SetLb(ivector &rv, const int &l) throw() { rv.l=l; rv.u=l+rv.size-1; return rv; }
	//! Sets the upper bound of the vector
	friend INLINE ivector &SetUb(ivector &rv, const int &u) throw() { rv.u=u; rv.l=u-rv.size+1; return rv; }
	//! Operator for accessing the single elements of the vector (read-only)
	INLINE interval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector 
	INLINE interval & operator [](const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	
	//! Operator for accessing the whole vector
	INLINE ivector & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE ivector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE ivector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	INLINE operator void*() throw();
//#else
//#endif
};


//! The Data Type ivector_slice
/*!
This data type represents a partial ivector.

\sa ivector
*/
class ivector_slice
{
	friend class ivector;
	friend class imatrix;
	friend class civector;
	friend class cimatrix;
	friend class l_ivector;
	friend class l_imatrix;
	private:
	interval *dat;
	int l,u,size;
	int start,end;

	public:
//#if(CXSC_INDEX_CHECK)	
#ifdef _CXSC_FRIEND_TPL
//------------------------- Templates -------------------------------------------
// Interval / Interval

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
 template <class VS,class S> friend 	 VS &_vsssetinf(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vsssetsup(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vssusetinf(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vssusetsup(VS &vs, const S &s) throw();

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

 template <class VS,class V,class E> friend 	 E _vsvimult(const VS & sl, const V &rv)
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
 template <class VS,class M,class S> friend 	 VS &_vsmimultassign(VS &v,const M &m)
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
 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsdiam(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsmid(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsinf(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vssup(const VS &sl) throw();
 template <class VS1,class VS2,class E> friend 	 E _vsvsimult(const VS1 & sl1, const VS2 &sl2)
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
	
	// Interval / Real
 template <class V,class MS,class E> friend 	 E _vmsimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	// complex
 template <class VS,class V,class E> friend 	 E _vsvcimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
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

	// cinterval-- vector-scalar
	// cinterval -- vector-vector
template <class VS1,class VS2,class E> friend 	 E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	// cinterval -- vector-matrix
	// l_real

 template <class VS1,class VS2,class E> friend 	 E _vsvslimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif

	// l_interval ----------
	// l_interval-- vector-scalar
	// l_interval -- vector-vector
template <class VS,class V,class E> friend 	 E _vsvlimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif


#endif
	
	//--------------------- Konstruktoren -----------------------------------
	//! Constructor of class ivector_slice
	explicit INLINE ivector_slice(ivector &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	//! Constructor of class ivector_slice
	explicit INLINE ivector_slice(ivector_slice &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	public: 
	//! Constructor of class ivector_slice
	INLINE ivector_slice(const ivector_slice &a) throw():dat(a.dat),l(a.l),u(a.u),size(a.size),start(a.start),end(a.end) { }
	public:
	// Interval
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const sivector &sl);
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const sivector_slice &sl);
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const srvector &sl);
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const srvector_slice &sl);

	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice &operator =(const imatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice & operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE ivector_slice &operator =(const rmatrix_subv &mv) throw();

	//--------------------- Standardfunktionen ------------------------------

	friend INLINE interval::interval(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_IVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const ivector_slice &sl) throw() { return sl.start; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const ivector_slice &sl) throw() { return sl.end; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const ivector_slice &sl) throw() { return sl.end-sl.start+1; }
	//! Operator for accessing the single elements of the vector (read-only)
	INLINE interval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector
	INLINE interval & operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the whole vector
	INLINE ivector_slice & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE ivector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE ivector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	//! Implementation of division and allocation operation
	INLINE ivector_slice &operator /=(const interval &r) throw();
	//! Implementation of division and allocation operation
	INLINE ivector_slice &operator /=(const real &r) throw();
	//! Implementation of multiplication and allocation operation
	INLINE ivector_slice &operator *=(const interval &r) throw();
	//! Implementation of multiplication and allocation operation
	INLINE ivector_slice &operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE ivector_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE ivector_slice &operator *=(const real &r) throw();
	//! Implementation of addition and allocation operation
	INLINE ivector_slice &operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE ivector_slice &operator +=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector_slice &operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector_slice &operator -=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector_slice &operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector_slice &operator |=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector_slice &operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector_slice &operator &=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	
	//! Implementation of addition and allocation operation
	INLINE ivector_slice &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE ivector_slice &operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector_slice &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector_slice &operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector_slice &operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector_slice &operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector_slice &operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector_slice &operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	INLINE operator void*() throw();
//#else
//#endif

	//! Implementation of assignment and addition operator	
        INLINE ivector_slice &operator+=(const srvector&);
	//! Implementation of assignment and addition operator	
        INLINE ivector_slice &operator+=(const sivector&);
	//! Implementation of assignment and addition operator	
        INLINE ivector_slice &operator+=(const srvector_slice&);
	//! Implementation of assignment and addition operator	
        INLINE ivector_slice &operator+=(const sivector_slice&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector_slice &operator-=(const srvector&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector_slice &operator-=(const sivector&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector_slice &operator-=(const srvector_slice&);
	//! Implementation of assignment and substraction operator	
        INLINE ivector_slice &operator-=(const sivector_slice&);
	//! Implementation of assignment and hull operator	
        INLINE ivector_slice &operator|=(const srvector&);
	//! Implementation of assignment and hull operator	
        INLINE ivector_slice &operator|=(const sivector&);
	//! Implementation of assignment and hull operator	
        INLINE ivector_slice &operator|=(const srvector_slice&);
	//! Implementation of assignment and hull operator	
        INLINE ivector_slice &operator|=(const sivector_slice&);
	//! Implementation of assignment and intersection operator	
        INLINE ivector_slice &operator&=(const sivector&);
	//! Implementation of assignment and intersection operator	
        INLINE ivector_slice &operator&=(const sivector_slice&);

};

//=======================================================================
//======================== Vector Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE ivector _ivector(const interval &r) throw();
//	INLINE ivector _ivector(const imatrix &m) throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE ivector _ivector(const imatrix_slice &sl) throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE ivector _ivector(const real &r) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE ivector _ivector(const rvector_slice &rs) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE ivector _ivector(const rvector &rs) throw();
//	INLINE ivector _ivector(const rmatrix &m) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE ivector _ivector(const rmatrix_slice &sl) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE ivector _ivector(const rmatrix_subv &rs) throw();

	//! Returns the vector with the new given infimum vector
	INLINE ivector &SetInf(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE ivector_slice &SetInf(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE ivector &SetInf(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE ivector_slice &SetInf(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE ivector &UncheckedSetInf(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE ivector_slice &UncheckedSetInf(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE ivector &UncheckedSetInf(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE ivector_slice &UncheckedSetInf(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given supremum vector
	INLINE ivector &SetSup(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE ivector_slice &SetSup(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE ivector &SetSup(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE ivector_slice &SetSup(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE ivector &UncheckedSetSup(ivector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE ivector_slice &UncheckedSetSup(ivector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE ivector &UncheckedSetSup(ivector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE ivector_slice &UncheckedSetSup(ivector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given supremum value
	INLINE ivector &SetSup(ivector &iv,const real &r) throw();
	//! Returns the vector with the new given infimum value
	INLINE ivector &SetInf(ivector &iv,const real &r) throw();
	//! Returns the vector with the new unchecked given supremum value
	INLINE ivector &UncheckedSetSup(ivector &iv,const real &r) throw();
	//! Returns the vector with the new unchecked given infimum value
	INLINE ivector &SetUncheckedInf(ivector &iv,const real &r) throw();
	//! Returns the vector with the new given supremum value
	INLINE ivector_slice &SetSup(ivector_slice &iv,const real &r) throw();
	//! Returns the vector with the new given infimum value
	INLINE ivector_slice &SetInf(ivector_slice &iv,const real &r) throw();
	//! Returns the vector with the new unchecked given supremum value
	INLINE ivector_slice &UncheckedSetSup(ivector_slice &iv,const real &r) throw();
	//! Returns the vector with the new unchecked given infimum value
	INLINE ivector_slice &SetUncheckedInf(ivector_slice &iv,const real &r) throw();

	//! Resizes the vector
	INLINE void Resize(ivector &rv) throw();
	//! Resizes the vector
	INLINE void Resize(ivector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<ivector>);
#else
	throw();
#endif
	//! Resizes the vector
	INLINE void Resize(ivector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<ivector>);
#else
	throw();
#endif
	
	//! Returns the absolute value of the vector
	INLINE ivector abs(const ivector &rv) throw();
	//! Returns the absolute value of the vector
	INLINE ivector abs(const ivector_slice &sl) throw();
	//! Returns the absolute minimum value of the vector
	INLINE rvector absmin(const ivector &rv) throw();
	//! Returns the absolute maximum value of the vector
	INLINE rvector absmax(const ivector_slice &sl) throw();
	//! Returns the absolute minimum value of the vector
	INLINE rvector absmin(const ivector &rv) throw();
	//! Returns the absolute maximum value of the vector
	INLINE rvector absmax(const ivector_slice &sl) throw();
        //! Returns the rounded diameter of the vector
	INLINE rvector diam(const ivector &v) throw();
        //! Returns the rounded diameter of the vector
	INLINE rvector diam(const ivector_slice &v) throw();
	//! Returns the rounded middle of the vector
	INLINE rvector mid(const ivector &v) throw();
	//! Returns the rounded middle of the vector
	INLINE rvector mid(const ivector_slice &v) throw();
	//! Returns the infimum of the vector
	INLINE rvector Inf(const ivector &v) throw();
	//! Returns the infimum of the vector
	INLINE rvector Inf(const ivector_slice &v) throw();
	//! Returns the supremum of the vector
	INLINE rvector Sup(const ivector &v) throw();
	//! Returns the supremum of the vector
	INLINE rvector Sup(const ivector_slice &v) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const ivector &rv) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const ivector_slice &sl) throw();

//======================= Vector / Scalar ===============================

//----------------------------- Interval ---------------------------

	//! Implementation of multiplication operation
	INLINE ivector operator *(const ivector &rv, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const ivector_slice &sl, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const interval &s, const ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const interval &s, const ivector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE ivector &operator *=(ivector &rv,const interval &r) throw();

	//! Implementation of division operation
	INLINE ivector operator /(const ivector &rv, const interval &s) throw();
	//! Implementation of division operation
	INLINE ivector operator /(const ivector_slice &sl, const interval &s) throw();
	//! Implementation of division and allocation operation
	INLINE ivector &operator /=(ivector &rv,const interval &r) throw();

//---------------------------- Real --------------------------------------

	//! Implementation of multiplication operation
	INLINE ivector operator *(const ivector &rv, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const ivector_slice &sl, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const real &s, const ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const real &s, const ivector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE ivector &operator *=(ivector &rv,const real &r) throw();

	//! Implementation of division operation
	INLINE ivector operator /(const ivector &rv, const real &s) throw();
	//! Implementation of division operation
	INLINE ivector operator /(const ivector_slice &sl, const real &s) throw();
	//! Implementation of division and allocation operation
	INLINE ivector &operator /=(ivector &rv,const real &r) throw();

	//! Implementation of multiplication operation
	INLINE ivector operator *(const rvector &rv, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const rvector_slice &sl, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const interval &s, const rvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE ivector operator *(const interval &s, const rvector_slice &sl) throw();

	//! Implementation of division operation
	INLINE ivector operator /(const rvector &rv, const interval &s) throw();
	//! Implementation of division operation
	INLINE ivector operator /(const rvector_slice &sl, const interval &s) throw();

//======================= Vector / Vector ===============================


	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s, const ivector &rv) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &o, const ivector_slice &sl) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, ivector &rv) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, ivector_slice &rv) throw();
	
//----------------------- Interval / Interval ---------------------------

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(idotprecision &dp, const ivector &);

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(idotprecision &dp, const rvector &);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector_slice & sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	
	//! Implementation of positive sign operation
	INLINE const ivector &operator +(const ivector &rv) throw();
	//! Implementation of positive sign operation
	INLINE ivector operator +(const ivector_slice &sl) throw();

	//! Implementation of addition operation
	INLINE ivector operator +(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE ivector & operator +=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE ivector &operator +=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of negative sign operation
	INLINE ivector operator -(const ivector &rv) throw();
	//! Implementation of negative sign operation
	INLINE ivector operator -(const ivector_slice &sl) throw();
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector & operator -=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector &operator -=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector & operator |=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector &operator |=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector & operator &=(ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector &operator &=(ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of standard equality operation
	INLINE bool operator ==(const ivector &rv1, const ivector &rv2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const ivector_slice &sl1, const ivector_slice &sl2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const ivector_slice &sl, const ivector &rv) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const ivector &rv, const ivector_slice &sl) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const ivector &rv1, const ivector &rv2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const ivector_slice &sl1, const ivector_slice &sl2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const ivector_slice &sl, const ivector &rv) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const ivector &rv, const ivector_slice &sl) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const ivector &rv1, const ivector &rv2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const ivector_slice &sl1, const ivector_slice &sl2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator < (const ivector_slice &sl, const ivector &rv) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator < (const ivector &rv, const ivector_slice &sl) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const ivector &rv1, const ivector &rv2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const ivector_slice &sl1, const ivector_slice &sl2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const ivector_slice &sl, const ivector &rv) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const ivector &rv, const ivector_slice &sl) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const ivector &rv1, const ivector &rv2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const ivector_slice &sl1, const ivector_slice &sl2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const ivector_slice &sl, const ivector &rv) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const ivector &rv, const ivector_slice &sl) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const ivector &rv1, const ivector &rv2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const ivector_slice &sl1, const ivector_slice &sl2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const ivector_slice &sl, const ivector &rv) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const ivector &rv, const ivector_slice &sl) throw();

//-------------------------------- Interval / Real --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rvector_slice & sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp,const ivector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rvector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rvector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp,const ivector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const imatrix_subv & rv1, const rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(idotprecision &dp, const rvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector_slice & sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const ivector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector_slice & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp,const ivector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const ivector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE interval operator *(const rvector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const rvector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE interval operator *(const ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE ivector operator +(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE ivector operator +(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE ivector operator +(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE ivector & operator +=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE ivector &operator +=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE ivector operator -(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE ivector operator -(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE ivector & operator -=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE ivector &operator -=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE ivector operator |(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector & operator |=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE ivector &operator |=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE ivector operator &(const rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE ivector operator &(const ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector & operator &=(ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE ivector &operator &=(ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>);
#else
	throw();
#endif

//! Checks if first argument is part of second argument
extern int      in         ( const ivector&, const ivector& );
//! Checks if first argument is part of second argument
extern int      in         ( int,      ivector& );
//! Performs an epsilon inflation
extern ivector  Blow       ( const ivector&, real );
//! Checks arguments for disjointness
extern int      Disjoint   ( ivector&, ivector& );
//! Checks if vector is zero vector
extern int      Zero       ( ivector& );
//! Returns the rounded middle of the vector
extern rvector  mid        ( ivector& );
//! Computes the relative diameter \f$ d_{rel}((x)) \f$
extern real     MaxRelDiam ( const ivector& );
//! Computes the relative diameter \f$ d_{rel}((x)) \f$
extern real     MaxRelDiam ( const ivector_slice& );
//! Checks if the diameter of the vector \f$ [ x ] \f$ is \f$ \le n \f$ ulps
extern int      UlpAcc     ( ivector&, int );
//! Doubles the dimension of the vector preserving existing components
extern void     DoubleSize ( ivector& );


} // namespace cxsc 

#ifdef _CXSC_INCL_INL
#include "vector.inl"
#include "ivector.inl"
#endif
		
#ifdef _CXSC_RMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "ivecrmat.inl"
# else
#  include "ivecrmat.hpp"
# endif
#endif
		
#ifdef _CXSC_CMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "iveccmat.inl"
# else
#  include "iveccmat.hpp"
# endif
#endif
		
#ifdef _CXSC_LRMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "iveclrmat.inl"
# else
#  include "iveclrmat.hpp"
# endif
#endif
		
#ifdef _CXSC_CVECTOR_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "iveccvec.inl"
# else
#  include "iveccvec.hpp"
# endif
#endif

#ifdef _CXSC_LRVECTOR_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "lrvecivec.inl"
# else
#  include "lrvecivec.hpp"
# endif
#endif

#ifdef CXSC_USE_BLAS
#define _CXSC_BLAS_IVECTOR
#include "cxsc_blas.inl"
#endif

#endif
