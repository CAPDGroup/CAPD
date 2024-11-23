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

/* CVS $Id: l_ivector.hpp,v 1.19 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_LIVECTOR_HPP_INCLUDED
#define _CXSC_LIVECTOR_HPP_INCLUDED

#include "xscclass.hpp"
#include "except.hpp"
#include "idot.hpp"
#include "l_interval.hpp" // used for declaration of Inf, Sup,...
//#include "cxscmatr.hpp"
#include "rvector.hpp"
#include "ivector.hpp"
#include "l_rvector.hpp"
#include "vector.hpp"


#include <iostream>

//#include "matrix.hpp" // hat hier eigentlich nichts zu suchen, sonst aber Internal Compiler Error #9

namespace cxsc {

class l_ivector_slice;

//! The Multiple-Precision Data Type l_ivector
/*!
The vectors of C-XSC are one dimensional arrays of the corresponding scalar base type. 

\sa l_rvector
*/
class l_ivector
{
	friend class l_ivector_slice;
	friend class l_imatrix;
	friend class l_imatrix_subv;
	private:
	l_interval *dat;
	int l,u,size;

	public:
//#if(CXSC_INDEX_CHECK)
#ifdef _CXSC_FRIEND_TPL
	//------------ Templates --------------------------------------------------
	// l_interval
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
 template <class V,class MV,class S> friend 	 S _vmvlimult(const V &rv1, const MV &rv2)
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
 template <class M,class V,class E> friend 	 E _mvlimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmlimult(const V &v,const M &m)
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
 template <class V,class M,class S> friend 	 V &_vmlimultassign(V &v,const M &m)
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
 template <class V,class MS,class E> friend 	 E _vmslimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class S> friend 	 V &_vmslimultassign(V &v,const MS &ms)
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

	// l_rvector x ivector ----------------
	// vector - matrix ------------
  /*   friend TINLINE l_ivector _mvslimult<imatrix,l_rvector_slice,l_ivector>(const imatrix &m,const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
  /*   friend TINLINE l_ivector _vsmlimult<l_rvector_slice,imatrix,l_ivector>(const l_rvector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */


#endif
	

	//------ Konstruktoren ----------------------------------------------------
	//! Constructor of class l_ivector
	INLINE l_ivector () throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const int &i) throw();
#ifdef OLD_CXSC
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const class index &i) throw(); // for backwards compatibility
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_WRONG_BOUNDARIES,ERROR_LIVECTOR_NO_MORE_MEMORY);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	INLINE l_ivector(const l_imatrix_subv &) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_interval &) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_imatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	INLINE l_ivector(const l_ivector_slice &rs) throw();
	//! Constructor of class l_ivector
	INLINE l_ivector(const l_ivector &v) throw();
	// Real
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const real &) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const rvector_slice &rs) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const rvector &v) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const rmatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const rmatrix_subv &) throw();
	
	// l_real
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_real &) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_rvector_slice &rs) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_rvector &v) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_rmatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const l_rmatrix_subv &) throw();
	
	// interval
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const interval &) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const ivector_slice &rs) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const ivector &v) throw();
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const imatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class l_ivector
	explicit INLINE l_ivector(const imatrix_subv &) throw();
	
	// l_interval
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_ivector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_ivector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_imatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_imatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const rvector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const rvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const rmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const rmatrix_subv &) throw();

	// l_real
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_rvector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_rvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_rmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const l_rmatrix_subv &) throw();

	// interval
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const ivector &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const ivector_slice &sl) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const imatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector &operator =(const imatrix_subv &) throw();

	//--------- Destruktor ----------------------------------------------------
	INLINE ~l_ivector() { delete [] dat; }

	//------ Standardfunktionen -----------------------------------------------
	
	friend INLINE l_interval::l_interval(const l_ivector &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LIVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const l_ivector &rv) throw() { return rv.l; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const l_ivector &rv) throw() { return rv.u; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const l_ivector &rv) throw() { return rv.size; }
	//! Sets the lower bound of the vector
	friend INLINE l_ivector & SetLb(l_ivector &rv, const int &l) throw() { rv.l=l; rv.u=l+rv.size-1; return rv;}
	//! Sets the upper bound of the vector
	friend INLINE l_ivector & SetUb(l_ivector &rv, const int &u) throw() { rv.u=u; rv.l=u-rv.size+1; return rv;}
	//! Operator for accessing the single elements of the vector
	INLINE l_interval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	//! Operator for accessing the whole vector
	INLINE l_ivector & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE l_ivector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE l_ivector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	INLINE operator void*() throw();
//#else
//#endif
};


//! The Multiple-Precision Data Type l_ivector_slice
/*!
This data type represents a partial ivector.

\sa l_ivector
*/
class l_ivector_slice
{
	friend class l_ivector;
	friend class l_imatrix;
	private:
	l_interval *dat;
	int l,u,size;
	int start,end;

	public:
//#if(CXSC_INDEX_CHECK)	
#ifdef _CXSC_FRIEND_TPL
//------------------------- Templates -------------------------------------------
// l_interval / l_interval

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

 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsdiam(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsmid(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsinf(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vssup(const VS &sl) throw();
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
 template <class VS,class V,class E> friend 	 E _vsvlimult(const VS & sl, const V &rv)
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
 template <class VS,class M,class S> friend 	 VS &_vsmlimultassign(VS &v,const M &m)
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
 template <class VS1,class VS2,class E> friend 	 E _vsvslimult(const VS1 & sl1, const VS2 &sl2)
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
	
	// l_interval / Real
 template <class V,class MS,class E> friend 	 E _vmslimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
	// l_interval / l_real
	 template <class DP,class V1,class V2> friend 	 void _vvaccu(DP &dp, const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	// l_real
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

	//--l_real -------- vector-scalar ------------
 template <class MV,class S,class E> friend 	 E _mvsmult(const MV &rv, const S &s) throw();
 template <class V,class S,class E> friend 	 E _vsmult(const V &rv, const S &s) throw();
 template <class V,class S,class E> friend 	 E _vsdiv(const V &rv, const S &s) throw();
 template <class V,class S> friend 	 V &_vsdivassign(V &rv,const S &r) throw();
 template <class V,class S> friend 	 V &_vsmultassign(V &rv,const S &r) throw();

	//--l_real--------- Vector-vector---------
 template <class V1,class V2,class E> friend 	 E _vvlimult(const V1 & rv1, const V2 &rv2)
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
	
	//-- l_real -------- Vector-matrix ----------
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
 template <class V,class MS,class S> friend 	 V &_vmslimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class M,class S> friend 	 V &_vmlimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

  /*	friend TINLINE l_ivector &_vsmassign<l_ivector_slice,imatrix,l_interval>(l_ivector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
	#endif */

	//-- interval -------- Vector-matrix ----------
  /*   friend TINLINE l_ivector _mvslimult<imatrix,l_ivector_slice,l_ivector>(const imatrix &m,const l_ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
  /*   friend TINLINE l_ivector _msvslimult<imatrix_slice,l_ivector_slice,l_ivector>(const imatrix_slice &ms,const l_ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
  /*   friend TINLINE l_ivector _vsmlimult<l_ivector_slice,imatrix,l_ivector>(const l_ivector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
  /*   friend TINLINE l_ivector _vsmslimult<l_ivector_slice,imatrix_slice,l_ivector>(const l_ivector_slice &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
  /*   friend TINLINE l_ivector &_vsmslimultassign<l_ivector_slice,imatrix_slice,l_interval>(l_ivector_slice &v,const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */

#endif
	
	//--------------------- Konstruktoren -----------------------------------
	//! Constructor of class l_ivector_slice
	explicit INLINE l_ivector_slice(l_ivector &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	//! Constructor of class l_ivector_slice
	explicit INLINE l_ivector_slice(l_ivector_slice &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	public: 
	//! Constructor of class l_ivector_slice
	INLINE l_ivector_slice(const l_ivector_slice &a) throw():dat(a.dat),l(a.l),u(a.u),size(a.size),start(a.start),end(a.end) { }
	public:
	// l_interval
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>,ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>,ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice &operator =(const l_imatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice &operator =(const rmatrix_subv &mv) throw();

	// l_real
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>,ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>,ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice &operator =(const l_rmatrix_subv &mv) throw();

	// interval
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice & operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_ivector_slice &operator =(const imatrix_subv &mv) throw();

	//--------------------- Standardfunktionen ------------------------------

	friend INLINE l_interval::l_interval(const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_LIVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const l_ivector_slice &sl) throw() { return sl.start; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const l_ivector_slice &sl) throw() { return sl.end; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const l_ivector_slice &sl) throw() { return sl.end-sl.start+1; }
	//! Operator for accessing the single elements of the vector
	INLINE l_interval & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	//! Operator for accessing the whole vector
	INLINE l_ivector_slice & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE l_ivector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE l_ivector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();
	
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const l_interval &r) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector_slice &operator /=(const l_interval &r) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const real &r) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector_slice &operator /=(const real &r) throw();
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const l_real &r) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector_slice &operator /=(const l_real &r) throw();
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const interval &r) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector_slice &operator /=(const interval &r) throw();
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector_slice &operator +=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector_slice &operator -=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector_slice &operator |=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector_slice &operator &=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector_slice &operator *=(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
//#else
//#endif
};

//=======================================================================
//======================== Vector Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const l_interval &r) throw();
//	INLINE l_ivector _l_ivector(const l_imatrix &m) throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE l_ivector _l_ivector(const l_imatrix_slice &sl) throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const real &r) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const rvector_slice &rs) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const rvector &rs) throw();
//	INLINE l_ivector _l_ivector(const rmatrix &m) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE l_ivector _l_ivector(const rmatrix_slice &sl) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_ivector _l_ivector(const rmatrix_subv &rs) throw();

	//! Returns the vector with the new given infimum vector
	INLINE l_ivector &SetInf(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE l_ivector_slice &SetInf(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE l_ivector &SetInf(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given infimum vector
	INLINE l_ivector_slice &SetInf(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE l_ivector &UncheckedSetInf(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE l_ivector_slice &UncheckedSetInf(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE l_ivector &UncheckedSetInf(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given infimum vector
	INLINE l_ivector_slice &UncheckedSetInf(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given supremum vector
	INLINE l_ivector &SetSup(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE l_ivector_slice &SetSup(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE l_ivector &SetSup(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given supremum vector
	INLINE l_ivector_slice &SetSup(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE l_ivector &UncheckedSetSup(l_ivector &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE l_ivector_slice &UncheckedSetSup(l_ivector_slice &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE l_ivector &UncheckedSetSup(l_ivector &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new unchecked given supremum vector
	INLINE l_ivector_slice &UncheckedSetSup(l_ivector_slice &iv,const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given supremum value
	INLINE l_ivector &SetSup(l_ivector &iv,const l_real &r) throw();
	//! Returns the vector with the new given infimum value
	INLINE l_ivector &SetInf(l_ivector &iv,const l_real &r) throw();
	//! Returns the vector with the new unchecked given supremum value
	INLINE l_ivector &UncheckedSetSup(l_ivector &iv,const l_real &r) throw();
	//! Returns the vector with the new unchecked given infimum value
	INLINE l_ivector &SetUncheckedInf(l_ivector &iv,const l_real &r) throw();

	//! Returns the vector with the new given supremum value
	INLINE l_ivector_slice &SetSup(l_ivector_slice &iv,const l_real &r) throw();
	//! Returns the vector with the new given infimum value
	INLINE l_ivector_slice &SetInf(l_ivector_slice &iv,const l_real &r) throw();
	//! Returns the vector with the new unchecked given supremum value
	INLINE l_ivector_slice &UncheckedSetSup(l_ivector_slice &iv,const l_real &r) throw();
	//! Returns the vector with the new unchecked given infimum value
	INLINE l_ivector_slice &SetUncheckedInf(l_ivector_slice &iv,const l_real &r) throw();

	//! Resizes the vector
	INLINE void Resize(l_ivector &rv) throw();
	//! Resizes the vector
	INLINE void Resize(l_ivector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<l_ivector>);
#else
	throw();
#endif
	//! Resizes the vector
	INLINE void Resize(l_ivector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<l_ivector>);
#else
	throw();
#endif

	//! Returns the absolute value of the vector
	INLINE l_ivector abs(const l_ivector &rv) throw();
	//! Returns the absolute value of the vector
	INLINE l_ivector abs(const l_ivector_slice &sl) throw();
        //! Returns the rounded diameter of the vector
	INLINE l_rvector diam(const l_ivector &v) throw();
        //! Returns the rounded diameter of the vector
	INLINE l_rvector diam(const l_ivector_slice &v) throw();
	//! Returns the rounded middle of the vector
	INLINE l_rvector mid(const l_ivector &v) throw();
	//! Returns the rounded middle of the vector
	INLINE l_rvector mid(const l_ivector_slice &v) throw();
	//! Returns the infimum of the vector
	INLINE l_rvector Inf(const l_ivector &v) throw();
	//! Returns the infimum of the vector
	INLINE l_rvector Inf(const l_ivector_slice &v) throw();
	//! Returns the supremum of the vector
	INLINE l_rvector Sup(const l_ivector &v) throw();
	//! Returns the supremum of the vector
	INLINE l_rvector Sup(const l_ivector_slice &v) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const l_ivector &rv) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const l_ivector_slice &sl) throw();

//======================= Vector / Scalar ===============================

//----------------------------- l_interval ---------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &rv, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &sl, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const l_ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const l_ivector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &rv,const l_interval &r) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector &rv, const l_interval &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector_slice &sl, const l_interval &s) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector &operator /=(l_ivector &rv,const l_interval &r) throw();

//---------------------------- Real --------------------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &rv, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &sl, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const real &s, const l_ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const real &s, const l_ivector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &rv,const real &r) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector &rv, const real &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector_slice &sl, const real &s) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector &operator /=(l_ivector &rv,const real &r) throw();

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const rvector &rv, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const rvector_slice &sl, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const rvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const rvector_slice &sl) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const rvector &rv, const l_interval &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const rvector_slice &sl, const l_interval &s) throw();

//---------------------------- Complex --------------------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &rv, const l_real &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &sl, const l_real &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_real &s, const l_ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_real &s, const l_ivector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &rv,const l_real &r) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector &rv, const l_real &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector_slice &sl, const l_real &s) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector &operator /=(l_ivector &rv,const l_real &r) throw();

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rvector &rv, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rvector_slice &sl, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const l_rvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const l_rvector_slice &sl) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const l_rvector &rv, const l_interval &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const l_rvector_slice &sl, const l_interval &s) throw();

//---------------------------- interval --------------------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &rv, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &sl, const interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const interval &s, const l_ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const interval &s, const l_ivector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &rv,const interval &r) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector &rv, const interval &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const l_ivector_slice &sl, const interval &s) throw();
	//! Implementation of division and allocation operation
	INLINE l_ivector &operator /=(l_ivector &rv,const interval &r) throw();

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector &rv, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector_slice &sl, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const ivector &rv) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const ivector_slice &sl) throw();

	//! Implementation of division operation
	INLINE l_ivector operator /(const ivector &rv, const l_interval &s) throw();
	//! Implementation of division operation
	INLINE l_ivector operator /(const ivector_slice &sl, const l_interval &s) throw();

//======================= Vector / Vector ===============================


	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s, const l_ivector &rv) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &o, const l_ivector_slice &sl) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, l_ivector &rv) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, l_ivector_slice &rv) throw();
	
//----------------------- l_interval / l_interval ---------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of positive sign operation
	INLINE const l_ivector &operator +(const l_ivector &rv) throw();
	//! Implementation of positive sign operation
	INLINE l_ivector operator +(const l_ivector_slice &sl) throw();

	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector & operator +=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector &operator +=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of negative sign operation
	INLINE l_ivector operator -(const l_ivector &rv) throw();
	//! Implementation of negative sign operation
	INLINE l_ivector operator -(const l_ivector_slice &sl) throw();
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector & operator -=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector &operator -=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector & operator |=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector &operator |=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector & operator &=(l_ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector &operator &=(l_ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_ivector &rv1, const l_ivector &rv2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_ivector_slice &sl1, const l_ivector_slice &sl2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_ivector_slice &sl, const l_ivector &rv) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_ivector &rv, const l_ivector_slice &sl) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_ivector &rv1, const l_ivector &rv2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_ivector_slice &sl1, const l_ivector_slice &sl2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_ivector_slice &sl, const l_ivector &rv) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_ivector &rv, const l_ivector_slice &sl) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_ivector &rv1, const l_ivector &rv2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_ivector_slice &sl1, const l_ivector_slice &sl2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator < (const l_ivector_slice &sl, const l_ivector &rv) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator < (const l_ivector &rv, const l_ivector_slice &sl) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_ivector &rv1, const l_ivector &rv2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_ivector_slice &sl1, const l_ivector_slice &sl2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_ivector_slice &sl, const l_ivector &rv) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_ivector &rv, const l_ivector_slice &sl) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_ivector &rv1, const l_ivector &rv2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_ivector_slice &sl1, const l_ivector_slice &sl2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_ivector_slice &sl, const l_ivector &rv) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_ivector &rv, const l_ivector_slice &sl) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_ivector &rv1, const l_ivector &rv2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_ivector_slice &sl1, const l_ivector_slice &sl2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_ivector_slice &sl, const l_ivector &rv) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_ivector &rv, const l_ivector_slice &sl) throw();

//-------------------------------- l_interval / Real --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const l_ivector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const l_ivector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_interval operator *(const rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE l_ivector operator +(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE l_ivector & operator +=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector &operator +=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE l_ivector & operator -=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector &operator -=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector & operator |=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector &operator |=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector & operator &=(l_ivector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector &operator &=(l_ivector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
//-------------------------------- l_interval / l_real --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const l_ivector_slice &sl,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const l_ivector &rv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_rvector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice & sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE l_ivector & operator +=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector &operator +=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE l_ivector & operator -=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector &operator -=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector & operator |=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector &operator |=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector & operator &=(l_ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector &operator &=(l_ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

//-------------------------------- l_interval / interval --------------------------------

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const l_ivector_slice &sl,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp,const l_ivector &rv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const ivector_slice & sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE l_ivector & operator +=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_ivector &operator +=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction and allocation operation
	INLINE l_ivector & operator -=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_ivector &operator -=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector & operator |=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_ivector &operator |=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector &rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector &rv, const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector_slice &sl, const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector_slice &sl1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_ivector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector & operator &=(l_ivector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_ivector &operator &=(l_ivector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

//------------- real x l_real ------------------------
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector_slice &sl,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector &rv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

//------------- l_real x l_real ------------------------
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

//-------------------------------- interval / l_real --------------------------------

	// multiplication in lrvecivec.hpp
	
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif


	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif


	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_ivector operator |(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const l_rvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector &rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector &rv, const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector_slice &sl, const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_ivector operator &(const ivector_slice &sl1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>);
#else
	throw();
#endif

} // namespace cxsc 

#ifdef _CXSC_INCL_INL
#include "vector.inl"
#include "l_ivector.inl"
#endif

#ifdef _CXSC_RMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "livecrmat.inl"
# else
#  include "livecrmat.hpp"
# endif
#endif

#ifdef _CXSC_LRMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "liveclrmat.inl"
# else
#  include "liveclrmat.hpp"
# endif
#endif

#ifdef _CXSC_IMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "livecimat.inl"
# else
#  include "livecimat.hpp"
# endif
#endif


#endif

