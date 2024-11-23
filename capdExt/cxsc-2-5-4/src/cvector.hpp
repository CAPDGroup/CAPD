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

/* CVS $Id: cvector.hpp,v 1.31 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CVECTOR_HPP_INCLUDED
#define _CXSC_CVECTOR_HPP_INCLUDED

#include "xscclass.hpp"
#include "except.hpp"
#include "cdot.hpp"
#include "cidot.hpp"
#include "complex.hpp" // used for declaration of Inf, Sup,...
//#include "cxscmatr.hpp"
#include "rvector.hpp"
#include "vector.hpp"


#include <iostream>

//#include "matrix.hpp" // hat hier eigentlich nichts zu suchen, sonst aber Internal Compiler Error #9

namespace cxsc {

class cvector_slice;
class scvector;
class scvector_slice;
class srvector;
class srvector_slice;

//! The Data Type cvector
/*!
The vectors of C-XSC are one dimensional arrays of the corresponding scalar base type.

\sa rvector
*/
class cvector
{
	friend class cvector_slice;
	friend class cmatrix;
	friend class cmatrix_subv;
	friend class civector;
	friend class cimatrix;
	private:
	complex *dat;
	int l,u,size;

	public:
        double* to_blas_array() const { return (double*)dat; }
//#if(CXSC_INDEX_CHECK)
#ifdef _CXSC_FRIEND_TPL
	//------------ Templates --------------------------------------------------
	// complex
template <class V,class MS,class S> friend void _vmsconstr(V &v,const MS &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>);
#else
	throw();
#endif
template <class V,class M,class S> friend void _vmconstr(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
 template <class V> friend 	void _vresize(V &rv) throw();
 template <class V,class S> friend 	void _vresize(V &rv, const int &len)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>);
#else
	throw();
#endif
 template <class V,class S> friend 	void _vresize(V &rv, const int &lb, const int &ub)
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
 template <class V> friend 	 V _vconj(const V &rv) throw();
 template <class VS,class E> friend 	 E _vsconj(const VS &sl) throw();
 template <class V,class E> friend 	 E _vabs(const V &rv) throw();
 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();
template <class MV,class V> friend  V _mvabs(const MV &mv) throw();
 template <class V,class E> friend 	 E _vim(const V &rv) throw();
 template <class V,class E> friend 	 E _vre(const V &rv) throw();
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
 template <class V,class S> friend 	 V &_vssetre(V &v, const S &s) throw();
 template <class V,class S> friend 	 V &_vssetim(V &v, const S &s) throw();

//-------- vector-vector -----------------------
 template <class DP,class V1,class V2> friend 	void _vvaccu(DP &dp, const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
 template <class DP,class VS,class V> friend 	void _vsvaccu(DP &dp, const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvcmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>);
#else
	throw();
#endif
 template <class VS,class V,class E> friend 	 E _vsvcmult(const VS & sl, const V &rv)
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
 template <class V,class MV,class S> friend 	 S _vmvcmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
 template <class V,class MV,class S> friend 	 S _vmvcimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
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
template <class DP,class V,class SV> friend 	void _vmvaccu(DP &dp, const V & rv1, const SV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
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
 template <class V> friend 	void *_vvoid(const V &rv) throw();
 template <class VS1,class VS2> friend 	 bool _vsvseq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsneq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsless(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS1,class VS2> friend 	 bool _vsvsleq(const VS1 &sl1, const VS2 &sl2) throw();
 template <class VS> friend 	 bool _vsnot(const VS &sl) throw();
 template <class VS> friend 	void *_vsvoid(const VS &sl) throw();
 template <class V> friend 	std::ostream &_vout(std::ostream &s, const V &rv) throw();
 template <class V> friend 	std::istream &_vin(std::istream &s, V &rv) throw();

	//------------- vector-matrix ---------------
template <class V,class MV2,class S> friend  V &_vmvassign(V &v,const MV2 &rv) throw();
 template <class M,class V,class E> friend 	 E _mvcmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class M,class V,class E> friend 	 E _mvcimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmcmult(const V &v,const M &m)
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
 template <class V,class M,class S> friend 	 V &_vmcmultassign(V &v,const M &m)
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
 template <class MS,class V,class E> friend 	 E _msvcmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class MS,class V,class E> friend 	 E _msvcimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmscmult(const V &v,const MS &ms)
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
 template <class V,class MS,class S> friend 	 V &_vmscmultassign(V &v,const MS &ms)
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
	//--Real -------- vector-scalar ------------
	//--Real--------- Vector-vector---------
	//-- Real -------- Vector-matrix ----------
	// interval -----------------
	// vector-scalar
	// vector-vector
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
 template <class VS1,class VS2,class E> friend 	 E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
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

	// vector-matrix
	// cinterval -----------------
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
template <class MV,class V> friend  MV &_mvvsetinf(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvsetsup(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvusetinf(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvusetsup(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
 template <class V,class E> friend 	 E _vmid(const V &rv) throw();
 template <class V,class E> friend 	 E _vinf(const V &rv) throw();
 template <class V,class E> friend 	 E _vsup(const V &rv) throw();
 template <class V,class E> friend 	 E _vdiam(const V &rv) throw();
 template <class VS,class E> friend 	 E _vsmid(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsinf(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vssup(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsdiam(const VS &sl) throw();
template <class MV,class V> friend  V _mvdiam(const MV &mv) throw();
template <class MV,class V> friend  V _mvmid(const MV &mv) throw();
template <class MV,class V> friend  V _mvinf(const MV &mv) throw();
template <class MV,class V> friend  V _mvsup(const MV &mv) throw();

	// vector-vector
 template <class V1,class V2> friend 	 V1 &_vvconvassign(V1 &rv1, const V2 &rv2)
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
 template <class VS,class V> friend 	 VS &_vsvconvassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvsectassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvconvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvsectassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif


#endif

	//------ Konstruktoren ----------------------------------------------------
	//! Constructor of class cvector
	cvector () throw();
	//! Constructor of class cvector
	explicit cvector(const int &i) throw();
#ifdef OLD_CXSC
	//! Constructor of class cvector
	explicit cvector(const class index &i) throw(); // for backwards compatibility
#endif
	//! Constructor of class cvector
	explicit cvector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_WRONG_BOUNDARIES,ERROR_CVECTOR_NO_MORE_MEMORY);
#else
	throw();
#endif
	//! Constructor of class cvector
	cvector(const cmatrix_subv &) throw();
	//! Constructor of class cvector
	explicit cvector(const complex& r) throw();
	//! Constructor of class cvector
	explicit cvector(const cmatrix& )
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class cvector
	explicit cvector(const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class cvector
	cvector(const cvector_slice &rs) throw();
	//! Constructor of class cvector
	cvector(const cvector &v) throw();
	//! Constructor of class cvector
	cvector(const scvector_slice &rs);
	//! Constructor of class cvector
	cvector(const scvector &v);
	// Real
	//! Constructor of class cvector
	explicit cvector(const srvector_slice &rs);
	//! Constructor of class cvector
	explicit cvector(const srvector &v);
	//! Constructor of class cvector
	explicit cvector(const real &) throw();
	//! Constructor of class cvector
	explicit cvector(const rvector_slice &rs) throw();
	//! Constructor of class cvector
	explicit cvector(const rvector &v) throw();
	//! Constructor of class cvector
	explicit cvector(const rmatrix &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class cvector
	explicit cvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Constructor of class cvector
	explicit cvector(const rmatrix_subv &) throw();
	
	// complex
	//! Implementation of standard assigning operator
	cvector &operator =(const cvector &rv) throw();
	//! Implementation of standard assigning operator
	cvector &operator =(const cvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	cvector &operator =(const scvector &rv);
	//! Implementation of standard assigning operator
	cvector &operator =(const scvector_slice &sl);
	//! Implementation of standard assigning operator
	cvector &operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	cvector &operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	cvector &operator =(const cmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	cvector &operator =(const cmatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	cvector &operator =(const rvector &rv) throw();
	//! Implementation of standard assigning operator
	cvector &operator =(const rvector_slice &sl) throw();
	//! Implementation of standard assigning operator
	cvector &operator =(const srvector &rv);
	//! Implementation of standard assigning operator
	cvector &operator =(const srvector_slice &sl);
	//! Implementation of standard assigning operator
	cvector &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	cvector &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	cvector &operator =(const rmatrix_slice &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	cvector &operator =(const rmatrix_subv &) throw();

        cvector& operator+=(const srvector&);
        cvector& operator+=(const scvector&);
        cvector& operator+=(const srvector_slice&);
        cvector& operator+=(const scvector_slice&);
        cvector& operator-=(const srvector&);
        cvector& operator-=(const scvector&);
        cvector& operator-=(const srvector_slice&);
        cvector& operator-=(const scvector_slice&);

        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE cvector operator()(const intvector& p);
        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE cvector operator()(const intmatrix& P);

	//--------- Destruktor ----------------------------------------------------
	INLINE ~cvector() { delete [] dat; }

	//------ Standardfunktionen -----------------------------------------------
	
	friend INLINE complex::complex(const cvector &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const cvector &rv) throw() { return rv.l; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const cvector &rv) throw() { return rv.u; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const cvector &rv) throw() { return rv.size; }
	//! Sets the lower bound of the vector
	friend INLINE cvector & SetLb(cvector &rv, const int &l) throw() { rv.l=l; rv.u=l+rv.size-1; return rv;}
	//! Sets the upper bound of the vector
	friend INLINE cvector & SetUb(cvector &rv, const int &u) throw() { rv.u=u; rv.l=u-rv.size+1; return rv;}
	//! Operator for accessing the single elements of the vector (read-only)
	INLINE complex & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector
	INLINE complex & operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the whole vector
	INLINE cvector & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE cvector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	cvector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	operator void*() throw();
//#else
//#endif
};


//! The Data Type cvector_slice
/*!
This data type represents a partial cvector.

\sa cvector
*/
class cvector_slice
{
	friend class cvector;
	friend class cmatrix;
	friend class civector;
	friend class cimatrix;
	private:
	complex *dat;
	int l,u,size;
	int start,end;

	public:
//#if(CXSC_INDEX_CHECK)	
#ifdef _CXSC_FRIEND_TPL
//------------------------- Templates -------------------------------------------
// complex / complex

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

 template <class DP,class VS,class V> friend 	void _vsvaccu(DP &dp, const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
 template <class DP,class VS1,class VS2> friend 	void _vsvsaccu(DP &dp, const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

 template <class VS,class S,class E> friend 	 E _vssdiv(const VS &sl, const S &s) throw();
 template <class VS,class S,class E> friend 	 E _vssmult(const VS &sl, const S &s) throw();

 template <class VS,class V,class E> friend 	 E _vsvcmult(const VS & sl, const V &rv)
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
 template <class VS,class V> friend 	 bool _vsveq(const VS &sl, const V &rv) throw();
 template <class VS,class V> friend 	 bool _vsvneq(const VS &sl, const V &rv) throw();
 template <class VS,class V> friend 	 bool _vsvless(const VS &sl, const V &rv) throw();
 template <class VS,class V> friend 	 bool _vsvleq(const VS &sl, const V &rv) throw();
 template <class V,class VS> friend 	 bool _vvsless(const V &rv, const VS &sl) throw();
 template <class V,class VS> friend 	 bool _vvsleq(const V &rv, const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsconj(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsabs(const VS &sl) throw();

 template <class VS1,class VS2,class E> friend 	 E _vsvscmult(const VS1 & sl1, const VS2 &sl2)
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
 template <class VS> friend 	void *_vsvoid(const VS &sl) throw();
 template <class V> friend 	std::ostream &_vsout(std::ostream &s, const V &rv) throw();
 template <class V> friend 	std::istream &_vsin(std::istream &s, V &rv) throw();
 template <class VS,class E> friend 	 E _vsim(const VS &sl) throw();
 template <class VS,class E> friend 	 E _vsre(const VS &sl) throw();
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
/*	friend TINLINE cvector_slice &_vsmvsetim(cvector_slice &,const
rmatrix_subv &) #if(CXSC_INDEX_CHECK)
throw(ERROR__OP_WITH_WRONG_DIM<cvector_slice>); #else 	throw();
#endif
	friend TINLINE cvector_slice &_vsmvsetre(cvector_slice &,const rmatrix_subv &)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector_slice>);
#else
	throw();
#endif  */     // 4.10.00 S.W.

 template <class VS,class S> friend 	 VS &_vsssetim(VS &vs, const S &s) throw();
 template <class VS,class S> friend 	 VS &_vsssetre(VS &vs, const S &s) throw();

 template <class VS,class M,class S> friend 	 VS &_vsmcmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
	
	// complex / Real
 template <class VS,class V> friend 	 VS &_vsvplusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class VS,class V> friend 	 VS &_vsvminusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmscmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
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
	// interval -----------
	// vector-vector -------
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
 template <class VS1,class VS2,class E> friend 	 E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class V1,class V2,class E> friend 	 E _vvplus(const V1 &rv1, const V2 &rv2)
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

	// vector-matrix -------
/*   friend TINLINE civector _mvscimult<imatrix,cvector_slice,civector>(const
imatrix &m,const cvector_slice &v) #if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif */  // 4.10.00. S.W.
/*   friend TINLINE civector _vsmcimult<cvector_slice,imatrix,civector>(const
cvector_slice &v,const imatrix &m) #if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
*/  // 4.10.00 S.W.

	// cinterval
	// cinterval -- vector-vector
 template <class V,class VS> friend 	 V &_vvsconvassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class V,class VS> friend 	 V &_vvssectassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvsconvassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
 template <class VS1,class VS2> friend 	 VS1 &_vsvssectassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>);
#else
	throw();
#endif
#endif
	
	//--------------------- Konstruktoren -----------------------------------
	//! Constructor of class cvector_slice
	explicit INLINE cvector_slice(cvector &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	//! Constructor of class cvector_slice
	explicit INLINE cvector_slice(cvector_slice &a, const int &lb, const int &ub) throw():dat(a.dat),l(a.l),u(a.u),size(ub-lb+1),start(lb),end(ub) { }
	public:
	//! Constructor of class cvector_slice
	INLINE cvector_slice(const cvector_slice &a) throw():dat(a.dat),l(a.l),u(a.u),size(a.size),start(a.start),end(a.end) { }
	public:
	// complex
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const scvector &sl);
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const scvector_slice &sl);

	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice &operator =(const cmatrix_subv &) throw();
	// Real
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const srvector &rv);
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const srvector_slice &rv);

	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	INLINE cvector_slice & operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice & operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cvector_slice &operator =(const rmatrix_subv &mv) throw();

	// cinterval --------
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

	// cinterval -- vector-vector
 template <class V,class MS,class E> friend 	 E _vmscimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif


	//--------------------- Standardfunktionen ------------------------------

	friend INLINE complex::complex(const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CVECTOR_USE_OF_UNINITIALIZED_OBJ);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const cvector_slice &sl) throw() { return sl.start; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const cvector_slice &sl) throw() { return sl.end; }
	//! Returns the dimension of the vector
        friend INLINE int VecLen(const cvector_slice &sl) throw() { return sl.end-sl.start+1; }
	//! Operator for accessing the single elements of the vector (read-only)
	INLINE complex & operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector
	INLINE complex & operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	
	//! Operator for accessing the whole vector
	INLINE cvector_slice & operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE cvector_slice operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE cvector_slice operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	//! Implementation of division and allocation operation
	INLINE cvector_slice &operator /=(const complex &r) throw();
	//! Implementation of division and allocation operation
	INLINE cvector_slice &operator /=(const real &r) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cvector_slice &operator *=(const complex &r) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cvector_slice &operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cvector_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cvector_slice &operator *=(const real &r) throw();
	//! Implementation of addition and allocation operation
	INLINE cvector_slice &operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cvector_slice &operator +=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector_slice &operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector_slice &operator -=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cvector_slice &operator |=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cvector_slice &operator |=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cvector_slice &operator &=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cvector_slice &operator &=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	
	//! Implementation of addition and allocation operation
	INLINE cvector_slice &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cvector_slice &operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector_slice &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector_slice &operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cvector_slice &operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cvector_slice &operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cvector_slice &operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cvector_slice &operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	INLINE operator void*() throw();

        cvector_slice& operator+=(const srvector&);
        cvector_slice& operator+=(const scvector&);
        cvector_slice& operator+=(const srvector_slice&);
        cvector_slice& operator+=(const scvector_slice&);
        cvector_slice& operator-=(const srvector&);
        cvector_slice& operator-=(const scvector&);
        cvector_slice& operator-=(const srvector_slice&);
        cvector_slice& operator-=(const scvector_slice&);

//#else
//#endif
};

//=======================================================================
//======================== Vector Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cvector _cvector(const complex &r) throw();
//	INLINE cvector _cvector(const cmatrix &m) throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE cvector _cvector(const cmatrix_slice &sl) throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cvector _cvector(const real &r) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cvector _cvector(const rvector_slice &rs) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cvector _cvector(const rvector &rs) throw();
//	INLINE cvector _cvector(const rmatrix &m) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
//	INLINE cvector _cvector(const rmatrix_slice &sl) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ);
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cvector _cvector(const rmatrix_subv &rs) throw();

	//! Returns the vector with the new given real part vector
	INLINE cvector &SetRe(cvector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given real part vector
	INLINE cvector_slice &SetRe(cvector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given real part vector
	INLINE cvector &SetRe(cvector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given real part vector
	INLINE cvector_slice &SetRe(cvector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the vector with the new given imaginary part vector
	INLINE cvector &SetIm(cvector &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given imaginary part vector
	INLINE cvector_slice &SetIm(cvector_slice &iv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given imaginary part vector
	INLINE cvector &SetIm(cvector &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the vector with the new given imaginary part vector
	INLINE cvector_slice &SetIm(cvector_slice &iv,const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Sets componentwise the real parts of the cvector
	INLINE cvector &SetRe(cvector &iv,const real &r) throw();
	//! Sets componentwise the imaginary parts of the cvector
	INLINE cvector &SetIm(cvector &iv,const real &r) throw();
	//! Sets componentwise the real parts of the cvector
	INLINE cvector_slice &SetRe(cvector_slice &iv,const real &r) throw();
	//! Sets componentwise the imaginary parts of the cvector
	INLINE cvector_slice &SetIm(cvector_slice &iv,const real &r) throw();

	//! Resizes the vector
	INLINE void Resize(cvector &rv) throw();
	//! Resizes the vector
	INLINE void Resize(cvector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<cvector>);
#else
	throw();
#endif
	//! Resizes the vector
	INLINE void Resize(cvector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<cvector>);
#else
	throw();
#endif
	
	//! Returns the conjugated cvector
	INLINE cvector conj(const cvector &rv) throw();
	//! Returns the conjugated cvector
	INLINE cvector conj(const cvector_slice &sl) throw();
	
	//! Returns the absolute value of the vector
	INLINE rvector abs(const cvector &rv) throw();
	//! Returns the absolute value of the vector
	INLINE rvector abs(const cvector_slice &sl) throw();
	//! Returns the imaginary part of the vector
	INLINE rvector Im(const cvector &v) throw();
	//! Returns the imaginary part of the vector
	INLINE rvector Im(const cvector_slice &v) throw();
	//! Returns the real part of the cvector
	INLINE rvector Re(const cvector &v) throw();
	//! Returns the real part of the cvector
	INLINE rvector Re(const cvector_slice &v) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const cvector &rv) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const cvector_slice &sl) throw();

//======================= Vector / Scalar ===============================

//----------------------------- complex ---------------------------

	//! Implementation of multiplication operation
	INLINE cvector operator *(const cvector &rv, const complex &s) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const cvector_slice &sl, const complex &s) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const complex &s, const cvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const complex &s, const cvector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cvector &operator *=(cvector &rv,const complex &r) throw();

	//! Implementation of division operation
	INLINE cvector operator /(const cvector &rv, const complex &s) throw();
	//! Implementation of division operation
	INLINE cvector operator /(const cvector_slice &sl, const complex &s) throw();
	//! Implementation of division and allocation operation
	INLINE cvector &operator /=(cvector &rv,const complex &r) throw();

//---------------------------- Real --------------------------------------

	//! Implementation of multiplication operation
	INLINE cvector operator *(const cvector &rv, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const cvector_slice &sl, const real &s) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const real &s, const cvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const real &s, const cvector_slice &sl) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cvector &operator *=(cvector &rv,const real &r) throw();

	//! Implementation of division operation
	INLINE cvector operator /(const cvector &rv, const real &s) throw();
	//! Implementation of division operation
	INLINE cvector operator /(const cvector_slice &sl, const real &s) throw();
	//! Implementation of division and allocation operation
	INLINE cvector &operator /=(cvector &rv,const real &r) throw();

	//! Implementation of multiplication operation
	INLINE cvector operator *(const rvector &rv, const complex &s) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const rvector_slice &sl, const complex &s) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const complex &s, const rvector &rv) throw();
	//! Implementation of multiplication operation
	INLINE cvector operator *(const complex &s, const rvector_slice &sl) throw();

	//! Implementation of division operation
	INLINE cvector operator /(const rvector &rv, const complex &s) throw();
	//! Implementation of division operation
	INLINE cvector operator /(const rvector_slice &sl, const complex &s) throw();

//======================= Vector / Vector ===============================


	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s, const cvector &rv) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &o, const cvector_slice &sl) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, cvector &rv) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s, cvector_slice &rv) throw();
	
//----------------------- complex / complex ---------------------------

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(cdotprecision &dp, const cvector &);

	//! The accurate sum of the elements of the vector added to the first argument
	void accumulate(cdotprecision &dp, const rvector &);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector & rv1, const cvector &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector_slice & sl, const cvector &rv);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector &rv, const cvector_slice &sl);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector & rv1, const cmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cmatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cmatrix_subv & rv1, const cvector &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cmatrix_subv & rv1, const cmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector_slice & sl1, const cvector_slice &sl2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif


	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	
	//! Implementation of positive sign operation
	INLINE const cvector &operator +(const cvector &rv) throw();
	//! Implementation of positive sign operation
	INLINE cvector operator +(const cvector_slice &sl) throw();

	//! Implementation of addition operation
	INLINE cvector operator +(const cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const cvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cvector & operator +=(cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cvector &operator +=(cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

	//! Implementation of negative sign operation
	INLINE cvector operator -(const cvector &rv) throw();
	//! Implementation of negative sign operation
	INLINE cvector operator -(const cvector_slice &sl) throw();
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector & operator -=(cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector &operator -=(cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

	//! Implementation of standard equality operation
	INLINE bool operator ==(const cvector &rv1, const cvector &rv2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const cvector_slice &sl1, const cvector_slice &sl2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const cvector_slice &sl, const cvector &rv) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const cvector &rv, const cvector_slice &sl) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cvector &rv1, const cvector &rv2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cvector_slice &sl1, const cvector_slice &sl2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cvector_slice &sl, const cvector &rv) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cvector &rv, const cvector_slice &sl) throw();
/*	INLINE bool operator <(const cvector &rv1, const cvector &rv2) throw();
	INLINE bool operator <(const cvector_slice &sl1, const cvector_slice &sl2) throw();
	INLINE bool operator < (const cvector_slice &sl, const cvector &rv) throw();
	INLINE bool operator < (const cvector &rv, const cvector_slice &sl) throw();
	INLINE bool operator <=(const cvector &rv1, const cvector &rv2) throw();
	INLINE bool operator <=(const cvector_slice &sl1, const cvector_slice &sl2) throw();
	INLINE bool operator <=(const cvector_slice &sl, const cvector &rv) throw();
	INLINE bool operator <=(const cvector &rv, const cvector_slice &sl) throw();
	INLINE bool operator >(const cvector &rv1, const cvector &rv2) throw();
	INLINE bool operator >(const cvector_slice &sl1, const cvector_slice &sl2) throw();
	INLINE bool operator >(const cvector_slice &sl, const cvector &rv) throw();
	INLINE bool operator >(const cvector &rv, const cvector_slice &sl) throw();
	INLINE bool operator >=(const cvector &rv1, const cvector &rv2) throw();
	INLINE bool operator >=(const cvector_slice &sl1, const cvector_slice &sl2) throw();
	INLINE bool operator >=(const cvector_slice &sl, const cvector &rv) throw();
	INLINE bool operator >=(const cvector &rv, const cvector_slice &sl) throw();
*/
//-------------------------------- complex / Real --------------------------------


	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector & rv1, const rvector &rv2);


	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rvector & rv1, const cvector &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rvector_slice & sl, const cvector &rv);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp,const cvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp,const cvector_slice &sl,const rvector &rv);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rvector &rv, const cvector_slice &sl);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rvector & rv1, const cmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rvector_slice & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rvector_slice & rv1, const cmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp,const cvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp,const cvector &rv,const rvector_slice &sl);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cmatrix_subv & rv1, const rvector &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cmatrix_subv & rv1, const rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cmatrix_subv & rv1, const rvector_slice &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cvector_slice & sl1, const rvector_slice &sl2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rvector_slice & sl1, const cvector_slice &sl2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const rmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const rmatrix_subv & rv1, const cmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cdotprecision &dp, const cmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument (without error bound)
	void accumulate_approx(cdotprecision &dp, const cmatrix_subv & rv1, const rmatrix_subv &rv2);

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rvector_slice & sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp,const cvector_slice &sl,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rvector_slice & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp,const cvector &rv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const rvector_slice &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice &, const ivector &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const ivector & sl1, const cvector_slice &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector &, const ivector &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const ivector &, const cvector &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const ivector_slice &, const cvector &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector &, const ivector_slice &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cvector_slice &, const ivector_slice &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const ivector_slice &, const cvector_slice &)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE complex operator *(const rvector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const rvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const rvector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE complex operator *(const cvector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	
	//! Implementation of addition operation
	INLINE cvector operator +(const rvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const rvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const rvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE cvector operator +(const cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cvector operator +(const cvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE cvector & operator +=(cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cvector &operator +=(cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE cvector operator -(const rvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const rvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const rvector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const rvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cvector operator -(const cvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector & operator -=(cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cvector &operator -=(cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>);
#else
	throw();
#endif

} // namespace cxsc 

#ifdef _CXSC_INCL_INL
#include "vector.inl"
#include "cvector.inl"
#endif
		
#ifdef _CXSC_RMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "cvecrmat.inl"
# else
#  include "cvecrmat.hpp"
# endif
#endif

#ifdef _CXSC_IMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "cvecimat.inl"
# else
#  include "cvecimat.hpp"
# endif
#endif

#ifdef _CXSC_IVECTOR_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "iveccvec.inl"
# else
#  include "iveccvec.hpp"
# endif
#endif

#ifdef CXSC_USE_BLAS
#define _CXSC_BLAS_CVECTOR
#include "cxsc_blas.inl"
#endif


#endif
