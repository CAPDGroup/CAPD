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

/* CVS $Id: l_rmatrix.hpp,v 1.20 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_LRMATRIX_HPP_INCLUDED
#define _CXSC_LRMATRIX_HPP_INCLUDED

#include "xscclass.hpp"
#include "dot.hpp"
#include "idot.hpp"
#include "l_rvector.hpp"
#include "except.hpp"
#include "matrix.hpp"
#include "rmatrix.hpp"

namespace cxsc {

class l_rmatrix;
class l_rmatrix_slice;

//! The Multiple-Precision Data Type l_rmatrix_subv
/*!
This Multiple-Precision Data Type provides one column or row of a matrix as a vector.
*/
class l_rmatrix_subv
{
	friend class l_rvector;
	friend class l_ivector;
	friend class l_rmatrix;
	friend class l_rmatrix_slice;
	private:
	l_real *dat;
	int lb,ub;
	int size,start,offset; // start=first element index 0..n-1
	
	public:
	//! Returns one row of the matrix as a vector
	friend INLINE l_rmatrix_subv Row(l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one column of the matrix as a vector
	friend INLINE l_rmatrix_subv Col(l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one row of the matrix as a vector
	friend INLINE l_rmatrix_subv Row(const l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one column of the matrix as a vector
	friend INLINE l_rmatrix_subv Col(const l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
//#if(CXSC_INDEX_CHECK)

template<typename S, typename T>
friend inline void addDot(const S &x, const T &y, dotprecision &val);


#ifdef _CXSC_FRIEND_TPL
	//----------------- Templates ---------------------------------------
template <class MV1,class MV2> friend  MV1 &_mvmvassign(MV1 &v,const MV2 &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif
template <class MV,class S> friend  MV &_mvsassign(MV &v,const  S &r) throw();
template <class MV,class V> friend  MV &_mvvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class V,class MV2,class S> friend  V &_vmvassign(V &v,const MV2 &rv) throw();
template <class MV,class V> friend  V _mvabs(const MV &mv) throw();
template <class DP,class V,class SV> friend 	 void _vmvaccu(DP &dp, const V & rv1, const SV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
template <class DP,class MV1,class MV2> friend 	 void _mvmvaccu(DP &dp, const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif


 template <class MV1,class MV2,class S> friend 	 S _mvmvlmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif
 template <class V,class MV,class S> friend 	 S _vmvlmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
 template <class MV,class S,class E> friend 	 E _mvsmult(const MV &rv, const S &s) throw();
 template <class MV1,class MV2,class E> friend 	 E _mvmvplus(const MV1 &rv1, const MV2 &rv2)
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
 template <class MV,class S,class E> friend 	 E _mvsdiv(const MV &rv, const S &s) throw();
template <class MV,class S> friend  MV &_mvsmultassign(MV &v,const S &r) throw();
template <class MV, class S> friend  MV &_mvsplusassign(MV &v,const S &r) throw();
template <class MV,class S> friend  MV &_mvsminusassign(MV &v,const S &r) throw();
template <class MV,class S> friend  MV &_mvsdivassign(MV &v,const S &r) throw();
template <class MV,class V> friend  MV &_mvvplusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvplusassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif
template <class MV,class V> friend  MV &_mvvminusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif
template <class V,class MV> friend  V &_vmvminusassign(V &rv,const MV &v)
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

#endif

	//----------------- Konstruktoren ----------------------------------

	//! Constructor of class l_rmatrix_subv
	explicit INLINE l_rmatrix_subv (l_real *d, const int &l, const int &u, const int &s, const int &st, const int &o) throw():dat(d),lb(l),ub(u),size(s),start(st),offset(o) { }
	public:
	//! Constructor of class l_rmatrix_subv
	INLINE l_rmatrix_subv(const l_rmatrix_subv &v) throw():dat(v.dat),lb(v.lb),ub(v.ub),size(v.size),start(v.start),offset(v.offset) { }
	public:

	//---------------------- Standardfunktionen ------------------------
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const l_rmatrix_subv &rv) throw() { return rv.lb; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const l_rmatrix_subv &rv) throw() { return rv.ub; }
	//! Operator for accessing the single elements of the vector
	INLINE l_real &operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	//! Operator for accessing the whole vector
	INLINE l_rmatrix_subv &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE l_rmatrix_subv operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE l_rmatrix_subv operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	

	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const l_rmatrix_subv &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_subv &operator *=(const l_real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_subv &operator +=(const l_real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_subv &operator -=(const l_real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_rmatrix_subv &operator /=(const l_real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_subv &operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_subv &operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_subv &operator -=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_subv &operator +=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// Real
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const rmatrix_subv &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_subv &operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_subv &operator *=(const real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_subv &operator +=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_subv &operator -=(const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_rmatrix_subv &operator /=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_subv &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_subv &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_subv &operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_subv &operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
//#else
//#endif	

};


//----------------------- Matrix -----------------------------------------------

class l_rmatrix_slice;

//! The Multiple-Precision Data Type l_rmatrix
/*!
\sa rmatrix
*/
class l_rmatrix
{
	friend class l_rmatrix_slice;
	friend class l_rmatrix_subv;
	friend class l_imatrix;
	private:
	l_real *dat;
	int lb1,ub1,lb2,ub2,xsize,ysize;

	public:
//#if(CXSC_INDEX_CHECK)
#ifdef _CXSC_FRIEND_TPL
	//----------------- Templates ---------------------------------------
template <class S,class M> friend  void _smconstr(S &s,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>,ERROR__USE_OF_UNINITIALIZED_OBJ<M>);
#else
	throw();
#endif
template <class V,class M,class S> friend  void _vmconstr(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
 template <class M1,class M2,class S> friend 	 M1 &_mmassign(M1 &m1,const M2 &m,S ms) throw();
 template <class M,class MS2,class S> friend 	 M &_mmsassign(M &m,const MS2 &ms) throw();
 template <class MS,class M> friend 	 MS &_msmassign(MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class M,class S> friend 	 M &_msassign(M &m,const S &r) throw();
template <class V,class M,class S> friend  V &_vmassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif
template <class M,class V,class S> friend  M &_mvassign(M &m,const V &v) throw();
 template <class M> friend 	 int _mlb(const M &m, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
 template <class M> friend 	 int _mub(const M &m, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
 template <class M> friend 	 M &_msetlb(M &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
 template <class M> friend 	 M &_msetub(M &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
 template <class M> friend 	 void _mresize(M &A) throw();
 template <class M,class S> friend 	 void _mresize(M &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<M>);
#else
	throw();
#endif
 template <class M,class S> friend 	 void _mresize(M &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<M>);
#else
	throw();
#endif
 template <class M,class E> friend 	 E _mabs(const M &m) throw();
 template <class MS,class E> friend 	 E _msabs(const MS &ms) throw();
	//-------------- matrix-matrix -------------
 template <class M1,class M2,class E> friend 	 E _mmplus(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmsplus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class M> friend 	 M _mminus(const M &m) throw();
 template <class MS,class E> friend 	 E _msminus(const MS &ms) throw();
 template <class M1,class M2,class E> friend 	 E _mmminus(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class M2> friend 	 M1 &_mmplusassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS> friend 	 M &_mmsplusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmplusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmsplus(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmsminus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmminus(const MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmsminus(const MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M1,class M2> friend 	 M1 &_mmminusassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS> friend 	 M &_mmsminusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmminusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class M1,class M2,class E> friend 	 E _mmlmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M1,class M2,class S> friend 	 M1 &_mmlmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmslmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmlmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS,class S> friend 	 M &_mmslmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmslmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M1,class M2,class E> friend 	 E _mmconv(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmsconv(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
	//-------- matrix-scalar -----------------
 template <class S,class M,class E> friend 	 E _smmult(const S &c, const M &m) throw();
 template <class M,class S> friend 	 M &_msmultassign(M &m,const S &c) throw();
 template <class S,class MS,class E> friend 	 E _smsmult(const S &c, const MS &ms) throw();
 template <class M,class S,class E> friend 	 E _msdiv(const M &m,const S &c) throw();
 template <class M,class S> friend 	 M &_msdivassign(M &m,const S &c) throw();
 template <class MS,class S,class E> friend 	 E _mssdiv(const MS &ms, const S &c) throw();
	//-------- matrix-vector ---------------------
 template <class M,class V,class E> friend 	 E _mvlmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class E> friend 	 E _vmlmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class V,class M,class S> friend 	 V &_vmlmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class VS,class M,class S> friend 	 VS &_vsmlmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
	
 template <class M> friend 	 void *_mvoid(const M &m) throw();
 template <class M> friend 	 bool _mnot(const M &m) throw();
 template <class MS> friend 	 void *_msvoid(const MS &ms) throw();
 template <class MS> friend 	 bool _msnot(const MS &ms) throw();
 template <class M1,class M2> friend 	 bool _mmeq(const M1 &m1,const M2 &m2) throw();
 template <class M1,class M2> friend 	 bool _mmneq(const M1 &m1,const M2 &m2) throw();
 template <class M1,class M2> friend 	 bool _mmless(const M1 &m1,const M2 &m2) throw();
 template <class M1,class M2> friend 	 bool _mmleq(const M1 &m1,const M2 &m2) throw();
 template <class M,class MS> friend 	 bool _mmseq(const M &m1,const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmsneq(const M &m1,const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmsless(const M &m1,const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmsleq(const M &m1,const MS &ms) throw();
 template <class MS,class M> friend 	 bool _msmless(const MS &ms,const M &m1) throw();
 template <class MS,class M> friend 	 bool _msmleq(const MS &ms,const M &m1) throw();
 template <class M> friend 	std::ostream &_mout(std::ostream &s,const M &r) throw();
 template <class M> friend 	std::istream &_min(std::istream &s,M &r) throw();

	//--- Real --------- matrix-vector ----------------------
 template <class MS,class V,class E> friend 	 E _msvlmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmslmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

 template <class M1,class M2,class E> friend 	 E _mmlimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

 template <class M,class MS,class E> friend 	 E _mmslimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

 template <class MS,class M,class E> friend 	 E _msmlimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

 template <class MS1,class MS2,class E> friend 	 E _msmslimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

 template <class MS1,class MS2,class E> friend 	 E _msmsconv(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M1,class M2,class E> friend 	 E _mmsect(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmssect(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

 template <class MS1,class MS2,class E> friend 	 E _msmssect(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	// matrix-scalar

	// matrix-vector

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
 template <class MS,class V,class E> friend 	 E _msvlimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	// l_interval ----------

 template <class M,class E> friend 	 E _mdiam(const M &m) throw();
 template <class M,class E> friend 	 E _mmid(const M &m) throw();
 template <class MS,class E> friend 	 E _msdiam(const MS &ms) throw();
 template <class MS,class E> friend 	 E _msmid(const MS &ms) throw();
 template <class M,class E> friend 	 E _minf(const M &m) throw();
 template <class MS,class E> friend 	 E _msinf(const MS &ms) throw();
 template <class M,class E> friend 	 E _msup(const M &m) throw();
 template <class MS,class E> friend 	 E _mssup(const MS &ms) throw();
 template <class M1,class M2> friend 	 M1 &_mmsetinf(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class M2> friend 	 M1 &_mmsetsup(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class M2> friend 	 M1 &_mmusetinf(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class M2> friend 	 M1 &_mmusetsup(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class MS1,class M2> friend 	 MS1 &_msmsetinf(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class M2> friend 	 MS1 &_msmsetsup(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class M2> friend 	 MS1 &_msmusetinf(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class M2> friend 	 MS1 &_msmusetsup(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif

 template <class M1,class M2,class S> friend 	 M1 &_mmlimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

 template <class M,class MS,class S> friend 	 M &_mmslimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

 template <class M1,class M2> friend 	 M1 &_mmconvassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS> friend 	 M &_mmsconvassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmconvassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

 template <class M1,class M2> friend 	 M1 &_mmsectassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS> friend 	 M &_mmssectassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmsectassign(MS &ms,const M &m1)
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
 template <class VS,class M,class S> friend 	 VS &_vsmlimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
	
#endif
	
	//--------------------------  Konstruktoren ----------------------------

// l_real
	//! Constructor of class l_rmatrix
	INLINE l_rmatrix(const l_rmatrix &rm) throw();
	//! Constructor of class l_rmatrix
	INLINE l_rmatrix(const l_rmatrix_slice &rm) throw();
	//! Constructor of class l_rmatrix
	INLINE l_rmatrix() throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const int &m1, const int &n1, const int &m2, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const l_rvector &v) throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const l_rvector_slice &v) throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const l_rmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const l_rmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const l_rvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const l_rvector_slice &v) throw();
//  real
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const real &) throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const rmatrix &rm) throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const rmatrix_slice &rm) throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const rvector &v) throw();
	//! Constructor of class l_rmatrix
	explicit INLINE l_rmatrix(const rvector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const rmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const rmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const rvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix &operator =(const rvector_slice &v) throw();

	//--------------------------- Destruktoren -----------------------------

	INLINE ~l_rmatrix() throw() { delete [] dat; }

	//------------------------- Standardfunktionen -------------------------

	//! Operator for accessing a single row of the matrix
	INLINE l_rmatrix_subv operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE l_rmatrix_subv operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing the whole matrix
	INLINE l_rmatrix &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE l_rmatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE l_rmatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();
//#else
//#endif
};

		
//! The Multiple-Precision Data Type l_rmatrix_slice
/*!
This data type represents a partial l_rmatrix.

\sa l_rmatrix
*/
class l_rmatrix_slice
{
	friend class l_rmatrix;
	friend class l_imatrix;
	private:
	l_real *dat;
	int offset1,offset2,mxsize,mysize;
	int start1,end1,start2,end2,sxsize,sysize;     // slice size

	public:
//#if(CXSC_INDEX_CHECK)
#ifdef _CXSC_FRIEND_TPL
	//----------------- Templates ---------------------------------------
template <class V,class MS,class S> friend  void _vmsconstr(V &v,const MS &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmassign(MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmsassign(MS1 &ms1,const MS2 &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class M,class MS2,class S> friend 	 M &_mmsassign(M &m,const MS2 &ms) throw();
 template <class MS,class S> friend 	 MS &_mssassign(MS &ms,const S &r) throw();

 template <class MS> friend 	 int _mslb(const MS &ms, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<MS>);
#else
	throw();
#endif
 template <class MS> friend 	 int _msub(const MS &ms, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<MS>);
#else
	throw();
#endif
 template <class MS,class E> friend 	 E _msabs(const MS &ms) throw();
	//-------- matrix-matrix --------------
 template <class MS,class E> friend 	 E _msminus(const MS &ms) throw();
 template <class M,class MS,class E> friend 	 E _mmsplus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmsplus(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS> friend 	 M &_mmsplusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmplusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmsplusassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmsminus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmminus(const MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmsminus(const MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS> friend 	 M &_mmsminusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS,class M> friend 	 MS &_msmminusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmsminusassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmslmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmlmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS,class S> friend 	 M &_mmslmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmslmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmsconv(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmsconv(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	//--------- matrix-vector --------------
 template <class MS,class V,class E> friend 	 E _msvlmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class E> friend 	 E _vmslmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
 template <class V,class MS,class S> friend 	 V &_vmslmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
	//--------- matrix-scalar --------------
 template <class S,class MS,class E> friend 	 E _smsmult(const S &c, const MS &ms) throw();
 template <class MS,class S> friend 	 MS &_mssmultassign(MS &ms,const S &c) throw();
 template <class MS,class S,class E> friend 	 E _mssdiv(const MS &ms, const S &c) throw();
 template <class MS,class S> friend 	 MS &_mssdivassign(MS &ms,const S &c) throw();
	
 template <class MS> friend 	 void *_msvoid(const MS &ms) throw();
 template <class MS> friend 	 bool _msnot(const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmseq(const M &m1,const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmsneq(const M &m1,const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmsless(const M &m1,const MS &ms) throw();
 template <class M,class MS> friend 	 bool _mmsleq(const M &m1,const MS &ms) throw();
 template <class MS,class M> friend 	 bool _msmless(const MS &ms,const M &m1) throw();
 template <class MS,class M> friend 	 bool _msmleq(const MS &ms,const M &m1) throw();
 template <class MS1,class MS2> friend 	 bool _msmseq(const MS1 &ms1,const MS2 &ms2) throw();
 template <class MS1,class MS2> friend 	 bool _msmsneq(const MS1 &ms1,const MS2 &ms2) throw();
 template <class MS1,class MS2> friend 	 bool _msmsless(const MS1 &ms1,const MS2 &ms2) throw();
 template <class MS1,class MS2> friend 	 bool _msmsleq(const MS1 &ms1,const MS2 &ms2) throw();
 template <class MS> friend 	std::ostream &_msout(std::ostream &s,const MS &r) throw();
 template <class MS> friend 	std::istream &_msin(std::istream &s,MS &r) throw();

	//--- Real ------------ matrix-matrix -----------

 template <class M,class MS,class E> friend 	 E _mmslimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmlimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

 template <class MS1,class MS2,class E> friend 	 E _msmslimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

 template <class M,class MS,class E> friend 	 E _mmssect(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmssect(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	// matrix-scalar
	// matrix-vector

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

	// l_interval --------------
 
 template <class M1,class MS2> friend 	 M1 &_mmssetinf(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class MS2> friend 	 M1 &_mmssetsup(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class MS2> friend 	 M1 &_mmsusetinf(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class MS2> friend 	 M1 &_mmsusetsup(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmssetinf(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmssetsup(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmsusetinf(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmsusetsup(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif


 template <class M,class MS,class S> friend 	 M &_mmslimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif


 template <class M,class MS> friend 	 M &_mmsconvassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmsconvassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif

 template <class M,class MS> friend 	 M &_mmssectassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmssectassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif

  /*   friend TINLINE l_ivector &_vsmslimultassign<l_ivector_slice,l_rmatrix_slice,l_interval>(l_ivector_slice &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
	#endif */
 template <class V,class MS,class S> friend 	 V &_vmslimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

#endif

	//--------------- Konstruktoren ----------------------------------------

	//! Constructor of class l_rmatrix_slice
	explicit INLINE l_rmatrix_slice(l_rmatrix &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(l1-a.lb1),offset2(l2-a.lb2),mxsize(a.xsize),mysize(a.ysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	//! Constructor of class l_rmatrix_slice
	explicit INLINE l_rmatrix_slice(l_rmatrix_slice &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(a.offset1+l1-a.start1),offset2(a.offset2+l2-a.start2),mxsize(a.mxsize),mysize(a.mysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	public: 
	//! Constructor of class l_rmatrix_slice
	INLINE l_rmatrix_slice(const l_rmatrix_slice &ms) throw():dat(ms.dat),offset1(ms.offset1),offset2(ms.offset2),mxsize(ms.mxsize),mysize(ms.mysize),start1(ms.start1),end1(ms.end1),start2(ms.start2),end2(ms.end2),sxsize(ms.sxsize),sysize(ms.sysize) { }
	public:

	//---------------- Standardfunktionen -----------------------------------

	friend INLINE l_rvector::l_rvector(const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	friend INLINE l_rmatrix::l_rmatrix(const l_rmatrix_slice &) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_rmatrix_slice &operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Operator for accessing a single row of the matrix
	INLINE l_rmatrix_subv operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE l_rmatrix_subv operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing the whole matrix
	INLINE l_rmatrix_slice &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE l_rmatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE l_rmatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_slice &operator *=(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_slice &operator *=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_slice &operator *=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_slice &operator +=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_slice &operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_slice &operator +=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix_slice &operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_slice &operator -=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_slice &operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_slice &operator -=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix_slice &operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_rmatrix_slice &operator |=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_rmatrix_slice &operator |=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_rmatrix_slice &operator |=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_rmatrix_slice &operator |=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_slice &operator *=(const l_real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix_slice &operator *=(const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_rmatrix_slice &operator /=(const l_real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_rmatrix_slice &operator /=(const real &c) throw();
	INLINE operator void*() throw();
//#else
//#endif
};

//================================================================
//====================== Subvector Functions =====================

//=======================Vector / Scalar =========================

	//! Implementation of division operation
	INLINE l_rvector operator /(const l_rmatrix_subv &rv, const l_real &s) throw();
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rmatrix_subv &rv, const l_real &s) throw();
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_real &s, const l_rmatrix_subv &rv) throw();

//======================== Vector / Vector ========================
	
	//! Returns the absolute value of the vector
	INLINE l_rvector abs(const l_rmatrix_subv &mv) throw();
	
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const l_rvector_slice & sl1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_real operator *(const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_real operator *(const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_real operator *(const l_rmatrix_subv &rv1,const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_real operator *(const l_rvector_slice &sl,const l_rmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_real operator *(const l_rmatrix_subv &mv,const l_rvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE l_rvector operator +(const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rvector operator +(const l_rmatrix_subv &rv1,const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rvector operator +(const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rvector operator +(const l_rvector_slice &sl,const l_rmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rvector operator +(const l_rmatrix_subv &mv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_rvector operator -(const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rvector operator -(const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rvector operator -(const l_rmatrix_subv &rv1,const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rvector operator -(const l_rvector_slice &sl,const l_rmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rvector operator -(const l_rmatrix_subv &mv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//====================================================================
//===================== Matrix Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_rmatrix _l_rmatrix(const l_rmatrix &rm) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_rmatrix _l_rmatrix(const l_rvector &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_rmatrix _l_rmatrix(const l_rvector_slice &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_rmatrix _l_rmatrix(const l_real &r) throw();

	//! Returns the lower bound index
	INLINE int Lb(const l_rmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const l_rmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the lower bound index
	INLINE int Lb(const l_rmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const l_rmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the lower bound index
	INLINE l_rmatrix &SetLb(l_rmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the upper bound index
	INLINE l_rmatrix &SetUb(l_rmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(l_rmatrix &A) throw();
	//! Resizes the matrix
	INLINE void Resize(l_rmatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(l_rmatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif

	//! Returns the absolute value of the matrix
	INLINE l_rmatrix abs(const l_rmatrix &m) throw();
	//! Returns the absolute value of the matrix
	INLINE l_rmatrix abs(const l_rmatrix_slice &ms) throw();

//===================== Matrix / Scalar ===============================

	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_real &c, const l_rmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_real &c, const l_rmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix &m,const l_real &c) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms,const l_real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix &operator *=(l_rmatrix &m,const l_real &c) throw();
	//! Implementation of division operation
	INLINE l_rmatrix operator /(const l_rmatrix &m,const l_real &c) throw();
	//! Implementation of division operation
	INLINE l_rmatrix operator /(const l_rmatrix_slice &ms, const l_real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_rmatrix &operator /=(l_rmatrix &m,const l_real &c) throw();
	
//------------ real - l_rmatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const real &c, const l_rmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const real &c, const l_rmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix &m,const real &c) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms,const real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix &operator *=(l_rmatrix &m,const real &c) throw();
	//! Implementation of division operation
	INLINE l_rmatrix operator /(const l_rmatrix &m,const real &c) throw();
	//! Implementation of division operation
	INLINE l_rmatrix operator /(const l_rmatrix_slice &ms, const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_rmatrix &operator /=(l_rmatrix &m,const real &c) throw();
//----------------- rmatrix - l_real ----------------

	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_real &c, const rmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_real &c, const rmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const rmatrix &m,const l_real &c) throw();
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const rmatrix_slice &ms,const l_real &c) throw();
	//! Implementation of division operation
	INLINE l_rmatrix operator /(const rmatrix &m,const l_real &c) throw();
	//! Implementation of division operation
	INLINE l_rmatrix operator /(const rmatrix_slice &ms, const l_real &c) throw();
	

//============================ Matrix / Vector ===================================

	
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rmatrix &m,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rmatrix_slice &ms,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rvector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rvector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rvector &operator *=(l_rvector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rvector &operator *=(l_rvector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rvector_slice &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rvector_slice &v,const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- real -------------------------------------

	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const rvector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const rvector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const rvector_slice &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rmatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rvector operator *(const l_rmatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//================ Matrix / Matrix ============================

	//! Implementation of positive sign operation
	INLINE const l_rmatrix &operator +(const l_rmatrix &m1) throw();
	//! Implementation of positive sign operation
	INLINE l_rmatrix operator +(const l_rmatrix_slice &ms) throw();
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of negative sign operation
	INLINE l_rmatrix operator -(const l_rmatrix &m) throw();
	//! Implementation of negative sign operation
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms) throw();
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	
	//---------- rmatrix-l_rmatrix ------------------
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const rmatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_rmatrix operator +(const l_rmatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const rmatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const rmatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const rmatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const rmatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const rmatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//============== Compare Operator ==========================

//-------------- Matrix - Matrix   -------------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_rmatrix &m1,const l_rmatrix &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_rmatrix &m1,const l_rmatrix &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_rmatrix &m1,const l_rmatrix &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_rmatrix &m1,const l_rmatrix &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_rmatrix &m1,const l_rmatrix &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_rmatrix &m1,const l_rmatrix &m2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw();

//---------------- Matrix - Matrix_slice ----------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw();

//=================== Not Operator =============================

	//! Implementation of standard negation operation
	INLINE bool operator !(const l_rmatrix &ms) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const l_rmatrix_slice &ms) throw();

//======================== Input / Output ========================

	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const l_rmatrix &r) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const l_rmatrix_slice &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,l_rmatrix &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,l_rmatrix_slice &r) throw();

        //! Returns the row dimension
        INLINE int RowLen ( const l_rmatrix& );
        //! Returns the column dimension
        INLINE int ColLen ( const l_rmatrix& );
         //! Returns the row dimension
        INLINE int RowLen ( const l_rmatrix_slice& );
        //! Returns the column dimension
        INLINE int ColLen ( const l_rmatrix_slice& );
       //! Returns the Identity matrix
        l_rmatrix  Id     ( const l_rmatrix& );
        //! Returns the transposed matrix
        l_rmatrix  transp ( const l_rmatrix& );
        //! Doubles the size of the matrix
        void       DoubleSize ( l_rmatrix& );
	
} // namespace cxsc 

#ifdef _CXSC_INCL_INL
#include "matrix.inl"
#include "l_rmatrix.inl"
#endif

#ifdef _CXSC_LIVECTOR_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "liveclrmat.inl"
# else
#  include "liveclrmat.hpp"
# endif
#endif

#ifdef _CXSC_IVECTOR_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "iveclrmat.inl"
# else
#  include "iveclrmat.hpp"
# endif
#endif

#ifdef _CXSC_IMATRIX_HPP_INCLUDED
# ifdef _CXSC_INCL_INL
#  include "lrmatimat.inl"
# else
#  include "lrmatimat.hpp"
# endif
#endif


#endif

