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

/* CVS $Id: intmatrix.hpp,v 1.20 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_INTMATRIX_HPP_INCLUDED
#define _CXSC_INTMATRIX_HPP_INCLUDED

#include "xscclass.hpp"

#include "dot.hpp"
#include "intvector.hpp"
#include "except.hpp"
#include "matrix.hpp"

namespace cxsc {

class intmatrix;
class intmatrix_slice;

//! The Data Type intmatrix_subv
/*!
This Data Type provides one column or row of a matrix as a vector.
*/
class intmatrix_subv
{
	friend class intvector;
	friend class intmatrix;
	friend class intmatrix_slice;
	private:
	int *dat;
	int lb,ub;
	int size,start,offset; // start=first element index 0..n-1
	
	public:
	//! Returns one row of the matrix as a vector
	friend INLINE intmatrix_subv Row(intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one column of the matrix as a vector
	friend INLINE intmatrix_subv Col(intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	friend INLINE intmatrix_subv Row(const intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one column of the matrix as a vector
	friend INLINE intmatrix_subv Col(const intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
//#if(CXSC_INDEX_CHECK)
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

#endif
	
	//----------------- Konstruktoren ----------------------------------

	//! Constructor of class intmatrix_subv
	explicit INLINE intmatrix_subv (int *d, const int &l, const int &u, const int &s, const int &st, const int &o) throw():dat(d),lb(l),ub(u),size(s),start(st),offset(o) { }
        public:
	//! Constructor of class intmatrix_subv
	INLINE intmatrix_subv(const intmatrix_subv &v) throw():dat(v.dat),lb(v.lb),ub(v.ub),size(v.size),start(v.start),offset(v.offset) { }
	public:

	//---------------------- Standardfunktionen ------------------------

	friend INLINE intvector::intvector(const intmatrix_subv &) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix_subv &operator =(const intmatrix_subv &rv) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix_subv &operator =(const int &r) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix_subv &operator =(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_subv &operator =(const intmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_subv &operator =(const intvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_subv &operator =(const intvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the lower bound of the vector
	friend INLINE int Lb(const intmatrix_subv &rv) throw() { return rv.lb; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const intmatrix_subv &rv) throw() { return rv.ub; }
	//! Operator for accessing the single elements of the vector
	INLINE int &operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	//! Operator for accessing the whole vector
	INLINE intmatrix_subv &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE intmatrix_subv operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE intmatrix_subv operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	//! Implementation of multiplication and allocation operation
	INLINE intmatrix_subv &operator *=(const int &c) throw();
	//! Implementation of addition and allocation operation
	INLINE intmatrix_subv &operator +=(const int &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix_subv &operator -=(const int &c) throw();
	//! Implementation of division and allocation operation
	INLINE intmatrix_subv &operator /=(const int &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix_subv &operator -=(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix_subv &operator -=(const intvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE intmatrix_subv &operator +=(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE intmatrix_subv &operator +=(const intvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
//#else
//#endif	

};


//! Returns one row of the matrix as a vector
INLINE intmatrix_subv Row(intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
//! Returns one column of the matrix as a vector
INLINE intmatrix_subv Col(intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif


//----------------------- Matrix -----------------------------------------------

class intmatrix_slice;
//! The Data Type intmatrix
/*!
\sa rmatrix
*/
class intmatrix
{
	friend class intmatrix_slice;
	friend class intmatrix_subv;
	private:
	int *dat;
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
	//------- matrix-matrix --------------
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
 template <class MS1,class MS2,class E> friend 	 E _msmsplus(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
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
	//-------- matrix-scalar ---------------------
 template <class S,class M,class E> friend 	 E _smmult(const S &c, const M &m) throw();
 template <class M,class S> friend 	 M &_msmultassign(M &m,const S &c) throw();
 template <class S,class MS,class E> friend 	 E _smsmult(const S &c, const MS &ms) throw();
 template <class M,class S,class E> friend 	 E _msdiv(const M &m,const S &c) throw();
 template <class M,class S> friend 	 M &_msdivassign(M &m,const S &c) throw();
 template <class MS,class S,class E> friend 	 E _mssdiv(const MS &ms, const S &c) throw();
	//--------- matrix-vector --------------------

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

#endif

	//--------------------------  Konstruktoren ----------------------------

	//! Constructor of class intmatrix
	INLINE intmatrix(const intmatrix &rm) throw();
	//! Constructor of class intmatrix
	INLINE intmatrix(const intmatrix_slice &rm) throw();
	//! Constructor of class intmatrix
	INLINE intmatrix() throw();
	//! Constructor of class intmatrix
	explicit INLINE intmatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class intmatrix
	explicit INLINE intmatrix(const int &m1, const int &n1, const int &m2, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class intmatrix
	explicit INLINE intmatrix(const intvector &v) throw();
	//! Constructor of class intmatrix
	explicit INLINE intmatrix(const intvector_slice &v) throw();
	//! Constructor of class intmatrix
	explicit INLINE intmatrix(const int &r) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix &operator =(const int &r) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix &operator =(const intmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix &operator =(const intmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix &operator =(const intvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix &operator =(const intvector_slice &v) throw();

	//--------------------------- Destruktoren -----------------------------

	INLINE ~intmatrix() throw() { delete [] dat; }

	//------------------------- Standardfunktionen -------------------------

	friend INLINE intvector::intvector(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Operator for accessing a single row of the matrix
	INLINE intmatrix_subv operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE intmatrix_subv operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing the whole matrix
	INLINE intmatrix &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE intmatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE intmatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();
//#else
//#endif
};

//! The Data Type intmatrix_slice
/*!
This data type represents a partial intmatrix.

\sa intmatrix
*/
class intmatrix_slice
{
	friend class intmatrix;
	private:
	int *dat;
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
	//------------ matrix-matrix --------------------
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
	//--------- matrix-scalar -----------------
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
	
#endif


	//--------------- Konstruktoren ----------------------------------------

	//! Constructor of class intmatrix_slice
	explicit INLINE intmatrix_slice(intmatrix &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(l1-a.lb1),offset2(l2-a.lb2),mxsize(a.xsize),mysize(a.ysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	//! Constructor of class intmatrix_slice
	explicit INLINE intmatrix_slice(intmatrix_slice &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(a.offset1+l1-a.start1),offset2(a.offset2+l2-a.start2),mxsize(a.mxsize),mysize(a.mysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	public: 
	//! Constructor of class intmatrix_slice
	INLINE intmatrix_slice(const intmatrix_slice &ms) throw():dat(ms.dat),offset1(ms.offset1),offset2(ms.offset2),mxsize(ms.mxsize),mysize(ms.mysize),start1(ms.start1),end1(ms.end1),start2(ms.start2),end2(ms.end2),sxsize(ms.sxsize),sysize(ms.sysize) { }
	public:

	//---------------- Standardfunktionen -----------------------------------

	friend  intvector::intvector(const intmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	friend INLINE intmatrix::intmatrix(const intmatrix_slice &) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix_slice &operator =(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_slice &operator =(const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_slice &operator =(const int &r) throw();
	//! Implementation of standard assigning operator
	INLINE intmatrix_slice &operator =(const intvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_slice &operator =(const intvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE intmatrix_slice &operator =(const intmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Operator for accessing a single row of the matrix
	INLINE intmatrix_subv operator [](const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE intmatrix_subv operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing the whole matrix
	INLINE intmatrix_slice &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE intmatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE intmatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE intmatrix_slice& operator +=(const intmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE intmatrix_slice& operator +=(const intmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix_slice& operator -=(const intmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix_slice& operator -=(const intmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE intmatrix_slice& operator *=(const int &c) throw();
	//! Implementation of division and allocation operation
	INLINE intmatrix_slice& operator /=(const int &c) throw();
	INLINE operator void*() throw();
//#else
//#endif
};

//================================================================
//====================== Subvector Functions =====================

//=======================Vector / Scalar =========================

	//! Implementation of division operation
	INLINE intvector operator /(const intmatrix_subv &rv, const int &s) throw();
	//! Implementation of multiplication operation
	INLINE intvector operator *(const intmatrix_subv &rv, const int &s) throw();
	//! Implementation of multiplication operation
	INLINE intvector operator *(const int &s, const intmatrix_subv &rv) throw();
	//! Returns the absolute value of the matrix
	INLINE intvector abs(const intmatrix_subv &mv) throw();

//======================== Vector / Vector ========================

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const intmatrix_subv & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const intvector & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp, const intmatrix_subv & rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp,const intvector_slice &sl,const intmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(dotprecision &dp,const intmatrix_subv &mv,const intvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE intvector operator +(const intmatrix_subv & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intvector operator +(const intmatrix_subv &rv1,const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intvector operator +(const intvector & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intvector operator +(const intvector_slice &sl,const intmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intvector operator +(const intmatrix_subv &mv,const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE intvector operator -(const intmatrix_subv & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intvector operator -(const intvector & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intvector operator -(const intmatrix_subv &rv1,const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intvector operator -(const intvector_slice &sl,const intmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intvector operator -(const intmatrix_subv &mv,const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//====================================================================
//===================== Matrix Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE intmatrix _intmatrix(const intmatrix &rm) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE intmatrix _intmatrix(const intvector &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE intmatrix _intmatrix(const intvector_slice &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE intmatrix _intmatrix(const int &r) throw();

	//! Returns the lower bound index
	INLINE int Lb(const intmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const intmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the lower bound index
	INLINE int Lb(const intmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const intmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the lower bound index
	INLINE intmatrix &SetLb(intmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the upper bound index
	INLINE intmatrix &SetUb(intmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(intmatrix &A) throw();
	//! Resizes the matrix
	INLINE void Resize(intmatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(intmatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif

	//! Returns the absolute value of the matrix
	INLINE intmatrix abs(const intmatrix &m) throw();
	//! Returns the absolute value of the matrix
	INLINE intmatrix abs(const intmatrix_slice &ms) throw();

//===================== Matrix / Scalar ===============================

	//! Implementation of multiplication operation
	INLINE intmatrix operator *(const int &c, const intmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE intmatrix operator *(const int &c, const intmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE intmatrix operator *(const intmatrix &m,const int &c) throw();
	//! Implementation of multiplication operation
	INLINE intmatrix operator *(const intmatrix_slice &ms,const int &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE intmatrix &operator *=(intmatrix &m,const int &c) throw();
	//! Implementation of division operation
	INLINE intmatrix operator /(const intmatrix &m,const int &c) throw();
	//! Implementation of division operation
	INLINE intmatrix operator /(const intmatrix_slice &ms, const int &c) throw();
	//! Implementation of division and allocation operation
	INLINE intmatrix &operator /=(intmatrix &m,const int &c) throw();


//================ Matrix / Matrix ============================

	//! Implementation of positive sign operation
	INLINE const intmatrix &operator +(const intmatrix &m1) throw();
	//! Implementation of positive sign operation
	INLINE intmatrix operator +(const intmatrix_slice &ms) throw();
	//! Implementation of addition operation
	INLINE intmatrix operator +(const intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intmatrix operator +(const intmatrix &m,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intmatrix operator +(const intmatrix_slice &ms,const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE intmatrix operator +(const intmatrix_slice &m1,const intmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE intmatrix &operator +=(intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE intmatrix &operator +=(intmatrix &m1,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of negative sign operation
	INLINE intmatrix operator -(const intmatrix &m) throw();
	//! Implementation of negative sign operation
	INLINE intmatrix operator -(const intmatrix_slice &ms) throw();
	//! Implementation of subtraction operation
	INLINE intmatrix operator -(const intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intmatrix operator -(const intmatrix &m,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intmatrix operator -(const intmatrix_slice &ms,const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE intmatrix operator -(const intmatrix_slice &ms1,const intmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix &operator -=(intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE intmatrix &operator -=(intmatrix &m1,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//============== Compare Operator ==========================

//-------------- Matrix - Matrix   -------------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const intmatrix &m1,const intmatrix &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const intmatrix &m1,const intmatrix &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const intmatrix &m1,const intmatrix &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const intmatrix &m1,const intmatrix &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const intmatrix &m1,const intmatrix &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const intmatrix &m1,const intmatrix &m2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const intmatrix &m1,const intmatrix_slice &ms) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const intmatrix &m1,const intmatrix_slice &ms) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const intmatrix &m1,const intmatrix_slice &ms) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const intmatrix &m1,const intmatrix_slice &ms) throw();
	//! Implementation of standard more-than operation
	INLINE bool operator >(const intmatrix &m1,const intmatrix_slice &ms) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const intmatrix &m1,const intmatrix_slice &ms) throw();

//---------------- Matrix - Matrix_slice ----------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const intmatrix_slice &m1,const intmatrix_slice &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const intmatrix_slice &m1,const intmatrix_slice &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const intmatrix_slice &m1,const intmatrix_slice &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const intmatrix_slice &m1,const intmatrix_slice &m2) throw();
	//! Implementation of standard more-than operation
	INLINE bool operator >(const intmatrix_slice &m1,const intmatrix_slice &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const intmatrix_slice &m1,const intmatrix_slice &m2) throw();

//=================== Not Operator =============================

	//! Implementation of standard negation operation
	INLINE bool operator !(const intmatrix &ms) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const intmatrix_slice &ms) throw();

//======================== Input / Output ========================

	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const intmatrix &r) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const intmatrix_slice &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,intmatrix &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,intmatrix_slice &r) throw();

        //! Returns the row dimension
        INLINE int      RowLen     ( const intmatrix& );
        //! Returns the column dimension
        INLINE int      ColLen     ( const intmatrix& );
        //! Doubles the size of the matrix
        //! Returns the row dimension
        INLINE int      RowLen     ( const intmatrix_slice& );
        //! Returns the column dimension
        INLINE int      ColLen     ( const intmatrix_slice& );
        //! Doubles the size of the matrix
        intmatrix    Id     ( const intmatrix& );
        //! Returns the transposed matrix
        intmatrix    transp ( const intmatrix& );
        //! Doubles the size of the matrix
        void     DoubleSize ( intmatrix& );

//=================== Permutation matrix/vector functions =======================

        INLINE intvector permvec(const intmatrix&);
        INLINE intmatrix permmat(const intvector&);
        INLINE intmatrix perminv(const intmatrix&);

} // namespace cxsc 

#ifdef _CXSC_INCL_INL
# include "matrix.inl"
# include "intmatrix.inl"
#endif

#endif

