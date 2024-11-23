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

/* CVS $Id: l_imatrix.hpp,v 1.22 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_LIMATRIX_HPP_INCLUDED
#define _CXSC_LIMATRIX_HPP_INCLUDED

#include "xscclass.hpp"
#include "idot.hpp"
#include "l_ivector.hpp"
#include "except.hpp"
#include "matrix.hpp"
#include "imatrix.hpp"
#include "l_rmatrix.hpp"

namespace cxsc {

class l_imatrix;
class l_imatrix_slice;

//! The Multiple-Precision Data Type l_imatrix_subv
/*!
This Data Type provides one column or row of a matrix as a vector.
*/
class l_imatrix_subv
{
	friend class l_ivector;
	friend class l_imatrix;
	friend class l_imatrix_slice;
	private:
	l_interval *dat;
	int lb,ub;
	int size,start,offset; // start=first element index 0..n-1
	
	public:
	//! Returns one row of the matrix as a vector
	friend INLINE l_imatrix_subv Row(l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one column of the matrix as a vector
	friend INLINE l_imatrix_subv Col(l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one row of the matrix as a vector
	friend INLINE l_imatrix_subv Row(const l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Returns one column of the matrix as a vector
	friend INLINE l_imatrix_subv Col(const l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
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
template <class MV,class V> friend  V _mvabs(const MV &mv) throw();
template <class MV,class V> friend  V _mvdiam(const MV &mv) throw();
template <class MV,class V> friend  V _mvmid(const MV &mv) throw();
template <class MV,class V> friend  V _mvinf(const MV &mv) throw();
template <class MV,class V> friend  V _mvsup(const MV &mv) throw();

 template <class MV,class S> friend 	 MV &_mvssetinf(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvssetsup(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvsusetinf(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvsusetsup(MV &mv, const S &s) throw();

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
 template <class MV1,class MV2,class S> friend 	 S _mvmvlimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif
 template <class V,class MV,class S> friend 	 S _vmvlimult(const V &rv1, const MV &rv2)
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


#endif

	//----------------- Konstruktoren ----------------------------------

	//! Constructor of class l_imatrix_subv
	explicit INLINE l_imatrix_subv (l_interval *d, const int &l, const int &u, const int &s, const int &st, const int &o) throw():dat(d),lb(l),ub(u),size(s),start(st),offset(o) { }
	public:
	//! Constructor of class l_imatrix_subv
	INLINE l_imatrix_subv(const l_imatrix_subv &v) throw():dat(v.dat),lb(v.lb),ub(v.ub),size(v.size),start(v.start),offset(v.offset) { }
	public:

	//---------------------- Standardfunktionen ------------------------

	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_imatrix_subv &rv) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// Real
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const rmatrix_subv &rv) throw();

	// l_real
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const l_rmatrix_subv &rv) throw();

	// interval
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_subv &operator =(const imatrix_subv &rv) throw();

	//! Returns the lower bound of the vector
	friend INLINE int Lb(const l_imatrix_subv &rv) throw() { return rv.lb; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const l_imatrix_subv &rv) throw() { return rv.ub; }
	//! Operator for accessing the single elements of the vector
	INLINE l_interval &operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif
	//! Operator for accessing the whole vector
	INLINE l_imatrix_subv &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE l_imatrix_subv operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE l_imatrix_subv operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_subv &operator *=(const l_interval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const l_interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const l_interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_subv &operator /=(const l_interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// real
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_subv &operator *=(const real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_subv &operator /=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// l_real
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_subv &operator *=(const l_real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const l_real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const l_real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_subv &operator /=(const l_real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// interval
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_subv &operator *=(const interval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_subv &operator /=(const interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_subv &operator -=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_subv &operator +=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_subv &operator |=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_subv &operator &=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
//#else
//#endif	

};


//----------------------- Matrix -----------------------------------------------

class l_imatrix_slice;

//! The Multiple-Precision Data Type l_imatrix
/*!
\sa l_rmatrix
*/
class l_imatrix
{
	friend class l_imatrix_slice;
	friend class l_imatrix_subv;
	private:
	l_interval *dat;
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
 template <class M,class E> friend 	 E _mdiam(const M &m) throw();
 template <class M,class E> friend 	 E _mmid(const M &m) throw();
 template <class M,class E> friend 	 E _minf(const M &m) throw();
 template <class M,class E> friend 	 E _msup(const M &m) throw();
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
 template <class M1,class M2,class E> friend 	 E _mmlimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M1,class M2,class S> friend 	 M1 &_mmlimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
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
 template <class M,class MS,class S> friend 	 M &_mmslimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmslimult(const MS1 &ms1, const MS2 &ms2)
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
 template <class MS1,class MS2,class E> friend 	 E _msmssect(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
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

	//--- Real --------- matrix-scalar ----------------------
 template <class M1,class MS2> friend 	 M1 &_mmssetre(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class MS2> friend 	 M1 &_mmssetim(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	//--- Real --------- matrix-vector ----------------------
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


#endif
	
	//--------------------------  Konstruktoren ----------------------------

// l_interval
	//! Constructor of class l_imatrix
	INLINE l_imatrix(const l_imatrix &rm) throw();
	//! Constructor of class l_imatrix
	INLINE l_imatrix(const l_imatrix_slice &rm) throw();
	//! Constructor of class l_imatrix
	INLINE l_imatrix() throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const int &m1, const int &n1, const int &m2, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_ivector &v) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_ivector_slice &v) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_imatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_imatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_ivector &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_ivector_slice &v) throw();
//  real
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const real &r) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const rmatrix &rm) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const rmatrix_slice &rm) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const rvector &v) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const rvector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const rmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const rmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const rvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const rvector_slice &v) throw();

//  l_real
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_real &r) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_rmatrix &rm) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_rmatrix_slice &rm) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_rvector &v) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const l_rvector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_rmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_rmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_rvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const l_rvector_slice &v) throw();

//  interval
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const interval &r) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const imatrix &rm) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const imatrix_slice &rm) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const ivector &v) throw();
	//! Constructor of class l_imatrix
	explicit INLINE l_imatrix(const ivector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const imatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const imatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const ivector &v) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix &operator =(const ivector_slice &v) throw();

	//--------------------------- Destruktoren -----------------------------

	INLINE ~l_imatrix() throw() { delete [] dat; }

	//------------------------- Standardfunktionen -------------------------

	//! Operator for accessing a single row of the matrix
	INLINE l_imatrix_subv operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE l_imatrix_subv operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing the whole matrix
	INLINE l_imatrix &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE l_imatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE l_imatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();
//#else
//#endif
};

	
//! The Multiple-Precision Data Type l_imatrix_slice
/*!
This data type represents a partial l_imatrix.

\sa l_imatrix
*/
class l_imatrix_slice
{
	friend class l_imatrix;
	private:
	l_interval *dat;
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
 template <class MS,class E> friend 	 E _msinf(const MS &ms) throw();
 template <class MS,class E> friend 	 E _mssup(const MS &ms) throw();
 template <class MS,class E> friend 	 E _msdiam(const MS &ms) throw();
 template <class MS,class E> friend 	 E _msmid(const MS &ms) throw();
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
 template <class M,class MS,class S> friend 	 M &_mmslimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmslimult(const MS1 &ms1, const MS2 &ms2)
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
 template <class MS1,class MS2> friend 	 MS1 &_msmsconvassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmsconv(const MS1 &m1,const MS2 &m2)
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
 template <class MS1,class MS2> friend 	 MS1 &_msmssectassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmssect(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	//--------- matrix-vector --------------
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


#endif

	//--------------- Konstruktoren ----------------------------------------

	//! Constructor of class l_imatrix_slice
	explicit INLINE l_imatrix_slice(l_imatrix &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(l1-a.lb1),offset2(l2-a.lb2),mxsize(a.xsize),mysize(a.ysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	//! Constructor of class l_imatrix_slice
	explicit INLINE l_imatrix_slice(l_imatrix_slice &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(a.offset1+l1-a.start1),offset2(a.offset2+l2-a.start2),mxsize(a.mxsize),mysize(a.mysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	public: 
	//! Constructor of class l_imatrix_slice
	INLINE l_imatrix_slice(const l_imatrix_slice &ms) throw():dat(ms.dat),offset1(ms.offset1),offset2(ms.offset2),mxsize(ms.mxsize),mysize(ms.mysize),start1(ms.start1),end1(ms.end1),start2(ms.start2),end2(ms.end2),sxsize(ms.sxsize),sysize(ms.sysize) { }
	public:

	//---------------- Standardfunktionen -----------------------------------

	friend INLINE l_ivector::l_ivector(const l_imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	friend INLINE l_imatrix::l_imatrix(const l_imatrix_slice &) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_imatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// real
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	// interval
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const imatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	// l_real
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_real &r) throw();
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE l_imatrix_slice &operator =(const l_rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Operator for accessing a single row of the matrix
	INLINE l_imatrix_subv operator [](const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE l_imatrix_subv operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing the whole matrix
	INLINE l_imatrix_slice &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE l_imatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE l_imatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();

	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const l_interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const l_interval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const l_interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_slice &operator /=(const l_interval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_slice &operator /=(const real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const l_real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const l_real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const l_real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_slice &operator /=(const l_real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const interval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix_slice &operator /=(const interval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix_slice &operator +=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix_slice &operator -=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix_slice &operator |=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix_slice &operator &=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix_slice &operator *=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//#else
//#endif
};

//================================================================
//====================== Subvector Functions =====================

//=======================Vector / Scalar =========================

	//! Implementation of division operation
	INLINE l_ivector operator /(const l_imatrix_subv &rv, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix_subv &rv, const l_interval &s) throw();
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_interval &s, const l_imatrix_subv &rv) throw();
	//! Returns the absolute value of the matrix
	INLINE l_ivector abs(const l_imatrix_subv &mv) throw();
	//! Returns the rounded diameter of the matrix
	INLINE l_rvector diam(const l_imatrix_subv &mv) throw();
	//! Returns the rounded middle of the matrix
	INLINE l_rvector mid(const l_imatrix_subv &mv) throw();
	//! Returns the infimum of the matrix
	INLINE l_rvector Inf(const l_imatrix_subv &mv) throw();
	//! Returns the supremum of the matrix
	INLINE l_rvector Sup(const l_imatrix_subv &mv) throw();
	//! Returns the matrix with the new given infimum value
	INLINE l_imatrix_subv &SetInf(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE l_imatrix_subv &SetSup(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE l_imatrix_subv &UncheckedSetInf(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE l_imatrix_subv &UncheckedSetSup(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the matrix with the new given supremum value
	INLINE l_imatrix_subv &SetSup(l_imatrix_subv &iv,const l_real &r) throw();
	//! Returns the matrix with the new given infimum value
	INLINE l_imatrix_subv &SetInf(l_imatrix_subv &iv,const l_real &r) throw();
	//! Returns the matrix with the new unchecked given supremum value
	INLINE l_imatrix_subv &UncheckedSetSup(l_imatrix_subv &iv,const l_real &r) throw();
	//! Returns the matrix with the new unchecked given infimum value
	INLINE l_imatrix_subv &SetUncheckedInf(l_imatrix_subv &iv,const l_real &r) throw();

//======================== Vector / Vector ========================

	
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
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
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_ivector_slice &sl,const l_imatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_interval operator *(const l_imatrix_subv &mv,const l_ivector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_ivector_slice &sl,const l_imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_ivector operator +(const l_imatrix_subv &mv,const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_ivector_slice &sl,const l_imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_ivector operator -(const l_imatrix_subv &mv,const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//  real

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rvector_slice &sl2)
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
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
// l_real

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rvector_slice &sl2)
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
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
// interval

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const ivector_slice &sl2)
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
	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//====================================================================
//===================== Matrix Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_imatrix _imatrix(const l_imatrix &rm) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_imatrix _imatrix(const l_ivector &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_imatrix _imatrix(const l_ivector_slice &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE l_imatrix _imatrix(const l_interval &r) throw();

	//! Returns the lower bound index
	INLINE int Lb(const l_imatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const l_imatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the lower bound index
	INLINE int Lb(const l_imatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const l_imatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the lower bound index
	INLINE l_imatrix &SetLb(l_imatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the upper bound index
	INLINE l_imatrix &SetUb(l_imatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(l_imatrix &A) throw();
	//! Resizes the matrix
	INLINE void Resize(l_imatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(l_imatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif

	//! Returns the absolute value of the matrix
	INLINE l_imatrix abs(const l_imatrix &m) throw();
	//! Returns the absolute value of the matrix
	INLINE l_imatrix abs(const l_imatrix_slice &ms) throw();
	//! Returns the rounded diameter of the matrix
	INLINE l_rmatrix diam(const l_imatrix &m) throw();
	//! Returns the rounded diameter of the matrix
	INLINE l_rmatrix diam(const l_imatrix_slice &m) throw();
	//! Returns the rounded middle of the matrix
	INLINE l_rmatrix mid(const l_imatrix &m) throw();
	//! Returns the rounded middle of the matrix
	INLINE l_rmatrix mid(const l_imatrix_slice &m) throw();
	//! Returns the infimum of the matrix
	INLINE l_rmatrix Inf(const l_imatrix &m) throw();
	//! Returns the supremum of the matrix
	INLINE l_rmatrix Sup(const l_imatrix &m) throw();
	//! Returns the infimum of the matrix
	INLINE l_rmatrix Inf(const l_imatrix_slice &m) throw();
	//! Returns the supremum of the matrix
	INLINE l_rmatrix Sup(const l_imatrix_slice &m) throw();
	//! Returns the matrix with the new given infimum value
	INLINE l_imatrix &SetInf(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given infimum value
	INLINE l_imatrix_slice &SetInf(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given infimum value
	INLINE l_imatrix &SetInf(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given infimum value
	INLINE l_imatrix_slice &SetInf(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE l_imatrix &SetSup(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE l_imatrix_slice &SetSup(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE l_imatrix &SetSup(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE l_imatrix_slice &SetSup(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE l_imatrix &UncheckedSetInf(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE l_imatrix_slice &UncheckedSetInf(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE l_imatrix &UncheckedSetInf(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE l_imatrix_slice &UncheckedSetInf(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE l_imatrix &UncheckedSetSup(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE l_imatrix_slice &UncheckedSetSup(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE l_imatrix &UncheckedSetSup(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE l_imatrix_slice &UncheckedSetSup(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//===================== Matrix / Scalar ===============================

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const l_imatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const l_imatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m,const l_interval &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const l_interval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix &m,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const l_interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix &operator /=(l_imatrix &m,const l_interval &c) throw();
	
//------------ real - l_imatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const real &c, const l_imatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const real &c, const l_imatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m,const real &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m,const real &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix &m,const real &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix &operator /=(l_imatrix &m,const real &c) throw();
//----------------- rmatrix - l_interval ----------------

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const rmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const rmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const rmatrix &m,const l_interval &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const rmatrix_slice &ms,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const rmatrix &m,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const rmatrix_slice &ms, const l_interval &c) throw();
	
//------------ l_real - l_imatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_real &c, const l_imatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_real &c, const l_imatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m,const l_real &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const l_real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m,const l_real &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix &m,const l_real &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const l_real &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix &operator /=(l_imatrix &m,const l_real &c) throw();
//----------------- l_rmatrix - l_interval ----------------

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const l_rmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const l_rmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix &m,const l_interval &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_rmatrix &m,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_rmatrix_slice &ms, const l_interval &c) throw();
	
//------------ interval - l_imatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const interval &c, const l_imatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const interval &c, const l_imatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m,const interval &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const interval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m,const interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix &m,const interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE l_imatrix &operator /=(l_imatrix &m,const interval &c) throw();
//----------------- imatrix - l_interval ----------------

	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const imatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_interval &c, const imatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix &m,const l_interval &c) throw();
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix_slice &ms,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const imatrix &m,const l_interval &c) throw();
	//! Implementation of division operation
	INLINE l_imatrix operator /(const imatrix_slice &ms, const l_interval &c) throw();
	

//============================ Matrix / Vector ===================================


	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix &m,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_ivector &operator *=(l_ivector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_ivector_slice &v,const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- real -------------------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const rvector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const rvector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const rvector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- l_real -------------------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rvector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rvector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_rvector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix &m,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- interval -------------------------------------

	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const ivector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//================ Matrix / Matrix ============================

	//! Implementation of positive sign operation
	INLINE const l_imatrix &operator +(const l_imatrix &m1) throw();
	//! Implementation of positive sign operation
	INLINE l_imatrix operator +(const l_imatrix_slice &ms) throw();
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of negative sign operation
	INLINE l_imatrix operator -(const l_imatrix &m) throw();
	//! Implementation of negative sign operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms) throw();
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- rmatrix-l_imatrix ------------------
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const rmatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const rmatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const rmatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const rmatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const rmatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- l_rmatrix-l_imatrix ------------------
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- imatrix-l_imatrix ------------------
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE l_imatrix &operator +=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE l_imatrix &operator -=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE l_imatrix &operator *=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE l_imatrix &operator |=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE l_imatrix &operator &=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- l_rmatrix-imatrix ------------------
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const l_rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE l_imatrix operator +(const imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE l_imatrix operator -(const imatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE l_imatrix operator *(const imatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const l_rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE l_imatrix operator &(const imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//------------- real x l_real ------------------------
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix &rv1, const l_rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &rv1, const rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &rv, const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix_slice &sl,const l_rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl, const rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix &rv,const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl1, const rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const rmatrix_slice &sl1, const l_rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	

//------------- l_real x l_real ------------------------
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &rv1, const l_rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &rv1, const l_rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &rv, const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl,const l_rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl, const l_rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix &rv,const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl1, const l_rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl1, const l_rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>);
#else
	throw();
#endif
	

//============== Compare Operator ==========================

//-------------- Matrix - Matrix   -------------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_imatrix &m1,const l_imatrix &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_imatrix &m1,const l_imatrix &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_imatrix &m1,const l_imatrix &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_imatrix &m1,const l_imatrix &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_imatrix &m1,const l_imatrix &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_imatrix &m1,const l_imatrix &m2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_imatrix &m1,const l_imatrix_slice &ms) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_imatrix &m1,const l_imatrix_slice &ms) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_imatrix &m1,const l_imatrix_slice &ms) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_imatrix &m1,const l_imatrix_slice &ms) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_imatrix &m1,const l_imatrix_slice &ms) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_imatrix &m1,const l_imatrix_slice &ms) throw();

//---------------- Matrix - Matrix_slice ----------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw();

//=================== Not Operator =============================

	//! Implementation of standard negation operation
	INLINE bool operator !(const l_imatrix &ms) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const l_imatrix_slice &ms) throw();

//======================== Input / Output ========================

	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const l_imatrix &r) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const l_imatrix_slice &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,l_imatrix &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,l_imatrix_slice &r) throw();

        //! Returns the Identity matrix
        l_imatrix  Id         ( const l_imatrix& );
        //! Returns the transposed matrix
        l_imatrix  transp     ( const l_imatrix& );
        //! Computes the relative diameter \f$ d_{rel}((x)) \f$
        l_real     MaxRelDiam ( const l_imatrix_subv& );
        //! Returns the row dimension
        INLINE int RowLen   ( const l_imatrix& );
        //! Returns the column dimension
        INLINE int ColLen   ( const l_imatrix& );
        //! Returns the row dimension
        INLINE int RowLen   ( const l_imatrix_slice& );
        //! Returns the column dimension
        INLINE int ColLen   ( const l_imatrix_slice& );
        //! Doubles the size of the matrix
        void     DoubleSize ( l_imatrix& );

} // namespace cxsc 

#ifdef _CXSC_INCL_INL
#include "matrix.inl"
#include "l_imatrix.inl"
#endif


#endif

