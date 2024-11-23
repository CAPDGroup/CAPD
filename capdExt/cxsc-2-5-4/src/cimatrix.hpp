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

/* CVS $Id: cimatrix.hpp,v 1.41 2014/01/30 17:23:43 cxsc Exp $ */

#ifndef _CXSC_CIMATRIX_HPP_INCLUDED
#define _CXSC_CIMATRIX_HPP_INCLUDED

#include "xscclass.hpp"
#include "cidot.hpp"
#include "civector.hpp"
#include "except.hpp"
#include "matrix.hpp"
#include "imatrix.hpp"
#include "cmatrix.hpp"


namespace cxsc {

class cimatrix;       // forward declaration
class cimatrix_slice; // forward declaration
class srmatrix;
class srmatrix_slice;
class srmatrix_subv;
class simatrix;
class simatrix_slice;
class simatrix_subv;
class scmatrix;
class scmatrix_slice;
class scmatrix_subv;
class scimatrix;
class scimatrix_slice;
class scimatrix_subv;



// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- class cimatrix_subv (declaration)                                 ----
// ----                                                                   ----
// ---------------------------------------------------------------------------

//! The Data Type cimatrix_subv
/*!
This Data Type provides one column or row of a matrix as a vector.
*/
class cimatrix_subv
{
   friend class civector;
   friend class cimatrix;
   friend class cimatrix_slice;

   private:
      cinterval *dat;
      int lb,ub;
      int size,start,offset; // start=first element index 0..n-1
	
   public:
      //! Returns one row of the matrix as a vector
     friend INLINE cimatrix_subv Row(cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
     throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
     throw();
#endif

      //! Returns one column of the matrix as a vector
      friend INLINE cimatrix_subv Col(cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
      throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
      throw();
#endif
      //! Returns one row of the matrix as a vector
     friend INLINE cimatrix_subv Row(const cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
     throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
     throw();
#endif

      //! Returns one column of the matrix as a vector
      friend INLINE cimatrix_subv Col(const cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
      throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
      throw();
#endif


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
template <class MV,class V> friend  V _mvim(const MV &mv) throw();
template <class MV,class V> friend  V _mvre(const MV &mv) throw();
template <class MV,class V> friend  V _mvdiam(const MV &mv) throw();
template <class MV,class V> friend  V _mvmid(const MV &mv) throw();
template <class MV,class V> friend  V _mvinf(const MV &mv) throw();
template <class MV,class V> friend  V _mvsup(const MV &mv) throw();

 template <class MV,class S> friend 	 MV &_mvssetinf(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvssetsup(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvsusetinf(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvsusetsup(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvssetim(MV &mv, const S &s) throw();
 template <class MV,class S> friend 	 MV &_mvssetre(MV &mv, const S &s) throw();
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
 template <class MV1,class MV2,class S> friend 	 S _mvmvcimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif
 template <class V,class MV,class S> friend 	 S _vmvcimult(const V &rv1, const MV &rv2)
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

	// Real

	// complex

	// interval

#endif

	//----------------- Konstruktoren ----------------------------------

	//! Constructor of class cimatrix_subv
	explicit INLINE cimatrix_subv (cinterval *d, const int &l, const int &u, const int &s, const int &st, const int &o) throw():dat(d),lb(l),ub(u),size(s),start(st),offset(o) { }
        public:
	//! Constructor of class cimatrix_subv
	INLINE cimatrix_subv(const cimatrix_subv &v) throw():dat(v.dat),lb(v.lb),ub(v.ub),size(v.size),start(v.start),offset(v.offset) { }
	public:

	//---------------------- Standardfunktionen ------------------------

	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const scimatrix_subv &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const srmatrix_subv &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const simatrix_subv &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const scmatrix_subv &rv);

	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const scivector &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const srvector &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const scvector &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const sivector &rv);

	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const scivector_slice &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const srvector_slice &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const scvector_slice &rv);
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const sivector_slice &rv);

	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const cimatrix_subv &rv) throw();
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const cinterval &r) throw();
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	cimatrix_subv &operator =(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const civector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// Real
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const rmatrix_subv &rv) throw();

	// complex
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const cvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const cmatrix_subv &rv) throw();

	// interval
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_subv &operator =(const imatrix_subv &rv) throw();

	//! Returns the lower bound of the vector
	friend INLINE int Lb(const cimatrix_subv &rv) throw() { return rv.lb; }
	//! Returns the upper bound of the vector
	friend INLINE int Ub(const cimatrix_subv &rv) throw() { return rv.ub; }
	//! Returns the size of the vector
	friend INLINE int VecLen(const cimatrix_subv &rv) throw() { return rv.size; }

	//! Operator for accessing the single elements of the vector (read-only)
	INLINE cinterval &operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the single elements of the vector
	INLINE cinterval &operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC);
#else
	throw();
#endif

	//! Operator for accessing the whole vector
	INLINE cimatrix_subv &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the vector
	INLINE cimatrix_subv operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the vector
	INLINE cimatrix_subv operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_subv &operator *=(const cinterval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const cinterval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const cinterval &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_subv &operator /=(const cinterval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const scivector &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const scivector &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const scivector &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const scivector &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const scivector_slice &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const scivector_slice &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const scivector_slice &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const scivector_slice &rv);
        //! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const scimatrix_subv &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const scimatrix_subv &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const scimatrix_subv &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const scimatrix_subv &rv);

	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// real
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_subv &operator *=(const real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_subv &operator /=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const srvector &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const srvector &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const srvector &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const srvector_slice &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const srvector_slice &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const srvector_slice &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const srmatrix_subv &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const srmatrix_subv &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const srmatrix_subv &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// complex
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_subv &operator *=(const complex &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const complex &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const complex &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_subv &operator /=(const complex &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const scvector &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const scvector &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const scvector &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const scvector_slice &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const scvector_slice &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const scvector_slice &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const scmatrix_subv &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const scmatrix_subv &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const scmatrix_subv &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const scmatrix_subv &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// interval
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_subv &operator *=(const interval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_subv &operator /=(const interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const sivector &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const sivector &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const sivector &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const sivector &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const sivector_slice &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const sivector_slice &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const sivector_slice &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const sivector_slice &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const simatrix_subv &rv);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const simatrix_subv &rv);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_subv &operator &=(const simatrix_subv &rv);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_subv &operator |=(const simatrix_subv &rv);
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_subv &operator -=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_subv &operator +=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_subv &operator |=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_subv &operator &=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
//#else
//#endif	

};


//! Returns one row of the matrix as a vector
INLINE cimatrix_subv Row(cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
//! Returns one column of the matrix as a vector
INLINE cimatrix_subv Col(cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif




class cimatrix_slice;

// ---------------------------------------------------------------------------
// ----                                                                   ----
// ---- class cimatrix (declaration)                                      ----
// ----                                                                   ----
// ---------------------------------------------------------------------------


//! The Data Type cimatrix
/*!
\sa rmatrix
*/
class cimatrix
{
	friend class cimatrix_slice;
	friend class cimatrix_subv;
	private:
	cinterval *dat;
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
 template <class M,class E> friend 	 E _mdiam(const M &m) throw();
 template <class M,class E> friend 	 E _mmid(const M &m) throw();
 template <class M,class E> friend 	 E _mre(const M &m) throw();
 template <class M,class E> friend 	 E _mim(const M &m) throw();
	friend INLINE rmatrix SupRe(const cimatrix &v) throw();
	friend INLINE rmatrix SupIm(const cimatrix &v) throw();
	friend INLINE rmatrix InfRe(const cimatrix &v) throw();
	friend INLINE rmatrix InfIm(const cimatrix &v) throw();
 template <class M1,class M2> friend 	 M1 &_mmsetre(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M1,class M2> friend 	 M1 &_mmsetim(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
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
 template <class M1,class M2,class E> friend 	 E _mmcimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M1,class M2,class S> friend 	 M1 &_mmcimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
 template <class M,class MS,class E> friend 	 E _mmscimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmcimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS,class S> friend 	 M &_mmscimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmscimult(const MS1 &ms1, const MS2 &ms2)
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
 template <class V,class M,class S> friend 	 V &_vmcimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class VS,class M,class S> friend 	 VS &_vsmcimultassign(VS &v,const M &m)
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

	// Real

	//--- Real --------- matrix-matrix ----------------------

	//--- Real --------- matrix-scalar ----------------------

	//--- Real --------- matrix-vector ----------------------
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

	// interval

	//--- interval --------- matrix-matrix ----------------------

	//--- interval --------- matrix-scalar ----------------------


	//--- interval --------- matrix-vector ----------------------

	// complex

	//--- complex --------- matrix-matrix ----------------------

	//--- complex --------- matrix-scalar ----------------------

	//--- complex --------- matrix-vector ----------------------

	// --- complex x real ----------------
	// -- complex x interval ----------------------
	// ---- complex x interval --- scalar--------


	// ---- complex x interval --- vector --------
	// ---- complex x interval --- matrix ------------

	// complex x complex --------------------


#endif
	
	//--------------------------  Konstruktoren ----------------------------

// cinterval
	//! Constructor of class cimatrix
	INLINE cimatrix(const cimatrix &rm) throw();
	//! Constructor of class cimatrix
	INLINE cimatrix(const cimatrix_slice &rm) throw();
	//! Constructor of class cimatrix
	INLINE cimatrix(const scimatrix &rm);
	//! Constructor of class cimatrix
	INLINE cimatrix(const scimatrix_slice &rm);
	//! Constructor of class cimatrix
	INLINE cimatrix() throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const int &m1, const int &n1, const int &m2, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const civector &v) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const civector_slice &v) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const cinterval &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cinterval &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cimatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cimatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const scimatrix &m);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const scimatrix_slice &ms);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const civector &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const civector_slice &v) throw();
//  real
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const real &r) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const rmatrix &rm) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const rmatrix_slice &rm) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const srmatrix &rm);
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const srmatrix_slice &rm);
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const rvector &v) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const rvector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const rmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const rmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const srmatrix &m);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const srmatrix_slice &ms);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const rvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const rvector_slice &v) throw();

//  complex
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const complex &r) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const cmatrix &rm) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const cmatrix_slice &rm) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const scmatrix &rm);
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const scmatrix_slice &rm);
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const cvector &v) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const cvector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cmatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cmatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const scmatrix &m);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const scmatrix_slice &ms);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cvector &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const cvector_slice &v) throw();

//  interval
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const interval &r) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const imatrix &rm) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const imatrix_slice &rm) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const simatrix &rm);
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const simatrix_slice &rm);
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const ivector &v) throw();
	//! Constructor of class cimatrix
	explicit INLINE cimatrix(const ivector_slice &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const imatrix &m) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const imatrix_slice &ms) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const simatrix &m);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const simatrix_slice &ms);
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const ivector &v) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix &operator =(const ivector_slice &v) throw();

	//--------------------------- Destruktoren -----------------------------

	INLINE ~cimatrix() throw() { delete [] dat; }

	//--------------------------- Operatoren -----------------------------
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const scimatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const scimatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const scimatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const scimatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const scimatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const scimatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const scimatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const scimatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const scimatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const scimatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const srmatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const srmatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const srmatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const srmatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const srmatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const srmatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const srmatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const srmatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const srmatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const srmatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const scmatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const scmatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const scmatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const scmatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const scmatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const scmatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const scmatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const scmatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const scmatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const scmatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const simatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(const simatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const simatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix &operator -=(const simatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const simatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix &operator |=(const simatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const simatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix &operator &=(const simatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const simatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(const simatrix_slice &m1);

        //! Computes permutation of matrix according to permutation vectors, C=PAQ
        INLINE cimatrix operator()(const intvector& p, const intvector& q);
        //! Computes permutation of matrix according to permutation matrices, C=PAQ
        INLINE cimatrix operator()(const intmatrix& P, const intmatrix& Q);
        //! Computes permutation of matrix according to permutation vector, C=PA
        INLINE cimatrix operator()(const intvector& p);
        //! Computes permutation of matrix according to permutation matrix, C=PAQ
        INLINE cimatrix operator()(const intmatrix& P);

	//------------------------- Standardfunktionen -------------------------

	//! Operator for accessing a single row of the matrix
	INLINE cimatrix_subv operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE cimatrix_subv operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single row of the matrix
	INLINE cimatrix_subv operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE cimatrix_subv operator [](const cxscmatrix_column &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif

	//! Operator for accessing the whole matrix
	INLINE cimatrix &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE cimatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE cimatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();
//#else
//#endif
};

		
//! The Data Type cimatrix_slice
/*!
This data type represents a partial cimatrix.

\sa cimatrix
*/
class cimatrix_slice
{
	friend class cimatrix;
	private:
	cinterval *dat;
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
 template <class MS,class E> friend 	 E _msre(const MS &ms) throw();
 template <class MS,class E> friend 	 E _msim(const MS &ms) throw();
	friend INLINE rmatrix SupRe(const cimatrix_slice &v) throw();
	friend INLINE rmatrix SupIm(const cimatrix_slice &v) throw();
	friend INLINE rmatrix InfRe(const cimatrix_slice &v) throw();
	friend INLINE rmatrix InfIm(const cimatrix_slice &v) throw();
 template <class MS1,class M2> friend 	 MS1 &_msmsetre(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class M2> friend 	 MS1 &_msmsetim(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmssetre(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
 template <class MS1,class MS2> friend 	 MS1 &_msmssetim(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
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
 template <class M,class MS,class E> friend 	 E _mmscimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class MS,class M,class E> friend 	 E _msmcimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
 template <class M,class MS,class S> friend 	 M &_mmscimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
 template <class MS1,class MS2,class E> friend 	 E _msmscimult(const MS1 &ms1, const MS2 &ms2)
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

	// Real

	//--- Real ------------ matrix-scalar -----------
	
	
	//--- Real ------------ matrix-vector -----------

	//--- Real ------------ matrix-matrix -----------

	// interval

	//--- interval ------------ matrix-scalar -----------
	
	
	//--- interval ------------ matrix-vector -----------

	//--- interval ------------ matrix-matrix -----------

	// complex

	//--- complex ------------ matrix-scalar -----------
	
	
	//--- complex ------------ matrix-vector -----------


	//--- complex ------------ matrix-matrix -----------


#endif

	//--------------- Konstruktoren ----------------------------------------

	//! Constructor of class cimatrix_slice
	explicit INLINE cimatrix_slice(cimatrix &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(l1-a.lb1),offset2(l2-a.lb2),mxsize(a.xsize),mysize(a.ysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	//! Constructor of class cimatrix_slice
	explicit INLINE cimatrix_slice(cimatrix_slice &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),offset1(a.offset1+l1-a.start1),offset2(a.offset2+l2-a.start2),mxsize(a.mxsize),mysize(a.mysize),start1(l1),end1(u1),start2(l2),end2(u2),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	public: 
	//! Constructor of class cimatrix_slice
	INLINE cimatrix_slice(const cimatrix_slice &ms) throw():dat(ms.dat),offset1(ms.offset1),offset2(ms.offset2),mxsize(ms.mxsize),mysize(ms.mysize),start1(ms.start1),end1(ms.end1),start2(ms.start2),end2(ms.end2),sxsize(ms.sxsize),sysize(ms.sysize) { }
	public:

	//---------------- Standardfunktionen -----------------------------------

	friend INLINE civector::civector(const cimatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ);
#else
	throw();
#endif
	friend INLINE cimatrix::cimatrix(const cimatrix_slice &) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cinterval &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const scimatrix &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const scimatrix_subv &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const scimatrix_slice &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const civector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cimatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	// real
        //! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const srmatrix &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const srmatrix_subv &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const srmatrix_slice &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const real &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	// interval
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const simatrix_subv &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const simatrix &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const simatrix_slice &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const interval &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const imatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	// complex
        //! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const scmatrix &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const scmatrix_subv &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const scmatrix_slice &v);
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const complex &r) throw();
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of standard assigning operator
	INLINE cimatrix_slice &operator =(const cmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Operator for accessing a single row of the matrix
	INLINE cimatrix_subv operator [](const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE cimatrix_subv operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif

	//! Operator for accessing a single row of the matrix
	INLINE cimatrix_subv operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif
	//! Operator for accessing a single column of the matrix
	INLINE cimatrix_subv operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT);
#else
	throw();
#endif

	//! Operator for accessing the whole matrix
	INLINE cimatrix_slice &operator ()() throw() { return *this; }
	//! Operator for accessing a part of the matrix
	INLINE cimatrix_slice operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	//! Operator for accessing a part of the matrix
	INLINE cimatrix_slice operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG);
#else
	throw();
#endif
	INLINE operator void*() throw();

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const cinterval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const cinterval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const cinterval &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_slice &operator /=(const cinterval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const scimatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const scimatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const scimatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const scimatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const scimatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const scimatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const scimatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const scimatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const scimatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const scimatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const real &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_slice &operator /=(const real &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const srmatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const srmatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const srmatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const srmatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const srmatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const srmatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const srmatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const srmatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const srmatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const srmatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const complex &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const complex &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const complex &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_slice &operator /=(const complex &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const scmatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const scmatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const scmatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const scmatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const scmatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const scmatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const scmatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const scmatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const scmatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const scmatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const interval &c) throw();
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const interval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix_slice &operator /=(const interval &c) throw();
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const simatrix &m1);
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const simatrix_slice &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const simatrix &m1);
	//! Implementation of substraction and allocation operation
	INLINE cimatrix_slice &operator -=(const simatrix_slice &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const simatrix &m1);
	//! Implementation of hull and allocation operation
	INLINE cimatrix_slice &operator |=(const simatrix_slice &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const simatrix &m1);
	//! Implementation of intersection and allocation operation
	INLINE cimatrix_slice &operator &=(const simatrix_slice &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const simatrix &m1);
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const simatrix_slice &m1);

	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix_slice &operator +=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix_slice &operator -=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix_slice &operator |=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix_slice &operator &=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix_slice &operator *=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
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
	INLINE civector operator /(const cimatrix_subv &rv, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix_subv &rv, const cinterval &s) throw();
	//! Implementation of multiplication operation
	INLINE civector operator *(const cinterval &s, const cimatrix_subv &rv) throw();
	//! Returns the absolute value of the matrix
	INLINE ivector abs(const cimatrix_subv &mv) throw();
	//! Returns the diameter of the matrix
	INLINE cvector diam(const cimatrix_subv &mv) throw();
	//! Returns the middle of the matrix
	INLINE cvector mid(const cimatrix_subv &mv) throw();
	//! Returns the infimum of the matrix
	INLINE cvector Inf(const cimatrix_subv &mv) throw();
	//! Returns the supremum of the matrix
	INLINE cvector Sup(const cimatrix_subv &mv) throw();
	//! Returns the imaginary part of the matrix
	INLINE ivector Im(const cimatrix_subv &mv) throw();
	//! Returns the real part of the matrix
	INLINE ivector Re(const cimatrix_subv &mv) throw();
	//! Returns the supremum of real part of the matrix
	INLINE rmatrix SupRe(const cimatrix &v) throw();
	//! Returns the supremum of imaginary part of the matrix
	INLINE rmatrix SupIm(const cimatrix &v) throw();
	//! Returns the infimum of real part of the matrix
	INLINE rmatrix InfRe(const cimatrix &v) throw();
	//! Returns the infimum of imaginary part of the matrix
	INLINE rmatrix InfIm(const cimatrix &v) throw();
	//! Returns the supremum of real part of the matrix
	INLINE rmatrix SupRe(const cimatrix_slice &v) throw();
	//! Returns the supremum of imaginary part of the matrix
	INLINE rmatrix SupIm(const cimatrix_slice &v) throw();
	//! Returns the infimum of real part of the matrix
	INLINE rmatrix InfRe(const cimatrix_slice &v) throw();
	//! Returns the infimum of imaginary part of the matrix
	INLINE rmatrix InfIm(const cimatrix_slice &v) throw();
	//! Returns the matrix with the new given infimum value
	INLINE cimatrix_subv &SetInf(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE cimatrix_subv &SetSup(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE cimatrix_subv &UncheckedSetInf(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE cimatrix_subv &UncheckedSetSup(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the imaginary parts of the matrix
	INLINE cimatrix_subv &SetIm(cimatrix_subv &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the real parts of the matrix
	INLINE cimatrix_subv &SetRe(cimatrix_subv &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the matrix with the new given supremum value
	INLINE cimatrix_subv &SetSup(cimatrix_subv &iv,const complex &r) throw();
	//! Returns the matrix with the new given infimum value
	INLINE cimatrix_subv &SetInf(cimatrix_subv &iv,const complex &r) throw();
	//! Returns the matrix with the new unchecked given supremum value
	INLINE cimatrix_subv &UncheckedSetSup(cimatrix_subv &iv,const complex &r) throw();
	//! Returns the matrix with the new unchecked given infimum value
	INLINE cimatrix_subv &SetUncheckedInf(cimatrix_subv &iv,const complex &r) throw();
	//! Sets componentwise the real parts of the matrix
	INLINE cimatrix_subv &SetRe(cimatrix_subv &iv,const interval &r) throw();
	//! Sets componentwise the imaginary parts of the matrix
	INLINE cimatrix_subv &SetIm(cimatrix_subv &iv,const interval &r) throw();

//======================== Vector / Vector ========================

	
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cimatrix_subv &rv2)
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
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const civector_slice & sl1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const civector_slice &sl,const cimatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cinterval operator *(const cimatrix_subv &mv,const civector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of addition operation
	INLINE civector operator +(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const civector_slice &sl,const cimatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE civector operator +(const cimatrix_subv &mv,const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE civector operator -(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const civector_slice &sl,const cimatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE civector operator -(const cimatrix_subv &mv,const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//  real

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const rvector_slice &sl2)
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
	 void accumulate(cidotprecision &dp, const rvector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rmatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const rvector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
// complex

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const cvector_slice &sl2)
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
	 void accumulate(cidotprecision &dp, const cvector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cvector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
// interval

	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const cimatrix_subv & rv1, const ivector_slice &sl2)
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
	 void accumulate(cidotprecision &dp, const ivector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	 void accumulate(cidotprecision &dp, const ivector_slice & sl1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const cmatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! The accurate scalar product of the last two arguments added to the value of the first argument
	void accumulate(cidotprecision &dp, const imatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//====================================================================
//===================== Matrix Functions =============================

	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cimatrix _imatrix(const cimatrix &rm) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cimatrix _imatrix(const civector &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cimatrix _imatrix(const civector_slice &v) throw();
	//! Deprecated typecast, which only exist for the reason of compatibility with older versions of C-XSC
	INLINE cimatrix _imatrix(const cinterval &r) throw();

	//! Returns the lower bound index
	INLINE int Lb(const cimatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const cimatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the lower bound index
	INLINE int Lb(const cimatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Returns the upper bound index
	INLINE int Ub(const cimatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the lower bound index
	INLINE cimatrix &SetLb(cimatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Sets the upper bound index
	INLINE cimatrix &SetUb(cimatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(cimatrix &A) throw();
	//! Resizes the matrix
	INLINE void Resize(cimatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif
	//! Resizes the matrix
	INLINE void Resize(cimatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES);
#else
	throw();
#endif

	//! Returns the absolute value of the matrix
	INLINE imatrix abs(const cimatrix &m) throw();
	//! Returns the absolute value of the matrix
	INLINE imatrix abs(const cimatrix_slice &ms) throw();
	//! Returns the rounded diameter of the matrix
	INLINE cmatrix diam(const cimatrix &m) throw();
	//! Returns the rounded diameter of the matrix
	INLINE cmatrix diam(const cimatrix_slice &m) throw();
	//! Returns the rounded middle of the matrix
	INLINE cmatrix mid(const cimatrix &m) throw();
	//! Returns the rounded middle of the matrix
	INLINE cmatrix mid(const cimatrix_slice &m) throw();
	//! Returns the infimum of the matrix
	INLINE cmatrix Inf(const cimatrix &m) throw();
	//! Returns the supremum of the matrix
	INLINE cmatrix Sup(const cimatrix &m) throw();
	//! Returns the infimum of the matrix
	INLINE cmatrix Inf(const cimatrix_slice &m) throw();
	//! Returns the supremum of the matrix
	INLINE cmatrix Sup(const cimatrix_slice &m) throw();
	//! Returns the matrix with the new given infimum value
	INLINE cimatrix &SetInf(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given infimum value
	INLINE cimatrix_slice &SetInf(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given infimum value
	INLINE cimatrix &SetInf(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given infimum value
	INLINE cimatrix_slice &SetInf(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE cimatrix &SetSup(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE cimatrix_slice &SetSup(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE cimatrix &SetSup(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new given supremum value
	INLINE cimatrix_slice &SetSup(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE cimatrix &UncheckedSetInf(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE cimatrix_slice &UncheckedSetInf(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE cimatrix &UncheckedSetInf(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given infimum value
	INLINE cimatrix_slice &UncheckedSetInf(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE cimatrix &UncheckedSetSup(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE cimatrix_slice &UncheckedSetSup(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE cimatrix &UncheckedSetSup(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the matrix with the new unchecked given supremum value
	INLINE cimatrix_slice &UncheckedSetSup(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Returns the imaginary part of the matrix
	INLINE imatrix Im(const cimatrix &m) throw();
	//! Returns the real part of the matrix
	INLINE imatrix Re(const cimatrix &m) throw();
	//! Returns the imaginary part of the matrix
	INLINE imatrix Im(const cimatrix_slice &m) throw();
	//! Returns the real part of the matrix
	INLINE imatrix Re(const cimatrix_slice &m) throw();
	//! Sets componentwise the imaginary parts of the matrix
	INLINE cimatrix &SetIm(cimatrix &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the imaginary parts of the matrix
	INLINE cimatrix_slice &SetIm(cimatrix_slice &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the imaginary parts of the matrix
	INLINE cimatrix &SetIm(cimatrix &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the imaginary parts of the matrix
	INLINE cimatrix_slice &SetIm(cimatrix_slice &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the real parts of the matrix
	INLINE cimatrix &SetRe(cimatrix &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the real parts of the matrix
	INLINE cimatrix_slice &SetRe(cimatrix_slice &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the real parts of the matrix
	INLINE cimatrix &SetRe(cimatrix &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Sets componentwise the real parts of the matrix
	INLINE cimatrix_slice &SetRe(cimatrix_slice &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

//===================== Matrix / Scalar ===============================

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const cimatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const cimatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m,const cinterval &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms,const cinterval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix &m,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix_slice &ms, const cinterval &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix &operator /=(cimatrix &m,const cinterval &c) throw();
	
//------------ real - cimatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const real &c, const cimatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const real &c, const cimatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m,const real &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms,const real &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m,const real &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix &m,const real &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix_slice &ms, const real &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix &operator /=(cimatrix &m,const real &c) throw();
//----------------- rmatrix - cinterval ----------------

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const rmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const rmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const rmatrix &m,const cinterval &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const rmatrix_slice &ms,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const rmatrix &m,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const rmatrix_slice &ms, const cinterval &c) throw();
	
//------------ complex - cimatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const complex &c, const cimatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const complex &c, const cimatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m,const complex &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms,const complex &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m,const complex &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix &m,const complex &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix_slice &ms, const complex &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix &operator /=(cimatrix &m,const complex &c) throw();
//----------------- cmatrix - cinterval ----------------

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const cmatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const cmatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix &m,const cinterval &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix_slice &ms,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cmatrix &m,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cmatrix_slice &ms, const cinterval &c) throw();
	
//------------ interval - cimatrix -----------------------------------------------

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const interval &c, const cimatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const interval &c, const cimatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m,const interval &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms,const interval &c) throw();
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m,const interval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix &m,const interval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const cimatrix_slice &ms, const interval &c) throw();
	//! Implementation of division and allocation operation
	INLINE cimatrix &operator /=(cimatrix &m,const interval &c) throw();
//----------------- imatrix - cinterval ----------------

	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const imatrix &m) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cinterval &c, const imatrix_slice &ms) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix &m,const cinterval &c) throw();
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix_slice &ms,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const imatrix &m,const cinterval &c) throw();
	//! Implementation of division operation
	INLINE cimatrix operator /(const imatrix_slice &ms, const cinterval &c) throw();
	

//============================ Matrix / Vector ===================================


	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix &m,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix_slice &ms,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector &operator *=(civector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE civector &operator *=(civector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	//! Implementation of multiplication operation
	INLINE civector operator *(const civector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const civector_slice &v,const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- real -------------------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const rvector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const rvector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const rvector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- complex -------------------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const cvector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const cvector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const cvector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix &m,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix_slice &ms,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//----------------- interval -------------------------------------

	//! Implementation of multiplication operation
	INLINE civector operator *(const ivector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const ivector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const ivector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE civector operator *(const cimatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	

//================ Matrix / Matrix ============================

	//! Implementation of positive sign operation
	INLINE const cimatrix &operator +(const cimatrix &m1) throw();
	//! Implementation of positive sign operation
	INLINE cimatrix operator +(const cimatrix_slice &ms) throw();
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of negative sign operation
	INLINE cimatrix operator -(const cimatrix &m) throw();
	//! Implementation of negative sign operation
	INLINE cimatrix operator -(const cimatrix_slice &ms) throw();
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- rmatrix-cimatrix ------------------
	//! Implementation of addition operation
	INLINE cimatrix operator +(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const rmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const rmatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const rmatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const rmatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const rmatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const rmatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const rmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- cmatrix-cimatrix ------------------
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms, const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- imatrix-cimatrix ------------------
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cimatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition and allocation operation
	INLINE cimatrix &operator +=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction and allocation operation
	INLINE cimatrix &operator -=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication and allocation operation
	INLINE cimatrix &operator *=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cimatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the convex hull of the arguments to the first argument
	INLINE cimatrix &operator |=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cimatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Allocates the intersection of the arguments to the first argument
	INLINE cimatrix &operator &=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//---------- cmatrix-imatrix ------------------
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const cmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of addition operation
	INLINE cimatrix operator +(const imatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const cmatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of subtraction operation
	INLINE cimatrix operator -(const imatrix_slice &ms1,const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix &m1, const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix &m1, const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix_slice &ms, const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const cmatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Implementation of multiplication operation
	INLINE cimatrix operator *(const imatrix_slice &ms1, const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const imatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const cmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	//! Returns the intersection of the arguments
	INLINE cimatrix operator &(const imatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM);
#else
	throw();
#endif
	
//------------- real x complex ------------------------
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix &rv1, const cmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &rv1, const rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &rv, const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix_slice &sl,const cmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &sl, const rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix &rv,const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &sl1, const rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const rmatrix_slice &sl1, const cmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	

//------------- complex x complex ------------------------
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &rv1, const cmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &rv1, const cmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &rv, const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &sl,const cmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &sl, const cmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix &rv,const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &sl1, const cmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	//! Returns the convex hull of the arguments
	INLINE cimatrix operator |(const cmatrix_slice &sl1, const cmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>);
#else
	throw();
#endif
	

//============== Compare Operator ==========================

//-------------- Matrix - Matrix   -------------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const cimatrix &m1,const cimatrix &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cimatrix &m1,const cimatrix &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const cimatrix &m1,const cimatrix &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const cimatrix &m1,const cimatrix &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const cimatrix &m1,const cimatrix &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const cimatrix &m1,const cimatrix &m2) throw();
	//! Implementation of standard equality operation
	INLINE bool operator ==(const cimatrix &m1,const cimatrix_slice &ms) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cimatrix &m1,const cimatrix_slice &ms) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const cimatrix &m1,const cimatrix_slice &ms) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const cimatrix &m1,const cimatrix_slice &ms) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const cimatrix &m1,const cimatrix_slice &ms) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const cimatrix &m1,const cimatrix_slice &ms) throw();

//---------------- Matrix - Matrix_slice ----------------------

	//! Implementation of standard equality operation
	INLINE bool operator ==(const cimatrix_slice &m1,const cimatrix_slice &m2) throw();
	//! Implementation of standard negated equality operation
	INLINE bool operator !=(const cimatrix_slice &m1,const cimatrix_slice &m2) throw();
	//! Implementation of standard less-than operation
	INLINE bool operator <(const cimatrix_slice &m1,const cimatrix_slice &m2) throw();
	//! Implementation of standard less-or-equal-than operation
	INLINE bool operator <=(const cimatrix_slice &m1,const cimatrix_slice &m2) throw();
	//! Implementation of standard greater-than operation
	INLINE bool operator >(const cimatrix_slice &m1,const cimatrix_slice &m2) throw();
	//! Implementation of standard greater-or-equal-than operation
	INLINE bool operator >=(const cimatrix_slice &m1,const cimatrix_slice &m2) throw();

//=================== Not Operator =============================

	//! Implementation of standard negation operation
	INLINE bool operator !(const cimatrix &ms) throw();
	//! Implementation of standard negation operation
	INLINE bool operator !(const cimatrix_slice &ms) throw();

//======================== Input / Output ========================

	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const cimatrix &r) throw();
	//! Implementation of standard output method
	INLINE std::ostream &operator <<(std::ostream &s,const cimatrix_slice &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,cimatrix &r) throw();
	//! Implementation of standard input method
	INLINE std::istream &operator >>(std::istream &s,cimatrix_slice &r) throw();

        //! Returns Ostrowski's comparison matrix
        rmatrix  CompMat    ( const cimatrix& );
        //! Returns the Identity matrix
        cimatrix  Id         ( const cimatrix& );
        //! Returns the transposed matrix
        cimatrix  transp     ( const cimatrix& );
        //! Returns the row dimension
        INLINE int      RowLen     ( const cimatrix& );
        //! Returns the column dimension
        INLINE int      ColLen     ( const cimatrix& );
        //! Returns the row dimension
        INLINE int      RowLen     ( const cimatrix_slice& );
        //! Returns the column dimension
        INLINE int      ColLen     ( const cimatrix_slice& );
        //! Doubles the size of the matrix
        void     DoubleSize ( cimatrix& );

} // namespace cxsc 


#ifdef _CXSC_INCL_INL
#include "matrix.inl"
#include "cimatrix.inl"
#endif


#ifdef CXSC_USE_BLAS
#define _CXSC_BLAS_CIMATRIX
#include "cxsc_blas.inl"
#endif

#endif
