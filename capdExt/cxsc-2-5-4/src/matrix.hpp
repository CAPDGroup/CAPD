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

/* CVS $Id: matrix.hpp,v 1.25 2014/01/30 17:23:47 cxsc Exp $ */

#ifndef _CXSC_MATRIX_HPP_INCLUDED
#define _CXSC_MATRIX_HPP_INCLUDED

#include "except.hpp"

namespace cxsc {

class cxscmatrix_column
{
	private:
	int c;
	public:
	inline cxscmatrix_column(const int &i) throw():c(i) { }
	inline int col() const throw() { return c; }
};

#define ROW 1
#define COL 2
inline int Row(const int &i) throw() { return i; }
inline cxscmatrix_column Col(const int &i) throw() { return cxscmatrix_column(i); }

} // namespace cxsc

#include <iostream>
#include "vector.inl"

namespace cxsc {

class l_real;
class l_interval;
	
template <class S,class M>
 void _smconstr(S &s,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>,ERROR__USE_OF_UNINITIALIZED_OBJ<M>);
#else
	throw();
#endif

template <class V,class M,class S>
 void _vmconstr(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif

template <class V,class MS,class S>
 void _vmsconstr(V &v,const MS &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<MS>);
#else
	throw();
#endif


	
	template <class M>
	std::ostream &_mout(std::ostream &s,const M &r) throw();

	template <class M>
	std::istream &_min(std::istream &s,M &r) throw();

	template <class MS>
	std::ostream &_msout(std::ostream &s,const MS &r) throw();

	template <class MS>
	std::istream &_msin(std::istream &s,MS &r) throw();

	template <class M1,class M2,class S>
	 M1 &_mmassign(M1 &m1,const M2 &m,S m2) throw();

	template <class M,class MS2,class S>
	 M &_mmsassign(M &m,const MS2 &ms) throw();

	template <class MS,class M>
	 MS &_msmassign(MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class MS1,class MS2>
	 MS1 &_msmsassign(MS1 &ms1,const MS2 &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif

	template <class M,class S>
	 M &_msassign(M &m,const S &r) throw();

	template <class MS,class S>
	 MS &_mssassign(MS &ms,const S &r) throw();

template <class V,class M,class S>
 V &_vmassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__TYPE_CAST_OF_THICK_OBJ<M>);
#else
	throw();
#endif

// result depends on the matrix (left arg), as defined in the reference to C-XSC
template <class M,class V,class S>
 M &_mvassign(M &m,const V &v) throw();

	template <class M>
	 int _mlb(const M &m, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
	
	template <class M>
	 int _mub(const M &m, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
	
	template <class MS>
	 int _mslb(const MS &ms, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<MS>);
#else
	throw();
#endif
	
	template <class MS>
	 int _msub(const MS &ms, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<MS>);
#else
	throw();
#endif
	
	template <class M>
	 M &_msetlb(M &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
	
	template <class M>
	 M &_msetub(M &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_ROW_OR_COL<M>);
#else
	throw();
#endif
	
	template <class M>
	 void _mresize(M &A) throw();

	template <class M,class S>
	 void _mresize(M &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<M>);
#else
	throw();
#endif

	template <class M,class S>
	 void _mresize(M &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<M>);
#else
	throw();
#endif
	
	template <class M,class E>
	 E _mabs(const M &m) throw();

	template <class MS,class E>
	 E _msabs(const MS &ms) throw();
	
	template <class M,class E>
	 E _mdiam(const M &m) throw();

	template <class MS,class E>
	 E _msdiam(const MS &ms) throw();
	
	template <class M,class E>
	 E _mmid(const M &m) throw();

	template <class MS,class E>
	 E _msmid(const MS &ms) throw();
	
	template <class M,class E>
	 E _minf(const M &m) throw();

	template <class MS,class E>
	 E _msinf(const MS &ms) throw();
	
	template <class M,class E>
	 E _msup(const M &m) throw();

	template <class MS,class E>
	 E _mssup(const MS &ms) throw();
	
	template <class M1,class M2>
	 M1 &_mmsetinf(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class MS1,class M2>
	 MS1 &_msmsetinf(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class MS2>
	 M1 &_mmssetinf(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class MS1,class MS2>
	 MS1 &_msmssetinf(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class M2>
	 M1 &_mmsetsup(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class MS1,class M2>
	 MS1 &_msmsetsup(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class MS2>
	 M1 &_mmssetsup(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class MS1,class MS2>
	 MS1 &_msmssetsup(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class M2>
	 M1 &_mmusetinf(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class MS1,class M2>
	 MS1 &_msmusetinf(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class MS2>
	 M1 &_mmsusetinf(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class MS1,class MS2>
	 MS1 &_msmsusetinf(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class M2>
	 M1 &_mmusetsup(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class MS1,class M2>
	 MS1 &_msmusetsup(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class MS2>
	 M1 &_mmsusetsup(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class MS1,class MS2>
	 MS1 &_msmsusetsup(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class M2>
	 M1 &_mmsetre(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class MS1,class M2>
	 MS1 &_msmsetre(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class MS2>
	 M1 &_mmssetre(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class MS1,class MS2>
	 MS1 &_msmssetre(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class M2>
	 M1 &_mmsetim(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class MS1,class M2>
	 MS1 &_msmsetim(MS1 &ms1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M1,class MS2>
	 M1 &_mmssetim(M1 &m1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class MS1,class MS2>
	 MS1 &_msmssetim(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
	
	template <class M,class E>
	 E _mim(const M &m) throw();

	template <class MS,class E>
	 E _msim(const MS &ms) throw();
	
	template <class M,class E>
	 E _mre(const M &m) throw();

	template <class MS,class E>
	 E _msre(const MS &ms) throw();
	
	template <class M1,class M2,class E>
	 E _mmsect(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class M,class MS,class E>
	 E _mmssect(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmconv(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class M,class MS,class E>
	 E _mmsconv(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmplus(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
	
	template <class M,class MS,class E>
	 E _mmsplus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M>
	 M _mminus(const M &m) throw();
	
	template <class MS,class E>
	 E _msminus(const MS &ms) throw();
	
	template <class M1,class M2,class E>
	 E _mmminus(const M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif


	template <class M1,class M2>
	 M1 &_mmconvassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
		
	template <class M,class MS>
	 M &_mmsconvassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
		
	template <class MS,class M>
	 MS &_msmconvassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
		
	template <class MS1,class MS2>
	 MS1 &_msmsconvassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
		
	template <class M1,class M2>
	 M1 &_mmsectassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
		
	template <class M,class MS>
	 M &_mmssectassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
		
	template <class MS,class M>
	 MS &_msmsectassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
		
	template <class MS1,class MS2>
	 MS1 &_msmssectassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
		
	template <class M1,class M2>
	 M1 &_mmplusassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
		
	template <class M,class MS>
	 M &_mmsplusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
		
	template <class MS,class M>
	 MS &_msmplusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
		
	template <class MS1,class MS2>
	 MS1 &_msmsplusassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
		
	template <class M,class MS,class E>
	 E _mmsminus(const M &m,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS,class M,class E>
	 E _msmminus(const MS &ms,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS1,class MS2,class E>
	 E _msmsminus(const MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class M1,class M2>
	 M1 &_mmminusassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif
		
	template <class M,class MS>
	 M &_mmsminusassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif
		
	template <class MS,class M>
	 MS &_msmminusassign(MS &ms,const M &m1)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif
		
	template <class MS1,class MS2>
	 MS1 &_msmsminusassign(MS1 &ms1,const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS1>);
#else
	throw();
#endif
		
	template <class M1,class M2,class E>
	 E _mmmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmlmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmlimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmcmult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M1,class M2,class E>
	 E _mmcimult(const M1 &m1, const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M1,class M2,class S>
	 M1 &_mmmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class M1,class M2,class S>
	 M1 &_mmimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class M1,class M2,class S>
	 M1 &_mmlmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class M1,class M2,class S>
	 M1 &_mmlimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class M1,class M2,class S>
	 M1 &_mmcmultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class M1,class M2,class S>
	 M1 &_mmcimultassign(M1 &m1,const M2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M1>);
#else
	throw();
#endif

	template <class M,class MS,class E>
	 E _mmsmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M,class MS,class E>
	 E _mmsimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M,class MS,class E>
	 E _mmslmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M,class MS,class E>
	 E _mmslimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M,class MS,class E>
	 E _mmscmult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M,class MS,class E>
	 E _mmscimult(const M &m1, const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS,class M,class E>
	 E _msmmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS,class M,class E>
	 E _msmimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS,class M,class E>
	 E _msmlmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS,class M,class E>
	 E _msmlimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS,class M,class E>
	 E _msmcmult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS,class M,class E>
	 E _msmcimult(const MS &ms, const M &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M,class MS,class S>
	 M &_mmsmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class MS,class S>
	 M &_mmsimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class MS,class S>
	 M &_mmslmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class MS,class S>
	 M &_mmslimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class MS,class S>
	 M &_mmscmultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class MS,class S>
	 M &_mmscimultassign(M &m1,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class MS1,class MS2,class E>
	 E _msmsmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS1,class MS2,class E>
	 E _msmsimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS1,class MS2,class E>
	 E _msmslmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS1,class MS2,class E>
	 E _msmslimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS1,class MS2,class E>
	 E _msmscmult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class MS1,class MS2,class E>
	 E _msmscimult(const MS1 &ms1, const MS2 &ms2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif
	
	template <class S,class M,class E>
	 E _smmult(const S &c, const M &m) throw();

	template <class M,class S>
	 M &_msmultassign(M &m,const S &c) throw();
	
	template <class S,class MS,class E>
	 E _smsmult(const S &c, const MS &ms) throw();

	template <class MS,class S>
	 MS &_mssmultassign(MS &ms,const S &c) throw();

	template <class M,class V,class E>
	 E _mvmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class V,class E>
	 E _mvimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class V,class E>
	 E _mvlmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class V,class E>
	 E _mvlimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class V,class E>
	 E _mvcmult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class M,class V,class E>
	 E _mvcimult(const M &m,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class E>
	 E _vmmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class E>
	 E _vmimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class E>
	 E _vmcmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class E>
	 E _vmcimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class E>
	 E _vmlmult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class E>
	 E _vmlimult(const V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class S>
	 V &_vmmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class S>
	 V &_vmimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class S>
	 V &_vmlmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class S>
	 V &_vmlimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class S>
	 V &_vmcmultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class V,class M,class S>
	 V &_vmcimultassign(V &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class VS,class M,class S>
	 VS &_vsmmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class VS,class M,class S>
	 VS &_vsmimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class VS,class M,class S>
	 VS &_vsmlmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class VS,class M,class S>
	 VS &_vsmlimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class VS,class M,class S>
	 VS &_vsmcmultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class VS,class M,class S>
	 VS &_vsmcimultassign(VS &v,const M &m)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<M>);
#else
	throw();
#endif

	template <class MS,class V,class E>
	 E _msvmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class MS,class V,class E>
	 E _msvimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class MS,class V,class E>
	 E _msvlmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class MS,class V,class E>
	 E _msvlimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class MS,class V,class E>
	 E _msvcmult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class MS,class V,class E>
	 E _msvcimult(const MS &ms,const V &v)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class E>
	 E _vmsmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class E>
	 E _vmsimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class E>
	 E _vmscmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class E>
	 E _vmscimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class E>
	 E _vmslmult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class E>
	 E _vmslimult(const V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class S>
	 V &_vmsmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class S>
	 V &_vmsimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class S>
	 V &_vmslmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class S>
	 V &_vmslimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class S>
	 V &_vmscmultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class V,class MS,class S>
	 V &_vmscimultassign(V &v,const MS &ms)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MS>);
#else
	throw();
#endif

	template <class M,class S,class E>
	 E _msdiv(const M &m,const S &c) throw();

	template <class M,class S>
	 M &_msdivassign(M &m,const S &c) throw();

	template <class MS,class S,class E>
	 E _mssdiv(const MS &ms, const S &c) throw();

	template <class MS,class S>
	 MS &_mssdivassign(MS &ms,const S &c) throw();
	
	template <class MS1,class MS2,class E>
	 E _msmsconv(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS1,class MS2,class E>
	 E _msmssect(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MS1,class MS2,class E>
	 E _msmsplus(const MS1 &m1,const MS2 &m2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class M>
	 void *_mvoid(const M &m) throw();

	template <class M>
	 bool _mnot(const M &m) throw();

	template <class MS>
	 void *_msvoid(const MS &ms) throw();

	template <class MS>
	 bool _msnot(const MS &ms) throw();

	template <class M1,class M2>
	 bool _mmeq(const M1 &m1,const M2 &m2) throw();

	template <class M1,class M2>
	 bool _mmneq(const M1 &m1,const M2 &m2) throw();

	template <class M1,class M2>
	 bool _mmless(const M1 &m1,const M2 &m2) throw();

	template <class M1,class M2>
	 bool _mmleq(const M1 &m1,const M2 &m2) throw();

	template <class M,class MS>
	 bool _mmseq(const M &m1,const MS &ms) throw();

	template <class M,class MS>
	 bool _mmsneq(const M &m1,const MS &ms) throw();

	template <class M,class MS>
	 bool _mmsless(const M &m1,const MS &ms) throw();

	template <class M,class MS>
	 bool _mmsleq(const M &m1,const MS &ms) throw();

	template <class MS,class M>
	 bool _msmless(const MS &ms,const M &m1) throw();

	template <class MS,class M>
	 bool _msmleq(const MS &ms,const M &m1) throw();

	template <class MS1,class MS2>
	 bool _msmseq(const MS1 &ms1,const MS2 &ms2) throw();

	template <class MS1,class MS2>
	 bool _msmsneq(const MS1 &ms1,const MS2 &ms2) throw();

	template <class MS1,class MS2>
	 bool _msmsless(const MS1 &ms1,const MS2 &ms2) throw();

	template <class MS1,class MS2>
	 bool _msmsleq(const MS1 &ms1,const MS2 &ms2) throw();


//------------------- matrix_subv -----------------------------


template <class V,class MV2,class S>
 V &_vmvassign(V &v,const MV2 &rv) throw();

template <class MV1,class MV2>
 MV1 &_mvmvassign(MV1 &v,const MV2 &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

template <class MV,class S>
 MV &_mvsassign(MV &v,const  S &r) throw();

template <class MV,class V>
 MV &_mvvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class MV,class V>
 V _mvabs(const MV &mv) throw();

template <class MV,class V>
 V _mvdiam(const MV &mv) throw();

template <class MV,class V>
 V _mvmid(const MV &mv) throw();

template <class MV,class V>
 V _mvinf(const MV &mv) throw();

template <class MV,class V>
 V _mvsup(const MV &mv) throw();

template <class MV,class V>
 V _mvim(const MV &mv) throw();

template <class MV,class V>
 V _mvre(const MV &mv) throw();

template <class DP,class V,class SV>
	 void _vmvaccu(DP &dp, const V & rv1, const SV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

template <class DP,class MV1,class MV2>
	 void _mvmvaccu(DP &dp, const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM);
#else
	throw();
#endif

	template <class MV1,class MV2,class S>
	 S _mvmvmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

	template <class MV1,class MV2,class S>
	 S _mvmvimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

	template <class MV1,class MV2,class S>
	 S _mvmvlmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

	template <class MV1,class MV2,class S>
	 S _mvmvlimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

	template <class MV1,class MV2,class S>
	 S _mvmvcmult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

	template <class MV1,class MV2,class S>
	 S _mvmvcimult(const MV1 & rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV1>);
#else
	throw();
#endif

	template <class V,class MV,class S>
	 S _vmvmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

	template <class V,class MV,class S>
	 S _vmvimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

	template <class V,class MV,class S>
	 S _vmvlmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

	template <class V,class MV,class S>
	 S _vmvlimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

	template <class V,class MV,class S>
	 S _vmvcmult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

	template <class V,class MV,class S>
	 S _vmvcimult(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

	template <class MV,class S,class E>
	 E _mvsmult(const MV &rv, const S &s) throw();


	template <class MV1,class MV2,class E>
	 E _mvmvconv(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV,class V,class E>
	 E _mvvconv(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV1,class MV2,class E>
	 E _mvmvsect(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV,class V,class E>
	 E _mvvsect(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV1,class MV2,class E>
	 E _mvmvplus(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV1,class MV2,class E>
	 E _mvmvminus(const MV1 &rv1, const MV2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV,class V,class E>
	 E _mvvplus(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV,class V,class E>
	 E _mvvminus(const MV &rv1, const V &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class V,class MV,class E>
	 E _vmvminus(const V &rv1, const MV &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>);
#else
	throw();
#endif

	template <class MV,class S,class E>
	 E _mvsdiv(const MV &rv, const S &s) throw();

template <class MV,class S>
 MV &_mvsmultassign(MV &v,const S &r) throw();

template <class MV, class S>
 MV &_mvsplusassign(MV &v,const S &r) throw();
	
template <class MV,class S>
 MV &_mvsminusassign(MV &v,const S &r) throw();
	
template <class MV,class S>
 MV &_mvsdivassign(MV &v,const S &r) throw();
	
template <class MV,class V>
 MV &_mvvconvassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvconvassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvsectassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvsectassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvplusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvplusassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvminusassign(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvminusassign(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

	template <class MV,class S>
	 MV &_mvsusetsup(MV &mv, const S &s) throw();

	template <class MV,class S>
	 MV &_mvsusetinf(MV &mv, const S &s) throw();

	template <class MV,class S>
	 MV &_mvssetinf(MV &mv, const S &s) throw();

	template <class MV,class S>
	 MV &_mvssetsup(MV &mv, const S &s) throw();

	template <class MV,class S>
	 MV &_mvssetre(MV &mv, const S &s) throw();

	template <class MV,class S>
	 MV &_mvssetim(MV &mv, const S &s) throw();

template <class MV,class V>
 MV &_mvvsetinf(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvsetinf(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvusetinf(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvusetinf(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvsetsup(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvsetsup(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvusetsup(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvusetsup(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvsetim(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvsetim(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif

template <class MV,class V>
 MV &_mvvsetre(MV &v,const V &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<MV>);
#else
	throw();
#endif

template <class V,class MV>
 V &_vmvsetre(V &rv,const MV &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<V>);
#else
	throw();
#endif





/*
template <class S,class S,class U>
rmatrix<U> _ms_mult(const rmatrix<S> &m,const S &s,const U &u=0) throw();

template <class S,class S,class R>
cxscvector<R> _mv_mult(const rmatrix<S> &m,const rvector &v,R) throw();

template <class S,class S, class R>
rvector _vm_mult(const rvector &v,const rmatrix<S> &m,R) throw();
*/

} // namespace cxsc 

#endif // _CXSC_MATRIX_HPP_INCLUDED
