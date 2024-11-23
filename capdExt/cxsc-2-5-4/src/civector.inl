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

/* CVS $Id: civector.inl,v 1.28 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CIVECTOR_INL_INCLUDED
#define _CXSC_CIVECTOR_INL_INCLUDED

namespace cxsc {

	INLINE civector::civector () throw():dat(NULL),l(1),u(0),size(0)
	{
	}

	INLINE civector::civector(const int &i) throw():l(1),u(i),size(i)
	{
		dat=new cinterval[i];
	}

#ifdef OLD_CXSC
	INLINE civector::civector(const class index &i) throw():l(1),u(i._int()),size(i._int())
	{
		dat=new cinterval[i._int()];
	}
#endif

	INLINE civector::civector(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_WRONG_BOUNDARIES,ERROR_CIVECTOR_NO_MORE_MEMORY):l(i1),u(i2),size(i2-i1+1)
#else
	throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_CIVECTOR_WRONG_BOUNDARIES("civector::civector(const int &i1,const int &i2)"));
#endif
		dat=new cinterval[size];
	}

	INLINE civector::civector(const civector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new cinterval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE civector::civector(const cvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new cinterval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE civector::civector(const ivector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new cinterval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE civector::civector(const rvector_slice &rs) throw():l(rs.start),u(rs.end),size(rs.end-rs.start+1)
	{
		dat=new cinterval[size];
		for(int i=0, j=l-rs.l;i<size;i++,j++)
			dat[i]=rs.dat[j];
	}

	INLINE civector::civector(const civector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE civector::civector(const cinterval &r) throw():l(1),u(1),size(1)
	{
		dat=new cinterval[1];
		*dat=r;
	}
	
	INLINE civector::civector(const cvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE civector::civector(const ivector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE civector::civector(const rvector &v) throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new cinterval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	INLINE civector::civector(const real &r) throw():l(1),u(1),size(1)
	{
		dat=new cinterval[1];
		*dat=r;
	}
	
	INLINE civector::civector(const interval &r) throw():l(1),u(1),size(1)
	{
		dat=new cinterval[1];
		*dat=r;
	}
	
	INLINE civector::civector(const complex &r) throw():l(1),u(1),size(1)
	{
		dat=new cinterval[1];
		*dat=r;
	}

	INLINE cinterval & civector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC("cinterval & civector::operator [](const int &i)const "));
#endif
		return dat[i-l];
	}

	INLINE cinterval & civector::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC("cinterval & civector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}

	INLINE cinterval & civector_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC("cinterval & civector_slice::operator [](const int &i) const"));
#endif
		return dat[i-l];
	}


	INLINE cinterval & civector_slice::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<start||i>end) cxscthrow(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC("cinterval & civector_slice::operator [](const int &i)"));
#endif
		return dat[i-l];
	}


	
	INLINE civector_slice civector::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<l||i>u) cxscthrow(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG("civector_slice civector::operator ()(const int &i)"));
#endif
		return civector_slice(*this,1,i);
	}
	
   INLINE civector_slice civector::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<l||i2>u) cxscthrow(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG("civector_slice civector::operator ()(const int &i1,const int &i2)"));
#endif
		return civector_slice(*this,i1,i2);
	}
	
	INLINE civector_slice civector_slice::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(1<start||i>end) cxscthrow(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG("civector_slice civector_slice::operator ()(const int &i)"));
#endif
		return civector_slice(*this,1,i);
	}
	
   INLINE civector_slice civector_slice::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1<start||i2>end) cxscthrow(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG("civector_slice civector_slice::operator ()(const int &i1,const int &i2)"));
#endif
		return civector_slice(*this,i1,i2);
	}
	
	INLINE cinterval::cinterval(const civector &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size>1) cxscthrow(ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ("cinterval::cinterval(const civector &rv)"));
		else if(rv.size<1) cxscthrow(ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ("cinterval::cinterval(const civector &rv)"));
#endif
		*this=rv.dat[0];
	}
	
	INLINE cinterval::cinterval(const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ,ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size>1) cxscthrow(ERROR_CIVECTOR_TYPE_CAST_OF_THICK_OBJ("cinterval::cinterval(const civector_slice &sl)"));
		else if(sl.size<1) cxscthrow(ERROR_CIVECTOR_USE_OF_UNINITIALIZED_OBJ("cinterval::cinterval(const civector_slice &sl)"));
#endif
		*this=sl.dat[sl.start-sl.l];
	}
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const interval &)
	*/
	INLINE civector _civector(const cinterval &r) throw() { return civector(r); }
	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const real &)
	*/
	INLINE civector _civector(const real &r) throw() { return civector(r); }	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const rvector_slice &rs)
	*/
	INLINE civector _civector(const rvector_slice &rs) throw() { return civector(rs); }	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const rvector &v)
	*/
	INLINE civector _civector(const rvector &rs) throw() { return civector(rs); }
		
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const complex &)
	*/
	INLINE civector _civector(const complex &r) throw() { return civector(r); }	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const cvector_slice &rs)
	*/
	INLINE civector _civector(const cvector_slice &rs) throw() { return civector(rs); }	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const cvector &v)
	*/
	INLINE civector _civector(const cvector &rs) throw() { return civector(rs); }
		
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const interval &)
	*/
	INLINE civector _civector(const interval &r) throw() { return civector(r); }	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const ivector_slice &rs)
	*/
	INLINE civector _civector(const ivector_slice &rs) throw() { return civector(rs); }	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::civector::civector(const ivector &v)
	*/
	INLINE civector _civector(const ivector &rs) throw() { return civector(rs); }
	
	INLINE civector &civector::operator =(const civector &rv) throw() { return _vvassign<civector,civector,cinterval>(*this,rv); }
	INLINE civector &civector::operator =(const cinterval &r) throw() { return _vsassign<civector,cinterval>(*this,r); }
	INLINE civector::operator void*() throw() { return _vvoid(*this); }
	INLINE civector_slice & civector_slice::operator =(const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }
	INLINE civector_slice & civector_slice::operator =(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE civector_slice & civector_slice::operator =(const cinterval &r) throw() { return _vssassign<civector_slice,cinterval>(*this,r); }
	INLINE civector_slice & civector_slice::operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>,ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,civector(m)); }
	INLINE civector_slice::operator void*() throw() { return _vsvoid(*this); }
	
//=======================================================================
//======================== Vector Functions =============================


	INLINE civector &SetInf(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetinf(iv,rv); }
	INLINE civector_slice &SetInf(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetinf(iv,rv); }
	INLINE civector &SetInf(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetinf(iv,rv); }
	INLINE civector_slice &SetInf(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetinf(iv,rv); }
	INLINE civector &UncheckedSetInf(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvusetinf(iv,rv); }
	INLINE civector_slice &UncheckedSetInf(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvusetinf(iv,rv); }
	INLINE civector &UncheckedSetInf(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsusetinf(iv,rv); }
	INLINE civector_slice &UncheckedSetInf(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsusetinf(iv,rv); }

	INLINE civector &SetSup(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetsup(iv,rv); }
	INLINE civector_slice &SetSup(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetsup(iv,rv); }
	INLINE civector &SetSup(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetsup(iv,rv); }
	INLINE civector_slice &SetSup(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetsup(iv,rv); }
	INLINE civector &UncheckedSetSup(civector &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvusetsup(iv,rv); }
	INLINE civector_slice &UncheckedSetSup(civector_slice &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvusetsup(iv,rv); }
	INLINE civector &UncheckedSetSup(civector &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsusetsup(iv,rv); }
	INLINE civector_slice &UncheckedSetSup(civector_slice &iv,const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsusetsup(iv,rv); }

	INLINE civector &SetRe(civector &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetre(iv,rv); }
	INLINE civector_slice &SetRe(civector_slice &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetre(iv,rv); }
	INLINE civector &SetRe(civector &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetre(iv,rv); }
	INLINE civector_slice &SetRe(civector_slice &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetre(iv,rv); }

	INLINE civector &SetIm(civector &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvsetim(iv,rv); }
	INLINE civector_slice &SetIm(civector_slice &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvsetim(iv,rv); }
	INLINE civector &SetIm(civector &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vvssetim(iv,rv); }
	INLINE civector_slice &SetIm(civector_slice &iv,const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsvssetim(iv,rv); }

	INLINE civector &SetSup(civector &iv,const complex &r) throw() { return _vssetsup(iv,r); }
	INLINE civector &SetInf(civector &iv,const complex &r) throw() { return _vssetinf(iv,r); }
	INLINE civector &UncheckedSetSup(civector &iv,const complex &r) throw() { return _vsusetsup(iv,r); }
	INLINE civector &SetUncheckedInf(civector &iv,const complex &r) throw() { return _vsusetinf(iv,r); }
	INLINE civector &SetRe(civector &iv,const interval &r) throw() { return _vssetre(iv,r); }
	INLINE civector &SetIm(civector &iv,const interval &r) throw() { return _vssetim(iv,r); }

	INLINE civector_slice &SetSup(civector_slice &iv,const complex &r) throw() { return _vsssetsup(iv,r); }
	INLINE civector_slice &SetInf(civector_slice &iv,const complex &r) throw() { return _vsssetinf(iv,r); }
	INLINE civector_slice &UncheckedSetSup(civector_slice &iv,const complex &r) throw() { return _vssusetsup(iv,r); }
	INLINE civector_slice &SetUncheckedInf(civector_slice &iv,const complex &r) throw() { return _vssusetinf(iv,r); }
	INLINE civector_slice &SetRe(civector_slice &iv,const interval &r) throw() { return _vsssetre(iv,r); }
	INLINE civector_slice &SetIm(civector_slice &iv,const interval &r) throw() { return _vsssetim(iv,r); }

	INLINE ivector Re(const civector &v) throw()
	{
		ivector erg(v.l,v.u);
		
		for(int i=0;i<v.size;i++)
			erg[i+v.l]=Re(v.dat[i]);

		return erg;
	}

	INLINE ivector Re(const civector_slice &sl) throw()
	{
		ivector erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg[i+sl.start]=Re(sl.dat[j]);

		return erg;
	}

	INLINE rvector InfRe(const civector &v) throw()
	{
		rvector erg(v.l,v.u);
		
		for(int i=0;i<v.size;i++)
			erg.dat[i]=InfRe(v.dat[i]);

		return erg;
	}

	INLINE rvector InfRe(const civector_slice &sl) throw()
	{
		rvector erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=InfRe(sl.dat[j]);

		return erg;
	}

	INLINE ivector Im(const civector &v) throw()
	{
		ivector erg(v.l,v.u);
		
		for(int i=0;i<v.size;i++)
			erg[i+v.l]=Im(v.dat[i]);

		return erg;
	}

	INLINE ivector Im(const civector_slice &sl) throw()
	{
		ivector erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg[i+sl.start]=Im(sl.dat[j]);

		return erg;
	}

	INLINE rvector InfIm(const civector &v) throw()
	{
		rvector erg(v.l,v.u);
		
		for(int i=0;i<v.size;i++)
			erg.dat[i]=InfIm(v.dat[i]);

		return erg;
	}

	INLINE rvector InfIm(const civector_slice &sl) throw()
	{
		rvector erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=InfIm(sl.dat[j]);

		return erg;
	}


	INLINE rvector SupIm(const civector &v) throw()
	{
		rvector erg(v.l,v.u);
		
		for(int i=0;i<v.size;i++)
			erg.dat[i]=SupIm(v.dat[i]);

		return erg;
	}

	INLINE rvector SupIm(const civector_slice &sl) throw()
	{
		rvector erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=SupIm(sl.dat[j]);

		return erg;
	}



	INLINE rvector SupRe(const civector &v) throw()
	{
		rvector erg(v.l,v.u);
		
		for(int i=0;i<v.size;i++)
			erg.dat[i]=SupRe(v.dat[i]);

		return erg;
	}

	INLINE rvector SupRe(const civector_slice &sl) throw()
	{
		rvector erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=SupRe(sl.dat[j]);

		return erg;
	}

	INLINE void Resize(civector &rv) throw() { _vresize(rv); } 
	INLINE void Resize(civector &rv, const int &len)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<civector>)
#else
	throw()
#endif
	{ _vresize<class civector,class cinterval>(rv,len); }
	INLINE void Resize(civector &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__WRONG_BOUNDARIES<civector>)
#else
	throw()
#endif
	{ _vresize<class civector,class cinterval>(rv,lb,ub); }
	
	INLINE civector conj(const civector &rv) throw() { return _vconj<civector>(rv); }
	INLINE civector conj(const civector_slice &sl) throw() { return _vsconj<civector_slice,civector>(sl); }
	
	INLINE ivector abs(const civector &rv) throw() { return _vabs<civector,ivector>(rv); }
	INLINE ivector abs(const civector_slice &sl) throw() { return _vsabs<civector_slice,ivector>(sl); }
	INLINE cvector diam(const civector &v) throw() { return _vdiam<civector,cvector>(v); }
	INLINE cvector diam(const civector_slice &v) throw() { return _vsdiam<civector_slice,cvector>(v); }
	INLINE cvector mid(const civector &v) throw() { return _vmid<civector,cvector>(v); }
	INLINE cvector mid(const civector_slice &v) throw() { return _vsmid<civector_slice,cvector>(v); }
	INLINE cvector Inf(const civector &v) throw() { return _vinf<civector,cvector>(v); }
	INLINE cvector Inf(const civector_slice &v) throw() { return _vsinf<civector_slice,cvector>(v); }
	INLINE cvector Sup(const civector &v) throw() { return _vsup<civector,cvector>(v); }
	INLINE cvector Sup(const civector_slice &v) throw() { return _vssup<civector_slice,cvector>(v); }
	INLINE bool operator !(const civector &rv) throw() { return _vnot(rv); }
	INLINE bool operator !(const civector_slice &sl) throw() { return _vsnot(sl); }

//======================= Vector / Scalar ===============================

//----------------------------- cinterval ---------------------------

	INLINE civector operator *(const civector &rv, const cinterval &s) throw() { return _vsmult<civector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const civector_slice &sl, const cinterval &s) throw() { return _vssmult<civector_slice,cinterval,civector>(sl,s); }
	INLINE civector operator *(const cinterval &s, const civector &rv) throw() { return _vsmult<civector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cinterval &s, const civector_slice &sl) throw() { return _vssmult<civector_slice,cinterval,civector>(sl,s); }
	INLINE civector &operator *=(civector &rv,const cinterval &r) throw() { return _vsmultassign(rv,r); }
	INLINE civector_slice &civector_slice::operator *=(const cinterval &r) throw() { return _vssmultassign(*this,r); }

	INLINE civector operator /(const civector &rv, const cinterval &s) throw() { return _vsdiv<civector,cinterval,civector>(rv,s); }
	INLINE civector operator /(const civector_slice &sl, const cinterval &s) throw() { return _vssdiv<civector_slice,cinterval,civector>(sl,s); }
	INLINE civector &operator /=(civector &rv,const cinterval &r) throw() { return _vsdivassign(rv,r); }
	INLINE civector_slice &civector_slice::operator /=(const cinterval &r) throw() { return _vssdivassign(*this,r); }

//---------------------------- Real --------------------------------------

	INLINE civector operator *(const civector &rv, const real &s) throw() { return _vsmult<civector,real,civector>(rv,s); }
	INLINE civector operator *(const civector_slice &sl, const real &s) throw() { return _vssmult<civector_slice,real,civector>(sl,s); }
	INLINE civector operator *(const real &s, const civector &rv) throw() { return _vsmult<civector,real,civector>(rv,s); }
	INLINE civector operator *(const real &s, const civector_slice &sl) throw() { return _vssmult<civector_slice,real,civector>(sl,s); }
	INLINE civector &operator *=(civector &rv,const real &r) throw() { return _vsmultassign(rv,r); }
	INLINE civector_slice &civector_slice::operator *=(const real &r) throw() { return _vssmultassign(*this,r); }

	INLINE civector operator /(const civector &rv, const real &s) throw() { return _vsdiv<civector,real,civector>(rv,s); }
	INLINE civector operator /(const civector_slice &sl, const real &s) throw() { return _vssdiv<civector_slice,real,civector>(sl,s); }
	INLINE civector &operator /=(civector &rv,const real &r) throw() { return _vsdivassign(rv,r); }
	INLINE civector_slice &civector_slice::operator /=(const real &r) throw() { return _vssdivassign(*this,r); }

	INLINE civector operator *(const rvector &rv, const cinterval &s) throw() { return _vsmult<rvector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const rvector_slice &sl, const cinterval &s) throw() { return _vssmult<rvector_slice,cinterval,civector>(sl,s); }
	INLINE civector operator *(const cinterval &s, const rvector &rv) throw() { return _vsmult<rvector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cinterval &s, const rvector_slice &sl) throw() { return _vssmult<rvector_slice,cinterval,civector>(sl,s); }

	INLINE civector operator /(const rvector &rv, const cinterval &s) throw() { return _vsdiv<rvector,cinterval,civector>(rv,s); }
	INLINE civector operator /(const rvector_slice &sl, const cinterval &s) throw() { return _vssdiv<rvector_slice,cinterval,civector>(sl,s); }

//---------------------------- complex --------------------------------------

	INLINE civector operator *(const civector &rv, const complex &s) throw() { return _vsmult<civector,complex,civector>(rv,s); }
	INLINE civector operator *(const civector_slice &sl, const complex &s) throw() { return _vssmult<civector_slice,complex,civector>(sl,s); }
	INLINE civector operator *(const complex &s, const civector &rv) throw() { return _vsmult<civector,complex,civector>(rv,s); }
	INLINE civector operator *(const complex &s, const civector_slice &sl) throw() { return _vssmult<civector_slice,complex,civector>(sl,s); }
	INLINE civector &operator *=(civector &rv,const complex &r) throw() { return _vsmultassign(rv,r); }
	INLINE civector_slice &civector_slice::operator *=(const complex &r) throw() { return _vssmultassign(*this,r); }

	INLINE civector operator /(const civector &rv, const complex &s) throw() { return _vsdiv<civector,complex,civector>(rv,s); }
	INLINE civector operator /(const civector_slice &sl, const complex &s) throw() { return _vssdiv<civector_slice,complex,civector>(sl,s); }
	INLINE civector &operator /=(civector &rv,const complex &r) throw() { return _vsdivassign(rv,r); }
	INLINE civector_slice &civector_slice::operator /=(const complex &r) throw() { return _vssdivassign(*this,r); }

	INLINE civector operator *(const cvector &rv, const cinterval &s) throw() { return _vsmult<cvector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cvector_slice &sl, const cinterval &s) throw() { return _vssmult<cvector_slice,cinterval,civector>(sl,s); }
	INLINE civector operator *(const cinterval &s, const cvector &rv) throw() { return _vsmult<cvector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cinterval &s, const cvector_slice &sl) throw() { return _vssmult<cvector_slice,cinterval,civector>(sl,s); }

	INLINE civector operator /(const cvector &rv, const cinterval &s) throw() { return _vsdiv<cvector,cinterval,civector>(rv,s); }
	INLINE civector operator /(const cvector_slice &sl, const cinterval &s) throw() { return _vssdiv<cvector_slice,cinterval,civector>(sl,s); }

//---------------------------- interval --------------------------------------

	INLINE civector operator *(const civector &rv, const interval &s) throw() { return _vsmult<civector,interval,civector>(rv,s); }
	INLINE civector operator *(const civector_slice &sl, const interval &s) throw() { return _vssmult<civector_slice,interval,civector>(sl,s); }
	INLINE civector operator *(const interval &s, const civector &rv) throw() { return _vsmult<civector,interval,civector>(rv,s); }
	INLINE civector operator *(const interval &s, const civector_slice &sl) throw() { return _vssmult<civector_slice,interval,civector>(sl,s); }
	INLINE civector &operator *=(civector &rv,const interval &r) throw() { return _vsmultassign(rv,r); }
	INLINE civector_slice &civector_slice::operator *=(const interval &r) throw() { return _vssmultassign(*this,r); }

	INLINE civector operator /(const civector &rv, const interval &s) throw() { return _vsdiv<civector,interval,civector>(rv,s); }
	INLINE civector operator /(const civector_slice &sl, const interval &s) throw() { return _vssdiv<civector_slice,interval,civector>(sl,s); }
	INLINE civector &operator /=(civector &rv,const interval &r) throw() { return _vsdivassign(rv,r); }
	INLINE civector_slice &civector_slice::operator /=(const interval &r) throw() { return _vssdivassign(*this,r); }

	INLINE civector operator *(const ivector &rv, const cinterval &s) throw() { return _vsmult<ivector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const ivector_slice &sl, const cinterval &s) throw() { return _vssmult<ivector_slice,cinterval,civector>(sl,s); }
	INLINE civector operator *(const cinterval &s, const ivector &rv) throw() { return _vsmult<ivector,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cinterval &s, const ivector_slice &sl) throw() { return _vssmult<ivector_slice,cinterval,civector>(sl,s); }

	INLINE civector operator /(const ivector &rv, const cinterval &s) throw() { return _vsdiv<ivector,cinterval,civector>(rv,s); }
	INLINE civector operator /(const ivector_slice &sl, const cinterval &s) throw() { return _vssdiv<ivector_slice,cinterval,civector>(sl,s); }

//======================= Vector / Vector ===============================


	INLINE std::ostream &operator <<(std::ostream &s, const civector &rv) throw() { return _vout(s,rv); }
	INLINE std::ostream &operator <<(std::ostream &o, const civector_slice &sl) throw() { return _vsout(o,sl); }
	INLINE std::istream &operator >>(std::istream &s, civector &rv) throw() { return _vin(s,rv); }
	INLINE std::istream &operator >>(std::istream &s, civector_slice &rv) throw() { return _vsin(s,rv); }
	
//----------------------- cinterval / cinterval ---------------------------
	INLINE civector & civector::operator =(const civector_slice &sl) throw() { return _vvsassign<civector,civector_slice,cinterval>(*this,sl); }


	INLINE cinterval operator *(const civector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<civector,civector,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<civector_slice,civector_slice,cinterval>(sl1,sl2); }
	
	INLINE const civector &operator +(const civector &rv) throw() { return rv; }
	INLINE civector operator +(const civector_slice &sl) throw() { return sl; }
	INLINE civector operator +(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<civector,civector,civector>(rv1,rv2); }
	INLINE civector operator +(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<civector_slice,civector_slice,civector>(sl1,sl2); }
	INLINE civector & operator +=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE civector &operator +=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator +=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator +=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE civector operator -(const civector &rv) throw() { return _vminus(rv); }
	INLINE civector operator -(const civector_slice &sl) throw() { return _vsminus<civector_slice,civector>(sl); }
	INLINE civector operator -(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<civector,civector,civector>(rv1,rv2); }
	INLINE civector operator -(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator -(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<civector_slice,civector,civector>(sl,rv); }
	INLINE civector operator -(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<civector_slice,civector_slice,civector>(sl1,sl2); }
	INLINE civector & operator -=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE civector &operator -=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator -=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator -=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE civector operator |(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<civector,civector,civector>(rv1,rv2); }
	INLINE civector operator |(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<civector_slice,civector_slice,civector>(sl1,sl2); }
	INLINE civector & operator |=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE civector &operator |=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator |=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator |=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE civector operator &(const civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<civector,civector,civector>(rv1,rv2); }
	INLINE civector operator &(const civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<civector_slice,civector_slice,civector>(sl1,sl2); }
	INLINE civector & operator &=(civector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE civector &operator &=(civector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator &=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator &=(const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssectassign(*this,sl2); }

	INLINE bool operator ==(const civector &rv1, const civector &rv2) throw() { return _vveq(rv1,rv2); }
	INLINE bool operator ==(const civector_slice &sl1, const civector_slice &sl2) throw() { return _vsvseq(sl1,sl2); }
	INLINE bool operator ==(const civector_slice &sl, const civector &rv) throw() { return _vsveq(sl,rv); }
	INLINE bool operator ==(const civector &rv, const civector_slice &sl) throw() { return _vsveq(sl,rv); }
	INLINE bool operator !=(const civector &rv1, const civector &rv2) throw() { return _vvneq(rv1,rv2); }
	INLINE bool operator !=(const civector_slice &sl1, const civector_slice &sl2) throw() { return _vsvsneq(sl1,sl2); }
	INLINE bool operator !=(const civector_slice &sl, const civector &rv) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator !=(const civector &rv, const civector_slice &sl) throw() { return _vsvneq(sl,rv); }
	INLINE bool operator <(const civector &rv1, const civector &rv2) throw() { return _vvless(rv1,rv2); }
	INLINE bool operator <(const civector_slice &sl1, const civector_slice &sl2) throw() { return _vsvsless(sl1,sl2); }
	INLINE bool operator < (const civector_slice &sl, const civector &rv) throw() { return _vsvless(sl,rv); }
	INLINE bool operator < (const civector &rv, const civector_slice &sl) throw() { return _vvsless(rv,sl); }
	INLINE bool operator <=(const civector &rv1, const civector &rv2) throw() { return _vvleq(rv1,rv2); }
	INLINE bool operator <=(const civector_slice &sl1, const civector_slice &sl2) throw() { return _vsvsleq(sl1,sl2); }
	INLINE bool operator <=(const civector_slice &sl, const civector &rv) throw() { return _vsvleq(sl,rv); }
	INLINE bool operator <=(const civector &rv, const civector_slice &sl) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >(const civector &rv1, const civector &rv2) throw() { return _vvless(rv2,rv1); }
	INLINE bool operator >(const civector_slice &sl1, const civector_slice &sl2) throw() { return _vsvsless(sl2,sl1); }
	INLINE bool operator >(const civector_slice &sl, const civector &rv) throw() { return _vvsless(rv,sl); }
	INLINE bool operator >(const civector &rv, const civector_slice &sl) throw() { return _vsvless(sl,rv); }
	INLINE bool operator >=(const civector &rv1, const civector &rv2) throw() { return _vvleq(rv2,rv1); }
	INLINE bool operator >=(const civector_slice &sl1, const civector_slice &sl2) throw() { return _vsvsleq(sl2,sl1); }
	INLINE bool operator >=(const civector_slice &sl, const civector &rv) throw() { return _vvsleq(rv,sl); }
	INLINE bool operator >=(const civector &rv, const civector_slice &sl) throw() { return _vsvleq(sl,rv); }

//-------------------------------- cinterval / Real --------------------------------

	INLINE civector &civector::operator =(const rvector &rv) throw() { return _vvassign<civector,rvector,cinterval>(*this,rv); }
	INLINE civector &civector::operator =(const real &r) throw() { return _vsassign<civector,real>(*this,r); }
	INLINE civector & civector::operator =(const rvector_slice &sl) throw() { return _vvsassign<civector,rvector_slice,cinterval>(*this,sl); }
	INLINE civector_slice &civector_slice::operator =(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator =(const real &r) throw() { return _vssassign<civector_slice,real>(*this,r); }
	INLINE civector_slice & civector_slice::operator =(const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }


	INLINE cinterval operator *(const rvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<rvector,civector,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<rvector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,rvector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const rvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<rvector_slice,civector_slice,cinterval>(sl1,sl2); }
	
	INLINE cinterval operator *(const civector & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<rvector,civector,cinterval>(rv2,rv1); }
	INLINE cinterval operator *(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,rvector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<rvector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector_slice & sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<rvector_slice,civector_slice,cinterval>(sl2,sl1); }
	
	INLINE civector operator +(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,civector,civector>(rv1,rv2); }
	INLINE civector operator +(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator +(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator +(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<rvector,civector,civector>(rv2,rv1); }
	INLINE civector operator +(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<rvector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator +=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE civector &operator +=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator +=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE civector operator -(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<rvector,civector,civector>(rv1,rv2); }
	INLINE civector operator -(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator -(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<rvector_slice,civector,civector>(sl,rv); }
	INLINE civector operator -(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<rvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator -(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<civector,rvector,civector>(rv1,rv2); }
	INLINE civector operator -(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator -(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<civector_slice,rvector,civector>(sl,rv); }
	INLINE civector operator -(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<civector_slice,rvector_slice,civector>(sl1,sl2); }

	INLINE civector & operator -=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE civector &operator -=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator -=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE civector operator |(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,civector,civector>(rv1,rv2); }
	INLINE civector operator |(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator |(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,civector,civector>(rv2,rv1); }
	INLINE civector operator |(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator |=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE civector &operator |=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator |=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE civector operator &(const rvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<rvector,civector,civector>(rv1,rv2); }
	INLINE civector operator &(const rvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const rvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator &(const rvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<rvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator &(const civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<rvector,civector,civector>(rv2,rv1); }
	INLINE civector operator &(const civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<rvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<rvector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator &=(civector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE civector &operator &=(civector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator &=(const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssectassign(*this,sl2); }

//-------------------------------- cinterval / complex --------------------------------

	INLINE civector &civector::operator =(const cvector &rv) throw() { return _vvassign<civector,cvector,cinterval>(*this,rv); }
	INLINE civector &civector::operator =(const complex &r) throw() { return _vsassign<civector,complex>(*this,r); }
	INLINE civector & civector::operator =(const cvector_slice &sl) throw() { return _vvsassign<civector,cvector_slice,cinterval>(*this,sl); }
	INLINE civector_slice &civector_slice::operator =(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator =(const complex &r) throw() { return _vssassign<civector_slice,complex>(*this,r); }
	INLINE civector_slice & civector_slice::operator =(const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }


	INLINE cinterval operator *(const cvector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<cvector,civector,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<cvector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,cvector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const cvector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<cvector_slice,civector_slice,cinterval>(sl1,sl2); }
	
	INLINE cinterval operator *(const civector & rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<cvector,civector,cinterval>(rv2,rv1); }
	INLINE cinterval operator *(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,cvector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<cvector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector_slice & sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<cvector_slice,civector_slice,cinterval>(sl2,sl1); }
	
	INLINE civector operator +(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<cvector,civector,civector>(rv1,rv2); }
	INLINE civector operator +(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator +(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<cvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator +(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<cvector,civector,civector>(rv2,rv1); }
	INLINE civector operator +(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<cvector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator +=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE civector &operator +=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator +=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE civector operator -(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<cvector,civector,civector>(rv1,rv2); }
	INLINE civector operator -(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator -(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<cvector_slice,civector,civector>(sl,rv); }
	INLINE civector operator -(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<cvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator -(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<civector,cvector,civector>(rv1,rv2); }
	INLINE civector operator -(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator -(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<civector_slice,cvector,civector>(sl,rv); }
	INLINE civector operator -(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<civector_slice,cvector_slice,civector>(sl1,sl2); }

	INLINE civector & operator -=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE civector &operator -=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator -=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE civector operator |(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<cvector,civector,civector>(rv1,rv2); }
	INLINE civector operator |(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<cvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator |(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<cvector,civector,civector>(rv2,rv1); }
	INLINE civector operator |(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<cvector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator |=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE civector &operator |=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator |=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator |=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE civector operator &(const cvector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<cvector,civector,civector>(rv1,rv2); }
	INLINE civector operator &(const cvector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const cvector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator &(const cvector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<cvector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator &(const civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<cvector,civector,civector>(rv2,rv1); }
	INLINE civector operator &(const civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<cvector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<cvector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator &=(civector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE civector &operator &=(civector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator &=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator &=(const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssectassign(*this,sl2); }

//-------------------------------- cinterval / interval --------------------------------

	INLINE civector &civector::operator =(const ivector &rv) throw() { return _vvassign<civector,ivector,cinterval>(*this,rv); }
	INLINE civector &civector::operator =(const interval &r) throw() { return _vsassign<civector,interval>(*this,r); }
	INLINE civector & civector::operator =(const ivector_slice &sl) throw() { return _vvsassign<civector,ivector_slice,cinterval>(*this,sl); }
	INLINE civector_slice &civector_slice::operator =(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator =(const interval &r) throw() { return _vssassign<civector_slice,interval>(*this,r); }
	INLINE civector_slice & civector_slice::operator =(const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsassign(*this,sl); }


	INLINE cinterval operator *(const ivector & rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<ivector,civector,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<ivector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,ivector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const ivector_slice & sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<ivector_slice,civector_slice,cinterval>(sl1,sl2); }
	
	INLINE cinterval operator *(const civector & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvcimult<ivector,civector,cinterval>(rv2,rv1); }
	INLINE cinterval operator *(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<civector_slice,ivector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvcimult<ivector_slice,civector,cinterval>(sl,rv); }
	INLINE cinterval operator *(const civector_slice & sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvscimult<ivector_slice,civector_slice,cinterval>(sl2,sl1); }
	
	INLINE civector operator +(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<ivector,civector,civector>(rv1,rv2); }
	INLINE civector operator +(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator +(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<ivector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator +(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<ivector,civector,civector>(rv2,rv1); }
	INLINE civector operator +(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator +(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<ivector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator +=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplusassign(rv1,rv2); }
	INLINE civector &operator +=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvplusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator +=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplusassign(*this,sl2); }

	INLINE civector operator -(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<ivector,civector,civector>(rv1,rv2); }
	INLINE civector operator -(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator -(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<ivector_slice,civector,civector>(sl,rv); }
	INLINE civector operator -(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<ivector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator -(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<civector,ivector,civector>(rv1,rv2); }
	INLINE civector operator -(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator -(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<civector_slice,ivector,civector>(sl,rv); }
	INLINE civector operator -(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<civector_slice,ivector_slice,civector>(sl1,sl2); }

	INLINE civector & operator -=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminusassign(rv1,rv2); }
	INLINE civector &operator -=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminusassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminusassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator -=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminusassign(*this,sl2); }

	INLINE civector operator |(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<ivector,civector,civector>(rv1,rv2); }
	INLINE civector operator |(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator |(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<ivector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator |(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<ivector,civector,civector>(rv2,rv1); }
	INLINE civector operator |(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator |(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<ivector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator |=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconvassign(rv1,rv2); }
	INLINE civector &operator |=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconvassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvconvassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator |=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconvassign(*this,sl2); }

	INLINE civector operator &(const ivector &rv1, const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<ivector,civector,civector>(rv1,rv2); }
	INLINE civector operator &(const ivector &rv, const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const ivector_slice &sl, const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator &(const ivector_slice &sl1, const civector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<ivector_slice,civector_slice,civector>(sl1,sl2); }

	INLINE civector operator &(const civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<ivector,civector,civector>(rv2,rv1); }
	INLINE civector operator &(const civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<civector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,civector_slice,civector>(rv,sl); }
	INLINE civector operator &(const civector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<ivector_slice,civector_slice,civector>(sl2,sl1); }

	INLINE civector & operator &=(civector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsectassign(rv1,rv2); }
	INLINE civector &operator &=(civector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssectassign(rv,sl); }
	INLINE civector_slice &civector_slice::operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsectassign(*this,rv); }
	INLINE civector_slice &civector_slice::operator &=(const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssectassign(*this,sl2); }

//------------- real x complex ------------------------
	INLINE civector operator |(const rvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,cvector,civector>(rv1,rv2); }
	INLINE civector operator |(const cvector &rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<rvector,cvector,civector>(rv2,rv1); }
	INLINE civector operator |(const cvector &rv, const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const rvector_slice &sl,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,rvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl, const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const rvector &rv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<rvector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,cvector_slice,civector>(sl2,sl1); }
	INLINE civector operator |(const rvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<rvector_slice,cvector_slice,civector>(sl1,sl2); }

//------------- complex x complex ------------------------
	INLINE civector operator |(const cvector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<cvector,cvector,civector>(rv2,rv1); }
	INLINE civector operator |(const cvector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<cvector_slice,cvector_slice,civector>(sl1,sl2); }

//-------------------------------- interval / complex --------------------------------

	INLINE civector operator +(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<cvector,ivector,civector>(rv1,rv2); }
	INLINE civector operator +(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator +(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator +(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<cvector_slice,ivector_slice,civector>(sl1,sl2); }

	INLINE civector operator +(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvplus<cvector,ivector,civector>(rv2,rv1); }
	INLINE civector operator +(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator +(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsplus<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator +(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsplus<cvector_slice,ivector_slice,civector>(sl2,sl1); }

	INLINE civector operator -(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<cvector,ivector,civector>(rv1,rv2); }
	INLINE civector operator -(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator -(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<cvector_slice,ivector,civector>(sl,rv); }
	INLINE civector operator -(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<cvector_slice,ivector_slice,civector>(sl1,sl2); }

	INLINE civector operator -(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvminus<ivector,cvector,civector>(rv1,rv2); }
	INLINE civector operator -(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsminus<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator -(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvminus<ivector_slice,cvector,civector>(sl,rv); }
	INLINE civector operator -(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsminus<ivector_slice,cvector_slice,civector>(sl1,sl2); }

	INLINE civector operator |(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<cvector,ivector,civector>(rv1,rv2); }
	INLINE civector operator |(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<cvector_slice,ivector_slice,civector>(sl1,sl2); }

	INLINE civector operator |(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvconv<cvector,ivector,civector>(rv2,rv1); }
	INLINE civector operator |(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator |(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsconv<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator |(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvsconv<cvector_slice,ivector_slice,civector>(sl2,sl1); }

	INLINE civector operator &(const cvector &rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<cvector,ivector,civector>(rv1,rv2); }
	INLINE civector operator &(const cvector &rv, const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator &(const cvector_slice &sl, const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator &(const cvector_slice &sl1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<cvector_slice,ivector_slice,civector>(sl1,sl2); }

	INLINE civector operator &(const ivector &rv1, const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvsect<cvector,ivector,civector>(rv2,rv1); }
	INLINE civector operator &(const ivector &rv, const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<ivector,cvector_slice,civector>(rv,sl); }
	INLINE civector operator &(const ivector_slice &sl, const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vvssect<cvector,ivector_slice,civector>(rv,sl); }
	INLINE civector operator &(const ivector_slice &sl1, const cvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>)
#else
	throw()
#endif
	{ return _vsvssect<cvector_slice,ivector_slice,civector>(sl2,sl1); }


        //! Computes permutation of vector according to permutation vector, C=Px
        INLINE civector civector::operator()(const intvector& p) {
          civector x(*this);
          for(int i=0 ; i<VecLen(x) ; i++)
              x[i+Lb(x)] = (*this)[p[i+Lb(p)]+Lb(*this)];
          return x;

        }

	INLINE bool in(const civector& v1, const civector& v2) {
          int n = VecLen(v1);
          bool ret = true;
          for(int i=0 ; i<n && ret ; i++) {
             ret = in(v1[Lb(v1)+i], v2[Lb(v2)+i]);
          }
          return ret;
        }


} // namespace cxsc

#endif

