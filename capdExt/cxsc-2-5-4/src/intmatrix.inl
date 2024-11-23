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

/* CVS $Id: intmatrix.inl,v 1.19 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_INTMATRIX_INL_INCLUDED
#define _CXSC_INTMATRIX_INL_INCLUDED

namespace cxsc {

INLINE intmatrix::intmatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE intmatrix::intmatrix(const int &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new int[1];
	*dat=r;
}

INLINE intmatrix::intmatrix(const intmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new int[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE intmatrix::intmatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_INTMATRIX_WRONG_BOUNDARIES("intmatrix::intmatrix(const int &m, const int &n)"));
#endif
	dat=new int[m*n];
}

INLINE intmatrix::intmatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_INTMATRIX_WRONG_BOUNDARIES("intmatrix::intmatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new int[xsize*ysize];
}

INLINE intvector::intvector(const intmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new int[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE intmatrix::intmatrix(const intvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new int[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE intmatrix::intmatrix(const intvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new int[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}
	
	INLINE int &intmatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_INTVECTOR_ELEMENT_NOT_IN_VEC("int &intmatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE intmatrix::intmatrix(const intmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new int[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE intmatrix_subv Row(intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE intmatrix_subv Col(intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE intmatrix_subv Row(const intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE intmatrix_subv Col(const intmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
		
	INLINE intmatrix_subv intmatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT("intmatrix_subv intmatrix::operator [](const int &i)"));
#endif
		return intmatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE intmatrix_subv intmatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT("intmatrix_subv intmatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return intmatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}
	
	INLINE intmatrix_slice intmatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG("intmatrix_slice intmatrix::operator ()(const int &m, const int &n)"));
#endif
		return intmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE intmatrix_slice intmatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG("intmatrix_slice intmatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return intmatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE intmatrix_subv intmatrix_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT("intmatrix_subv intmatrix_slice::operator [](const int &i)"));
#endif
		return intmatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE intmatrix_subv intmatrix_slice::operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_INTMATRIX_ROW_OR_COL_NOT_IN_MAT("intmatrix_subv intmatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return intmatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}
	
	INLINE intmatrix_slice intmatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG("intmatrix_slice intmatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return intmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE intmatrix_slice intmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_INTMATRIX_SUB_ARRAY_TOO_BIG("intmatrix_slice intmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return intmatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE intmatrix_subv intmatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG("intmatrix_subv intmatrix_subv::operator ()(const int &i)"));
#endif
	return intmatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE intmatrix_subv intmatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_INTVECTOR_SUB_ARRAY_TOO_BIG("intmatrix_subv intmatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return intmatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}

	INLINE intmatrix_subv &intmatrix_subv::operator =(const intmatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE intmatrix_subv &intmatrix_subv::operator =(const int &r) throw() { return _mvsassign(*this,r); }
	INLINE intmatrix_subv &intmatrix_subv::operator =(const intvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE intmatrix_subv &intmatrix_subv::operator =(const intvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,intvector(v)); }
	INLINE intmatrix &intmatrix::operator =(const int &r) throw() { return _msassign(*this,r); }
	INLINE intmatrix &intmatrix::operator =(const intmatrix &m) throw() { return _mmassign<intmatrix,intmatrix,int>(*this,m,0); }
	INLINE intmatrix &intmatrix::operator =(const intvector &v) throw() { return _mvassign<intmatrix,intvector,int>(*this,v); }
	INLINE intmatrix &intmatrix::operator =(const intvector_slice &v) throw() { return _mvassign<intmatrix,intvector,int>(*this,intvector(v)); }
	INLINE intmatrix::operator void*() throw() { return _mvoid(*this); }
	INLINE intmatrix_slice &intmatrix_slice::operator =(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE intmatrix_slice &intmatrix_slice::operator =(const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE intmatrix_slice &intmatrix_slice::operator =(const int &r) throw() { return _mssassign(*this,r); }
	INLINE intmatrix_slice &intmatrix_slice::operator =(const intvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,intmatrix(v)); }
	INLINE intmatrix_slice &intmatrix_slice::operator =(const intvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,intmatrix(intvector(v))); }
	INLINE intmatrix_slice &intmatrix_slice::operator =(const intmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,intmatrix(intvector(v))); }
	INLINE intmatrix_slice &intmatrix_slice::operator +=(const intmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE intmatrix_slice &intmatrix_slice::operator +=(const intmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE intmatrix_slice &intmatrix_slice::operator -=(const intmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE intmatrix_slice &intmatrix_slice::operator -=(const intmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE intmatrix_slice &intmatrix_slice::operator *=(const int &c) throw() { return _mssmultassign(*this,c); }
	INLINE intmatrix_slice &intmatrix_slice::operator /=(const int &c) throw() { return _mssdivassign(*this,c); }
	INLINE intmatrix_slice::operator void*() throw() { return _msvoid(*this); }
	INLINE intvector operator /(const intmatrix_subv &rv, const int &s) throw() { return _mvsdiv<intmatrix_subv,int,intvector>(rv,s); }
	INLINE intvector operator *(const intmatrix_subv &rv, const int &s) throw() { return _mvsmult<intmatrix_subv,int,intvector>(rv,s); }
	INLINE intvector operator *(const int &s, const intmatrix_subv &rv) throw() { return _mvsmult<intmatrix_subv,int,intvector>(rv,s); }
	INLINE intmatrix_subv &intmatrix_subv::operator *=(const int &c) throw() { return _mvsmultassign(*this,c); }
	INLINE intmatrix_subv &intmatrix_subv::operator +=(const int &c) throw() { return _mvsplusassign(*this,c); }
	INLINE intmatrix_subv &intmatrix_subv::operator -=(const int &c) throw() { return _mvsminusassign(*this,c); }
	INLINE intmatrix_subv &intmatrix_subv::operator /=(const int &c) throw() { return _mvsdivassign(*this,c); }
	INLINE intvector abs(const intmatrix_subv &mv) throw() { return _mvabs<intmatrix_subv,intvector>(mv); }
	INLINE intvector &intvector::operator =(const intmatrix_subv &mv) throw() { return _vmvassign<intvector,intmatrix_subv,int>(*this,mv); }
	INLINE intvector_slice &intvector_slice::operator =(const intmatrix_subv &mv) throw() { return _vsvassign(*this,intvector(mv)); }

	INLINE void accumulate(dotprecision &dp, const intmatrix_subv & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const intvector & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const intmatrix_subv & rv1, const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(dotprecision &dp,const intvector_slice &sl,const intmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,intvector(sl),sv); }
	INLINE void accumulate(dotprecision &dp,const intmatrix_subv &mv,const intvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,intvector(vs),mv); }
	INLINE intvector operator +(const intmatrix_subv & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<intmatrix_subv,intmatrix_subv,intvector>(rv1,rv2); }
	INLINE intvector operator +(const intmatrix_subv &rv1,const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<intmatrix_subv,intvector,intvector>(rv1,rv2); }
	INLINE intvector operator +(const intvector & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<intmatrix_subv,intvector,intvector>(rv2,rv1); }
	INLINE intvector operator +(const intvector_slice &sl,const intmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<intmatrix_subv,intvector,intvector>(mv,intvector(sl)); }
	INLINE intvector operator +(const intmatrix_subv &mv,const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<intmatrix_subv,intvector,intvector>(mv,intvector(sl)); }
	INLINE intmatrix_subv &intmatrix_subv::operator +=(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE intmatrix_subv &intmatrix_subv::operator +=(const intvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,intvector(rv)); }
	INLINE intvector operator -(const intmatrix_subv & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<intmatrix_subv,intmatrix_subv,intvector>(rv1,rv2); }
	INLINE intvector operator -(const intvector & rv1, const intmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<intvector,intmatrix_subv,intvector>(rv1,rv2); }
	INLINE intvector operator -(const intmatrix_subv &rv1,const intvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<intmatrix_subv,intvector,intvector>(rv1,rv2); }
	INLINE intvector operator -(const intvector_slice &sl,const intmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<intvector,intmatrix_subv,intvector>(intvector(sl),mv); }
	INLINE intvector operator -(const intmatrix_subv &mv,const intvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<intmatrix_subv,intvector,intvector>(mv,intvector(sl)); }
	INLINE intmatrix_subv &intmatrix_subv::operator -=(const intvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE intmatrix_subv &intmatrix_subv::operator -=(const intvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,intvector(rv)); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::intmatrix::intmatrix(const intmatrix &rm)
	*/
	INLINE intmatrix _intmatrix(const intmatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::intmatrix::intmatrix(const intvector &v)
	*/
	INLINE intmatrix _intmatrix(const intvector &v) throw() { return intmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::intmatrix::intmatrix(const intvector_slice &v)
	*/
	INLINE intmatrix _intmatrix(const intvector_slice &v) throw() { return intmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::intmatrix::intmatrix(const int &r)
	*/
	INLINE intmatrix _intmatrix(const int &r) throw() { return intmatrix(r); }
	INLINE intmatrix &intmatrix::operator =(const intmatrix_slice &ms) throw() { return _mmsassign<intmatrix,intmatrix_slice,int>(*this,ms); }
	INLINE int Lb(const intmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const intmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const intmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const intmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE intmatrix &SetLb(intmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE intmatrix &SetUb(intmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }
	
	INLINE int RowLen ( const intmatrix& A ) // Length of the rows of a integer matrix
        { return Ub(A,2)-Lb(A,2)+1; }            //---------------------------------------

        INLINE int ColLen ( const intmatrix& A ) // Length of the columns of an integer matrix
        { return Ub(A,1)-Lb(A,1)+1; }            //-------------------------------------------

	INLINE int RowLen ( const intmatrix_slice& A ) // Length of the rows of a integer matrix
        { return Ub(A,2)-Lb(A,2)+1; }                  //---------------------------------------

        INLINE int ColLen ( const intmatrix_slice& A ) // Length of the columns of an integer matrix
        { return Ub(A,1)-Lb(A,1)+1; }                  //-------------------------------------------
	
	INLINE void Resize(intmatrix &A) throw() { _mresize(A); }
	INLINE void Resize(intmatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<intmatrix,int>(A,m,n); }
	INLINE void Resize(intmatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<intmatrix,int>(A,m1,m2,n1,n2); }
	INLINE intmatrix abs(const intmatrix &m) throw() { return _mabs<intmatrix,intmatrix>(m); }
	INLINE intmatrix abs(const intmatrix_slice &ms) throw() { return _msabs<intmatrix_slice,intmatrix>(ms); }
	INLINE intmatrix operator *(const int &c, const intmatrix &m) throw() { return _smmult<int,intmatrix,intmatrix>(c,m); }
	INLINE intmatrix operator *(const int &c, const intmatrix_slice &ms) throw() { return _smsmult<int,intmatrix_slice,intmatrix>(c,ms); }
	INLINE intmatrix operator *(const intmatrix &m,const int &c) throw() { return _smmult<int,intmatrix,intmatrix>(c,m); }
	INLINE intmatrix operator *(const intmatrix_slice &ms,const int &c) throw() { return _smsmult<int,intmatrix_slice,intmatrix>(c,ms); }
	INLINE intmatrix &operator *=(intmatrix &m,const int &c) throw() { return _msmultassign(m,c); }
	INLINE intmatrix operator /(const intmatrix &m,const int &c) throw() { return _msdiv<intmatrix,int,intmatrix>(m,c); }
	INLINE intmatrix operator /(const intmatrix_slice &ms, const int &c) throw() { return _mssdiv<intmatrix_slice,int,intmatrix>(ms,c); }
	INLINE intmatrix &operator /=(intmatrix &m,const int &c) throw() { return _msdivassign(m,c); }
	INLINE intvector::intvector(const intmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<intvector,intmatrix,int>(*this,sl); }
	INLINE intvector::intvector(const intmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<intvector,intmatrix_slice,int>(*this,sl); }
	INLINE intvector &intvector::operator =(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<intvector,intmatrix,int>(*this,m); }
	INLINE intvector &intvector::operator =(const intmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<intvector,intmatrix,int>(*this,intmatrix(m)); }
	INLINE intvector_slice &intvector_slice::operator =(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>,ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,intvector(m)); }
	INLINE intvector_slice & intvector_slice::operator =(const intmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<intvector>,ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,intvector(intmatrix(m))); }
	INLINE intmatrix_subv &intmatrix_subv::operator =(const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,intvector(m)); }
	INLINE intmatrix_subv &intmatrix_subv::operator =(const intmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,intvector(intmatrix(m))); }

	INLINE const intmatrix &operator +(const intmatrix &m) throw() { return m; }
	INLINE intmatrix operator +(const intmatrix_slice &m) throw() { return intmatrix(m); }
	INLINE intmatrix operator +(const intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<intmatrix,intmatrix,intmatrix>(m1,m2); }
	INLINE intmatrix operator +(const intmatrix &m,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<intmatrix,intmatrix_slice,intmatrix>(m,ms); }
	INLINE intmatrix operator +(const intmatrix_slice &ms,const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<intmatrix,intmatrix_slice,intmatrix>(m,ms); }
	INLINE intmatrix operator +(const intmatrix_slice &m1,const intmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<intmatrix_slice,intmatrix_slice,intmatrix>(m1,m2); }
	INLINE intmatrix &operator +=(intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE intmatrix &operator +=(intmatrix &m1,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE intmatrix operator -(const intmatrix &m) throw() { return _mminus(m); }
	INLINE intmatrix operator -(const intmatrix_slice &m) throw() { return _msminus<intmatrix_slice,intmatrix>(m); }
	INLINE intmatrix operator -(const intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<intmatrix,intmatrix,intmatrix>(m1,m2); }
	INLINE intmatrix operator -(const intmatrix &m,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<intmatrix,intmatrix_slice,intmatrix>(m,ms); }
	INLINE intmatrix operator -(const intmatrix_slice &ms,const intmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<intmatrix_slice,intmatrix,intmatrix>(ms,m); }
	INLINE intmatrix operator -(const intmatrix_slice &ms1,const intmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<intmatrix_slice,intmatrix_slice,intmatrix>(ms1,ms2); }
	INLINE intmatrix &operator -=(intmatrix &m1,const intmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE intmatrix &operator -=(intmatrix &m1,const intmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_INTMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE bool operator ==(const intmatrix &m1,const intmatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const intmatrix &m1,const intmatrix &m2) throw() { return _mmneq(m1,m2); }
	INLINE bool operator <(const intmatrix &m1,const intmatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const intmatrix &m1,const intmatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const intmatrix &m1,const intmatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const intmatrix &m1,const intmatrix &m2) throw() { return _mmleq(m2,m1); }
	INLINE bool operator ==(const intmatrix &m1,const intmatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const intmatrix &m1,const intmatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
	INLINE bool operator <(const intmatrix &m1,const intmatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const intmatrix &m1,const intmatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const intmatrix &m1,const intmatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const intmatrix &m1,const intmatrix_slice &ms) throw() { return _msmleq(ms,m1); }
	INLINE bool operator ==(const intmatrix_slice &m1,const intmatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const intmatrix_slice &m1,const intmatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
	INLINE bool operator <(const intmatrix_slice &m1,const intmatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const intmatrix_slice &m1,const intmatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const intmatrix_slice &m1,const intmatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const intmatrix_slice &m1,const intmatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
	INLINE bool operator !(const intmatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const intmatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const intmatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const intmatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,intmatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,intmatrix_slice &r) throw() { return _msin(s,r); }

        INLINE intvector permvec(const intmatrix& A) {
          intvector p(RowLen(A));
          SetLb(p,0);
          for(int i=0 ; i<ColLen(A) ; i++)
            for(int j=0 ; j<RowLen(A) ; j++)
              if(A[i+Lb(A,1)][j+Lb(A,2)] != 0) {
                p[i] = j;
                j = RowLen(A);
              }
          return p;
        }

        INLINE intmatrix permmat(const intvector& x) {
          intmatrix A(0,VecLen(x)-1,0,VecLen(x)-1);
          for(int i=0 ; i<VecLen(x) ; i++)
              A[i][x[i+Lb(x)]] = 1;
          return A;
        }

        INLINE intmatrix perminv(const intmatrix& A) {
          return transp(A);
        }

} // namespace cxsc

#endif
