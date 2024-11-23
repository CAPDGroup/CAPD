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

/* CVS $Id: rmatrix.inl,v 1.29 2014/01/30 17:23:48 cxsc Exp $ */

#ifndef _CXSC_RMATRIX_INL_INCLUDED
#define _CXSC_RMATRIX_INL_INCLUDED

#include "intmatrix.hpp"

namespace cxsc {

INLINE rmatrix::rmatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE rmatrix::rmatrix(const real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new real[1];
	*dat=r;
}

INLINE rmatrix::rmatrix(const rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new real[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE rmatrix::rmatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_RMATRIX_WRONG_BOUNDARIES("rmatrix::rmatrix(const int &m, const int &n)"));
#endif
	dat=new real[m*n];
}

INLINE rmatrix::rmatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_RMATRIX_WRONG_BOUNDARIES("rmatrix::rmatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new real[xsize*ysize];
}

INLINE rvector::rvector(const rmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new real[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE rmatrix::rmatrix(const rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new real[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE rmatrix::rmatrix(const rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new real[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}


	INLINE rmatrix::rmatrix(const rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new real[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

INLINE rmatrix::rmatrix(const intmatrix& I) : lb1(Lb(I,1)),ub1(Ub(I,1)),lb2(Lb(I,2)),ub2(Ub(I,2)),xsize(RowLen(I)),ysize(ColLen(I)) {
	dat=new real[xsize*ysize];
	for(int i=0 ; i<ysize ; i++)
          for(int j=0 ; j<xsize ; j++)
            dat[i*xsize+j] = I[i+lb1][j+lb2];
}

	INLINE rmatrix_subv Row(rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE rmatrix_subv Col(rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
	INLINE rmatrix_subv Row(const rmatrix &m,const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE rmatrix_subv Col(const rmatrix &m,const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE real& rmatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC("real &rmatrix_subv::operator [](const int &i) const"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE real& rmatrix_subv::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_RVECTOR_ELEMENT_NOT_IN_VEC("real &rmatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}


		
	INLINE rmatrix_subv rmatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT("rmatrix_subv rmatrix::operator [](const int &i)"));
#endif
		return rmatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE rmatrix_subv rmatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT("rmatrix_subv rmatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return rmatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}
	
	INLINE rmatrix_slice rmatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG("rmatrix_slice rmatrix::operator ()(const int &m, const int &n)"));
#endif
		return rmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE rmatrix_slice rmatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG("rmatrix_slice rmatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return rmatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE rmatrix_subv rmatrix_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT("rmatrix_subv rmatrix_slice::operator [](const int &i)"));
#endif
		return rmatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE rmatrix_subv rmatrix_slice::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_RMATRIX_ROW_OR_COL_NOT_IN_MAT("rmatrix_subv rmatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return rmatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}
	
	INLINE rmatrix_slice rmatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG("rmatrix_slice rmatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return rmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE rmatrix_slice rmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_RMATRIX_SUB_ARRAY_TOO_BIG("rmatrix_slice rmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return rmatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE rmatrix_subv rmatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG("rmatrix_subv rmatrix_subv::operator ()(const int &i)"));
#endif
	return rmatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE rmatrix_subv rmatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_RVECTOR_SUB_ARRAY_TOO_BIG("rmatrix_subv rmatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return rmatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}

// the following is generated from .hpp

	INLINE rmatrix_subv &rmatrix_subv::operator =(const rmatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE rmatrix_subv &rmatrix_subv::operator =(const real &r) throw() { return _mvsassign(*this,r); }
	INLINE rmatrix_subv &rmatrix_subv::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE rmatrix_subv &rmatrix_subv::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(v)); }
	INLINE rmatrix &rmatrix::operator =(const real &r) throw() { return _msassign(*this,r); }
	INLINE rmatrix &rmatrix::operator =(const rmatrix &m) throw() { return _mmassign<rmatrix,rmatrix,real>(*this,m,0); }
	INLINE rmatrix &rmatrix::operator =(const rvector &v) throw() { return _mvassign<rmatrix,rvector,real>(*this,v); }
	INLINE rmatrix &rmatrix::operator =(const rvector_slice &v) throw() { return _mvassign<rmatrix,rvector,real>(*this,rvector(v)); }
	INLINE rmatrix::operator void*() throw() { return _mvoid(*this); }
	INLINE rmatrix_slice &rmatrix_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE rmatrix_slice &rmatrix_slice::operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE rmatrix_slice &rmatrix_slice::operator =(const real &r) throw() { return _mssassign(*this,r); }
	INLINE rmatrix_slice &rmatrix_slice::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(v)); }
	INLINE rmatrix_slice &rmatrix_slice::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(rvector(v))); }
	INLINE rmatrix_slice &rmatrix_slice::operator =(const rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(rvector(v))); }
	INLINE rmatrix_slice &rmatrix_slice::operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE rmatrix_slice &rmatrix_slice::operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE rmatrix_slice &rmatrix_slice::operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE rmatrix_slice &rmatrix_slice::operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE rmatrix_slice &rmatrix_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE rmatrix_slice &rmatrix_slice::operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE rmatrix_slice &rmatrix_slice::operator *=(const real &c) throw() { return _mssmultassign(*this,c); }
	INLINE rmatrix_slice &rmatrix_slice::operator /=(const real &c) throw() { return _mssdivassign(*this,c); }
	INLINE rmatrix_slice::operator void*() throw() { return _msvoid(*this); }
	INLINE rvector operator /(const rmatrix_subv &rv, const real &s) throw() { return _mvsdiv<rmatrix_subv,real,rvector>(rv,s); }
	INLINE rvector operator *(const rmatrix_subv &rv, const real &s) throw() { return _mvsmult<rmatrix_subv,real,rvector>(rv,s); }
	INLINE rvector operator *(const real &s, const rmatrix_subv &rv) throw() { return _mvsmult<rmatrix_subv,real,rvector>(rv,s); }
	INLINE rmatrix_subv &rmatrix_subv::operator *=(const real &c) throw() { return _mvsmultassign(*this,c); }
	INLINE rmatrix_subv &rmatrix_subv::operator +=(const real &c) throw() { return _mvsplusassign(*this,c); }
	INLINE rmatrix_subv &rmatrix_subv::operator -=(const real &c) throw() { return _mvsminusassign(*this,c); }
	INLINE rmatrix_subv &rmatrix_subv::operator /=(const real &c) throw() { return _mvsdivassign(*this,c); }
	INLINE rvector abs(const rmatrix_subv &mv) throw() { return _mvabs<rmatrix_subv,rvector>(mv); }
	INLINE rvector &rvector::operator =(const rmatrix_subv &mv) throw() { return _vmvassign<rvector,rmatrix_subv,real>(*this,mv); }
	INLINE rvector_slice &rvector_slice::operator =(const rmatrix_subv &mv) throw() { return _vsvassign(*this,rvector(mv)); }


	INLINE real operator *(const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvmult<rmatrix_subv,rmatrix_subv,real>(rv1,rv2); }
	INLINE real operator *(const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvmult<rvector,rmatrix_subv,real>(rv1,rv2); }
	INLINE real operator *(const rmatrix_subv &rv1,const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvmult<rvector,rmatrix_subv,real>(rv2,rv1); }
	INLINE real operator *(const rvector_slice &sl,const rmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvmult<rvector,rmatrix_subv,real>(rvector(sl),sv); }
	INLINE real operator *(const rmatrix_subv &mv,const rvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvmult<rvector,rmatrix_subv,real>(rvector(vs),mv); }
	INLINE rvector operator +(const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<rmatrix_subv,rmatrix_subv,rvector>(rv1,rv2); }
	INLINE rvector operator +(const rmatrix_subv &rv1,const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<rmatrix_subv,rvector,rvector>(rv1,rv2); }
	INLINE rvector operator +(const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<rmatrix_subv,rvector,rvector>(rv2,rv1); }
	INLINE rvector operator +(const rvector_slice &sl,const rmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<rmatrix_subv,rvector,rvector>(mv,rvector(sl)); }
	INLINE rvector operator +(const rmatrix_subv &mv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<rmatrix_subv,rvector,rvector>(mv,rvector(sl)); }
	INLINE rmatrix_subv &rmatrix_subv::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE rmatrix_subv &rmatrix_subv::operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rvector(rv)); }
	INLINE rvector operator -(const rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<rmatrix_subv,rmatrix_subv,rvector>(rv1,rv2); }
	INLINE rvector operator -(const rvector & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<rvector,rmatrix_subv,rvector>(rv1,rv2); }
	INLINE rvector operator -(const rmatrix_subv &rv1,const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<rmatrix_subv,rvector,rvector>(rv1,rv2); }
	INLINE rvector operator -(const rvector_slice &sl,const rmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<rvector,rmatrix_subv,rvector>(rvector(sl),mv); }
	INLINE rvector operator -(const rmatrix_subv &mv,const rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<rmatrix_subv,rvector,rvector>(mv,rvector(sl)); }
	INLINE rmatrix_subv &rmatrix_subv::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE rmatrix_subv &rmatrix_subv::operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rvector(rv)); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::rmatrix::rmatrix(const rmatrix &rm)
	*/
	INLINE rmatrix _rmatrix(const rmatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::rmatrix::rmatrix(const rvector &v)
	*/
	INLINE rmatrix _rmatrix(const rvector &v) throw() { return rmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::rmatrix::rmatrix(const rvector_slice &v)
	*/
	INLINE rmatrix _rmatrix(const rvector_slice &v) throw() { return rmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::rmatrix::rmatrix(const real &r)
	*/
	INLINE rmatrix _rmatrix(const real &r) throw() { return rmatrix(r); }
	INLINE rmatrix &rmatrix::operator =(const rmatrix_slice &ms) throw() { return _mmsassign<rmatrix,rmatrix_slice,real>(*this,ms); }
	INLINE int Lb(const rmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const rmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const rmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const rmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE rmatrix &SetLb(rmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE rmatrix &SetUb(rmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }


        INLINE int RowLen ( const rmatrix& A ) // Length of the rows of a real matrix
        { return Ub(A,2)-Lb(A,2)+1; }          //------------------------------------
      
        INLINE int ColLen ( const rmatrix& A ) // Length of the columns of a real matrix
        { return Ub(A,1)-Lb(A,1)+1; }          //---------------------------------------

        INLINE int RowLen ( const rmatrix_slice& A ) // Length of the rows of a real matrix
        { return Ub(A,2)-Lb(A,2)+1; }                //------------------------------------
      
        INLINE int ColLen ( const rmatrix_slice& A ) // Length of the columns of a real matrix
        { return Ub(A,1)-Lb(A,1)+1; }                //---------------------------------------

	INLINE void Resize(rmatrix &A) throw() { _mresize(A); }
	INLINE void Resize(rmatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<rmatrix,real>(A,m,n); }
	INLINE void Resize(rmatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<rmatrix,real>(A,m1,m2,n1,n2); }
	INLINE rmatrix abs(const rmatrix &m) throw() { return _mabs<rmatrix,rmatrix>(m); }
	INLINE rmatrix abs(const rmatrix_slice &ms) throw() { return _msabs<rmatrix_slice,rmatrix>(ms); }
	INLINE real::real(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{ _smconstr(*this,m); }
//	INLINE real real::_real(const rmatrix &m) throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_RMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE rmatrix operator *(const real &c, const rmatrix &m) throw() { return _smmult<real,rmatrix,rmatrix>(c,m); }
	INLINE rmatrix operator *(const real &c, const rmatrix_slice &ms) throw() { return _smsmult<real,rmatrix_slice,rmatrix>(c,ms); }
	INLINE rmatrix operator *(const rmatrix &m,const real &c) throw() { return _smmult<real,rmatrix,rmatrix>(c,m); }
	INLINE rmatrix operator *(const rmatrix_slice &ms,const real &c) throw() { return _smsmult<real,rmatrix_slice,rmatrix>(c,ms); }
	INLINE rmatrix &operator *=(rmatrix &m,const real &c) throw() { return _msmultassign(m,c); }
	INLINE rmatrix operator /(const rmatrix &m,const real &c) throw() { return _msdiv<rmatrix,real,rmatrix>(m,c); }
	INLINE rmatrix operator /(const rmatrix_slice &ms, const real &c) throw() { return _mssdiv<rmatrix_slice,real,rmatrix>(ms,c); }
	INLINE rmatrix &operator /=(rmatrix &m,const real &c) throw() { return _msdivassign(m,c); }
	INLINE rvector::rvector(const rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<rvector,rmatrix,real>(*this,sl); }
	INLINE rvector::rvector(const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<rvector,rmatrix_slice,real>(*this,sl); }
	INLINE rvector &rvector::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<rvector,rmatrix,real>(*this,m); }
	INLINE rvector &rvector::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<rvector,rmatrix,real>(*this,rmatrix(m)); }
	INLINE rvector_slice &rvector_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,rvector(m)); }
	INLINE rvector_slice & rvector_slice::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<rvector>,ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,rvector(rmatrix(m))); }
	INLINE rmatrix_subv &rmatrix_subv::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(m)); }
	INLINE rmatrix_subv &rmatrix_subv::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(rmatrix(m))); }
	INLINE rvector operator *(const rmatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmult<rmatrix,rvector,rvector>(m,v); };
	INLINE rvector operator *(const rmatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvmult<rmatrix_slice,rvector,rvector>(ms,v); }
	INLINE rvector operator *(const rvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmmult<rvector,rmatrix,rvector>(v,m); }
	INLINE rvector operator *(const rvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsmult<rvector,rmatrix_slice,rvector>(v,ms); }
	INLINE rvector &operator *=(rvector &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmmultassign<rvector,rmatrix,real>(v,m); }
	INLINE rvector &operator *=(rvector &v,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsmultassign<rvector,rmatrix_slice,real>(v,ms); }
	INLINE rvector_slice &rvector_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmmultassign<rvector_slice,rmatrix,real>(*this,m); }
	INLINE rvector operator *(const rvector_slice &v,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmmult<rvector,rmatrix,rvector>(rvector(v),m); }

	INLINE const rmatrix &operator +(const rmatrix &m) throw() { return m; }
	INLINE rmatrix operator +(const rmatrix_slice &m) throw() { return rmatrix(m); }
	INLINE rmatrix operator +(const rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,rmatrix,rmatrix>(m1,m2); }
	INLINE rmatrix operator +(const rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,rmatrix_slice,rmatrix>(m,ms); }
	INLINE rmatrix operator +(const rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,rmatrix_slice,rmatrix>(m,ms); }
	INLINE rmatrix operator +(const rmatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,rmatrix_slice,rmatrix>(m1,m2); }
	INLINE rmatrix &operator +=(rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE rmatrix &operator +=(rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE rmatrix operator -(const rmatrix &m) throw() { return _mminus(m); }
	INLINE rmatrix operator -(const rmatrix_slice &m) throw() { return _msminus<rmatrix_slice,rmatrix>(m); }
	INLINE rmatrix operator -(const rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<rmatrix,rmatrix,rmatrix>(m1,m2); }
	INLINE rmatrix operator -(const rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<rmatrix,rmatrix_slice,rmatrix>(m,ms); }
	INLINE rmatrix operator -(const rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<rmatrix_slice,rmatrix,rmatrix>(ms,m); }
	INLINE rmatrix operator -(const rmatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<rmatrix_slice,rmatrix_slice,rmatrix>(ms1,ms2); }
	INLINE rmatrix &operator -=(rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE rmatrix &operator -=(rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE rmatrix operator *(const rmatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmmult<rmatrix,rmatrix,rmatrix>(m1,m2); }
	INLINE rmatrix operator *(const rmatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsmult<rmatrix,rmatrix_slice,rmatrix>(m1,ms); }
	INLINE rmatrix operator *(const rmatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmmult<rmatrix_slice,rmatrix,rmatrix>(ms,m1); }
	INLINE rmatrix operator *(const rmatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsmult<rmatrix_slice,rmatrix_slice,rmatrix>(ms1,ms2); }
	INLINE rmatrix &operator *=(rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmmultassign<rmatrix,rmatrix,real>(m1,m2); }
	INLINE rmatrix &operator *=(rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsmultassign<rmatrix,rmatrix_slice,real>(m1,ms); }
	INLINE bool operator ==(const rmatrix &m1,const rmatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const rmatrix &m1,const rmatrix &m2) throw() { return _mmneq(m1,m2); }
	INLINE bool operator <(const rmatrix &m1,const rmatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const rmatrix &m1,const rmatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const rmatrix &m1,const rmatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const rmatrix &m1,const rmatrix &m2) throw() { return _mmleq(m2,m1); }
	INLINE bool operator ==(const rmatrix &m1,const rmatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const rmatrix &m1,const rmatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
	INLINE bool operator <(const rmatrix &m1,const rmatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const rmatrix &m1,const rmatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const rmatrix &m1,const rmatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const rmatrix &m1,const rmatrix_slice &ms) throw() { return _msmleq(ms,m1); }
	INLINE bool operator ==(const rmatrix_slice &m1,const rmatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const rmatrix_slice &m1,const rmatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
	INLINE bool operator <(const rmatrix_slice &m1,const rmatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const rmatrix_slice &m1,const rmatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const rmatrix_slice &m1,const rmatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const rmatrix_slice &m1,const rmatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
	INLINE bool operator !(const rmatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const rmatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const rmatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const rmatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,rmatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,rmatrix_slice &r) throw() { return _msin(s,r); }

        //! Computes permutation of matrix according to permutation vectors, C=PAQ
        INLINE rmatrix rmatrix::operator()(const intvector& p, const intvector& q) {
          rmatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
            for(int j=0 ; j<RowLen(A) ; j++)
              A[i+Lb(A,1)][j+Lb(A,2)] = (*this)[p[i+Lb(p)]+Lb(A,1)][q[j+Lb(q)]+Lb(A,2)];
          return A;
        }

        //! Computes permutation of matrix according to permutation vector, C=PA
        INLINE rmatrix rmatrix::operator()(const intvector& p) {
          rmatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
              A[i+Lb(A,1)] = (*this)[p[i+Lb(p)]+Lb(A,1)];
          return A;
        }

        //! Computes permutation of matrix according to permutation matrix, C=PA
	INLINE rmatrix rmatrix::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

        //! Computes permutation of matrix according to permutation matrices, C=PAQ
        INLINE rmatrix rmatrix::operator()(const intmatrix& P, const intmatrix& Q) {
          intvector p = permvec(P);
          intvector q = perminv(permvec(Q));
          return (*this)(p,q);
        }

        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE rvector rvector::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

} // namespace cxsc

#endif

