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

/* CVS $Id: imatrix.inl,v 1.30 2014/01/30 17:23:45 cxsc Exp $ */

#ifndef _CXSC_IMATRIX_INL_INCLUDED
#define _CXSC_IMATRIX_INL_INCLUDED

namespace cxsc {

INLINE imatrix::imatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE imatrix::imatrix(const interval &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new interval[1];
	*dat=r;
}

INLINE imatrix::imatrix(const real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new interval[1];
	*dat=r;
}

INLINE imatrix::imatrix(const imatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new interval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE imatrix::imatrix(const rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new interval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE imatrix::imatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_IMATRIX_WRONG_BOUNDARIES("imatrix::imatrix(const int &m, const int &n)"));
#endif
	dat=new interval[m*n];
}

INLINE imatrix::imatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_IMATRIX_WRONG_BOUNDARIES("imatrix::imatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new interval[xsize*ysize];
}

INLINE ivector::ivector(const imatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new interval[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE imatrix::imatrix(const ivector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new interval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE imatrix::imatrix(const rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new interval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE imatrix::imatrix(const ivector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new interval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE imatrix::imatrix(const rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new interval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}



	INLINE imatrix::imatrix(const imatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new interval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE imatrix::imatrix(const rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new interval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE imatrix_subv Row(imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE imatrix_subv Col(imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE imatrix_subv Row(const imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE imatrix_subv Col(const imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE interval& imatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC("interval &imatrix_subv::operator [](const int &i) const"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE interval& imatrix_subv::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC("interval &imatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}

		
	INLINE imatrix_subv imatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT("imatrix_subv imatrix::operator [](const int &i)"));
#endif
		return imatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE imatrix_subv imatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT("imatrix_subv imatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return imatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}
	
	INLINE imatrix_slice imatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG("imatrix_slice imatrix::operator ()(const int &m, const int &n)"));
#endif
		return imatrix_slice(*this,1,m,1,n);
	}
	
	INLINE imatrix_slice imatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG("imatrix_slice imatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return imatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE imatrix_subv imatrix_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT("imatrix_subv imatrix_slice::operator [](const int &i)"));
#endif
		return imatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE imatrix_subv imatrix_slice::operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT("imatrix_subv imatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return imatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}

	INLINE imatrix_subv imatrix_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT("imatrix_subv imatrix_slice::operator [](const int &i)"));
#endif
		return imatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE imatrix_subv imatrix_slice::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_IMATRIX_ROW_OR_COL_NOT_IN_MAT("imatrix_subv imatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return imatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}

	
	INLINE imatrix_slice imatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG("imatrix_slice imatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return imatrix_slice(*this,1,m,1,n);
	}
	
	INLINE imatrix_slice imatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_IMATRIX_SUB_ARRAY_TOO_BIG("imatrix_slice imatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return imatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE imatrix_subv imatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG("imatrix_subv imatrix_subv::operator ()(const int &i)"));
#endif
	return imatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE imatrix_subv imatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_IVECTOR_SUB_ARRAY_TOO_BIG("imatrix_subv imatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return imatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}

	INLINE imatrix_subv &imatrix_subv::operator =(const imatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator =(const interval &r) throw() { return _mvsassign(*this,r); }
	INLINE imatrix_subv &imatrix_subv::operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE imatrix_subv &imatrix_subv::operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(v)); }
	INLINE imatrix_subv &imatrix_subv::operator =(const rmatrix_subv &rv) throw() { return _mvvassign(*this,rvector(rv)); }
	INLINE imatrix_subv &imatrix_subv::operator =(const real &r) throw() { return _mvsassign(*this,r); }
	INLINE imatrix_subv &imatrix_subv::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE imatrix_subv &imatrix_subv::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(v)); }
	INLINE imatrix_subv &imatrix_subv::operator *=(const real &c) throw() { return _mvsmultassign(*this,c); }
	INLINE imatrix_subv &imatrix_subv::operator +=(const real &c) throw() { return _mvsplusassign(*this,c); }
	INLINE imatrix_subv &imatrix_subv::operator -=(const real &c) throw() { return _mvsminusassign(*this,c); }
	INLINE imatrix_subv &imatrix_subv::operator /=(const real &c) throw() { return _mvsdivassign(*this,c); }
	INLINE imatrix &imatrix::operator =(const interval &r) throw() { return _msassign(*this,r); }
	INLINE imatrix &imatrix::operator =(const imatrix &m) throw() { return _mmassign<imatrix,imatrix,interval>(*this,m,interval(0,0)); }
	INLINE imatrix &imatrix::operator =(const imatrix_slice &ms) throw() { return _mmsassign<imatrix,imatrix_slice,interval>(*this,ms); }
	INLINE imatrix &imatrix::operator =(const ivector &v) throw() { return _mvassign<imatrix,ivector,interval>(*this,v); }
	INLINE imatrix &imatrix::operator =(const ivector_slice &v) throw() { return _mvassign<imatrix,ivector,interval>(*this,ivector(v)); }
	INLINE imatrix &imatrix::operator =(const real &r) throw() { return _msassign(*this,interval(r)); }
	INLINE imatrix &imatrix::operator =(const rmatrix &m) throw() { return _mmassign<imatrix,rmatrix,interval>(*this,m,interval(0,0)); }
	INLINE imatrix &imatrix::operator =(const rmatrix_slice &ms) throw() { return _mmsassign<imatrix,rmatrix_slice,interval>(*this,ms); }
	INLINE imatrix &imatrix::operator =(const rvector &v) throw() { return _mvassign<imatrix,rvector,interval>(*this,v); }
	INLINE imatrix &imatrix::operator =(const rvector_slice &v) throw() { return _mvassign<imatrix,rvector,interval>(*this,rvector(v)); }
	INLINE imatrix::operator void*() throw() { return _mvoid(*this); }

	INLINE imatrix_slice &imatrix_slice::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE imatrix_slice &imatrix_slice::operator =(const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE imatrix_slice &imatrix_slice::operator =(const interval &r) throw() { return _mssassign(*this,r); }
	INLINE imatrix_slice &imatrix_slice::operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,imatrix(v)); }
	INLINE imatrix_slice &imatrix_slice::operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,imatrix(ivector(v))); }
	INLINE imatrix_slice &imatrix_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE imatrix_slice &imatrix_slice::operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE imatrix_slice &imatrix_slice::operator =(const real &r) throw() { return _mssassign(*this,r); }
	INLINE imatrix_slice &imatrix_slice::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(v)); }
	INLINE imatrix_slice &imatrix_slice::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(rvector(v))); }
	INLINE imatrix_slice::operator void*() throw() { return _msvoid(*this); }
	
	INLINE ivector operator /(const imatrix_subv &rv, const interval &s) throw() { return _mvsdiv<imatrix_subv,interval,ivector>(rv,s); }
	INLINE ivector operator *(const imatrix_subv &rv, const interval &s) throw() { return _mvsmult<imatrix_subv,interval,ivector>(rv,s); }
	INLINE ivector operator *(const interval &s, const imatrix_subv &rv) throw() { return _mvsmult<imatrix_subv,interval,ivector>(rv,s); }
	INLINE imatrix_subv &imatrix_subv::operator *=(const interval &c) throw() { return _mvsmultassign(*this,c); }
	INLINE imatrix_subv &imatrix_subv::operator +=(const interval &c) throw() { return _mvsplusassign(*this,c); }
	INLINE imatrix_subv &imatrix_subv::operator -=(const interval &c) throw() { return _mvsminusassign(*this,c); }
	INLINE imatrix_subv &imatrix_subv::operator /=(const interval &c) throw() { return _mvsdivassign(*this,c); }
	INLINE ivector abs(const imatrix_subv &mv) throw() { return _mvabs<imatrix_subv,ivector>(mv); }
	INLINE rvector absmin(const imatrix_subv &mv) throw() { 
          rvector x(Lb(mv),Ub(mv));
          for(int i=Lb(mv) ; i<=Ub(mv) ; i++)
            x[i] = AbsMin(mv[i]);
          return x;
        }
	INLINE rvector absmax(const imatrix_subv &mv) throw() { 
          rvector x(Lb(mv),Ub(mv));
          for(int i=Lb(mv) ; i<=Ub(mv) ; i++)
            x[i] = AbsMax(mv[i]);
          return x;
        }
	INLINE rvector diam(const imatrix_subv &mv) throw() { return _mvdiam<imatrix_subv,rvector>(mv); }
	INLINE rvector mid(const imatrix_subv &mv) throw() { return _mvmid<imatrix_subv,rvector>(mv); }
	INLINE rvector Inf(const imatrix_subv &mv) throw() { return _mvinf<imatrix_subv,rvector>(mv); }
	INLINE rvector Sup(const imatrix_subv &mv) throw() { return _mvsup<imatrix_subv,rvector>(mv); }
	INLINE imatrix_subv &SetInf(imatrix_subv &mv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetinf(mv,rv); }
	INLINE imatrix_subv &SetSup(imatrix_subv &mv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetsup(mv,rv); }
	INLINE imatrix_subv &UncheckedSetInf(imatrix_subv &mv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvusetinf(mv,rv); }
	INLINE imatrix_subv &UncheckedSetSup(imatrix_subv &mv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvusetsup(mv,rv); }
	INLINE imatrix_subv &SetSup(imatrix_subv &iv,const real &r) throw() { return _mvssetsup(iv,r); }
	INLINE imatrix_subv &SetInf(imatrix_subv &iv,const real &r) throw() { return _mvssetinf(iv,r); }
	INLINE imatrix_subv &UncheckedSetSup(imatrix_subv &iv,const real &r) throw() { return _mvsusetsup(iv,r); }
	INLINE imatrix_subv &SetUncheckedInf(imatrix_subv &iv,const real &r) throw() { return _mvsusetinf(iv,r); }
	INLINE ivector &ivector::operator =(const imatrix_subv &mv) throw() { return _vmvassign<ivector,imatrix_subv,interval>(*this,mv); }
	INLINE ivector_slice &ivector_slice::operator =(const imatrix_subv &mv) throw() { return _vsvassign(*this,ivector(mv)); }


	INLINE interval operator *(const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvimult<imatrix_subv,imatrix_subv,interval>(rv1,rv2); }
	INLINE interval operator *(const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvimult<ivector,imatrix_subv,interval>(rv1,rv2); }
	INLINE interval operator *(const imatrix_subv &rv1,const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvimult<ivector,imatrix_subv,interval>(rv2,rv1); }
	INLINE interval operator *(const ivector_slice &sl,const imatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvimult<ivector,imatrix_subv,interval>(ivector(sl),sv); }
	INLINE interval operator *(const imatrix_subv &mv,const ivector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvimult<ivector,imatrix_subv,interval>(ivector(vs),mv); }
	INLINE ivector operator +(const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<imatrix_subv,imatrix_subv,ivector>(rv1,rv2); }
	INLINE ivector operator +(const imatrix_subv &rv1,const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<imatrix_subv,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator +(const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<imatrix_subv,ivector,ivector>(rv2,rv1); }
	INLINE ivector operator +(const ivector_slice &sl,const imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE ivector operator +(const imatrix_subv &mv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE imatrix_subv &imatrix_subv::operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator +=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,ivector(rv)); }
	INLINE ivector operator -(const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<imatrix_subv,imatrix_subv,ivector>(rv1,rv2); }
	INLINE ivector operator -(const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<ivector,imatrix_subv,ivector>(rv1,rv2); }
	INLINE ivector operator -(const imatrix_subv &rv1,const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<imatrix_subv,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator -(const ivector_slice &sl,const imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<ivector,imatrix_subv,ivector>(ivector(sl),mv); }
	INLINE ivector operator -(const imatrix_subv &mv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE imatrix_subv &imatrix_subv::operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator -=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,ivector(rv)); }
	INLINE ivector operator |(const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvconv<imatrix_subv,imatrix_subv,ivector>(rv1,rv2); }
	INLINE ivector operator |(const imatrix_subv &rv1,const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<imatrix_subv,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator |(const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<imatrix_subv,ivector,ivector>(rv2,rv1); }
	INLINE ivector operator |(const ivector_slice &sl,const imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE ivector operator |(const imatrix_subv &mv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE imatrix_subv &imatrix_subv::operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator |=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,ivector(rv)); }
	INLINE ivector operator &(const imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvsect<imatrix_subv,imatrix_subv,ivector>(rv1,rv2); }
	INLINE ivector operator &(const imatrix_subv &rv1,const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<imatrix_subv,ivector,ivector>(rv1,rv2); }
	INLINE ivector operator &(const ivector & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<imatrix_subv,ivector,ivector>(rv2,rv1); }
	INLINE ivector operator &(const ivector_slice &sl,const imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE ivector operator &(const imatrix_subv &mv,const ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<imatrix_subv,ivector,ivector>(mv,ivector(sl)); }
	INLINE imatrix_subv &imatrix_subv::operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator &=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,ivector(rv)); }


	INLINE imatrix_subv &imatrix_subv::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rvector(rv)); }
	INLINE imatrix_subv &imatrix_subv::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rvector(rv)); }
	INLINE imatrix_subv &imatrix_subv::operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator |=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rvector(rv)); }
	INLINE imatrix_subv &imatrix_subv::operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE imatrix_subv &imatrix_subv::operator &=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rvector(rv)); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::imatrix::imatrix(const imatrix &rm)
	*/
	INLINE imatrix _imatrix(const imatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::imatrix::imatrix(const ivector &v)
	*/
	INLINE imatrix _imatrix(const ivector &v) throw() { return imatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::imatrix::imatrix(const ivector_slice &v)
	*/
	INLINE imatrix _imatrix(const ivector_slice &v) throw() { return imatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::imatrix::imatrix(const interval &r)
	*/
	INLINE imatrix _imatrix(const interval &r) throw() { return imatrix(r); }
	INLINE int Lb(const imatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const imatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const imatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const imatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE imatrix &SetLb(imatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE imatrix &SetUb(imatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }
	
	INLINE int RowLen ( const imatrix& A )  // Length of the rows of a interval matrix
        { return Ub(A,2)-Lb(A,2)+1; }     //----------------------------------------

        INLINE int ColLen ( const imatrix& A )  // Length of the columns of a interval matrix
        { return Ub(A,1)-Lb(A,1)+1; }     //-------------------------------------------

	INLINE int RowLen ( const imatrix_slice& A )  // Length of the rows of a interval matrix
        { return Ub(A,2)-Lb(A,2)+1; }           //----------------------------------------

        INLINE int ColLen ( const imatrix_slice& A )  // Length of the columns of a interval matrix
        { return Ub(A,1)-Lb(A,1)+1; }           //-------------------------------------------
	
	INLINE void Resize(imatrix &A) throw() { _mresize(A);}
	INLINE void Resize(imatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<imatrix,interval>(A,m,n); }
	INLINE void Resize(imatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<imatrix,interval>(A,m1,m2,n1,n2); }
	INLINE imatrix abs(const imatrix &m) throw() { return _mabs<imatrix,imatrix>(m); }
	INLINE rmatrix absmin(const imatrix &m) throw() { 
           rmatrix A(Lb(m,1),Ub(m,1),Lb(m,2),Ub(m,2));
           for(int i=Lb(m,1) ; i<=Ub(m,1) ; i++)
             for(int j=Lb(m,2) ; j<=Ub(m,2) ; j++)
               A[i][j] = AbsMin(m[i][j]);
           return A;
        }
	INLINE rmatrix absmax(const imatrix &m) throw() { 
           rmatrix A(Lb(m,1),Ub(m,1),Lb(m,2),Ub(m,2));
           for(int i=Lb(m,1) ; i<=Ub(m,1) ; i++)
             for(int j=Lb(m,2) ; j<=Ub(m,2) ; j++)
               A[i][j] = AbsMax(m[i][j]);
           return A;
        }
	INLINE imatrix abs(const imatrix_slice &ms) throw() { return _msabs<imatrix_slice,imatrix>(ms); }
	INLINE rmatrix absmin(const imatrix_slice &m) throw() { 
           rmatrix A(Lb(m,1),Ub(m,1),Lb(m,2),Ub(m,2));
           for(int i=Lb(m,1) ; i<=Ub(m,1) ; i++)
             for(int j=Lb(m,2) ; j<=Ub(m,2) ; j++)
               A[i][j] = AbsMin(m[i][j]);
           return A;
        }
	INLINE rmatrix absmax(const imatrix_slice &m) throw() { 
           rmatrix A(Lb(m,1),Ub(m,1),Lb(m,2),Ub(m,2));
           for(int i=Lb(m,1) ; i<=Ub(m,1) ; i++)
             for(int j=Lb(m,2) ; j<=Ub(m,2) ; j++)
               A[i][j] = AbsMax(m[i][j]);
           return A;
        }
	INLINE rmatrix diam(const imatrix &m) throw() { return _mdiam<imatrix,rmatrix>(m); }
	INLINE rmatrix diam(const imatrix_slice &ms) throw() { return _msdiam<imatrix_slice,rmatrix>(ms); }
	INLINE rmatrix mid(const imatrix &m) throw() { return _mmid<imatrix,rmatrix>(m); }
	INLINE rmatrix mid(const imatrix_slice &ms) throw() { return _msmid<imatrix_slice,rmatrix>(ms); }
	INLINE rmatrix Inf(const imatrix &m) throw() { return _minf<imatrix,rmatrix>(m); }
	INLINE rmatrix Sup(const imatrix &m) throw() { return _msup<imatrix,rmatrix>(m); }
	INLINE rmatrix Inf(const imatrix_slice &m) throw() { return _msinf<imatrix_slice,rmatrix>(m); }
	INLINE rmatrix Sup(const imatrix_slice &m) throw() { return _mssup<imatrix_slice,rmatrix>(m); }
	INLINE imatrix &SetInf(imatrix &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetinf<imatrix,rmatrix>(cm,rm); }
	INLINE imatrix_slice &SetInf(imatrix_slice &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetinf<imatrix_slice,rmatrix>(cm,rm); }
	INLINE imatrix &SetInf(imatrix &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetinf<imatrix,rmatrix_slice>(cm,rm); }
	INLINE imatrix_slice &SetInf(imatrix_slice &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetinf<imatrix_slice,rmatrix_slice>(cm,rm); }
	INLINE imatrix &SetSup(imatrix &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetsup<imatrix,rmatrix>(cm,rm); }
	INLINE imatrix_slice &SetSup(imatrix_slice &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetsup<imatrix_slice,rmatrix>(cm,rm); }
	INLINE imatrix &SetSup(imatrix &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetsup<imatrix,rmatrix_slice>(cm,rm); }
	INLINE imatrix_slice &SetSup(imatrix_slice &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetsup<imatrix_slice,rmatrix_slice>(cm,rm); }
	INLINE imatrix &UncheckedSetInf(imatrix &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmusetinf<imatrix,rmatrix>(cm,rm); }
	INLINE imatrix_slice &UncheckedSetInf(imatrix_slice &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmusetinf<imatrix_slice,rmatrix>(cm,rm); }
	INLINE imatrix &UncheckedSetInf(imatrix &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsusetinf<imatrix,rmatrix_slice>(cm,rm); }
	INLINE imatrix_slice &UncheckedSetInf(imatrix_slice &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsusetinf<imatrix_slice,rmatrix_slice>(cm,rm); }
	INLINE imatrix &UncheckedSetSup(imatrix &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmusetsup<imatrix,rmatrix>(cm,rm); }
	INLINE imatrix_slice &UncheckedSetSup(imatrix_slice &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmusetsup<imatrix_slice,rmatrix>(cm,rm); }
	INLINE imatrix &UncheckedSetSup(imatrix &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsusetsup<imatrix,rmatrix_slice>(cm,rm); }
	INLINE imatrix_slice &UncheckedSetSup(imatrix_slice &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsusetsup<imatrix_slice,rmatrix_slice>(cm,rm); }
	INLINE interval::interval(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_IMATRIX_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{ _smconstr(*this,m); }
//	INLINE interval interval::_interval(const imatrix &m) throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_IMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE imatrix operator *(const interval &c, const imatrix &m) throw() { return _smmult<interval,imatrix,imatrix>(c,m); }
	INLINE imatrix operator *(const interval &c, const imatrix_slice &ms) throw() { return _smsmult<interval,imatrix_slice,imatrix>(c,ms); }
	INLINE imatrix operator *(const imatrix &m,const interval &c) throw() { return _smmult<interval,imatrix,imatrix>(c,m); }
	INLINE imatrix operator *(const imatrix_slice &ms,const interval &c) throw() { return _smsmult<interval,imatrix_slice,imatrix>(c,ms); }
	INLINE imatrix &operator *=(imatrix &m,const interval &c) throw() { return _msmultassign(m,c); }
	INLINE imatrix_slice &imatrix_slice::operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE imatrix_slice &imatrix_slice::operator *=(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE imatrix_slice &imatrix_slice::operator *=(const interval &c) throw() { return _mssmultassign(*this,c); }
	INLINE imatrix operator /(const imatrix &m,const interval &c) throw() { return _msdiv<imatrix,interval,imatrix>(m,c); }
	INLINE imatrix operator /(const imatrix_slice &ms, const interval &c) throw() { return _mssdiv<imatrix_slice,interval,imatrix>(ms,c); }
	INLINE imatrix &operator /=(imatrix &m,const interval &c) throw() { return _msdivassign(m,c); }
	INLINE imatrix_slice &imatrix_slice::operator /=(const interval &c) throw() { return _mssdivassign(*this,c); }
	INLINE imatrix operator *(const real &c, const imatrix &m) throw() { return _smmult<real,imatrix,imatrix>(c,m); }
	INLINE imatrix operator *(const real &c, const imatrix_slice &ms) throw() { return _smsmult<real,imatrix_slice,imatrix>(c,ms); }
	INLINE imatrix operator *(const imatrix &m,const real &c) throw() { return _smmult<real,imatrix,imatrix>(c,m); }
	INLINE imatrix operator *(const imatrix_slice &ms,const real &c) throw() { return _smsmult<real,imatrix_slice,imatrix>(c,ms); }
	INLINE imatrix &operator *=(imatrix &m,const real &c) throw() { return _msmultassign(m,c); }
	INLINE imatrix_slice &imatrix_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE imatrix_slice &imatrix_slice::operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE imatrix_slice &imatrix_slice::operator *=(const real &c) throw() { return _mssmultassign(*this,c); }
	INLINE imatrix operator /(const imatrix &m,const real &c) throw() { return _msdiv<imatrix,real,imatrix>(m,c); }
	INLINE imatrix operator /(const imatrix_slice &ms, const real &c) throw() { return _mssdiv<imatrix_slice,real,imatrix>(ms,c); }
	INLINE imatrix &operator /=(imatrix &m,const real &c) throw() { return _msdivassign(m,c); }
	INLINE imatrix_slice &imatrix_slice::operator /=(const real &c) throw() { return _mssdivassign(*this,c); }
//	INLINE interval::interval(const rmatrix &m) throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_IMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); }
//	INLINE interval interval::_interval(const imatrix &m) throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_IMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE imatrix operator *(const interval &c, const rmatrix &m) throw() { return _smmult<interval,rmatrix,imatrix>(c,m); }
	INLINE imatrix operator *(const interval &c, const rmatrix_slice &ms) throw() { return _smsmult<interval,rmatrix_slice,imatrix>(c,ms); }
	INLINE imatrix operator *(const rmatrix &m,const interval &c) throw() { return _smmult<interval,rmatrix,imatrix>(c,m); }
	INLINE imatrix operator *(const rmatrix_slice &ms,const interval &c) throw() { return _smsmult<interval,rmatrix_slice,imatrix>(c,ms); }
	INLINE imatrix operator /(const rmatrix &m,const interval &c) throw() { return _msdiv<rmatrix,interval,imatrix>(m,c); }
	INLINE imatrix operator /(const rmatrix_slice &ms, const interval &c) throw() { return _mssdiv<rmatrix_slice,interval,imatrix>(ms,c); }
	INLINE ivector::ivector(const imatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<ivector,imatrix,interval>(*this,sl); }
	INLINE ivector::ivector(const imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<ivector,imatrix_slice,interval>(*this,sl); }
	INLINE ivector &ivector::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<ivector,imatrix,interval>(*this,m); }
	INLINE ivector &ivector::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<ivector,imatrix,interval>(*this,imatrix(m)); }
	INLINE ivector_slice & ivector_slice::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<ivector>,ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,ivector(imatrix(m))); }
	INLINE imatrix_subv &imatrix_subv::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(m)); }
	INLINE imatrix_subv &imatrix_subv::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(imatrix(m))); }
	INLINE ivector operator *(const imatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvimult<imatrix,ivector,ivector>(m,v); }
	INLINE ivector operator *(const imatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvimult<imatrix_slice,ivector,ivector>(ms,v); }
	INLINE ivector operator *(const ivector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimult<ivector,imatrix,ivector>(v,m); }
	INLINE ivector operator *(const ivector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsimult<ivector,imatrix_slice,ivector>(v,ms); }
	INLINE ivector &operator *=(ivector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimultassign<ivector,imatrix,interval>(v,m); }
	INLINE ivector &operator *=(ivector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsimultassign<ivector,imatrix_slice,interval>(v,ms); }
	INLINE ivector_slice &ivector_slice::operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmimultassign<ivector_slice,imatrix,interval>(*this,m); }
	INLINE ivector operator *(const ivector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimult<ivector,imatrix,ivector>(ivector(v),m); }
	INLINE ivector operator *(const ivector_slice &v,const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsimult<ivector,imatrix_slice,ivector>(ivector(v),m); }
	INLINE imatrix_subv &imatrix_subv::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(m)); }
	INLINE imatrix_subv &imatrix_subv::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(rmatrix(m))); }
	INLINE ivector operator *(const rvector &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimult<rvector,imatrix,ivector>(v,m); }
	INLINE ivector operator *(const rvector &v,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmsimult<rvector,imatrix_slice,ivector>(v,ms); }
	INLINE ivector operator *(const rvector_slice &v,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmimult<ivector,imatrix,ivector>(ivector(v),m); }
	INLINE ivector operator *(const imatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvimult<imatrix,rvector,ivector>(m,v); }
	INLINE ivector operator *(const imatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvimult<imatrix_slice,rvector,ivector>(ms,v); }


	INLINE const imatrix &operator +(const imatrix &m1) throw() { return m1; }
	INLINE imatrix operator +(const imatrix_slice &ms) throw() { return imatrix(ms); }
	INLINE imatrix operator +(const imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<imatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator +(const imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator +(const imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator +(const imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<imatrix_slice,imatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix &operator +=(imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE imatrix &operator +=(imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator +=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator +=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE imatrix operator -(const imatrix &m) throw() { return _mminus(m); }
	INLINE imatrix operator -(const imatrix_slice &ms) throw() { return _msminus<imatrix_slice,imatrix>(ms); }
	INLINE imatrix operator -(const imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<imatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator -(const imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator -(const imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<imatrix_slice,imatrix,imatrix>(ms,m); }
	INLINE imatrix operator -(const imatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<imatrix_slice,imatrix_slice,imatrix>(ms1,ms2); }
	INLINE imatrix &operator -=(imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE imatrix &operator -=(imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator -=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator -=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE imatrix operator *(const imatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmimult<imatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator *(const imatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsimult<imatrix,imatrix_slice,imatrix>(m1,ms); }
	INLINE imatrix operator *(const imatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmimult<imatrix_slice,imatrix,imatrix>(ms,m1); }
	INLINE imatrix operator *(const imatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsimult<imatrix_slice,imatrix_slice,imatrix>(ms1,ms2); }
	INLINE imatrix &operator *=(imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmimultassign<imatrix,imatrix,interval>(m1,m2); }
	INLINE imatrix &operator *=(imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsimultassign<imatrix,imatrix_slice,interval>(m1,ms); }
	INLINE imatrix operator |(const imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<imatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator |(const imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<imatrix_slice,imatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix &operator |=(imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE imatrix &operator |=(imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator |=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator |=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE imatrix operator &(const imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<imatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator &(const imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator &(const imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator &(const imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<imatrix_slice,imatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix &operator &=(imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE imatrix &operator &=(imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator &=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator &=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }
	INLINE imatrix operator +(const rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator +(const imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,imatrix,imatrix>(m2,m1); }
	INLINE imatrix operator +(const rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator +(const imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator +(const rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator +(const imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator +(const rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,imatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix operator +(const imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,imatrix_slice,imatrix>(m2,m1); }
	INLINE imatrix &operator +=(imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE imatrix &operator +=(imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE imatrix operator -(const rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<rmatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator -(const imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<imatrix,rmatrix,imatrix>(m1,m2); }
	INLINE imatrix operator -(const rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator -(const imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator -(const rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<rmatrix_slice,imatrix,imatrix>(ms,m); }
	INLINE imatrix operator -(const imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<imatrix_slice,rmatrix,imatrix>(ms,m); }
	INLINE imatrix operator -(const rmatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<rmatrix_slice,imatrix_slice,imatrix>(ms1,ms2); }
	INLINE imatrix operator -(const imatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<imatrix_slice,rmatrix_slice,imatrix>(ms1,ms2); }
	INLINE imatrix &operator -=(imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE imatrix &operator -=(imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE imatrix operator *(const rmatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmimult<rmatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator *(const imatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmimult<imatrix,rmatrix,imatrix>(m1,m2); }
	INLINE imatrix operator *(const rmatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsimult<rmatrix,imatrix_slice,imatrix>(m1,ms); }
	INLINE imatrix operator *(const imatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsimult<imatrix,rmatrix_slice,imatrix>(m1,ms); }
	INLINE imatrix operator *(const rmatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmimult<rmatrix_slice,imatrix,imatrix>(ms,m1); }
	INLINE imatrix operator *(const imatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmimult<imatrix_slice,rmatrix,imatrix>(ms,m1); }
	INLINE imatrix operator *(const rmatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsimult<rmatrix_slice,imatrix_slice,imatrix>(ms1,ms2); }
	INLINE imatrix operator *(const imatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsimult<imatrix_slice,rmatrix_slice,imatrix>(ms1,ms2); }
	INLINE imatrix &operator *=(imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmimultassign<imatrix,rmatrix,interval>(m1,m2); }
	INLINE imatrix &operator *=(imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsimultassign<imatrix,rmatrix_slice,interval>(m1,ms); }
	INLINE imatrix operator |(const rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,rmatrix,imatrix>(m1,m2); }
	INLINE imatrix operator |(const rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const rmatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_RMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,rmatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix operator |(const rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator |(const imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,imatrix,imatrix>(m2,m1); }
	INLINE imatrix operator |(const rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator |(const rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,imatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix operator |(const imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,imatrix_slice,imatrix>(m2,m1); }
	INLINE imatrix &operator |=(imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE imatrix &operator |=(imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator |=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator |=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE imatrix operator &(const rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<rmatrix,imatrix,imatrix>(m1,m2); }
	INLINE imatrix operator &(const imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<rmatrix,imatrix,imatrix>(m2,m1); }
	INLINE imatrix operator &(const rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator &(const imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator &(const rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,rmatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator &(const imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<rmatrix,imatrix_slice,imatrix>(m,ms); }
	INLINE imatrix operator &(const rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<rmatrix_slice,imatrix_slice,imatrix>(m1,m2); }
	INLINE imatrix operator &(const imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<rmatrix_slice,imatrix_slice,imatrix>(m2,m1); }
	INLINE imatrix &operator &=(imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE imatrix &operator &=(imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE imatrix_slice &imatrix_slice::operator &=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE imatrix_slice &imatrix_slice::operator &=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_IMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }
	INLINE bool operator ==(const imatrix &m1,const imatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const imatrix &m1,const imatrix &m2) throw() { return _mmneq(m1,m2); }
	INLINE bool operator <(const imatrix &m1,const imatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const imatrix &m1,const imatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const imatrix &m1,const imatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const imatrix &m1,const imatrix &m2) throw() { return _mmleq(m2,m1); }
	INLINE bool operator ==(const imatrix &m1,const imatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const imatrix &m1,const imatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
	INLINE bool operator <(const imatrix &m1,const imatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const imatrix &m1,const imatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const imatrix &m1,const imatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const imatrix &m1,const imatrix_slice &ms) throw() { return _msmleq(ms,m1); }
	INLINE bool operator ==(const imatrix_slice &m1,const imatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const imatrix_slice &m1,const imatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
	INLINE bool operator <(const imatrix_slice &m1,const imatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const imatrix_slice &m1,const imatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const imatrix_slice &m1,const imatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const imatrix_slice &m1,const imatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
	INLINE bool operator !(const imatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const imatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const imatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const imatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,imatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,imatrix_slice &r) throw() { return _msin(s,r); }

        //! Computes permutation of matrix according to permutation vectors, C=PAQ
        INLINE imatrix imatrix::operator()(const intvector& p, const intvector& q) {
          imatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
            for(int j=0 ; j<RowLen(A) ; j++)
              A[i+Lb(A,1)][j+Lb(A,2)] = (*this)[p[i+Lb(p)]+Lb(A,1)][q[j+Lb(q)]+Lb(A,2)];
          return A;
        }

        //! Computes permutation of matrix according to permutation vector, C=PA
        INLINE imatrix imatrix::operator()(const intvector& p) {
          imatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
              A[i+Lb(A,1)] = (*this)[p[i+Lb(p)]+Lb(A,1)];
          return A;
        }

        //! Computes permutation of matrix according to permutation matrix, C=PA
	INLINE imatrix imatrix::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

        //! Computes permutation of matrix according to permutation matrices, C=PAQ
        INLINE imatrix imatrix::operator()(const intmatrix& P, const intmatrix& Q) {
          intvector p = permvec(P);
          intvector q = perminv(permvec(Q));
          return (*this)(p,q);
        }

        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE ivector ivector::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

} // namespace cxsc

#endif

