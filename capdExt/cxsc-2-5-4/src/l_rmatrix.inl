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

/* CVS $Id: l_rmatrix.inl,v 1.21 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_LRMATRIX_INL_INCLUDED
#define _CXSC_LRMATRIX_INL_INCLUDED

namespace cxsc {

INLINE l_rmatrix::l_rmatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE l_rmatrix::l_rmatrix(const l_real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new l_real[1];
	*dat=r;
}

INLINE l_rmatrix::l_rmatrix(const real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new l_real[1];
	*dat=r;
}

INLINE l_rmatrix::l_rmatrix(const l_rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new l_real[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE l_rmatrix::l_rmatrix(const rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new l_real[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE l_rmatrix::l_rmatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_LRMATRIX_WRONG_BOUNDARIES("l_rmatrix::l_rmatrix(const int &m, const int &n)"));
#endif
	dat=new l_real[m*n];
}

INLINE l_rmatrix::l_rmatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_LRMATRIX_WRONG_BOUNDARIES("l_rmatrix::l_rmatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new l_real[xsize*ysize];
}

INLINE l_rvector::l_rvector(const l_rmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new l_real[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE l_rvector_slice & l_rvector_slice::operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>,ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
throw()
#endif
{ return _vsvassign(*this,l_rvector(m)); }

INLINE l_rvector _l_rvector(const rmatrix_subv &rs) throw() { return l_rvector(rs); }

INLINE l_rmatrix::l_rmatrix(const l_rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_real[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE l_rmatrix::l_rmatrix(const rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_real[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE l_rmatrix::l_rmatrix(const l_rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_real[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE l_rmatrix::l_rmatrix(const rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_real[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

	INLINE l_real &l_rmatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_LRVECTOR_ELEMENT_NOT_IN_VEC("l_real &l_rmatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE l_rmatrix::l_rmatrix(const l_rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new l_real[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE l_rmatrix::l_rmatrix(const rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new l_real[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE l_rmatrix_subv Row(l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE l_rmatrix_subv Col(l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
		
	INLINE l_rmatrix_subv Row(const l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE l_rmatrix_subv Col(const l_rmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
		
	INLINE l_rmatrix_subv l_rmatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT("l_rmatrix_subv l_rmatrix::operator [](const int &i)"));
#endif
		return l_rmatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE l_rmatrix_subv l_rmatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT("l_rmatrix_subv l_rmatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return l_rmatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}
	
	INLINE l_rmatrix_slice l_rmatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG("l_rmatrix_slice l_rmatrix::operator ()(const int &m, const int &n)"));
#endif
		return l_rmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE l_rmatrix_slice l_rmatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG("l_rmatrix_slice l_rmatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return l_rmatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE l_rmatrix_subv l_rmatrix_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT("l_rmatrix_subv l_rmatrix_slice::operator [](const int &i)"));
#endif
		return l_rmatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE l_rmatrix_subv l_rmatrix_slice::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_LRMATRIX_ROW_OR_COL_NOT_IN_MAT("l_rmatrix_subv l_rmatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return l_rmatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}
	
	INLINE l_rmatrix_slice l_rmatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG("l_rmatrix_slice l_rmatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return l_rmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE l_rmatrix_slice l_rmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_LRMATRIX_SUB_ARRAY_TOO_BIG("l_rmatrix_slice l_rmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return l_rmatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE l_rmatrix_subv l_rmatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG("l_rmatrix_subv l_rmatrix_subv::operator ()(const int &i)"));
#endif
	return l_rmatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE l_rmatrix_subv l_rmatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_LRVECTOR_SUB_ARRAY_TOO_BIG("l_rmatrix_subv l_rmatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return l_rmatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}



	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const l_rmatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const l_real &r) throw() { return _mvsassign(*this,r); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_rvector(v)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const rmatrix_subv &rv) throw() { return _mvvassign(*this,rvector(rv)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const real &r) throw() { return _mvsassign(*this,r); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_rvector(v)); }
	INLINE l_rmatrix &l_rmatrix::operator =(const l_real &r) throw() { return _msassign(*this,r); }
	INLINE l_rmatrix &l_rmatrix::operator =(const l_rmatrix &m) throw() { return _mmassign<l_rmatrix,l_rmatrix,l_real>(*this,m, l_real(0)); }
	INLINE l_rmatrix &l_rmatrix::operator =(const l_rmatrix_slice &ms) throw() { return _mmsassign<l_rmatrix,l_rmatrix_slice,l_real>(*this,ms); }
	INLINE l_rmatrix &l_rmatrix::operator =(const l_rvector &v) throw() { return _mvassign<l_rmatrix,l_rvector,l_real>(*this,v); }
	INLINE l_rmatrix &l_rmatrix::operator =(const l_rvector_slice &v) throw() { return _mvassign<l_rmatrix,l_rvector,l_real>(*this,l_rvector(v)); }
	INLINE l_rmatrix &l_rmatrix::operator =(const real &r) throw() { return _msassign(*this,l_real(r)); }
	INLINE l_rmatrix &l_rmatrix::operator =(const rmatrix &m) throw() { return _mmassign<l_rmatrix,rmatrix,l_real>(*this,m, l_real(0)); }
	INLINE l_rmatrix &l_rmatrix::operator =(const rmatrix_slice &ms) throw() { return _mmsassign<l_rmatrix,rmatrix_slice,l_real>(*this,ms); }
	INLINE l_rmatrix &l_rmatrix::operator =(const rvector &v) throw() { return _mvassign<l_rmatrix,rvector,l_real>(*this,v); }
	INLINE l_rmatrix &l_rmatrix::operator =(const rvector_slice &v) throw() { return _mvassign<l_rmatrix,rvector,l_real>(*this,rvector(v)); }
	
	INLINE l_rmatrix::operator void*() throw() { return _mvoid(*this); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const l_real &r) throw() { return _mssassign(*this,r); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_rmatrix(v)); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_rmatrix(l_rvector(v))); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const real &r) throw() { return _mssassign(*this,r); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(v)); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(rvector(v))); }
	INLINE l_rmatrix_slice::operator void*() throw() { return _msvoid(*this); }
	INLINE l_rvector operator /(const l_rmatrix_subv &rv, const l_real &s) throw() { return _mvsdiv<l_rmatrix_subv,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const l_rmatrix_subv &rv, const l_real &s) throw() { return _mvsmult<l_rmatrix_subv,l_real,l_rvector>(rv,s); }
	INLINE l_rvector operator *(const l_real &s, const l_rmatrix_subv &rv) throw() { return _mvsmult<l_rmatrix_subv,l_real,l_rvector>(rv,s); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator *=(const l_real &c) throw() { return _mvsmultassign(*this,c); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator +=(const l_real &c) throw() { return _mvsplusassign(*this,c); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator -=(const l_real &c) throw() { return _mvsminusassign(*this,c); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator /=(const l_real &c) throw() { return _mvsdivassign(*this,c); }
	INLINE l_rvector abs(const l_rmatrix_subv &mv) throw() { return _mvabs<l_rmatrix_subv,l_rvector>(mv); }
	INLINE l_rvector &l_rvector::operator =(const l_rmatrix_subv &mv) throw() { return _vmvassign<l_rvector,l_rmatrix_subv,l_real>(*this,mv); }
	INLINE l_rvector_slice &l_rvector_slice::operator =(const l_rmatrix_subv &mv) throw() { return _vsvassign(*this,l_rvector(mv)); }
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(dotprecision &dp, const l_rvector_slice & sl1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_rvector(sl1),rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_rvector(sl2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_rvector(sl1),rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_rvector(sl2),rv1); }
	INLINE l_real operator *(const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvlmult<l_rmatrix_subv,l_rmatrix_subv,l_real>(rv1,rv2); }
	INLINE l_real operator *(const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlmult<l_rvector,l_rmatrix_subv,l_real>(rv1,rv2); }
	INLINE l_real operator *(const l_rmatrix_subv &rv1,const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlmult<l_rvector,l_rmatrix_subv,l_real>(rv2,rv1); }
	INLINE l_real operator *(const l_rvector_slice &sl,const l_rmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlmult<l_rvector,l_rmatrix_subv,l_real>(l_rvector(sl),sv); }
	INLINE l_real operator *(const l_rmatrix_subv &mv,const l_rvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlmult<l_rvector,l_rmatrix_subv,l_real>(l_rvector(vs),mv); }
	INLINE l_rvector operator +(const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<l_rmatrix_subv,l_rmatrix_subv,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator +(const l_rmatrix_subv &rv1,const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_rmatrix_subv,l_rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator +(const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_rmatrix_subv,l_rvector,l_rvector>(rv2,rv1); }
	INLINE l_rvector operator +(const l_rvector_slice &sl,const l_rmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_rmatrix_subv,l_rvector,l_rvector>(mv,l_rvector(sl)); }
	INLINE l_rvector operator +(const l_rmatrix_subv &mv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_rmatrix_subv,l_rvector,l_rvector>(mv,l_rvector(sl)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator +=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,l_rvector(rv)); }
	INLINE l_rvector operator -(const l_rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<l_rmatrix_subv,l_rmatrix_subv,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator -(const l_rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<l_rvector,l_rmatrix_subv,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator -(const l_rmatrix_subv &rv1,const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<l_rmatrix_subv,l_rvector,l_rvector>(rv1,rv2); }
	INLINE l_rvector operator -(const l_rvector_slice &sl,const l_rmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<l_rvector,l_rmatrix_subv,l_rvector>(l_rvector(sl),mv); }
	INLINE l_rvector operator -(const l_rmatrix_subv &mv,const l_rvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<l_rmatrix_subv,l_rvector,l_rvector>(mv,l_rvector(sl)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator -=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,l_rvector(rv)); }
//  real
	INLINE void accumulate(dotprecision &dp, const rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(dotprecision &dp, const rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(dotprecision &dp, const rvector_slice & sl1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rvector(sl1),rv2); }
	INLINE void accumulate(dotprecision &dp, const l_rmatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rvector(sl2),rv1); }

	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rvector(sl1),rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rvector(sl2),rv1); }

	INLINE l_rmatrix_subv &l_rmatrix_subv::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,l_rvector(rv)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rvector(rv)); }

// l_rmatrix x l_rmatrix	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rmatrix::l_rmatrix(const l_rmatrix &rm)
	*/
	INLINE l_rmatrix _l_rmatrix(const l_rmatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rmatrix::l_rmatrix(const l_rvector &v)
	*/
	INLINE l_rmatrix _l_rmatrix(const l_rvector &v) throw() { return l_rmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rmatrix::l_rmatrix(const l_rvector_slice &v)
	*/
	INLINE l_rmatrix _l_rmatrix(const l_rvector_slice &v) throw() { return l_rmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::l_rmatrix::l_rmatrix(const l_real &r)
	*/
	INLINE l_rmatrix _l_rmatrix(const l_real &r) throw() { return l_rmatrix(r); }
	INLINE int Lb(const l_rmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const l_rmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const l_rmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const l_rmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE l_rmatrix &SetLb(l_rmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE l_rmatrix &SetUb(l_rmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }
	
        INLINE int RowLen ( const l_rmatrix& A ) // Length of the rows of a matrix
        { return Ub(A,2)-Lb(A,2)+1; }            //-------------------------------
      
        INLINE int ColLen ( const l_rmatrix& A ) // Length of the columns of a matrix
        { return Ub(A,1)-Lb(A,1)+1; }            //----------------------------------

        INLINE int RowLen ( const l_rmatrix_slice& A ) // Length of the rows of a matrix
        { return Ub(A,2)-Lb(A,2)+1; }                  //-------------------------------
      
        INLINE int ColLen ( const l_rmatrix_slice& A ) // Length of the columns of a matrix
        { return Ub(A,1)-Lb(A,1)+1; }                  //----------------------------------
	
	INLINE void Resize(l_rmatrix &A) throw() { _mresize(A);}
	INLINE void Resize(l_rmatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<l_rmatrix,l_real>(A,m,n); }
	INLINE void Resize(l_rmatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<l_rmatrix,l_real>(A,m1,m2,n1,n2); }
	INLINE l_rmatrix abs(const l_rmatrix &m) throw() { return _mabs<l_rmatrix,l_rmatrix>(m); }
	INLINE l_rmatrix abs(const l_rmatrix_slice &ms) throw() { return _msabs<l_rmatrix_slice,l_rmatrix>(ms); }
	INLINE l_real::l_real(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{ _smconstr(*this,m); }
//	INLINE l_real l_real::_l_real(const l_rmatrix &m) throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE l_rmatrix operator *(const l_real &c, const l_rmatrix &m) throw() { return _smmult<l_real,l_rmatrix,l_rmatrix>(c,m); }
	INLINE l_rmatrix operator *(const l_real &c, const l_rmatrix_slice &ms) throw() { return _smsmult<l_real,l_rmatrix_slice,l_rmatrix>(c,ms); }
	INLINE l_rmatrix operator *(const l_rmatrix &m,const l_real &c) throw() { return _smmult<l_real,l_rmatrix,l_rmatrix>(c,m); }
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms,const l_real &c) throw() { return _smsmult<l_real,l_rmatrix_slice,l_rmatrix>(c,ms); }
	INLINE l_rmatrix &operator *=(l_rmatrix &m,const l_real &c) throw() { return _msmultassign(m,c); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator *=(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator *=(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator *=(const l_real &c) throw() { return _mssmultassign(*this,c); }
	INLINE l_rmatrix operator /(const l_rmatrix &m,const l_real &c) throw() { return _msdiv<l_rmatrix,l_real,l_rmatrix>(m,c); }
	INLINE l_rmatrix operator /(const l_rmatrix_slice &ms, const l_real &c) throw() { return _mssdiv<l_rmatrix_slice,l_real,l_rmatrix>(ms,c); }
	INLINE l_rmatrix &operator /=(l_rmatrix &m,const l_real &c) throw() { return _msdivassign(m,c); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator /=(const l_real &c) throw() { return _mssdivassign(*this,c); }
	INLINE l_rmatrix operator *(const real &c, const l_rmatrix &m) throw() { return _smmult<real,l_rmatrix,l_rmatrix>(c,m); }
	INLINE l_rmatrix operator *(const real &c, const l_rmatrix_slice &ms) throw() { return _smsmult<real,l_rmatrix_slice,l_rmatrix>(c,ms); }
	INLINE l_rmatrix operator *(const l_rmatrix &m,const real &c) throw() { return _smmult<real,l_rmatrix,l_rmatrix>(c,m); }
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms,const real &c) throw() { return _smsmult<real,l_rmatrix_slice,l_rmatrix>(c,ms); }
	INLINE l_rmatrix &operator *=(l_rmatrix &m,const real &c) throw() { return _msmultassign(m,c); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator *=(const real &c) throw() { return _mssmultassign(*this,c); }
	INLINE l_rmatrix operator /(const l_rmatrix &m,const real &c) throw() { return _msdiv<l_rmatrix,real,l_rmatrix>(m,c); }
	INLINE l_rmatrix operator /(const l_rmatrix_slice &ms, const real &c) throw() { return _mssdiv<l_rmatrix_slice,real,l_rmatrix>(ms,c); }
	INLINE l_rmatrix &operator /=(l_rmatrix &m,const real &c) throw() { return _msdivassign(m,c); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator /=(const real &c) throw() { return _mssdivassign(*this,c); }
//	INLINE l_real::l_real(const rmatrix &m) throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); }
//	INLINE l_real l_real::_l_real(const l_rmatrix &m) throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LRMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE l_rmatrix operator *(const l_real &c, const rmatrix &m) throw() { return _smmult<l_real,rmatrix,l_rmatrix>(c,m); }
	INLINE l_rmatrix operator *(const l_real &c, const rmatrix_slice &ms) throw() { return _smsmult<l_real,rmatrix_slice,l_rmatrix>(c,ms); }
	INLINE l_rmatrix operator *(const rmatrix &m,const l_real &c) throw() { return _smmult<l_real,rmatrix,l_rmatrix>(c,m); }
	INLINE l_rmatrix operator *(const rmatrix_slice &ms,const l_real &c) throw() { return _smsmult<l_real,rmatrix_slice,l_rmatrix>(c,ms); }
	INLINE l_rmatrix operator /(const rmatrix &m,const l_real &c) throw() { return _msdiv<rmatrix,l_real,l_rmatrix>(m,c); }
	INLINE l_rmatrix operator /(const rmatrix_slice &ms, const l_real &c) throw() { return _mssdiv<rmatrix_slice,l_real,l_rmatrix>(ms,c); }
	INLINE l_rvector::l_rvector(const l_rmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<l_rvector,l_rmatrix,l_real>(*this,sl); }
	INLINE l_rvector::l_rvector(const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<l_rvector,l_rmatrix_slice,l_real>(*this,sl); }
	INLINE l_rvector &l_rvector::operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_rvector,l_rmatrix,l_real>(*this,m); }
	INLINE l_rvector &l_rvector::operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_rvector,l_rmatrix,l_real>(*this,l_rmatrix(m)); }
	INLINE l_rvector_slice & l_rvector_slice::operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_rvector>,ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,l_rvector(l_rmatrix(m))); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_rvector(m)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_rvector(l_rmatrix(m))); }
	INLINE l_rvector operator *(const l_rmatrix &m,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlmult<l_rmatrix,l_rvector,l_rvector>(m,v); }
	INLINE l_rvector operator *(const l_rmatrix_slice &ms,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlmult<l_rmatrix_slice,l_rvector,l_rvector>(ms,v); }
	INLINE l_rvector operator *(const l_rvector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmult<l_rvector,l_rmatrix,l_rvector>(v,m); }
	INLINE l_rvector operator *(const l_rvector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslmult<l_rvector,l_rmatrix_slice,l_rvector>(v,ms); }
	INLINE l_rvector &operator *=(l_rvector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmultassign<l_rvector,l_rmatrix,l_real>(v,m); }
	INLINE l_rvector &operator *=(l_rvector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslmultassign<l_rvector,l_rmatrix_slice,l_real>(v,ms); }
	INLINE l_rvector_slice &l_rvector_slice::operator *=(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmlmultassign<l_rvector_slice,l_rmatrix,l_real>(*this,m); }
	INLINE l_rvector operator *(const l_rvector_slice &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmult<l_rvector,l_rmatrix,l_rvector>(l_rvector(v),m); }
	INLINE l_rvector operator *(const l_rvector_slice &v,const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslmult<l_rvector,l_rmatrix_slice,l_rvector>(l_rvector(v),m); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(m)); }
	INLINE l_rmatrix_subv &l_rmatrix_subv::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(rmatrix(m))); }
	INLINE l_rvector operator *(const rvector &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmult<rvector,l_rmatrix,l_rvector>(v,m); }
	INLINE l_rvector operator *(const rvector &v,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslmult<rvector,l_rmatrix_slice,l_rvector>(v,ms); }
	INLINE l_rvector operator *(const rvector_slice &v,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlmult<l_rvector,l_rmatrix,l_rvector>(l_rvector(v),m); }
	INLINE l_rvector operator *(const l_rmatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlmult<l_rmatrix,rvector,l_rvector>(m,v); }
	INLINE l_rvector operator *(const l_rmatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlmult<l_rmatrix_slice,rvector,l_rvector>(ms,v); }

	INLINE const l_rmatrix &operator +(const l_rmatrix &m) throw() { return m; }
	INLINE l_rmatrix operator +(const l_rmatrix_slice &m) throw() { return l_rmatrix(m); }
	INLINE l_rmatrix operator +(const l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<l_rmatrix,l_rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator +(const l_rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,l_rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator +(const l_rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,l_rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator +(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<l_rmatrix_slice,l_rmatrix_slice,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator +=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator +=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE l_rmatrix operator -(const l_rmatrix &m) throw() { return _mminus(m); }
	INLINE l_rmatrix operator -(const l_rmatrix_slice &m) throw() { return _msminus<l_rmatrix_slice,l_rmatrix>(m); }
	INLINE l_rmatrix operator -(const l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_rmatrix,l_rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator -(const l_rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_rmatrix,l_rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_rmatrix_slice,l_rmatrix,l_rmatrix>(ms,m); }
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_rmatrix_slice,l_rmatrix_slice,l_rmatrix>(ms1,ms2); }
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator -=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator -=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlmult<l_rmatrix,l_rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslmult<l_rmatrix,l_rmatrix_slice,l_rmatrix>(m1,ms); }
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlmult<l_rmatrix_slice,l_rmatrix,l_rmatrix>(ms,m1); }
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslmult<l_rmatrix_slice,l_rmatrix_slice,l_rmatrix>(ms1,ms2); }
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlmultassign<l_rmatrix,l_rmatrix,l_real>(m1,m2); }
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslmultassign<l_rmatrix,l_rmatrix_slice,l_real>(m1,ms); }
	INLINE l_rmatrix operator +(const rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,l_rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator +(const l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,l_rmatrix,l_rmatrix>(m2,m1); }
	INLINE l_rmatrix operator +(const rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,l_rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator +(const l_rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator +(const rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator +(const l_rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,l_rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator +(const rmatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,l_rmatrix_slice,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator +(const l_rmatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,l_rmatrix_slice,l_rmatrix>(m2,m1); }
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE l_rmatrix &operator +=(l_rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE l_rmatrix operator -(const rmatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<rmatrix,l_rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator -(const l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_rmatrix,rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator -(const rmatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<rmatrix,l_rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator -(const l_rmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_rmatrix,rmatrix_slice,l_rmatrix>(m,ms); }
	INLINE l_rmatrix operator -(const rmatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<rmatrix_slice,l_rmatrix,l_rmatrix>(ms,m); }
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_rmatrix_slice,rmatrix,l_rmatrix>(ms,m); }
	INLINE l_rmatrix operator -(const rmatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<rmatrix_slice,l_rmatrix_slice,l_rmatrix>(ms1,ms2); }
	INLINE l_rmatrix operator -(const l_rmatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_rmatrix_slice,rmatrix_slice,l_rmatrix>(ms1,ms2); }
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE l_rmatrix &operator -=(l_rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE l_rmatrix_slice &l_rmatrix_slice::operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE l_rmatrix operator *(const rmatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlmult<rmatrix,l_rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlmult<l_rmatrix,rmatrix,l_rmatrix>(m1,m2); }
	INLINE l_rmatrix operator *(const rmatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslmult<rmatrix,l_rmatrix_slice,l_rmatrix>(m1,ms); }
	INLINE l_rmatrix operator *(const l_rmatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslmult<l_rmatrix,rmatrix_slice,l_rmatrix>(m1,ms); }
	INLINE l_rmatrix operator *(const rmatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlmult<rmatrix_slice,l_rmatrix,l_rmatrix>(ms,m1); }
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlmult<l_rmatrix_slice,rmatrix,l_rmatrix>(ms,m1); }
	INLINE l_rmatrix operator *(const rmatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslmult<rmatrix_slice,l_rmatrix_slice,l_rmatrix>(ms1,ms2); }
	INLINE l_rmatrix operator *(const l_rmatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslmult<l_rmatrix_slice,rmatrix_slice,l_rmatrix>(ms1,ms2); }
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlmultassign<l_rmatrix,rmatrix,l_real>(m1,m2); }
	INLINE l_rmatrix &operator *=(l_rmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LRMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslmultassign<l_rmatrix,rmatrix_slice,l_real>(m1,ms); }
	INLINE bool operator ==(const l_rmatrix &m1,const l_rmatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const l_rmatrix &m1,const l_rmatrix &m2) throw() { return _mmneq(m1,m2); }
	INLINE bool operator <(const l_rmatrix &m1,const l_rmatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const l_rmatrix &m1,const l_rmatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const l_rmatrix &m1,const l_rmatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const l_rmatrix &m1,const l_rmatrix &m2) throw() { return _mmleq(m2,m1); }
	INLINE bool operator ==(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
	INLINE bool operator <(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const l_rmatrix &m1,const l_rmatrix_slice &ms) throw() { return _msmleq(ms,m1); }
	INLINE bool operator ==(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
	INLINE bool operator <(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const l_rmatrix_slice &m1,const l_rmatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
	INLINE bool operator !(const l_rmatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const l_rmatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const l_rmatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const l_rmatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,l_rmatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,l_rmatrix_slice &r) throw() { return _msin(s,r); }

} // namespace cxsc

#endif

