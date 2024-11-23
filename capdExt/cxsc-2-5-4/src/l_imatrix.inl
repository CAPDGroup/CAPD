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

/* CVS $Id: l_imatrix.inl,v 1.20 2014/01/30 17:23:46 cxsc Exp $ */

#ifndef _CXSC_LIMATRIX_INL_INCLUDED
#define _CXSC_LIMATRIX_INL_INCLUDED

namespace cxsc {

INLINE l_imatrix::l_imatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE l_imatrix::l_imatrix(const l_interval &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new l_interval[1];
	*dat=r;
}

INLINE l_imatrix::l_imatrix(const real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new l_interval[1];
	*dat=r;
}

INLINE l_imatrix::l_imatrix(const l_real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new l_interval[1];
	*dat=r;
}

INLINE l_imatrix::l_imatrix(const interval &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new l_interval[1];
	*dat=r;
}

INLINE l_imatrix::l_imatrix(const rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new l_interval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE l_imatrix::l_imatrix(const l_rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new l_interval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE l_imatrix::l_imatrix(const imatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new l_interval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE l_imatrix::l_imatrix(const l_imatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new l_interval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE l_imatrix::l_imatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_LIMATRIX_WRONG_BOUNDARIES("l_imatrix::l_imatrix(const int &m, const int &n)"));
#endif
	dat=new l_interval[m*n];
}

INLINE l_imatrix::l_imatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_LIMATRIX_WRONG_BOUNDARIES("l_imatrix::l_imatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new l_interval[xsize*ysize];
}

INLINE l_ivector::l_ivector(const l_imatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new l_interval[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE l_ivector_slice & l_ivector_slice::operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>,ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
throw()
#endif
{ return _vsvassign(*this,l_ivector(m)); }


INLINE l_imatrix::l_imatrix(const l_ivector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE l_imatrix::l_imatrix(const rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE l_imatrix::l_imatrix(const l_rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE l_imatrix::l_imatrix(const ivector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE l_imatrix::l_imatrix(const l_ivector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE l_imatrix::l_imatrix(const rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE l_imatrix::l_imatrix(const ivector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE l_imatrix::l_imatrix(const l_rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new l_interval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

	INLINE l_interval &l_imatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_LIVECTOR_ELEMENT_NOT_IN_VEC("l_interval &l_imatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE l_imatrix::l_imatrix(const l_imatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new l_interval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE l_imatrix::l_imatrix(const rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new l_interval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}
	INLINE l_imatrix::l_imatrix(const l_rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new l_interval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}
	INLINE l_imatrix::l_imatrix(const imatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new l_interval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE l_imatrix_subv Row(l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE l_imatrix_subv Col(l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
		
	INLINE l_imatrix_subv Row(const l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE l_imatrix_subv Col(const l_imatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
		
	INLINE l_imatrix_subv l_imatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT("l_imatrix_subv l_imatrix::operator [](const int &i)"));
#endif
		return l_imatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE l_imatrix_subv l_imatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT("l_imatrix_subv l_imatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return l_imatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}
	
	INLINE l_imatrix_slice l_imatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG("l_imatrix_slice l_imatrix::operator ()(const int &m, const int &n)"));
#endif
		return l_imatrix_slice(*this,1,m,1,n);
	}
	
	INLINE l_imatrix_slice l_imatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG("l_imatrix_slice l_imatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return l_imatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE l_imatrix_subv l_imatrix_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT("l_imatrix_subv l_imatrix_slice::operator [](const int &i)"));
#endif
		return l_imatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE l_imatrix_subv l_imatrix_slice::operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_LIMATRIX_ROW_OR_COL_NOT_IN_MAT("l_imatrix_subv l_imatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return l_imatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}
	
	INLINE l_imatrix_slice l_imatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG("l_imatrix_slice l_imatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return l_imatrix_slice(*this,1,m,1,n);
	}
	
	INLINE l_imatrix_slice l_imatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_LIMATRIX_SUB_ARRAY_TOO_BIG("l_imatrix_slice l_imatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return l_imatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE l_imatrix_subv l_imatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG("l_imatrix_subv l_imatrix_subv::operator ()(const int &i)"));
#endif
	return l_imatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE l_imatrix_subv l_imatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_LIVECTOR_SUB_ARRAY_TOO_BIG("l_imatrix_subv l_imatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return l_imatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}



	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_imatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_interval &r) throw() { return _mvsassign(*this,r); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_ivector(v)); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const rmatrix_subv &rv) throw() { return _mvvassign(*this,rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const real &r) throw() { return _mvsassign(*this,r); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_ivector(v)); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_rmatrix_subv &rv) throw() { return _mvvassign(*this,l_rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_real &r) throw() { return _mvsassign(*this,r); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_ivector(v)); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const imatrix_subv &rv) throw() { return _mvvassign(*this,ivector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const interval &r) throw() { return _mvsassign(*this,r); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_ivector(v)); }

	INLINE l_imatrix &l_imatrix::operator =(const l_interval &r) throw() { return _msassign(*this,r); }
	INLINE l_imatrix &l_imatrix::operator =(const l_imatrix &m) throw() { return _mmassign<l_imatrix,l_imatrix,l_interval>(*this,m, l_interval()); }
	INLINE l_imatrix &l_imatrix::operator =(const l_imatrix_slice &ms) throw() { return _mmsassign<l_imatrix,l_imatrix_slice,l_interval>(*this,ms); }
	INLINE l_imatrix &l_imatrix::operator =(const l_ivector &v) throw() { return _mvassign<l_imatrix,l_ivector,l_interval>(*this,v); }
	INLINE l_imatrix &l_imatrix::operator =(const l_ivector_slice &v) throw() { return _mvassign<l_imatrix,l_ivector,l_interval>(*this,l_ivector(v)); }
	
	INLINE l_imatrix &l_imatrix::operator =(const real &r) throw() { return _msassign(*this,l_interval(r)); }
	INLINE l_imatrix &l_imatrix::operator =(const rmatrix &m) throw() { return _mmassign<l_imatrix,rmatrix,l_interval>(*this,m, l_interval()); }
	INLINE l_imatrix &l_imatrix::operator =(const rmatrix_slice &ms) throw() { return _mmsassign<l_imatrix,rmatrix_slice,l_interval>(*this,ms); }
	INLINE l_imatrix &l_imatrix::operator =(const rvector &v) throw() { return _mvassign<l_imatrix,rvector,l_interval>(*this,v); }
	INLINE l_imatrix &l_imatrix::operator =(const rvector_slice &v) throw() { return _mvassign<l_imatrix,rvector,l_interval>(*this,rvector(v)); }
	
	INLINE l_imatrix &l_imatrix::operator =(const l_real &r) throw() { return _msassign(*this,l_interval(r)); }
	INLINE l_imatrix &l_imatrix::operator =(const l_rmatrix &m) throw() { return _mmassign<l_imatrix,l_rmatrix,l_interval>(*this,m, l_interval()); }
	INLINE l_imatrix &l_imatrix::operator =(const l_rmatrix_slice &ms) throw() { return _mmsassign<l_imatrix,l_rmatrix_slice,l_interval>(*this,ms); }
	INLINE l_imatrix &l_imatrix::operator =(const l_rvector &v) throw() { return _mvassign<l_imatrix,l_rvector,l_interval>(*this,v); }
	INLINE l_imatrix &l_imatrix::operator =(const l_rvector_slice &v) throw() { return _mvassign<l_imatrix,l_rvector,l_interval>(*this,l_rvector(v)); }
	
	INLINE l_imatrix &l_imatrix::operator =(const interval &r) throw() { return _msassign(*this,l_interval(r)); }
	INLINE l_imatrix &l_imatrix::operator =(const imatrix &m) throw() { return _mmassign<l_imatrix,imatrix,l_interval>(*this,m, l_interval()); }
	INLINE l_imatrix &l_imatrix::operator =(const imatrix_slice &ms) throw() { return _mmsassign<l_imatrix,imatrix_slice,l_interval>(*this,ms); }
	INLINE l_imatrix &l_imatrix::operator =(const ivector &v) throw() { return _mvassign<l_imatrix,ivector,l_interval>(*this,v); }
	INLINE l_imatrix &l_imatrix::operator =(const ivector_slice &v) throw() { return _mvassign<l_imatrix,ivector,l_interval>(*this,ivector(v)); }
	
	INLINE l_imatrix::operator void*() throw() { return _mvoid(*this); }

	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_interval &r) throw() { return _mssassign(*this,r); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(v)); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_imatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }

	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const real &r) throw() { return _mssassign(*this,r); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(v)); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }

	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_real &r) throw() { return _mssassign(*this,r); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(v)); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const l_rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }

	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const interval &r) throw() { return _mssassign(*this,r); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(v)); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator =(const imatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,l_imatrix(l_ivector(v))); }

	INLINE l_imatrix_slice::operator void*() throw() { return _msvoid(*this); }
	INLINE l_ivector operator /(const l_imatrix_subv &rv, const l_interval &s) throw() { return _mvsdiv<l_imatrix_subv,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_imatrix_subv &rv, const l_interval &s) throw() { return _mvsmult<l_imatrix_subv,l_interval,l_ivector>(rv,s); }
	INLINE l_ivector operator *(const l_interval &s, const l_imatrix_subv &rv) throw() { return _mvsmult<l_imatrix_subv,l_interval,l_ivector>(rv,s); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator *=(const l_interval &c) throw() { return _mvsmultassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const l_interval &c) throw() { return _mvsplusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const l_interval &c) throw() { return _mvsminusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator /=(const l_interval &c) throw() { return _mvsdivassign(*this,c); }
	INLINE l_ivector abs(const l_imatrix_subv &mv) throw() { return _mvabs<l_imatrix_subv,l_ivector>(mv); }
	INLINE l_rvector diam(const l_imatrix_subv &mv) throw() { return _mvdiam<l_imatrix_subv,l_rvector>(mv); }
	INLINE l_rvector mid(const l_imatrix_subv &mv) throw() { return _mvmid<l_imatrix_subv,l_rvector>(mv); }
	INLINE l_rvector Inf(const l_imatrix_subv &mv) throw() { return _mvinf<l_imatrix_subv,l_rvector>(mv); }
	INLINE l_rvector Sup(const l_imatrix_subv &mv) throw() { return _mvsup<l_imatrix_subv,l_rvector>(mv); }

	INLINE l_imatrix_subv &SetInf(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetinf(iv,rv); }
	INLINE l_imatrix_subv &SetSup(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetsup(iv,rv); }
	INLINE l_imatrix_subv &UncheckedSetInf(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvusetinf(iv,rv); }
	INLINE l_imatrix_subv &UncheckedSetSup(l_imatrix_subv &iv,const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvusetsup(iv,rv); }
	INLINE l_imatrix_subv &SetSup(l_imatrix_subv &iv,const l_real &r) throw() { return _mvssetsup(iv,r); }
	INLINE l_imatrix_subv &SetInf(l_imatrix_subv &iv,const l_real &r) throw() { return _mvssetinf(iv,r); }
	INLINE l_imatrix_subv &UncheckedSetSup(l_imatrix_subv &iv,const l_real &r) throw() { return _mvsusetsup(iv,r); }
	INLINE l_imatrix_subv &SetUncheckedInf(l_imatrix_subv &iv,const l_real &r) throw() { return _mvsusetinf(iv,r); }
	INLINE l_ivector &l_ivector::operator =(const l_imatrix_subv &mv) throw() { return _vmvassign<l_ivector,l_imatrix_subv,l_interval>(*this,mv); }
	INLINE l_ivector_slice &l_ivector_slice::operator =(const l_imatrix_subv &mv) throw() { return _vsvassign(*this,l_ivector(mv)); }
	
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_ivector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_ivector(sl1),rv2); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_ivector(sl2),rv1); }
	
	INLINE l_interval operator *(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvlimult<l_imatrix_subv,l_imatrix_subv,l_interval>(rv1,rv2); }
	INLINE l_interval operator *(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlimult<l_ivector,l_imatrix_subv,l_interval>(rv1,rv2); }
	INLINE l_interval operator *(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlimult<l_ivector,l_imatrix_subv,l_interval>(rv2,rv1); }
	INLINE l_interval operator *(const l_ivector_slice &sl,const l_imatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlimult<l_ivector,l_imatrix_subv,l_interval>(l_ivector(sl),sv); }
	INLINE l_interval operator *(const l_imatrix_subv &mv,const l_ivector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvlimult<l_ivector,l_imatrix_subv,l_interval>(l_ivector(vs),mv); }
	INLINE l_ivector operator +(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<l_imatrix_subv,l_imatrix_subv,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_imatrix_subv,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator +(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_imatrix_subv,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator +(const l_ivector_slice &sl,const l_imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_ivector operator +(const l_imatrix_subv &mv,const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,l_ivector(rv)); }
	INLINE l_ivector operator -(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<l_imatrix_subv,l_imatrix_subv,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<l_ivector,l_imatrix_subv,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<l_imatrix_subv,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator -(const l_ivector_slice &sl,const l_imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<l_ivector,l_imatrix_subv,l_ivector>(l_ivector(sl),mv); }
	INLINE l_ivector operator -(const l_imatrix_subv &mv,const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,l_ivector(rv)); }
	INLINE l_ivector operator |(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvconv<l_imatrix_subv,l_imatrix_subv,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<l_imatrix_subv,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator |(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<l_imatrix_subv,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator |(const l_ivector_slice &sl,const l_imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_ivector operator |(const l_imatrix_subv &mv,const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,l_ivector(rv)); }
	INLINE l_ivector operator &(const l_imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvsect<l_imatrix_subv,l_imatrix_subv,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const l_imatrix_subv &rv1,const l_ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<l_imatrix_subv,l_ivector,l_ivector>(rv1,rv2); }
	INLINE l_ivector operator &(const l_ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<l_imatrix_subv,l_ivector,l_ivector>(rv2,rv1); }
	INLINE l_ivector operator &(const l_ivector_slice &sl,const l_imatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_ivector operator &(const l_imatrix_subv &mv,const l_ivector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<l_imatrix_subv,l_ivector,l_ivector>(mv,l_ivector(sl)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const l_ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const l_ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,l_ivector(rv)); }

// real
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rvector(sl2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const rmatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const rvector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rvector(sl1),rv2); }
	INLINE void accumulate(idotprecision &dp, const rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
// l_real
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rvector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_rvector(sl2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const l_rvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_rmatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rvector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,l_rvector(sl1),rv2); }
	INLINE void accumulate(idotprecision &dp, const l_rvector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
// interval
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const ivector_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,ivector(sl2),rv1); }
	INLINE void accumulate(idotprecision &dp, const l_imatrix_subv & rv1, const ivector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv2,rv1); }
	INLINE void accumulate(idotprecision &dp, const imatrix_subv & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _mvmvaccu(dp,rv1,rv2); }
	INLINE void accumulate(idotprecision &dp, const ivector_slice & sl1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,ivector(sl1),rv2); }
	INLINE void accumulate(idotprecision &dp, const ivector & rv1, const l_imatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ _vmvaccu(dp,rv1,rv2); }
	


//  matrix x matrix
	/*!
	\deprecated use standard contructors for typecasting

	\sa l_imatrix(const l_imatrix &rm)
	*/
	INLINE l_imatrix _imatrix(const l_imatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa l_imatrix(const l_ivector &v)
	*/
	INLINE l_imatrix _imatrix(const l_ivector &v) throw() { return l_imatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa l_imatrix(const l_ivector_slice &v)
	*/
	INLINE l_imatrix _imatrix(const l_ivector_slice &v) throw() { return l_imatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa l_imatrix(const l_interval &r)
	*/
	INLINE l_imatrix _imatrix(const l_interval &r) throw() { return l_imatrix(r); }
	INLINE int Lb(const l_imatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const l_imatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const l_imatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const l_imatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE l_imatrix &SetLb(l_imatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE l_imatrix &SetUb(l_imatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }
	
	
	INLINE int RowLen ( const l_imatrix& A ) // Length of the rows of a l_interval matrix
        { return Ub(A,2)-Lb(A,2)+1; }            //------------------------------------------

        INLINE int ColLen ( const l_imatrix& A ) // Length of the columns of a l_interval matrix
        { return Ub(A,1)-Lb(A,1)+1; }            //---------------------------------------------
	
	INLINE int RowLen ( const l_imatrix_slice& A ) // Length of the rows of a l_interval matrix
        { return Ub(A,2)-Lb(A,2)+1; }                  //------------------------------------------

        INLINE int ColLen ( const l_imatrix_slice& A ) // Length of the columns of a l_interval matrix
        { return Ub(A,1)-Lb(A,1)+1; }                  //---------------------------------------------
	
	INLINE void Resize(l_imatrix &A) throw() { _mresize(A);}
	INLINE void Resize(l_imatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<l_imatrix,l_interval>(A,m,n); }
	INLINE void Resize(l_imatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<l_imatrix,l_interval>(A,m1,m2,n1,n2); }
	INLINE l_imatrix abs(const l_imatrix &m) throw() { return _mabs<l_imatrix,l_imatrix>(m); }
	INLINE l_imatrix abs(const l_imatrix_slice &ms) throw() { return _msabs<l_imatrix_slice,l_imatrix>(ms); }
	INLINE l_rmatrix diam(const l_imatrix &m) throw() { return _mdiam<l_imatrix,l_rmatrix>(m); }
	INLINE l_rmatrix diam(const l_imatrix_slice &m) throw() { return _msdiam<l_imatrix_slice,l_rmatrix>(m); }
	INLINE l_rmatrix mid(const l_imatrix &m) throw() { return _mmid<l_imatrix,l_rmatrix>(m); }
	INLINE l_rmatrix mid(const l_imatrix_slice &m) throw() { return _msmid<l_imatrix_slice,l_rmatrix>(m); }
	INLINE l_rmatrix Inf(const l_imatrix &m) throw() { return _minf<l_imatrix,l_rmatrix>(m); }
	INLINE l_rmatrix Sup(const l_imatrix &m) throw() { return _msup<l_imatrix,l_rmatrix>(m); }
	INLINE l_rmatrix Inf(const l_imatrix_slice &m) throw() { return _msinf<l_imatrix_slice,l_rmatrix>(m); }
	INLINE l_rmatrix Sup(const l_imatrix_slice &m) throw() { return _mssup<l_imatrix_slice,l_rmatrix>(m); }
	INLINE l_imatrix &SetInf(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetinf<l_imatrix,l_rmatrix>(cm,rm); }
	INLINE l_imatrix_slice &SetInf(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetinf<l_imatrix_slice,l_rmatrix>(cm,rm); }
	INLINE l_imatrix &SetInf(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetinf<l_imatrix,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix_slice &SetInf(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetinf<l_imatrix_slice,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix &SetSup(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetsup<l_imatrix,l_rmatrix>(cm,rm); }
	INLINE l_imatrix_slice &SetSup(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetsup<l_imatrix_slice,l_rmatrix>(cm,rm); }
	INLINE l_imatrix &SetSup(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetsup<l_imatrix,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix_slice &SetSup(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetsup<l_imatrix_slice,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix &UncheckedSetInf(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmusetinf<l_imatrix,l_rmatrix>(cm,rm); }
	INLINE l_imatrix_slice &UncheckedSetInf(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmusetinf<l_imatrix_slice,l_rmatrix>(cm,rm); }
	INLINE l_imatrix &UncheckedSetInf(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsusetinf<l_imatrix,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix_slice &UncheckedSetInf(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsusetinf<l_imatrix_slice,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix &UncheckedSetSup(l_imatrix &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmusetsup<l_imatrix,l_rmatrix>(cm,rm); }
	INLINE l_imatrix_slice &UncheckedSetSup(l_imatrix_slice &cm,const l_rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmusetsup<l_imatrix_slice,l_rmatrix>(cm,rm); }
	INLINE l_imatrix &UncheckedSetSup(l_imatrix &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsusetsup<l_imatrix,l_rmatrix_slice>(cm,rm); }
	INLINE l_imatrix_slice &UncheckedSetSup(l_imatrix_slice &cm,const l_rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsusetsup<l_imatrix_slice,l_rmatrix_slice>(cm,rm); }
	INLINE l_interval::l_interval(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LIMATRIX_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{ _smconstr(*this,m); }
//	INLINE l_interval l_interval::_interval(const l_imatrix &m) throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LIMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }

	INLINE l_imatrix_subv &l_imatrix_subv::operator *=(const real &c) throw() { return _mvsmultassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const real &c) throw() { return _mvsplusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const real &c) throw() { return _mvsminusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator /=(const real &c) throw() { return _mvsdivassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rvector(rv)); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator *=(const l_real &c) throw() { return _mvsmultassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const l_real &c) throw() { return _mvsplusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const l_real &c) throw() { return _mvsminusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator /=(const l_real &c) throw() { return _mvsdivassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,l_rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,l_rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,l_rvector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const l_rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const l_rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,l_rvector(rv)); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator *=(const interval &c) throw() { return _mvsmultassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const interval &c) throw() { return _mvsplusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const interval &c) throw() { return _mvsminusassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator /=(const interval &c) throw() { return _mvsdivassign(*this,c); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator +=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,ivector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator -=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,ivector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator |=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,ivector(rv)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator &=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,ivector(rv)); }


	INLINE l_imatrix operator *(const l_interval &c, const l_imatrix &m) throw() { return _smmult<l_interval,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_interval &c, const l_imatrix_slice &ms) throw() { return _smsmult<l_interval,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m,const l_interval &c) throw() { return _smmult<l_interval,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const l_interval &c) throw() { return _smsmult<l_interval,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix &operator *=(l_imatrix &m,const l_interval &c) throw() { return _msmultassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const l_interval &c) throw() { return _mssmultassign(*this,c); }
	INLINE l_imatrix operator /(const l_imatrix &m,const l_interval &c) throw() { return _msdiv<l_imatrix,l_interval,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const l_interval &c) throw() { return _mssdiv<l_imatrix_slice,l_interval,l_imatrix>(ms,c); }
	INLINE l_imatrix &operator /=(l_imatrix &m,const l_interval &c) throw() { return _msdivassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator /=(const l_interval &c) throw() { return _mssdivassign(*this,c); }
	
	INLINE l_imatrix operator *(const real &c, const l_imatrix &m) throw() { return _smmult<real,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const real &c, const l_imatrix_slice &ms) throw() { return _smsmult<real,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m,const real &c) throw() { return _smmult<real,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const real &c) throw() { return _smsmult<real,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix &operator *=(l_imatrix &m,const real &c) throw() { return _msmultassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const real &c) throw() { return _mssmultassign(*this,c); }
	INLINE l_imatrix operator /(const l_imatrix &m,const real &c) throw() { return _msdiv<l_imatrix,real,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const real &c) throw() { return _mssdiv<l_imatrix_slice,real,l_imatrix>(ms,c); }
	INLINE l_imatrix &operator /=(l_imatrix &m,const real &c) throw() { return _msdivassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator /=(const real &c) throw() { return _mssdivassign(*this,c); }

	INLINE l_imatrix operator *(const l_real &c, const l_imatrix &m) throw() { return _smmult<l_real,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_real &c, const l_imatrix_slice &ms) throw() { return _smsmult<l_real,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m,const l_real &c) throw() { return _smmult<l_real,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const l_real &c) throw() { return _smsmult<l_real,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix &operator *=(l_imatrix &m,const l_real &c) throw() { return _msmultassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const l_real &c) throw() { return _mssmultassign(*this,c); }
	INLINE l_imatrix operator /(const l_imatrix &m,const l_real &c) throw() { return _msdiv<l_imatrix,l_real,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const l_real &c) throw() { return _mssdiv<l_imatrix_slice,l_real,l_imatrix>(ms,c); }
	INLINE l_imatrix &operator /=(l_imatrix &m,const l_real &c) throw() { return _msdivassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator /=(const l_real &c) throw() { return _mssdivassign(*this,c); }

	INLINE l_imatrix operator *(const interval &c, const l_imatrix &m) throw() { return _smmult<interval,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const interval &c, const l_imatrix_slice &ms) throw() { return _smsmult<interval,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m,const interval &c) throw() { return _smmult<interval,l_imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms,const interval &c) throw() { return _smsmult<interval,l_imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix &operator *=(l_imatrix &m,const interval &c) throw() { return _msmultassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator *=(const interval &c) throw() { return _mssmultassign(*this,c); }
	INLINE l_imatrix operator /(const l_imatrix &m,const interval &c) throw() { return _msdiv<l_imatrix,interval,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const l_imatrix_slice &ms, const interval &c) throw() { return _mssdiv<l_imatrix_slice,interval,l_imatrix>(ms,c); }
	INLINE l_imatrix &operator /=(l_imatrix &m,const interval &c) throw() { return _msdivassign(m,c); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator /=(const interval &c) throw() { return _mssdivassign(*this,c); }

//	INLINE l_interval::l_interval(const rmatrix &m) throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LIMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); }
//	INLINE l_interval l_interval::_interval(const l_imatrix &m) throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_LIMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }

	INLINE l_imatrix operator *(const l_interval &c, const rmatrix &m) throw() { return _smmult<l_interval,rmatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_interval &c, const rmatrix_slice &ms) throw() { return _smsmult<l_interval,rmatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const rmatrix &m,const l_interval &c) throw() { return _smmult<l_interval,rmatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const rmatrix_slice &ms,const l_interval &c) throw() { return _smsmult<l_interval,rmatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator /(const rmatrix &m,const l_interval &c) throw() { return _msdiv<rmatrix,l_interval,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const rmatrix_slice &ms, const l_interval &c) throw() { return _mssdiv<rmatrix_slice,l_interval,l_imatrix>(ms,c); }

	INLINE l_imatrix operator *(const l_interval &c, const l_rmatrix &m) throw() { return _smmult<l_interval,l_rmatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_interval &c, const l_rmatrix_slice &ms) throw() { return _smsmult<l_interval,l_rmatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const l_rmatrix &m,const l_interval &c) throw() { return _smmult<l_interval,l_rmatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms,const l_interval &c) throw() { return _smsmult<l_interval,l_rmatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator /(const l_rmatrix &m,const l_interval &c) throw() { return _msdiv<l_rmatrix,l_interval,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const l_rmatrix_slice &ms, const l_interval &c) throw() { return _mssdiv<l_rmatrix_slice,l_interval,l_imatrix>(ms,c); }

	INLINE l_imatrix operator *(const l_interval &c, const imatrix &m) throw() { return _smmult<l_interval,imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const l_interval &c, const imatrix_slice &ms) throw() { return _smsmult<l_interval,imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator *(const imatrix &m,const l_interval &c) throw() { return _smmult<l_interval,imatrix,l_imatrix>(c,m); }
	INLINE l_imatrix operator *(const imatrix_slice &ms,const l_interval &c) throw() { return _smsmult<l_interval,imatrix_slice,l_imatrix>(c,ms); }
	INLINE l_imatrix operator /(const imatrix &m,const l_interval &c) throw() { return _msdiv<imatrix,l_interval,l_imatrix>(m,c); }
	INLINE l_imatrix operator /(const imatrix_slice &ms, const l_interval &c) throw() { return _mssdiv<imatrix_slice,l_interval,l_imatrix>(ms,c); }

	INLINE l_ivector::l_ivector(const l_imatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<l_ivector,l_imatrix,l_interval>(*this,sl); }
	INLINE l_ivector::l_ivector(const l_imatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<l_ivector,l_imatrix_slice,l_interval>(*this,sl); }
	INLINE l_ivector &l_ivector::operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_ivector,l_imatrix,l_interval>(*this,m); }
	INLINE l_ivector &l_ivector::operator =(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<l_ivector,l_imatrix,l_interval>(*this,l_imatrix(m)); }
	INLINE l_ivector_slice & l_ivector_slice::operator =(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_ivector>,ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,l_ivector(l_imatrix(m))); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_ivector(m)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_ivector(l_imatrix(m))); }
	INLINE l_ivector operator *(const l_imatrix &m,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlimult<l_imatrix,l_ivector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const l_ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlimult<l_imatrix_slice,l_ivector,l_ivector>(ms,v); }
	INLINE l_ivector operator *(const l_ivector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,l_imatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const l_ivector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<l_ivector,l_imatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector &operator *=(l_ivector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimultassign<l_ivector,l_imatrix,l_interval>(v,m); }
	INLINE l_ivector &operator *=(l_ivector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimultassign<l_ivector,l_imatrix_slice,l_interval>(v,ms); }
	INLINE l_ivector_slice &l_ivector_slice::operator *=(const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmlimultassign<l_ivector_slice,l_imatrix,l_interval>(*this,m); }
	INLINE l_ivector operator *(const l_ivector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,l_imatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector operator *(const l_ivector_slice &v,const l_imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<l_ivector,l_imatrix_slice,l_ivector>(l_ivector(v),m); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(m)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(rmatrix(m))); }
	INLINE l_ivector operator *(const rvector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<rvector,l_imatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const rvector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<rvector,l_imatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector operator *(const rvector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,l_imatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector operator *(const l_imatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlimult<l_imatrix,rvector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlimult<l_imatrix_slice,rvector,l_ivector>(ms,v); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_rvector(m)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const l_rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,l_rvector(l_rmatrix(m))); }
	INLINE l_ivector operator *(const l_rvector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_rvector,l_imatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const l_rvector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<l_rvector,l_imatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector operator *(const l_rvector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,l_imatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector operator *(const l_imatrix &m,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlimult<l_imatrix,l_rvector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const l_rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlimult<l_imatrix_slice,l_rvector,l_ivector>(ms,v); }

	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(m)); }
	INLINE l_imatrix_subv &l_imatrix_subv::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(imatrix(m))); }
	INLINE l_ivector operator *(const ivector &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<ivector,l_imatrix,l_ivector>(v,m); }
	INLINE l_ivector operator *(const ivector &v,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmslimult<ivector,l_imatrix_slice,l_ivector>(v,ms); }
	INLINE l_ivector operator *(const ivector_slice &v,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmlimult<l_ivector,l_imatrix,l_ivector>(l_ivector(v),m); }
	INLINE l_ivector operator *(const l_imatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvlimult<l_imatrix,ivector,l_ivector>(m,v); }
	INLINE l_ivector operator *(const l_imatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvlimult<l_imatrix_slice,ivector,l_ivector>(ms,v); }

	INLINE const l_imatrix &operator +(const l_imatrix &m) throw() { return m; }
	INLINE l_imatrix operator +(const l_imatrix_slice &m) throw() { return l_imatrix(m); }
	INLINE l_imatrix operator +(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<l_imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<l_imatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE l_imatrix operator -(const l_imatrix &m) throw() { return _mminus(m); }
	INLINE l_imatrix operator -(const l_imatrix_slice &m) throw() { return _msminus<l_imatrix_slice,l_imatrix>(m); }
	INLINE l_imatrix operator -(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_imatrix_slice,l_imatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_imatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<l_imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<l_imatrix,l_imatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<l_imatrix_slice,l_imatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<l_imatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimultassign<l_imatrix,l_imatrix,l_interval>(m1,m2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimultassign<l_imatrix,l_imatrix_slice,l_interval>(m1,ms); }
	INLINE l_imatrix operator |(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<l_imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<l_imatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE l_imatrix operator &(const l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<l_imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<l_imatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE l_imatrix operator +(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator +(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE l_imatrix operator -(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_imatrix,rmatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<rmatrix_slice,l_imatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_imatrix_slice,rmatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const rmatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<rmatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_imatrix_slice,rmatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE l_imatrix operator *(const rmatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<l_imatrix,rmatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const rmatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<rmatrix,l_imatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<l_imatrix,rmatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const rmatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<rmatrix_slice,l_imatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<l_imatrix_slice,rmatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const rmatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<rmatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<l_imatrix_slice,rmatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimultassign<l_imatrix,rmatrix,l_interval>(m1,m2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimultassign<l_imatrix,rmatrix_slice,l_interval>(m1,ms); }
	INLINE l_imatrix operator |(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator |(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE l_imatrix operator &(const rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<rmatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator &(const rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<rmatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<rmatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE l_imatrix operator +(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<l_rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<l_rmatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator +(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE l_imatrix operator -(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_imatrix,l_rmatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_rmatrix_slice,l_imatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_imatrix_slice,l_rmatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_imatrix_slice,l_rmatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE l_imatrix operator *(const l_rmatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<l_rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<l_imatrix,l_rmatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const l_rmatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<l_rmatrix,l_imatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<l_imatrix,l_rmatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<l_rmatrix_slice,l_imatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<l_imatrix_slice,l_rmatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<l_imatrix_slice,l_rmatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimultassign<l_imatrix,l_rmatrix,l_interval>(m1,m2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimultassign<l_imatrix,l_rmatrix_slice,l_interval>(m1,ms); }
	INLINE l_imatrix operator |(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<l_rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<l_rmatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator |(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE l_imatrix operator &(const l_rmatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<l_rmatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<l_rmatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator &(const l_rmatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_rmatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_rmatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_rmatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<l_rmatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE l_imatrix operator +(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<imatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator +(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<imatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const l_imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<imatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE l_imatrix &operator +=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator +=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE l_imatrix operator -(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_imatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<imatrix_slice,l_imatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_imatrix_slice,imatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const imatrix_slice &ms1,const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<imatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator -(const l_imatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_imatrix_slice,imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE l_imatrix &operator -=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator -=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE l_imatrix operator *(const imatrix &m1, const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<l_imatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const imatrix &m1, const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<imatrix,l_imatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const l_imatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<l_imatrix,imatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const imatrix_slice &ms, const l_imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<imatrix_slice,l_imatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<l_imatrix_slice,imatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const imatrix_slice &ms1, const l_imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<imatrix_slice,l_imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator *(const l_imatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<l_imatrix_slice,imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimultassign<l_imatrix,imatrix,l_interval>(m1,m2); }
	INLINE l_imatrix &operator *=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimultassign<l_imatrix,imatrix_slice,l_interval>(m1,ms); }
	INLINE l_imatrix operator |(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<imatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator |(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<imatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const l_imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<imatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE l_imatrix &operator |=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator |=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE l_imatrix operator &(const imatrix &m1,const l_imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<imatrix,l_imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<imatrix,l_imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator &(const imatrix &m,const l_imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const imatrix_slice &ms,const l_imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_imatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_imatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,l_imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const imatrix_slice &m1,const l_imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<imatrix_slice,l_imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const l_imatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<imatrix_slice,l_imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE l_imatrix &operator &=(l_imatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE l_imatrix_slice &l_imatrix_slice::operator &=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE l_imatrix operator +(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<l_rmatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<l_rmatrix,imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator +(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator +(const l_rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<l_rmatrix_slice,imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator +(const imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<l_rmatrix_slice,imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator -(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<l_rmatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<imatrix,l_rmatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator -(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<l_rmatrix_slice,imatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<imatrix_slice,l_rmatrix,l_imatrix>(ms,m); }
	INLINE l_imatrix operator -(const l_rmatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<l_rmatrix_slice,imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator -(const imatrix_slice &ms1,const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<imatrix_slice,l_rmatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator *(const l_rmatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<l_rmatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const imatrix &m1, const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmlimult<imatrix,l_rmatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator *(const l_rmatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<l_rmatrix,imatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const imatrix &m1, const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmslimult<imatrix,l_rmatrix_slice,l_imatrix>(m1,ms); }
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<l_rmatrix_slice,imatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const imatrix_slice &ms, const l_rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmlimult<imatrix_slice,l_rmatrix,l_imatrix>(ms,m1); }
	INLINE l_imatrix operator *(const l_rmatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<l_rmatrix_slice,imatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator *(const imatrix_slice &ms1, const l_rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmslimult<imatrix_slice,l_rmatrix_slice,l_imatrix>(ms1,ms2); }
	INLINE l_imatrix operator |(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<l_rmatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<l_rmatrix,imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator |(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<l_rmatrix_slice,imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator |(const imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<l_rmatrix_slice,imatrix_slice,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator &(const l_rmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<l_rmatrix,imatrix,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const imatrix &m1,const l_rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<l_rmatrix,imatrix,l_imatrix>(m2,m1); }
	INLINE l_imatrix operator &(const l_rmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const imatrix &m,const l_rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_rmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,l_rmatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const imatrix_slice &ms,const l_rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<l_rmatrix,imatrix_slice,l_imatrix>(m,ms); }
	INLINE l_imatrix operator &(const l_rmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<l_rmatrix_slice,imatrix_slice,l_imatrix>(m1,m2); }
	INLINE l_imatrix operator &(const imatrix_slice &m1,const l_rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_LIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<l_rmatrix_slice,imatrix_slice,l_imatrix>(m2,m1); }

//------------- real x l_real ------------------------
	INLINE l_imatrix operator |(const rmatrix &rv1, const l_rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,l_rmatrix,l_imatrix>(rv1,rv2); }
	INLINE l_imatrix operator |(const l_rmatrix &rv1, const rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,l_rmatrix,l_imatrix>(rv2,rv1); }
	INLINE l_imatrix operator |(const l_rmatrix &rv, const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,rmatrix_slice,l_imatrix>(rv,sl); }
	INLINE l_imatrix operator |(const rmatrix_slice &sl,const l_rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,rmatrix_slice,l_imatrix>(rv,sl); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl, const rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,l_rmatrix_slice,l_imatrix>(rv,sl); }
	INLINE l_imatrix operator |(const rmatrix &rv,const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,l_rmatrix_slice,l_imatrix>(rv,sl); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl1, const rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,l_rmatrix_slice,l_imatrix>(sl2,sl1); }
	INLINE l_imatrix operator |(const rmatrix_slice &sl1, const l_rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,l_rmatrix_slice,l_imatrix>(sl1,sl2); }
	
//------------- l_real x l_real ------------------------
	INLINE l_imatrix operator |(const l_rmatrix &rv1, const l_rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmconv<l_rmatrix,l_rmatrix,l_imatrix>(rv1,rv2); }
	INLINE l_imatrix operator |(const l_rmatrix &rv, const l_rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,l_rmatrix_slice,l_imatrix>(rv,sl); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl,const l_rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _mmsconv<l_rmatrix,l_rmatrix_slice,l_imatrix>(rv,sl); }
	INLINE l_imatrix operator |(const l_rmatrix_slice &sl1, const l_rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<l_imatrix>)
#else
	throw()
#endif
	{ return _msmsconv<l_rmatrix_slice,l_rmatrix_slice,l_imatrix>(sl1,sl2); }
	
	INLINE bool operator ==(const l_imatrix &m1,const l_imatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const l_imatrix &m1,const l_imatrix &m2) throw() { return _mmneq(m1,m2); }
	INLINE bool operator <(const l_imatrix &m1,const l_imatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const l_imatrix &m1,const l_imatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const l_imatrix &m1,const l_imatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const l_imatrix &m1,const l_imatrix &m2) throw() { return _mmleq(m2,m1); }
	INLINE bool operator ==(const l_imatrix &m1,const l_imatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const l_imatrix &m1,const l_imatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
	INLINE bool operator <(const l_imatrix &m1,const l_imatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const l_imatrix &m1,const l_imatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const l_imatrix &m1,const l_imatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const l_imatrix &m1,const l_imatrix_slice &ms) throw() { return _msmleq(ms,m1); }
	INLINE bool operator ==(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
	INLINE bool operator <(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const l_imatrix_slice &m1,const l_imatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
	INLINE bool operator !(const l_imatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const l_imatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const l_imatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const l_imatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,l_imatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,l_imatrix_slice &r) throw() { return _msin(s,r); }

} // namespace cxsc

#endif

