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

/* CVS $Id: cimatrix.inl,v 1.30 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CIMATRIX_INL_INCLUDED
#define _CXSC_CIMATRIX_INL_INCLUDED

namespace cxsc {

INLINE cimatrix::cimatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE cimatrix::cimatrix(const cinterval &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new cinterval[1];
	*dat=r;
}

INLINE cimatrix::cimatrix(const real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new cinterval[1];
	*dat=r;
}

INLINE cimatrix::cimatrix(const complex &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new cinterval[1];
	*dat=r;
}

INLINE cimatrix::cimatrix(const interval &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new cinterval[1];
	*dat=r;
}

INLINE cimatrix::cimatrix(const rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new cinterval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE cimatrix::cimatrix(const cmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new cinterval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE cimatrix::cimatrix(const imatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new cinterval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE cimatrix::cimatrix(const cimatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new cinterval[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE cimatrix::cimatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_CIMATRIX_WRONG_BOUNDARIES("cimatrix::cimatrix(const int &m, const int &n)"));
#endif
	dat=new cinterval[m*n];
}

INLINE cimatrix::cimatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_CIMATRIX_WRONG_BOUNDARIES("cimatrix::cimatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new cinterval[xsize*ysize];
}

INLINE civector::civector(const cimatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new cinterval[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE cimatrix::cimatrix(const civector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE cimatrix::cimatrix(const rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE cimatrix::cimatrix(const cvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE cimatrix::cimatrix(const ivector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE cimatrix::cimatrix(const civector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE cimatrix::cimatrix(const rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE cimatrix::cimatrix(const ivector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE cimatrix::cimatrix(const cvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new cinterval[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}


	INLINE cimatrix::cimatrix(const cimatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new cinterval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE cimatrix::cimatrix(const rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new cinterval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}
	INLINE cimatrix::cimatrix(const cmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new cinterval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}
	INLINE cimatrix::cimatrix(const imatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new cinterval[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE cimatrix_subv Row(cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE cimatrix_subv Col(cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE cimatrix_subv Row(const cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE cimatrix_subv Col(const cimatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}
		
	INLINE cinterval& cimatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC("cinterval &cimatrix_subv::operator [](const int &i) const"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE cinterval& cimatrix_subv::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_CIVECTOR_ELEMENT_NOT_IN_VEC("cinterval &cimatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE cimatrix_subv cimatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix::operator [](const int &i)"));
#endif
		return cimatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE cimatrix_subv cimatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return cimatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}

	INLINE cimatrix_subv cimatrix::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix::operator [](const int &i)"));
#endif
		return cimatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE cimatrix_subv cimatrix::operator [](const cxscmatrix_column &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return cimatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}

	
	INLINE cimatrix_slice cimatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG("cimatrix_slice cimatrix::operator ()(const int &m, const int &n)"));
#endif
		return cimatrix_slice(*this,1,m,1,n);
	}
	
	INLINE cimatrix_slice cimatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG("cimatrix_slice cimatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return cimatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE cimatrix_subv cimatrix_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix_slice::operator [](const int &i)"));
#endif
		return cimatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE cimatrix_subv cimatrix_slice::operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return cimatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}

	INLINE cimatrix_subv cimatrix_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix_slice::operator [](const int &i)"));
#endif
		return cimatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE cimatrix_subv cimatrix_slice::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_CIMATRIX_ROW_OR_COL_NOT_IN_MAT("cimatrix_subv cimatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return cimatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}

	
	INLINE cimatrix_slice cimatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG("cimatrix_slice cimatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return cimatrix_slice(*this,1,m,1,n);
	}
	
	INLINE cimatrix_slice cimatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_CIMATRIX_SUB_ARRAY_TOO_BIG("cimatrix_slice cimatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return cimatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE cimatrix_subv cimatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG("cimatrix_subv cimatrix_subv::operator ()(const int &i)"));
#endif
	return cimatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE cimatrix_subv cimatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_CIVECTOR_SUB_ARRAY_TOO_BIG("cimatrix_subv cimatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return cimatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}



	INLINE cimatrix_subv &cimatrix_subv::operator =(const cimatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const cinterval &r) throw() { return _mvsassign(*this,r); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const civector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,civector(v)); }

	INLINE cimatrix_subv &cimatrix_subv::operator =(const rmatrix_subv &rv) throw() { return _mvvassign(*this,rvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const real &r) throw() { return _mvsassign(*this,r); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,civector(v)); }

	INLINE cimatrix_subv &cimatrix_subv::operator =(const cmatrix_subv &rv) throw() { return _mvvassign(*this,cvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const complex &r) throw() { return _mvsassign(*this,r); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const cvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,civector(v)); }

	INLINE cimatrix_subv &cimatrix_subv::operator =(const imatrix_subv &rv) throw() { return _mvvassign(*this,ivector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const interval &r) throw() { return _mvsassign(*this,r); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,civector(v)); }

	INLINE cimatrix &cimatrix::operator =(const cinterval &r) throw() { return _msassign(*this,r); }
	INLINE cimatrix &cimatrix::operator =(const cimatrix &m) throw() { return _mmassign<cimatrix,cimatrix,cinterval>(*this,m,cinterval()); }
	INLINE cimatrix &cimatrix::operator =(const cimatrix_slice &ms) throw() { return _mmsassign<cimatrix,cimatrix_slice,cinterval>(*this,ms); }
	INLINE cimatrix &cimatrix::operator =(const civector &v) throw() { return _mvassign<cimatrix,civector,cinterval>(*this,v); }
	INLINE cimatrix &cimatrix::operator =(const civector_slice &v) throw() { return _mvassign<cimatrix,civector,cinterval>(*this,civector(v)); }
	
	INLINE cimatrix &cimatrix::operator =(const real &r) throw() { return _msassign(*this,cinterval(r)); }
	INLINE cimatrix &cimatrix::operator =(const rmatrix &m) throw() { return _mmassign<cimatrix,rmatrix,cinterval>(*this,m,cinterval()); }
	INLINE cimatrix &cimatrix::operator =(const rmatrix_slice &ms) throw() { return _mmsassign<cimatrix,rmatrix_slice,cinterval>(*this,ms); }
	INLINE cimatrix &cimatrix::operator =(const rvector &v) throw() { return _mvassign<cimatrix,rvector,cinterval>(*this,v); }
	INLINE cimatrix &cimatrix::operator =(const rvector_slice &v) throw() { return _mvassign<cimatrix,rvector,cinterval>(*this,rvector(v)); }
	
	INLINE cimatrix &cimatrix::operator =(const complex &r) throw() { return _msassign(*this,cinterval(r)); }
	INLINE cimatrix &cimatrix::operator =(const cmatrix &m) throw() { return _mmassign<cimatrix,cmatrix,cinterval>(*this,m,cinterval()); }
	INLINE cimatrix &cimatrix::operator =(const cmatrix_slice &ms) throw() { return _mmsassign<cimatrix,cmatrix_slice,cinterval>(*this,ms); }
	INLINE cimatrix &cimatrix::operator =(const cvector &v) throw() { return _mvassign<cimatrix,cvector,cinterval>(*this,v); }
	INLINE cimatrix &cimatrix::operator =(const cvector_slice &v) throw() { return _mvassign<cimatrix,cvector,cinterval>(*this,cvector(v)); }
	
	INLINE cimatrix &cimatrix::operator =(const interval &r) throw() { return _msassign(*this,cinterval(r)); }
	INLINE cimatrix &cimatrix::operator =(const imatrix &m) throw() { return _mmassign<cimatrix,imatrix,cinterval>(*this,m,cinterval()); }
	INLINE cimatrix &cimatrix::operator =(const imatrix_slice &ms) throw() { return _mmsassign<cimatrix,imatrix_slice,cinterval>(*this,ms); }
	INLINE cimatrix &cimatrix::operator =(const ivector &v) throw() { return _mvassign<cimatrix,ivector,cinterval>(*this,v); }
	INLINE cimatrix &cimatrix::operator =(const ivector_slice &v) throw() { return _mvassign<cimatrix,ivector,cinterval>(*this,ivector(v)); }
	
	INLINE cimatrix::operator void*() throw() { return _mvoid(*this); }

	INLINE cimatrix_slice &cimatrix_slice::operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cinterval &r) throw() { return _mssassign(*this,r); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(v)); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const civector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cimatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }

	INLINE cimatrix_slice &cimatrix_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const real &r) throw() { return _mssassign(*this,r); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(v)); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const rmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }

	INLINE cimatrix_slice &cimatrix_slice::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const complex &r) throw() { return _mssassign(*this,r); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(v)); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const cmatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }

	INLINE cimatrix_slice &cimatrix_slice::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const interval &r) throw() { return _mssassign(*this,r); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(v)); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const ivector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }
	INLINE cimatrix_slice &cimatrix_slice::operator =(const imatrix_subv &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cimatrix(civector(v))); }

	INLINE cimatrix_slice::operator void*() throw() { return _msvoid(*this); }
	INLINE civector operator /(const cimatrix_subv &rv, const cinterval &s) throw() { return _mvsdiv<cimatrix_subv,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cimatrix_subv &rv, const cinterval &s) throw() { return _mvsmult<cimatrix_subv,cinterval,civector>(rv,s); }
	INLINE civector operator *(const cinterval &s, const cimatrix_subv &rv) throw() { return _mvsmult<cimatrix_subv,cinterval,civector>(rv,s); }
	INLINE cimatrix_subv &cimatrix_subv::operator *=(const cinterval &c) throw() { return _mvsmultassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const cinterval &c) throw() { return _mvsplusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const cinterval &c) throw() { return _mvsminusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator /=(const cinterval &c) throw() { return _mvsdivassign(*this,c); }
	INLINE ivector abs(const cimatrix_subv &mv) throw() { return _mvabs<cimatrix_subv,ivector>(mv); }
	INLINE cvector diam(const cimatrix_subv &mv) throw() { return _mvdiam<cimatrix_subv,cvector>(mv); }
	INLINE cvector mid(const cimatrix_subv &mv) throw() { return _mvmid<cimatrix_subv,cvector>(mv); }
	INLINE cvector Inf(const cimatrix_subv &mv) throw() { return _mvinf<cimatrix_subv,cvector>(mv); }
	INLINE cvector Sup(const cimatrix_subv &mv) throw() { return _mvsup<cimatrix_subv,cvector>(mv); }
	INLINE ivector Im(const cimatrix_subv &mv) throw() { return _mvim<cimatrix_subv,ivector>(mv); }
	INLINE ivector Re(const cimatrix_subv &mv) throw() { return _mvre<cimatrix_subv,ivector>(mv); }

	INLINE rmatrix InfRe(const cimatrix &m) throw()
	{
		rmatrix erg(m.lb1,m.lb2,m.ub1,m.ub2);
		
		for(int i=0;i<m.xsize*m.ysize;i++)
			erg.dat[i]=InfRe(m.dat[i]);

		return erg;
	}

	INLINE rmatrix InfRe(const cimatrix_slice &sl) throw()
	{
		rmatrix erg(sl.start1,sl.start2,sl.end1,sl.end2);
		
		for(int i=0;i<sl.sysize;i++)
		{
			for(int j=0;j<sl.sxsize;j++)
			{
				erg.dat[i*sl.sxsize+j]=InfRe(sl.dat[(i+sl.offset1)*sl.mxsize+j+sl.offset2]);
			}
		}

		return erg;
	}

	INLINE rmatrix SupRe(const cimatrix &m) throw()
	{
		rmatrix erg(m.lb1,m.lb2,m.ub1,m.ub2);
		
		for(int i=0;i<m.xsize*m.ysize;i++)
			erg.dat[i]=SupRe(m.dat[i]);

		return erg;
	}

	INLINE rmatrix SupRe(const cimatrix_slice &sl) throw()
	{
		rmatrix erg(sl.start1,sl.start2,sl.end1,sl.end2);
		
		for(int i=0;i<sl.sysize;i++)
		{
			for(int j=0;j<sl.sxsize;j++)
			{
				erg.dat[i*sl.sxsize+j]=SupRe(sl.dat[(i+sl.offset1)*sl.mxsize+j+sl.offset2]);
			}
		}

		return erg;
	}

	INLINE rmatrix InfIm(const cimatrix &m) throw()
	{
		rmatrix erg(m.lb1,m.lb2,m.ub1,m.ub2);
		
		for(int i=0;i<m.xsize*m.ysize;i++)
			erg.dat[i]=InfIm(m.dat[i]);

		return erg;
	}

	INLINE rmatrix InfIm(const cimatrix_slice &sl) throw()
	{
		rmatrix erg(sl.start1,sl.start2,sl.end1,sl.end2);
		
		for(int i=0;i<sl.sysize;i++)
		{
			for(int j=0;j<sl.sxsize;j++)
			{
				erg.dat[i*sl.sxsize+j]=InfIm(sl.dat[(i+sl.offset1)*sl.mxsize+j+sl.offset2]);
			}
		}

		return erg;
	}

	INLINE rmatrix SupIm(const cimatrix &m) throw()
	{
		rmatrix erg(m.lb1,m.lb2,m.ub1,m.ub2);
		
		for(int i=0;i<m.xsize*m.ysize;i++)
			erg.dat[i]=SupIm(m.dat[i]);

		return erg;
	}

	INLINE rmatrix SupIm(const cimatrix_slice &sl) throw()
	{
		rmatrix erg(sl.start1,sl.start2,sl.end1,sl.end2);
		
		for(int i=0;i<sl.sysize;i++)
		{
			for(int j=0;j<sl.sxsize;j++)
			{
				erg.dat[i*sl.sxsize+j]=SupIm(sl.dat[(i+sl.offset1)*sl.mxsize+j+sl.offset2]);
			}
		}

		return erg;
	}

	INLINE cimatrix_subv &SetInf(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetinf(iv,rv); }
	INLINE cimatrix_subv &SetSup(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetsup(iv,rv); }
	INLINE cimatrix_subv &UncheckedSetInf(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvusetinf(iv,rv); }
	INLINE cimatrix_subv &UncheckedSetSup(cimatrix_subv &iv,const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvusetsup(iv,rv); }
	INLINE cimatrix_subv &SetIm(cimatrix_subv &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetim(iv,rv); }
	INLINE cimatrix_subv &SetRe(cimatrix_subv &iv,const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetre(iv,rv); }
	INLINE cimatrix_subv &SetSup(cimatrix_subv &iv,const complex &r) throw() { return _mvssetsup(iv,r); }
	INLINE cimatrix_subv &SetInf(cimatrix_subv &iv,const complex &r) throw() { return _mvssetinf(iv,r); }
	INLINE cimatrix_subv &UncheckedSetSup(cimatrix_subv &iv,const complex &r) throw() { return _mvsusetsup(iv,r); }
	INLINE cimatrix_subv &SetUncheckedInf(cimatrix_subv &iv,const complex &r) throw() { return _mvsusetinf(iv,r); }
	INLINE cimatrix_subv &SetRe(cimatrix_subv &iv,const interval &r) throw() { return _mvssetre(iv,r); }
	INLINE cimatrix_subv &SetIm(cimatrix_subv &iv,const interval &r) throw() { return _mvssetim(iv,r); }
	INLINE civector &civector::operator =(const cimatrix_subv &mv) throw() { return _vmvassign<civector,cimatrix_subv,cinterval>(*this,mv); }
	INLINE civector_slice &civector_slice::operator =(const cimatrix_subv &mv) throw() { return _vsvassign(*this,civector(mv)); }
	
	INLINE cinterval operator *(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvcimult<cimatrix_subv,cimatrix_subv,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcimult<civector,cimatrix_subv,cinterval>(rv1,rv2); }
	INLINE cinterval operator *(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcimult<civector,cimatrix_subv,cinterval>(rv2,rv1); }
	INLINE cinterval operator *(const civector_slice &sl,const cimatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcimult<civector,cimatrix_subv,cinterval>(civector(sl),sv); }
	INLINE cinterval operator *(const cimatrix_subv &mv,const civector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcimult<civector,cimatrix_subv,cinterval>(civector(vs),mv); }
	INLINE civector operator +(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<cimatrix_subv,cimatrix_subv,civector>(rv1,rv2); }
	INLINE civector operator +(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cimatrix_subv,civector,civector>(rv1,rv2); }
	INLINE civector operator +(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cimatrix_subv,civector,civector>(rv2,rv1); }
	INLINE civector operator +(const civector_slice &sl,const cimatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE civector operator +(const cimatrix_subv &mv,const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,civector(rv)); }
	INLINE civector operator -(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<cimatrix_subv,cimatrix_subv,civector>(rv1,rv2); }
	INLINE civector operator -(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<civector,cimatrix_subv,civector>(rv1,rv2); }
	INLINE civector operator -(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<cimatrix_subv,civector,civector>(rv1,rv2); }
	INLINE civector operator -(const civector_slice &sl,const cimatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<civector,cimatrix_subv,civector>(civector(sl),mv); }
	INLINE civector operator -(const cimatrix_subv &mv,const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,civector(rv)); }
	INLINE civector operator |(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvconv<cimatrix_subv,cimatrix_subv,civector>(rv1,rv2); }
	INLINE civector operator |(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<cimatrix_subv,civector,civector>(rv1,rv2); }
	INLINE civector operator |(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<cimatrix_subv,civector,civector>(rv2,rv1); }
	INLINE civector operator |(const civector_slice &sl,const cimatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE civector operator |(const cimatrix_subv &mv,const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconv<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,civector(rv)); }
	INLINE civector operator &(const cimatrix_subv & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvsect<cimatrix_subv,cimatrix_subv,civector>(rv1,rv2); }
	INLINE civector operator &(const cimatrix_subv &rv1,const civector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<cimatrix_subv,civector,civector>(rv1,rv2); }
	INLINE civector operator &(const civector & rv1, const cimatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<cimatrix_subv,civector,civector>(rv2,rv1); }
	INLINE civector operator &(const civector_slice &sl,const cimatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE civector operator &(const cimatrix_subv &mv,const civector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsect<cimatrix_subv,civector,civector>(mv,civector(sl)); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const civector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const civector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,civector(rv)); }

// real


//  matrix x matrix
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cimatrix::cimatrix(const cimatrix &rm)
	*/
	INLINE cimatrix _imatrix(const cimatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cimatrix::cimatrix(const civector &v)
	*/
	INLINE cimatrix _imatrix(const civector &v) throw() { return cimatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cimatrix::cimatrix(const civector_slice &v)
	*/
	INLINE cimatrix _imatrix(const civector_slice &v) throw() { return cimatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cimatrix::cimatrix(const cinterval &r)
	*/
	INLINE cimatrix _imatrix(const cinterval &r) throw() { return cimatrix(r); }
	INLINE int Lb(const cimatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const cimatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const cimatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const cimatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE cimatrix &SetLb(cimatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE cimatrix &SetUb(cimatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }
	
        INLINE int RowLen ( const cimatrix& A ) // Length of the rows of a cinterval matrix
        { return Ub(A,2)-Lb(A,2)+1; }           //------------------------------------------

        INLINE int ColLen ( const cimatrix& A ) // Length of the columns of a cinterval matrix
        { return Ub(A,1)-Lb(A,1)+1; }           //---------------------------------------------
	
        INLINE int RowLen ( const cimatrix_slice& A ) // Length of the rows of a cinterval matrix
        { return Ub(A,2)-Lb(A,2)+1; }                 //------------------------------------------

        INLINE int ColLen ( const cimatrix_slice& A ) // Length of the columns of a cinterval matrix
        { return Ub(A,1)-Lb(A,1)+1; }                 //---------------------------------------------
	
	INLINE void Resize(cimatrix &A) throw() { _mresize(A);}
	INLINE void Resize(cimatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<cimatrix,cinterval>(A,m,n); }
	INLINE void Resize(cimatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<cimatrix,cinterval>(A,m1,m2,n1,n2); }
	INLINE imatrix abs(const cimatrix &m) throw() { return _mabs<cimatrix,imatrix>(m); }
	INLINE imatrix abs(const cimatrix_slice &ms) throw() { return _msabs<cimatrix_slice,imatrix>(ms); }
	INLINE cmatrix diam(const cimatrix &m) throw() { return _mdiam<cimatrix,cmatrix>(m); }
	INLINE cmatrix diam(const cimatrix_slice &m) throw() { return _msdiam<cimatrix_slice,cmatrix>(m); }
	INLINE cmatrix mid(const cimatrix &m) throw() { return _mmid<cimatrix,cmatrix>(m); }
	INLINE cmatrix mid(const cimatrix_slice &m) throw() { return _msmid<cimatrix_slice,cmatrix>(m); }
	INLINE cmatrix Inf(const cimatrix &m) throw() { return _minf<cimatrix,cmatrix>(m); }
	INLINE cmatrix Sup(const cimatrix &m) throw() { return _msup<cimatrix,cmatrix>(m); }
	INLINE cmatrix Inf(const cimatrix_slice &m) throw() { return _msinf<cimatrix_slice,cmatrix>(m); }
	INLINE cmatrix Sup(const cimatrix_slice &m) throw() { return _mssup<cimatrix_slice,cmatrix>(m); }
	INLINE cimatrix &SetInf(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetinf<cimatrix,cmatrix>(cm,rm); }
	INLINE cimatrix_slice &SetInf(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetinf<cimatrix_slice,cmatrix>(cm,rm); }
	INLINE cimatrix &SetInf(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetinf<cimatrix,cmatrix_slice>(cm,rm); }
	INLINE cimatrix_slice &SetInf(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetinf<cimatrix_slice,cmatrix_slice>(cm,rm); }
	INLINE cimatrix &SetSup(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetsup<cimatrix,cmatrix>(cm,rm); }
	INLINE cimatrix_slice &SetSup(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetsup<cimatrix_slice,cmatrix>(cm,rm); }
	INLINE cimatrix &SetSup(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetsup<cimatrix,cmatrix_slice>(cm,rm); }
	INLINE cimatrix_slice &SetSup(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetsup<cimatrix_slice,cmatrix_slice>(cm,rm); }
	INLINE cimatrix &UncheckedSetInf(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmusetinf<cimatrix,cmatrix>(cm,rm); }
	INLINE cimatrix_slice &UncheckedSetInf(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmusetinf<cimatrix_slice,cmatrix>(cm,rm); }
	INLINE cimatrix &UncheckedSetInf(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsusetinf<cimatrix,cmatrix_slice>(cm,rm); }
	INLINE cimatrix_slice &UncheckedSetInf(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsusetinf<cimatrix_slice,cmatrix_slice>(cm,rm); }
	INLINE cimatrix &UncheckedSetSup(cimatrix &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmusetsup<cimatrix,cmatrix>(cm,rm); }
	INLINE cimatrix_slice &UncheckedSetSup(cimatrix_slice &cm,const cmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmusetsup<cimatrix_slice,cmatrix>(cm,rm); }
	INLINE cimatrix &UncheckedSetSup(cimatrix &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsusetsup<cimatrix,cmatrix_slice>(cm,rm); }
	INLINE cimatrix_slice &UncheckedSetSup(cimatrix_slice &cm,const cmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsusetsup<cimatrix_slice,cmatrix_slice>(cm,rm); }
	INLINE imatrix Im(const cimatrix &m) throw() { return _mim<cimatrix,imatrix>(m); }
	INLINE imatrix Re(const cimatrix &m) throw() { return _mre<cimatrix,imatrix>(m); }
	INLINE imatrix Im(const cimatrix_slice &m) throw() { return _msim<cimatrix_slice,imatrix>(m); }
	INLINE imatrix Re(const cimatrix_slice &m) throw() { return _msre<cimatrix_slice,imatrix>(m); }
	INLINE cimatrix &SetIm(cimatrix &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetim<cimatrix,imatrix>(cm,rm); }
	INLINE cimatrix_slice &SetIm(cimatrix_slice &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetim<cimatrix_slice,imatrix>(cm,rm); }
	INLINE cimatrix &SetIm(cimatrix &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetim<cimatrix,imatrix_slice>(cm,rm); }
	INLINE cimatrix_slice &SetIm(cimatrix_slice &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetim<cimatrix_slice,imatrix_slice>(cm,rm); }
	INLINE cimatrix &SetRe(cimatrix &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetre<cimatrix,imatrix>(cm,rm); }
	INLINE cimatrix_slice &SetRe(cimatrix_slice &cm,const imatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetre<cimatrix_slice,imatrix>(cm,rm); }
	INLINE cimatrix &SetRe(cimatrix &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetre<cimatrix,imatrix_slice>(cm,rm); }
	INLINE cimatrix_slice &SetRe(cimatrix_slice &cm,const imatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetre<cimatrix_slice,imatrix_slice>(cm,rm); }
	INLINE cinterval::cinterval(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{ _smconstr(*this,m); }
//	INLINE cinterval cinterval::_interval(const cimatrix &m) throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }

	INLINE cimatrix_subv &cimatrix_subv::operator *=(const real &c) throw() { return _mvsmultassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const real &c) throw() { return _mvsplusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const real &c) throw() { return _mvsminusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator /=(const real &c) throw() { return _mvsdivassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rvector(rv)); }

	INLINE cimatrix_subv &cimatrix_subv::operator *=(const complex &c) throw() { return _mvsmultassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const complex &c) throw() { return _mvsplusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const complex &c) throw() { return _mvsminusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator /=(const complex &c) throw() { return _mvsdivassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,cvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,cvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,cvector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,cvector(rv)); }

	INLINE cimatrix_subv &cimatrix_subv::operator *=(const interval &c) throw() { return _mvsmultassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const interval &c) throw() { return _mvsplusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const interval &c) throw() { return _mvsminusassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator /=(const interval &c) throw() { return _mvsdivassign(*this,c); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator +=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,ivector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator -=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,ivector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator |=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvconvassign(*this,ivector(rv)); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const ivector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,rv); }
	INLINE cimatrix_subv &cimatrix_subv::operator &=(const ivector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsectassign(*this,ivector(rv)); }


	INLINE cimatrix operator *(const cinterval &c, const cimatrix &m) throw() { return _smmult<cinterval,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cinterval &c, const cimatrix_slice &ms) throw() { return _smsmult<cinterval,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const cimatrix &m,const cinterval &c) throw() { return _smmult<cinterval,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cimatrix_slice &ms,const cinterval &c) throw() { return _smsmult<cinterval,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix &operator *=(cimatrix &m,const cinterval &c) throw() { return _msmultassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const cinterval &c) throw() { return _mssmultassign(*this,c); }
	INLINE cimatrix operator /(const cimatrix &m,const cinterval &c) throw() { return _msdiv<cimatrix,cinterval,cimatrix>(m,c); }
	INLINE cimatrix operator /(const cimatrix_slice &ms, const cinterval &c) throw() { return _mssdiv<cimatrix_slice,cinterval,cimatrix>(ms,c); }
	INLINE cimatrix &operator /=(cimatrix &m,const cinterval &c) throw() { return _msdivassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator /=(const cinterval &c) throw() { return _mssdivassign(*this,c); }
	
	INLINE cimatrix operator *(const real &c, const cimatrix &m) throw() { return _smmult<real,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const real &c, const cimatrix_slice &ms) throw() { return _smsmult<real,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const cimatrix &m,const real &c) throw() { return _smmult<real,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cimatrix_slice &ms,const real &c) throw() { return _smsmult<real,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix &operator *=(cimatrix &m,const real &c) throw() { return _msmultassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const real &c) throw() { return _mssmultassign(*this,c); }
	INLINE cimatrix operator /(const cimatrix &m,const real &c) throw() { return _msdiv<cimatrix,real,cimatrix>(m,c); }
	INLINE cimatrix operator /(const cimatrix_slice &ms, const real &c) throw() { return _mssdiv<cimatrix_slice,real,cimatrix>(ms,c); }
	INLINE cimatrix &operator /=(cimatrix &m,const real &c) throw() { return _msdivassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator /=(const real &c) throw() { return _mssdivassign(*this,c); }

	INLINE cimatrix operator *(const complex &c, const cimatrix &m) throw() { return _smmult<complex,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const complex &c, const cimatrix_slice &ms) throw() { return _smsmult<complex,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const cimatrix &m,const complex &c) throw() { return _smmult<complex,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cimatrix_slice &ms,const complex &c) throw() { return _smsmult<complex,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix &operator *=(cimatrix &m,const complex &c) throw() { return _msmultassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const complex &c) throw() { return _mssmultassign(*this,c); }
	INLINE cimatrix operator /(const cimatrix &m,const complex &c) throw() { return _msdiv<cimatrix,complex,cimatrix>(m,c); }
	INLINE cimatrix operator /(const cimatrix_slice &ms, const complex &c) throw() { return _mssdiv<cimatrix_slice,complex,cimatrix>(ms,c); }
	INLINE cimatrix &operator /=(cimatrix &m,const complex &c) throw() { return _msdivassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator /=(const complex &c) throw() { return _mssdivassign(*this,c); }

	INLINE cimatrix operator *(const interval &c, const cimatrix &m) throw() { return _smmult<interval,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const interval &c, const cimatrix_slice &ms) throw() { return _smsmult<interval,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const cimatrix &m,const interval &c) throw() { return _smmult<interval,cimatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cimatrix_slice &ms,const interval &c) throw() { return _smsmult<interval,cimatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix &operator *=(cimatrix &m,const interval &c) throw() { return _msmultassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cimatrix_slice &cimatrix_slice::operator *=(const interval &c) throw() { return _mssmultassign(*this,c); }
	INLINE cimatrix operator /(const cimatrix &m,const interval &c) throw() { return _msdiv<cimatrix,interval,cimatrix>(m,c); }
	INLINE cimatrix operator /(const cimatrix_slice &ms, const interval &c) throw() { return _mssdiv<cimatrix_slice,interval,cimatrix>(ms,c); }
	INLINE cimatrix &operator /=(cimatrix &m,const interval &c) throw() { return _msdivassign(m,c); }
	INLINE cimatrix_slice &cimatrix_slice::operator /=(const interval &c) throw() { return _mssdivassign(*this,c); }

//	INLINE cinterval::cinterval(const rmatrix &m) throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); }
//	INLINE cinterval cinterval::_interval(const cimatrix &m) throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CIMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }

	INLINE cimatrix operator *(const cinterval &c, const rmatrix &m) throw() { return _smmult<cinterval,rmatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cinterval &c, const rmatrix_slice &ms) throw() { return _smsmult<cinterval,rmatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const rmatrix &m,const cinterval &c) throw() { return _smmult<cinterval,rmatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const rmatrix_slice &ms,const cinterval &c) throw() { return _smsmult<cinterval,rmatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator /(const rmatrix &m,const cinterval &c) throw() { return _msdiv<rmatrix,cinterval,cimatrix>(m,c); }
	INLINE cimatrix operator /(const rmatrix_slice &ms, const cinterval &c) throw() { return _mssdiv<rmatrix_slice,cinterval,cimatrix>(ms,c); }

	INLINE cimatrix operator *(const cinterval &c, const cmatrix &m) throw() { return _smmult<cinterval,cmatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cinterval &c, const cmatrix_slice &ms) throw() { return _smsmult<cinterval,cmatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const cmatrix &m,const cinterval &c) throw() { return _smmult<cinterval,cmatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cmatrix_slice &ms,const cinterval &c) throw() { return _smsmult<cinterval,cmatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator /(const cmatrix &m,const cinterval &c) throw() { return _msdiv<cmatrix,cinterval,cimatrix>(m,c); }
	INLINE cimatrix operator /(const cmatrix_slice &ms, const cinterval &c) throw() { return _mssdiv<cmatrix_slice,cinterval,cimatrix>(ms,c); }

	INLINE cimatrix operator *(const cinterval &c, const imatrix &m) throw() { return _smmult<cinterval,imatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const cinterval &c, const imatrix_slice &ms) throw() { return _smsmult<cinterval,imatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator *(const imatrix &m,const cinterval &c) throw() { return _smmult<cinterval,imatrix,cimatrix>(c,m); }
	INLINE cimatrix operator *(const imatrix_slice &ms,const cinterval &c) throw() { return _smsmult<cinterval,imatrix_slice,cimatrix>(c,ms); }
	INLINE cimatrix operator /(const imatrix &m,const cinterval &c) throw() { return _msdiv<imatrix,cinterval,cimatrix>(m,c); }
	INLINE cimatrix operator /(const imatrix_slice &ms, const cinterval &c) throw() { return _mssdiv<imatrix_slice,cinterval,cimatrix>(ms,c); }

	INLINE civector::civector(const cimatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<civector,cimatrix,cinterval>(*this,sl); }
	INLINE civector::civector(const cimatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<civector,cimatrix_slice,cinterval>(*this,sl); }
	INLINE civector &civector::operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,cimatrix,cinterval>(*this,m); }
	INLINE civector &civector::operator =(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<civector,cimatrix,cinterval>(*this,cimatrix(m)); }
	INLINE civector_slice & civector_slice::operator =(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<civector>,ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,civector(cimatrix(m))); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,civector(m)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,civector(cimatrix(m))); }
	INLINE civector operator *(const cimatrix &m,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<cimatrix,civector,civector>(m,v); }
	INLINE civector operator *(const cimatrix_slice &ms,const civector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<cimatrix_slice,civector,civector>(ms,v); }
	INLINE civector operator *(const civector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cimatrix,civector>(v,m); }
	INLINE civector operator *(const civector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<civector,cimatrix_slice,civector>(v,ms); }
	INLINE civector &operator *=(civector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimultassign<civector,cimatrix,cinterval>(v,m); }
	INLINE civector &operator *=(civector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimultassign<civector,cimatrix_slice,cinterval>(v,ms); }
	INLINE civector_slice &civector_slice::operator *=(const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmcimultassign<civector_slice,cimatrix,cinterval>(*this,m); }
	INLINE civector operator *(const civector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cimatrix,civector>(civector(v),m); }
	INLINE civector operator *(const civector_slice &v,const cimatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<civector,cimatrix_slice,civector>(civector(v),m); }

	INLINE cimatrix_subv &cimatrix_subv::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(m)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(rmatrix(m))); }
	INLINE civector operator *(const rvector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<rvector,cimatrix,civector>(v,m); }
	INLINE civector operator *(const rvector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<rvector,cimatrix_slice,civector>(v,ms); }
	INLINE civector operator *(const rvector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cimatrix,civector>(civector(v),m); }
	INLINE civector operator *(const cimatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<cimatrix,rvector,civector>(m,v); }
	INLINE civector operator *(const cimatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<cimatrix_slice,rvector,civector>(ms,v); }

	INLINE cimatrix_subv &cimatrix_subv::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,cvector(m)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,cvector(cmatrix(m))); }
	INLINE civector operator *(const cvector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<cvector,cimatrix,civector>(v,m); }
	INLINE civector operator *(const cvector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<cvector,cimatrix_slice,civector>(v,ms); }
	INLINE civector operator *(const cvector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cimatrix,civector>(civector(v),m); }
	INLINE civector operator *(const cimatrix &m,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<cimatrix,cvector,civector>(m,v); }
	INLINE civector operator *(const cimatrix_slice &ms,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<cimatrix_slice,cvector,civector>(ms,v); }

	INLINE cimatrix_subv &cimatrix_subv::operator =(const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(m)); }
	INLINE cimatrix_subv &cimatrix_subv::operator =(const imatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,ivector(imatrix(m))); }
	INLINE civector operator *(const ivector &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<ivector,cimatrix,civector>(v,m); }
	INLINE civector operator *(const ivector &v,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscimult<ivector,cimatrix_slice,civector>(v,ms); }
	INLINE civector operator *(const ivector_slice &v,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcimult<civector,cimatrix,civector>(civector(v),m); }
	INLINE civector operator *(const cimatrix &m,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcimult<cimatrix,ivector,civector>(m,v); }
	INLINE civector operator *(const cimatrix_slice &ms,const ivector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcimult<cimatrix_slice,ivector,civector>(ms,v); }

	INLINE const cimatrix &operator +(const cimatrix &m) throw() { return m; }
	INLINE cimatrix operator +(const cimatrix_slice &m) throw() { return cimatrix(m); }
	INLINE cimatrix operator +(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<cimatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<cimatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix &operator +=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE cimatrix &operator +=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE cimatrix operator -(const cimatrix &m) throw() { return _mminus(m); }
	INLINE cimatrix operator -(const cimatrix_slice &m) throw() { return _msminus<cimatrix_slice,cimatrix>(m); }
	INLINE cimatrix operator -(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cimatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cimatrix_slice,cimatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cimatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE cimatrix operator *(const cimatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<cimatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const cimatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<cimatrix,cimatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const cimatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<cimatrix_slice,cimatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<cimatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimultassign<cimatrix,cimatrix,cinterval>(m1,m2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimultassign<cimatrix,cimatrix_slice,cinterval>(m1,ms); }
	INLINE cimatrix operator |(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<cimatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<cimatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix &operator |=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE cimatrix &operator |=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE cimatrix operator &(const cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<cimatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<cimatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix &operator &=(cimatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE cimatrix &operator &=(cimatrix &m1,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE cimatrix operator +(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator +(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const rmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator +=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE cimatrix &operator +=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE cimatrix operator -(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<rmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cimatrix,rmatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<rmatrix_slice,cimatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cimatrix_slice,rmatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const rmatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<rmatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cimatrix_slice,rmatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE cimatrix operator *(const rmatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<rmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const cimatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<cimatrix,rmatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const rmatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<rmatrix,cimatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const cimatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<cimatrix,rmatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const rmatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<rmatrix_slice,cimatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const cimatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<cimatrix_slice,rmatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const rmatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<rmatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<cimatrix_slice,rmatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimultassign<cimatrix,rmatrix,cinterval>(m1,m2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimultassign<cimatrix,rmatrix_slice,cinterval>(m1,ms); }
	INLINE cimatrix operator |(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator |(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const rmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator |=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE cimatrix &operator |=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE cimatrix operator &(const rmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<rmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<rmatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator &(const rmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const rmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,rmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<rmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const rmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<rmatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<rmatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator &=(cimatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE cimatrix &operator &=(cimatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE cimatrix operator +(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<cmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<cmatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator +(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<cmatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<cmatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator +=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE cimatrix &operator +=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE cimatrix operator -(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cimatrix,cmatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cmatrix_slice,cimatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cimatrix_slice,cmatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const cmatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cmatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cimatrix_slice,cmatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE cimatrix operator *(const cmatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<cmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const cimatrix &m1, const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<cimatrix,cmatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const cmatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<cmatrix,cimatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const cimatrix &m1, const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<cimatrix,cmatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const cmatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<cmatrix_slice,cimatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const cimatrix_slice &ms, const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<cimatrix_slice,cmatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const cmatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<cmatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<cimatrix_slice,cmatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimultassign<cimatrix,cmatrix,cinterval>(m1,m2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimultassign<cimatrix,cmatrix_slice,cinterval>(m1,ms); }
	INLINE cimatrix operator |(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<cmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<cmatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator |(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<cmatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<cmatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator |=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE cimatrix &operator |=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE cimatrix operator &(const cmatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<cmatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<cmatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator &(const cmatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cmatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cmatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cmatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<cmatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<cmatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator &=(cimatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE cimatrix &operator &=(cimatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE cimatrix operator +(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<imatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<imatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator +(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const imatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<imatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const cimatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<imatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator +=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE cimatrix &operator +=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator +=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE cimatrix operator -(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<imatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cimatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<imatrix_slice,cimatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cimatrix_slice,imatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const imatrix_slice &ms1,const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<imatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator -(const cimatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cimatrix_slice,imatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE cimatrix &operator -=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator -=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE cimatrix operator *(const imatrix &m1, const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<imatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const cimatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<cimatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const imatrix &m1, const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<imatrix,cimatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const cimatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<cimatrix,imatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const imatrix_slice &ms, const cimatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<imatrix_slice,cimatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const cimatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<cimatrix_slice,imatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const imatrix_slice &ms1, const cimatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<imatrix_slice,cimatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator *(const cimatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<cimatrix_slice,imatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimultassign<cimatrix,imatrix,cinterval>(m1,m2); }
	INLINE cimatrix &operator *=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimultassign<cimatrix,imatrix_slice,cinterval>(m1,ms); }
	INLINE cimatrix operator |(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<imatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<imatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator |(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const imatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<imatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const cimatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<imatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator |=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconvassign(m1,m2); }
	INLINE cimatrix &operator |=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconvassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmconvassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator |=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconvassign(*this,ms2); }
	INLINE cimatrix operator &(const imatrix &m1,const cimatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<imatrix,cimatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<imatrix,cimatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator &(const imatrix &m,const cimatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const imatrix_slice &ms,const cimatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cimatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cimatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,cimatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const imatrix_slice &m1,const cimatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<imatrix_slice,cimatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const cimatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<imatrix_slice,cimatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix &operator &=(cimatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsectassign(m1,m2); }
	INLINE cimatrix &operator &=(cimatrix &m1,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssectassign(m1,ms); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsectassign(*this,m1); }
	INLINE cimatrix_slice &cimatrix_slice::operator &=(const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssectassign(*this,ms2); }

	INLINE cimatrix operator +(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<cmatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<cmatrix,imatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator +(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator +(const cmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<cmatrix_slice,imatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator +(const imatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<cmatrix_slice,imatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix operator -(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cmatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<imatrix,cmatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator -(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator -(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cmatrix_slice,imatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<imatrix_slice,cmatrix,cimatrix>(ms,m); }
	INLINE cimatrix operator -(const cmatrix_slice &ms1,const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cmatrix_slice,imatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator -(const imatrix_slice &ms1,const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<imatrix_slice,cmatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator *(const cmatrix &m1, const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<cmatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const imatrix &m1, const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcimult<imatrix,cmatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator *(const cmatrix &m1, const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<cmatrix,imatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const imatrix &m1, const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscimult<imatrix,cmatrix_slice,cimatrix>(m1,ms); }
	INLINE cimatrix operator *(const cmatrix_slice &ms, const imatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<cmatrix_slice,imatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const imatrix_slice &ms, const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcimult<imatrix_slice,cmatrix,cimatrix>(ms,m1); }
	INLINE cimatrix operator *(const cmatrix_slice &ms1, const imatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<cmatrix_slice,imatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator *(const imatrix_slice &ms1, const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscimult<imatrix_slice,cmatrix_slice,cimatrix>(ms1,ms2); }
	INLINE cimatrix operator |(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<cmatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmconv<cmatrix,imatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator |(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator |(const cmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<cmatrix_slice,imatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator |(const imatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsconv<cmatrix_slice,imatrix_slice,cimatrix>(m2,m1); }
	INLINE cimatrix operator &(const cmatrix &m1,const imatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<cmatrix,imatrix,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const imatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsect<cmatrix,imatrix,cimatrix>(m2,m1); }
	INLINE cimatrix operator &(const cmatrix &m,const imatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const imatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cmatrix_slice &ms,const imatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<imatrix,cmatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const imatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssect<cmatrix,imatrix_slice,cimatrix>(m,ms); }
	INLINE cimatrix operator &(const cmatrix_slice &m1,const imatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<cmatrix_slice,imatrix_slice,cimatrix>(m1,m2); }
	INLINE cimatrix operator &(const imatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CIMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssect<cmatrix_slice,imatrix_slice,cimatrix>(m2,m1); }

//------------- real x complex ------------------------
	INLINE cimatrix operator |(const rmatrix &rv1, const cmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,cmatrix,cimatrix>(rv1,rv2); }
	INLINE cimatrix operator |(const cmatrix &rv1, const rmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmconv<rmatrix,cmatrix,cimatrix>(rv2,rv1); }
	INLINE cimatrix operator |(const cmatrix &rv, const rmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,rmatrix_slice,cimatrix>(rv,sl); }
	INLINE cimatrix operator |(const rmatrix_slice &sl,const cmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,rmatrix_slice,cimatrix>(rv,sl); }
	INLINE cimatrix operator |(const cmatrix_slice &sl, const rmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,cmatrix_slice,cimatrix>(rv,sl); }
	INLINE cimatrix operator |(const rmatrix &rv,const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmsconv<rmatrix,cmatrix_slice,cimatrix>(rv,sl); }
	INLINE cimatrix operator |(const cmatrix_slice &sl1, const rmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,cmatrix_slice,cimatrix>(sl2,sl1); }
	INLINE cimatrix operator |(const rmatrix_slice &sl1, const cmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _msmsconv<rmatrix_slice,cmatrix_slice,cimatrix>(sl1,sl2); }
	
//------------- complex x complex ------------------------
	INLINE cimatrix operator |(const cmatrix &rv1, const cmatrix &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmconv<cmatrix,cmatrix,cimatrix>(rv1,rv2); }
	INLINE cimatrix operator |(const cmatrix &rv, const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,cmatrix_slice,cimatrix>(rv,sl); }
	INLINE cimatrix operator |(const cmatrix_slice &sl,const cmatrix &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _mmsconv<cmatrix,cmatrix_slice,cimatrix>(rv,sl); }
	INLINE cimatrix operator |(const cmatrix_slice &sl1, const cmatrix_slice &sl2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cimatrix>)
#else
	throw()
#endif
	{ return _msmsconv<cmatrix_slice,cmatrix_slice,cimatrix>(sl1,sl2); }
	
	INLINE bool operator ==(const cimatrix &m1,const cimatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const cimatrix &m1,const cimatrix &m2) throw() { return _mmneq(m1,m2); }
	INLINE bool operator <(const cimatrix &m1,const cimatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const cimatrix &m1,const cimatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const cimatrix &m1,const cimatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const cimatrix &m1,const cimatrix &m2) throw() { return _mmleq(m2,m1); }
	INLINE bool operator ==(const cimatrix &m1,const cimatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const cimatrix &m1,const cimatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
	INLINE bool operator <(const cimatrix &m1,const cimatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const cimatrix &m1,const cimatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const cimatrix &m1,const cimatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const cimatrix &m1,const cimatrix_slice &ms) throw() { return _msmleq(ms,m1); }
	INLINE bool operator ==(const cimatrix_slice &m1,const cimatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const cimatrix_slice &m1,const cimatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
	INLINE bool operator <(const cimatrix_slice &m1,const cimatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const cimatrix_slice &m1,const cimatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const cimatrix_slice &m1,const cimatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const cimatrix_slice &m1,const cimatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
	INLINE bool operator !(const cimatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const cimatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const cimatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const cimatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,cimatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,cimatrix_slice &r) throw() { return _msin(s,r); }

        //! Computes permutation of matrix according to permutation vectors, C=PAQ
        INLINE cimatrix cimatrix::operator()(const intvector& p, const intvector& q) {
          cimatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
            for(int j=0 ; j<RowLen(A) ; j++)
              A[i+Lb(A,1)][j+Lb(A,2)] = (*this)[p[i+Lb(p)]+Lb(A,1)][q[j+Lb(q)]+Lb(A,2)];
          return A;
        }

        //! Computes permutation of matrix according to permutation vector, C=PA
        INLINE cimatrix cimatrix::operator()(const intvector& p) {
          cimatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
              A[i+Lb(A,1)] = (*this)[p[i+Lb(p)]+Lb(A,1)];
          return A;
        }

        //! Computes permutation of matrix according to permutation matrix, C=PA
	INLINE cimatrix cimatrix::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

        //! Computes permutation of matrix according to permutation matrices, C=PAQ
        INLINE cimatrix cimatrix::operator()(const intmatrix& P, const intmatrix& Q) {
          intvector p = permvec(P);
          intvector q = perminv(permvec(Q));
          return (*this)(p,q);
        }

        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE civector civector::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

} // namespace cxsc

#endif
