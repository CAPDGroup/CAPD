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

/* CVS $Id: cmatrix.inl,v 1.28 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CMATRIX_INL_INCLUDED
#define _CXSC_CMATRIX_INL_INCLUDED

namespace cxsc {

INLINE cmatrix::cmatrix() throw():dat(NULL),lb1(1),ub1(0),lb2(1),ub2(0),xsize(0),ysize(0)
{
}

INLINE cmatrix::cmatrix(const complex &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new complex[1];
	*dat=r;
}

INLINE cmatrix::cmatrix(const real &r) throw():lb1(1),ub1(1),lb2(1),ub2(1),xsize(1),ysize(1)
{
	dat=new complex[1];
	*dat=r;
}

INLINE cmatrix::cmatrix(const cmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new complex[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE cmatrix::cmatrix(const rmatrix &rm) throw():lb1(rm.lb1),ub1(rm.ub1),lb2(rm.lb2),ub2(rm.ub2),xsize(rm.xsize),ysize(rm.ysize)
{
	dat=new complex[xsize*ysize];
	for(int i=0;i<xsize*ysize;i++)
		dat[i]=rm.dat[i];
}

INLINE cmatrix::cmatrix(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_BOUNDARIES):lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#else
	throw():lb1(1),ub1(m),lb2(1),ub2(n),xsize(n),ysize(m)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((n<0)||(m<0)) cxscthrow(ERROR_CMATRIX_WRONG_BOUNDARIES("cmatrix::cmatrix(const int &m, const int &n)"));
#endif
	dat=new complex[m*n];
}

INLINE cmatrix::cmatrix(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_BOUNDARIES):lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#else
	throw():lb1(m1),ub1(m2),lb2(n1),ub2(n2),xsize(n2-n1+1),ysize(m2-m1+1)
#endif
{
#if(CXSC_INDEX_CHECK)
	if((m2<m1)||(n2<n1)) cxscthrow(ERROR_CMATRIX_WRONG_BOUNDARIES("cmatrix::cmatrix(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
	dat=new complex[xsize*ysize];
}

INLINE cvector::cvector(const cmatrix_subv &v) throw():l(v.lb),u(v.ub),size(v.size)
{
	dat=new complex[size];
	for (int i=0, j=v.start;i<v.size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

INLINE cmatrix::cmatrix(const cvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new complex[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE cmatrix::cmatrix(const rvector &v) throw():lb1(v.l),ub1(v.u),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new complex[v.size];
	for(int i=0;i<v.size;i++)
		dat[i]=v.dat[i];
}

INLINE cmatrix::cmatrix(const cvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new complex[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}

INLINE cmatrix::cmatrix(const rvector_slice &v) throw():lb1(v.start),ub1(v.end),lb2(1),ub2(1),xsize(1),ysize(v.size)
{
	dat=new complex[v.size];
	for(int i=0,j=v.start-v.l;i<v.size;i++,j++)
		dat[i]=v.dat[j];
}


	INLINE cmatrix::cmatrix(const cmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new complex[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE cmatrix::cmatrix(const rmatrix_slice &sl) throw():lb1(sl.start1),ub1(sl.end1),lb2(sl.start2),ub2(sl.end2),xsize(sl.sxsize),ysize(sl.sysize)
	{
		int i,j;
		
		dat=new complex[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(sl.offset1+i)*sl.mxsize+sl.offset2+j];
			}
		}
	}

	INLINE cmatrix_subv Row(cmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE cmatrix_subv Col(cmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE cmatrix_subv Row(const cmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[i];
	}

	INLINE cmatrix_subv Col(const cmatrix &m,const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	
	{
		return m[Col(i)];
	}

	INLINE complex& cmatrix_subv::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC("complex &cmatrix_subv::operator [](const int &i) const"));
#endif
		return dat[start+((i-lb)*offset)];
	}

	INLINE complex& cmatrix_subv::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb)||(i>ub)) cxscthrow(ERROR_CVECTOR_ELEMENT_NOT_IN_VEC("complex &cmatrix_subv::operator [](const int &i)"));
#endif
		return dat[start+((i-lb)*offset)];
	}

		
	INLINE cmatrix_subv cmatrix::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix::operator [](const int &i)"));
#endif
		return cmatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE cmatrix_subv cmatrix::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return cmatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}

	INLINE cmatrix_subv cmatrix::operator [](const int &i) 
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<lb1)||(i>ub1)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix::operator [](const int &i)"));
#endif
		return cmatrix_subv(dat, lb2, ub2, xsize, xsize*(i-lb1),1);
	}
	
	INLINE cmatrix_subv cmatrix::operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<lb2)||(i.col()>ub2)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix::operator [](const cxscmatrix_column &i)"));
#endif
		return cmatrix_subv(dat, lb1, ub1, ysize, i.col()-lb2, xsize);
	}
	
	INLINE cmatrix_slice cmatrix::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m<1)||(n<1)||(m<lb1)||(n<lb2)||(m>ub1)||(n>ub2)) cxscthrow(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG("cmatrix_slice cmatrix::operator ()(const int &m, const int &n)"));
#endif
		return cmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE cmatrix_slice cmatrix::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
	if((m1<lb1)||(n1<lb2)||(m2>ub1)||(n2>ub2)) cxscthrow(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG("cmatrix_slice cmatrix::operator ()(const int &m1, const int &n1, const int &m2, const int &n2)"));
#endif
		return cmatrix_slice(*this,m1,m2,n1,n2);
	}

	INLINE cmatrix_subv cmatrix_slice::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix_slice::operator [](const int &i)"));
#endif
		return cmatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE cmatrix_subv cmatrix_slice::operator [](const cxscmatrix_column &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return cmatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}

	INLINE cmatrix_subv cmatrix_slice::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i<start1)||(i>end1)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix_slice::operator [](const int &i)"));
#endif
		return cmatrix_subv(dat, start2, end2, sxsize, mxsize*(i-start1+offset1)+offset2,1);
	}
	
	INLINE cmatrix_subv cmatrix_slice::operator [](const cxscmatrix_column &i) const
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((i.col()<start2)||(i.col()>end2)) cxscthrow(ERROR_CMATRIX_ROW_OR_COL_NOT_IN_MAT("cmatrix_subv cmatrix_slice::operator [](const cxscmatrix_column &i)"));
#endif
		return cmatrix_subv(dat, start1, end1, sysize, offset1*mxsize+i.col()-start2+offset2, mxsize);
	}
	
	INLINE cmatrix_slice cmatrix_slice::operator ()(const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m<1)||(n<1)||(m<start1)||(n<start2)||(m>end1)||(n>end2)) cxscthrow(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG("cmatrix_slice cmatrix_slice::operator ()(const int &m, const int &n)"));
#endif
		return cmatrix_slice(*this,1,m,1,n);
	}
	
	INLINE cmatrix_slice cmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if((m1<start1)||(n1<start2)||(m2>end1)||(n2>end2)) cxscthrow(ERROR_CMATRIX_SUB_ARRAY_TOO_BIG("cmatrix_slice cmatrix_slice::operator ()(const int &m1, const int &m2, const int &n1, const int &n2)"));
#endif
		return cmatrix_slice(*this,m1,m2,n1,n2);
	}

INLINE cmatrix_subv cmatrix_subv::operator ()(const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(1<lb||i>ub) cxscthrow(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG("cmatrix_subv cmatrix_subv::operator ()(const int &i)"));
#endif
	return cmatrix_subv(dat,1,i,i,start+(1-lb)*offset,offset);
}

INLINE cmatrix_subv cmatrix_subv::operator ()(const int &i1,const int &i2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG)
#else
	throw()
#endif
{
#if(CXSC_INDEX_CHECK)
	if(i1<lb||i2>ub) cxscthrow(ERROR_CVECTOR_SUB_ARRAY_TOO_BIG("cmatrix_subv cmatrix_subv::operator ()(const int &i1,const int &i2)"));
#endif
	return cmatrix_subv(dat,i1,i2,i2-i1+1,start+(i1-lb)*offset,offset);
}



	INLINE cmatrix_subv &cmatrix_subv::operator =(const cmatrix_subv &rv) throw() { return _mvmvassign(*this,rv); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const complex &r) throw() { return _mvsassign(*this,r); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const cvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,cvector(v)); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const rmatrix_subv &rv) throw() { return _mvvassign(*this,rvector(rv)); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const real &r) throw() { return _mvsassign(*this,r); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,v); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvassign(*this,cvector(v)); }
	INLINE cmatrix &cmatrix::operator =(const complex &r) throw() { return _msassign(*this,r); }
	INLINE cmatrix &cmatrix::operator =(const cmatrix &m) throw() { return _mmassign<cmatrix,cmatrix,complex>(*this,m, complex(0,0)); }
	INLINE cmatrix &cmatrix::operator =(const cmatrix_slice &ms) throw() { return _mmsassign<cmatrix,cmatrix_slice,complex>(*this,ms); }
	INLINE cmatrix &cmatrix::operator =(const cvector &v) throw() { return _mvassign<cmatrix,cvector,complex>(*this,v); }
	INLINE cmatrix &cmatrix::operator =(const cvector_slice &v) throw() { return _mvassign<cmatrix,cvector,complex>(*this,cvector(v)); }
	INLINE cmatrix &cmatrix::operator =(const real &r) throw() { return _msassign(*this,complex(r)); }
	INLINE cmatrix &cmatrix::operator =(const rmatrix &m) throw() { return _mmassign<cmatrix,rmatrix,complex>(*this,m,complex(0,0)); }
	INLINE cmatrix &cmatrix::operator =(const rmatrix_slice &ms) throw() { return _mmsassign<cmatrix,rmatrix_slice,complex>(*this,ms); }
	INLINE cmatrix &cmatrix::operator =(const rvector &v) throw() { return _mvassign<cmatrix,rvector,complex>(*this,v); }
	INLINE cmatrix &cmatrix::operator =(const rvector_slice &v) throw() { return _mvassign<cmatrix,rvector,complex>(*this,rvector(v)); }
	
	INLINE cmatrix::operator void*() throw() { return _mvoid(*this); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const complex &r) throw() { return _mssassign(*this,r); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cmatrix(v)); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const cvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,cmatrix(cvector(v))); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,m); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsassign(*this,ms); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const real &r) throw() { return _mssassign(*this,r); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(v)); }
	INLINE cmatrix_slice &cmatrix_slice::operator =(const rvector_slice &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmassign(*this,rmatrix(rvector(v))); }
	INLINE cmatrix_slice::operator void*() throw() { return _msvoid(*this); }
	INLINE cvector operator /(const cmatrix_subv &rv, const complex &s) throw() { return _mvsdiv<cmatrix_subv,complex,cvector>(rv,s); }
	INLINE cvector operator *(const cmatrix_subv &rv, const complex &s) throw() { return _mvsmult<cmatrix_subv,complex,cvector>(rv,s); }
	INLINE cvector operator *(const complex &s, const cmatrix_subv &rv) throw() { return _mvsmult<cmatrix_subv,complex,cvector>(rv,s); }
	INLINE cmatrix_subv &cmatrix_subv::operator *=(const complex &c) throw() { return _mvsmultassign(*this,c); }
	INLINE cmatrix_subv &cmatrix_subv::operator +=(const complex &c) throw() { return _mvsplusassign(*this,c); }
	INLINE cmatrix_subv &cmatrix_subv::operator -=(const complex &c) throw() { return _mvsminusassign(*this,c); }
	INLINE cmatrix_subv &cmatrix_subv::operator /=(const complex &c) throw() { return _mvsdivassign(*this,c); }
	INLINE rvector abs(const cmatrix_subv &mv) throw() { return _mvabs<cmatrix_subv,rvector>(mv); }
	INLINE rvector Im(const cmatrix_subv &mv) throw() { return _mvim<cmatrix_subv,rvector>(mv); }
	INLINE rvector Re(const cmatrix_subv &mv) throw() { return _mvre<cmatrix_subv,rvector>(mv); }
	INLINE cmatrix_subv &SetIm(cmatrix_subv &mv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetim(mv,rv); }
	INLINE cmatrix_subv &SetRe(cmatrix_subv &mv,const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvsetre(mv,rv); }
	INLINE cmatrix_subv &SetRe(cmatrix_subv &iv,const real &r) throw() { return _mvssetre(iv,r); }
	INLINE cmatrix_subv &SetIm(cmatrix_subv &iv,const real &r) throw() { return _mvssetim(iv,r); }
	INLINE cvector &cvector::operator =(const cmatrix_subv &mv) throw() { return _vmvassign<cvector,cmatrix_subv,complex>(*this,mv); }
	INLINE cvector_slice &cvector_slice::operator =(const cmatrix_subv &mv) throw() { return _vsvassign(*this,cvector(mv)); }

	INLINE complex operator *(const cmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvcmult<cmatrix_subv,cmatrix_subv,complex>(rv1,rv2); }
	INLINE complex operator *(const cvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcmult<cvector,cmatrix_subv,complex>(rv1,rv2); }
	INLINE complex operator *(const cmatrix_subv &rv1,const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcmult<cvector,cmatrix_subv,complex>(rv2,rv1); }
	INLINE complex operator *(const cvector_slice &sl,const cmatrix_subv &sv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcmult<cvector,cmatrix_subv,complex>(cvector(sl),sv); }
	INLINE complex operator *(const cmatrix_subv &mv,const cvector_slice &vs)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvcmult<cvector,cmatrix_subv,complex>(cvector(vs),mv); }
	INLINE cvector operator +(const cmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvplus<cmatrix_subv,cmatrix_subv,cvector>(rv1,rv2); }
	INLINE cvector operator +(const cmatrix_subv &rv1,const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cmatrix_subv,cvector,cvector>(rv1,rv2); }
	INLINE cvector operator +(const cvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cmatrix_subv,cvector,cvector>(rv2,rv1); }
	INLINE cvector operator +(const cvector_slice &sl,const cmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cmatrix_subv,cvector,cvector>(mv,cvector(sl)); }
	INLINE cvector operator +(const cmatrix_subv &mv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplus<cmatrix_subv,cvector,cvector>(mv,cvector(sl)); }
	INLINE cmatrix_subv &cmatrix_subv::operator +=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE cmatrix_subv &cmatrix_subv::operator +=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,cvector(rv)); }
	INLINE cvector operator -(const cmatrix_subv & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvmvminus<cmatrix_subv,cmatrix_subv,cvector>(rv1,rv2); }
	INLINE cvector operator -(const cvector & rv1, const cmatrix_subv &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<cvector,cmatrix_subv,cvector>(rv1,rv2); }
	INLINE cvector operator -(const cmatrix_subv &rv1,const cvector &rv2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<cmatrix_subv,cvector,cvector>(rv1,rv2); }
	INLINE cvector operator -(const cvector_slice &sl,const cmatrix_subv &mv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmvminus<cvector,cmatrix_subv,cvector>(cvector(sl),mv); }
	INLINE cvector operator -(const cmatrix_subv &mv,const cvector_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminus<cmatrix_subv,cvector,cvector>(mv,cvector(sl)); }
	INLINE cmatrix_subv &cmatrix_subv::operator -=(const cvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE cmatrix_subv &cmatrix_subv::operator -=(const cvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,cvector(rv)); }
//  real

	INLINE cmatrix_subv &cmatrix_subv::operator +=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,rv); }
	INLINE cmatrix_subv &cmatrix_subv::operator +=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvplusassign(*this,cvector(rv)); }
	INLINE cmatrix_subv &cmatrix_subv::operator -=(const rvector &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rv); }
	INLINE cmatrix_subv &cmatrix_subv::operator -=(const rvector_slice &rv)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CVECTOR_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvvminusassign(*this,rvector(rv)); }

// matrix x matrix	
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cmatrix::cmatrix(const cmatrix &rm)
	*/
	INLINE cmatrix _cmatrix(const cmatrix &rm) throw() { return rm; }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cmatrix::cmatrix(const cvector &v)
	*/
	INLINE cmatrix _cmatrix(const cvector &v) throw() { return cmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cmatrix::cmatrix(const cvector_slice &v)
	*/
	INLINE cmatrix _cmatrix(const cvector_slice &v) throw() { return cmatrix(v); }
	/*!
	\deprecated use standard contructors for typecasting

	\sa cxsc::cmatrix::cmatrix(const complex &r)
	*/
	INLINE cmatrix _cmatrix(const complex &r) throw() { return cmatrix(r); }
	INLINE int Lb(const cmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mlb(rm,i); }
	INLINE int Ub(const cmatrix &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mub(rm,i); }
	INLINE int Lb(const cmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _mslb(rm,i); }
	INLINE int Ub(const cmatrix_slice &rm, const int &i)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msub(rm,i); }
	INLINE cmatrix &SetLb(cmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetlb(m,i,j); }
	INLINE cmatrix &SetUb(cmatrix &m, const int &i,const int &j)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_ROW_OR_COL)
#else
	throw()
#endif
	{ return _msetub(m,i,j); }
	
        INLINE int RowLen ( const cmatrix& A ) // Length of the rows of a complex matrix
        { return Ub(A,2)-Lb(A,2)+1; }          //---------------------------------------

        INLINE int ColLen ( const cmatrix& A ) // Length of the columns of a complex matrix
        { return Ub(A,1)-Lb(A,1)+1; }          //------------------------------------------
	
        INLINE int RowLen ( const cmatrix_slice& A ) // Length of the rows of a complex matrix
        { return Ub(A,2)-Lb(A,2)+1; }                //---------------------------------------

        INLINE int ColLen ( const cmatrix_slice& A ) // Length of the columns of a complex matrix
        { return Ub(A,1)-Lb(A,1)+1; }                //------------------------------------------
	
	INLINE void Resize(cmatrix &A) throw() { _mresize(A);}
	INLINE void Resize(cmatrix &A,const int &m, const int &n)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<cmatrix,complex>(A,m,n); }
	INLINE void Resize(cmatrix &A,const int &m1, const int &m2,const int &n1,const int &n2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_WRONG_BOUNDARIES)
#else
	throw()
#endif
	{ _mresize<cmatrix,complex>(A,m1,m2,n1,n2); }
	INLINE rmatrix abs(const cmatrix &m) throw() { return _mabs<cmatrix,rmatrix>(m); }
	INLINE rmatrix abs(const cmatrix_slice &ms) throw() { return _msabs<cmatrix_slice,rmatrix>(ms); }
	INLINE rmatrix Im(const cmatrix &m) throw() { return _mim<cmatrix,rmatrix>(m); }
	INLINE rmatrix Re(const cmatrix &m) throw() { return _mre<cmatrix,rmatrix>(m); }
	INLINE rmatrix Im(const cmatrix_slice &m) throw() { return _msim<cmatrix_slice,rmatrix>(m); }
	INLINE rmatrix Re(const cmatrix_slice &m) throw() { return _msre<cmatrix_slice,rmatrix>(m); }
	INLINE cmatrix &SetIm(cmatrix &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetim<cmatrix,rmatrix>(cm,rm); }
	INLINE cmatrix_slice &SetIm(cmatrix_slice &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetim<cmatrix_slice,rmatrix>(cm,rm); }
	INLINE cmatrix &SetIm(cmatrix &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetim<cmatrix,rmatrix_slice>(cm,rm); }
	INLINE cmatrix_slice &SetIm(cmatrix_slice &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetim<cmatrix_slice,rmatrix_slice>(cm,rm); }
	INLINE cmatrix &SetRe(cmatrix &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsetre<cmatrix,rmatrix>(cm,rm); }
	INLINE cmatrix_slice &SetRe(cmatrix_slice &cm,const rmatrix &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsetre<cmatrix_slice,rmatrix>(cm,rm); }
	INLINE cmatrix &SetRe(cmatrix &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmssetre<cmatrix,rmatrix_slice>(cm,rm); }
	INLINE cmatrix_slice &SetRe(cmatrix_slice &cm,const rmatrix_slice &rm)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmssetre<cmatrix_slice,rmatrix_slice>(cm,rm); }
	INLINE complex::complex(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ)
#else
	throw()
#endif
	{ _smconstr(*this,m); }
//	INLINE complex complex::_complex(const cmatrix &m) throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE cmatrix operator *(const complex &c, const cmatrix &m) throw() { return _smmult<complex,cmatrix,cmatrix>(c,m); }
	INLINE cmatrix operator *(const complex &c, const cmatrix_slice &ms) throw() { return _smsmult<complex,cmatrix_slice,cmatrix>(c,ms); }
	INLINE cmatrix operator *(const cmatrix &m,const complex &c) throw() { return _smmult<complex,cmatrix,cmatrix>(c,m); }
	INLINE cmatrix operator *(const cmatrix_slice &ms,const complex &c) throw() { return _smsmult<complex,cmatrix_slice,cmatrix>(c,ms); }
	INLINE cmatrix &operator *=(cmatrix &m,const complex &c) throw() { return _msmultassign(m,c); }
	INLINE cmatrix_slice &cmatrix_slice::operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cmatrix_slice &cmatrix_slice::operator *=(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cmatrix_slice &cmatrix_slice::operator *=(const complex &c) throw() { return _mssmultassign(*this,c); }
	INLINE cmatrix operator /(const cmatrix &m,const complex &c) throw() { return _msdiv<cmatrix,complex,cmatrix>(m,c); }
	INLINE cmatrix operator /(const cmatrix_slice &ms, const complex &c) throw() { return _mssdiv<cmatrix_slice,complex,cmatrix>(ms,c); }
	INLINE cmatrix &operator /=(cmatrix &m,const complex &c) throw() { return _msdivassign(m,c); }
	INLINE cmatrix_slice &cmatrix_slice::operator /=(const complex &c) throw() { return _mssdivassign(*this,c); }
	INLINE cmatrix operator *(const real &c, const cmatrix &m) throw() { return _smmult<real,cmatrix,cmatrix>(c,m); }
	INLINE cmatrix operator *(const real &c, const cmatrix_slice &ms) throw() { return _smsmult<real,cmatrix_slice,cmatrix>(c,ms); }
	INLINE cmatrix operator *(const cmatrix &m,const real &c) throw() { return _smmult<real,cmatrix,cmatrix>(c,m); }
	INLINE cmatrix operator *(const cmatrix_slice &ms,const real &c) throw() { return _smsmult<real,cmatrix_slice,cmatrix>(c,ms); }
	INLINE cmatrix &operator *=(cmatrix &m,const real &c) throw() { return _msmultassign(m,c); }
	INLINE cmatrix_slice &cmatrix_slice::operator *=(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cmatrix_slice &cmatrix_slice::operator *=(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return (*this=*this*m); }
	INLINE cmatrix_slice &cmatrix_slice::operator *=(const real &c) throw() { return _mssmultassign(*this,c); }
	INLINE cmatrix operator /(const cmatrix &m,const real &c) throw() { return _msdiv<cmatrix,real,cmatrix>(m,c); }
	INLINE cmatrix operator /(const cmatrix_slice &ms, const real &c) throw() { return _mssdiv<cmatrix_slice,real,cmatrix>(ms,c); }
	INLINE cmatrix &operator /=(cmatrix &m,const real &c) throw() { return _msdivassign(m,c); }
	INLINE cmatrix_slice &cmatrix_slice::operator /=(const real &c) throw() { return _mssdivassign(*this,c); }
//	INLINE complex::complex(const rmatrix &m) throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); }
//	INLINE complex complex::_complex(const cmatrix &m) throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ,ERROR_CMATRIX_USE_OF_UNINITIALIZED_OBJ) { _smconstr(*this,m); return *this; }
	INLINE cmatrix operator *(const complex &c, const rmatrix &m) throw() { return _smmult<complex,rmatrix,cmatrix>(c,m); }
	INLINE cmatrix operator *(const complex &c, const rmatrix_slice &ms) throw() { return _smsmult<complex,rmatrix_slice,cmatrix>(c,ms); }
	INLINE cmatrix operator *(const rmatrix &m,const complex &c) throw() { return _smmult<complex,rmatrix,cmatrix>(c,m); }
	INLINE cmatrix operator *(const rmatrix_slice &ms,const complex &c) throw() { return _smsmult<complex,rmatrix_slice,cmatrix>(c,ms); }
	INLINE cmatrix operator /(const rmatrix &m,const complex &c) throw() { return _msdiv<rmatrix,complex,cmatrix>(m,c); }
	INLINE cmatrix operator /(const rmatrix_slice &ms, const complex &c) throw() { return _mssdiv<rmatrix_slice,complex,cmatrix>(ms,c); }
	INLINE cvector::cvector(const cmatrix &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmconstr<cvector,cmatrix,complex>(*this,sl); }
	INLINE cvector::cvector(const cmatrix_slice &sl)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ _vmsconstr<cvector,cmatrix_slice,complex>(*this,sl); }
	INLINE cvector &cvector::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<cvector,cmatrix,complex>(*this,m); }
	INLINE cvector &cvector::operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vmassign<cvector,cmatrix,complex>(*this,cmatrix(m)); }
	INLINE cvector_slice & cvector_slice::operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR__OP_WITH_WRONG_DIM<cvector>,ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _vsvassign(*this,cvector(cmatrix(m))); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,cvector(m)); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,cvector(cmatrix(m))); }
	INLINE cvector operator *(const cmatrix &m,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcmult<cmatrix,cvector,cvector>(m,v); }
	INLINE cvector operator *(const cmatrix_slice &ms,const cvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcmult<cmatrix_slice,cvector,cvector>(ms,v); }
	INLINE cvector operator *(const cvector &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmult<cvector,cmatrix,cvector>(v,m); }
	INLINE cvector operator *(const cvector &v,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscmult<cvector,cmatrix_slice,cvector>(v,ms); }
	INLINE cvector &operator *=(cvector &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmultassign<cvector,cmatrix,complex>(v,m); }
	INLINE cvector &operator *=(cvector &v,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscmultassign<cvector,cmatrix_slice,complex>(v,ms); }
	INLINE cvector_slice &cvector_slice::operator *=(const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vsmcmultassign<cvector_slice,cmatrix,complex>(*this,m); }
	INLINE cvector operator *(const cvector_slice &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmult<cvector,cmatrix,cvector>(cvector(v),m); }
	INLINE cvector operator *(const cvector_slice &v,const cmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscmult<cvector,cmatrix_slice,cvector>(cvector(v),m); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(m)); }
	INLINE cmatrix_subv &cmatrix_subv::operator =(const rmatrix_slice &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_TYPE_CAST_OF_THICK_OBJ)
#else
	throw()
#endif
	{ return _mvvassign(*this,rvector(rmatrix(m))); }
	INLINE cvector operator *(const rvector &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmult<rvector,cmatrix,cvector>(v,m); }
	INLINE cvector operator *(const rvector &v,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmscmult<rvector,cmatrix_slice,cvector>(v,ms); }
	INLINE cvector operator *(const rvector_slice &v,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _vmcmult<cvector,cmatrix,cvector>(cvector(v),m); }
	INLINE cvector operator *(const cmatrix &m,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mvcmult<cmatrix,rvector,cvector>(m,v); }
	INLINE cvector operator *(const cmatrix_slice &ms,const rvector &v)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msvcmult<cmatrix_slice,rvector,cvector>(ms,v); }

	INLINE const cmatrix &operator +(const cmatrix &m) throw() { return m; }
	INLINE cmatrix operator +(const cmatrix_slice &m) throw() { return cmatrix(m); }
	INLINE cmatrix operator +(const cmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<cmatrix,cmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator +(const cmatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,cmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator +(const cmatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,cmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator +(const cmatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<cmatrix_slice,cmatrix_slice,cmatrix>(m1,m2); }
	INLINE cmatrix &operator +=(cmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE cmatrix &operator +=(cmatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE cmatrix_slice &cmatrix_slice::operator +=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE cmatrix_slice &cmatrix_slice::operator +=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE cmatrix operator -(const cmatrix &m) throw() { return _mminus(m); }
	INLINE cmatrix operator -(const cmatrix_slice &m) throw() { return _msminus<cmatrix_slice,cmatrix>(m); }
	INLINE cmatrix operator -(const cmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cmatrix,cmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator -(const cmatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cmatrix,cmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator -(const cmatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cmatrix_slice,cmatrix,cmatrix>(ms,m); }
	INLINE cmatrix operator -(const cmatrix_slice &ms1,const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cmatrix_slice,cmatrix_slice,cmatrix>(ms1,ms2); }
	INLINE cmatrix &operator -=(cmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE cmatrix &operator -=(cmatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE cmatrix_slice &cmatrix_slice::operator -=(const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE cmatrix_slice &cmatrix_slice::operator -=(const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE cmatrix operator *(const cmatrix &m1, const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcmult<cmatrix,cmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator *(const cmatrix &m1, const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscmult<cmatrix,cmatrix_slice,cmatrix>(m1,ms); }
	INLINE cmatrix operator *(const cmatrix_slice &ms, const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcmult<cmatrix_slice,cmatrix,cmatrix>(ms,m1); }
	INLINE cmatrix operator *(const cmatrix_slice &ms1, const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscmult<cmatrix_slice,cmatrix_slice,cmatrix>(ms1,ms2); }
	INLINE cmatrix &operator *=(cmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcmultassign<cmatrix,cmatrix,complex>(m1,m2); }
	INLINE cmatrix &operator *=(cmatrix &m1,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscmultassign<cmatrix,cmatrix_slice,complex>(m1,ms); }
	INLINE cmatrix operator +(const rmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,cmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator +(const cmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplus<rmatrix,cmatrix,cmatrix>(m2,m1); }
	INLINE cmatrix operator +(const rmatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,cmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator +(const cmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,rmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator +(const rmatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<cmatrix,rmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator +(const cmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplus<rmatrix,cmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator +(const rmatrix_slice &m1,const cmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,cmatrix_slice,cmatrix>(m1,m2); }
	INLINE cmatrix operator +(const cmatrix_slice &m1,const rmatrix_slice &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplus<rmatrix_slice,cmatrix_slice,cmatrix>(m2,m1); }
	INLINE cmatrix &operator +=(cmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmplusassign(m1,m2); }
	INLINE cmatrix &operator +=(cmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsplusassign(m1,ms); }
	INLINE cmatrix_slice &cmatrix_slice::operator +=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmplusassign(*this,m1); }
	INLINE cmatrix_slice &cmatrix_slice::operator +=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsplusassign(*this,ms2); }
	INLINE cmatrix operator -(const rmatrix &m1,const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<rmatrix,cmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator -(const cmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminus<cmatrix,rmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator -(const rmatrix &m,const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<rmatrix,cmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator -(const cmatrix &m,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminus<cmatrix,rmatrix_slice,cmatrix>(m,ms); }
	INLINE cmatrix operator -(const rmatrix_slice &ms,const cmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<rmatrix_slice,cmatrix,cmatrix>(ms,m); }
	INLINE cmatrix operator -(const cmatrix_slice &ms,const rmatrix &m)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminus<cmatrix_slice,rmatrix,cmatrix>(ms,m); }
	INLINE cmatrix operator -(const rmatrix_slice &ms1,const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<rmatrix_slice,cmatrix_slice,cmatrix>(ms1,ms2); }
	INLINE cmatrix operator -(const cmatrix_slice &ms1,const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminus<cmatrix_slice,rmatrix_slice,cmatrix>(ms1,ms2); }
	INLINE cmatrix &operator -=(cmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmminusassign(m1,m2); }
	INLINE cmatrix &operator -=(cmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmsminusassign(m1,ms); }
	INLINE cmatrix_slice &cmatrix_slice::operator -=(const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmminusassign(*this,m1); }
	INLINE cmatrix_slice &cmatrix_slice::operator -=(const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmsminusassign(*this,ms2); }
	INLINE cmatrix operator *(const rmatrix &m1, const cmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcmult<rmatrix,cmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator *(const cmatrix &m1, const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcmult<cmatrix,rmatrix,cmatrix>(m1,m2); }
	INLINE cmatrix operator *(const rmatrix &m1, const cmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscmult<rmatrix,cmatrix_slice,cmatrix>(m1,ms); }
	INLINE cmatrix operator *(const cmatrix &m1, const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscmult<cmatrix,rmatrix_slice,cmatrix>(m1,ms); }
	INLINE cmatrix operator *(const rmatrix_slice &ms, const cmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcmult<rmatrix_slice,cmatrix,cmatrix>(ms,m1); }
	INLINE cmatrix operator *(const cmatrix_slice &ms, const rmatrix &m1)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmcmult<cmatrix_slice,rmatrix,cmatrix>(ms,m1); }
	INLINE cmatrix operator *(const rmatrix_slice &ms1, const cmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscmult<rmatrix_slice,cmatrix_slice,cmatrix>(ms1,ms2); }
	INLINE cmatrix operator *(const cmatrix_slice &ms1, const rmatrix_slice &ms2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _msmscmult<cmatrix_slice,rmatrix_slice,cmatrix>(ms1,ms2); }
	INLINE cmatrix &operator *=(cmatrix &m1,const rmatrix &m2)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmcmultassign<cmatrix,rmatrix,complex>(m1,m2); }
	INLINE cmatrix &operator *=(cmatrix &m1,const rmatrix_slice &ms)
#if(CXSC_INDEX_CHECK)
	throw(ERROR_CMATRIX_OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{ return _mmscmultassign<cmatrix,rmatrix_slice,complex>(m1,ms); }
	INLINE bool operator ==(const cmatrix &m1,const cmatrix &m2) throw() { return _mmeq(m1,m2); }
	INLINE bool operator !=(const cmatrix &m1,const cmatrix &m2) throw() { return _mmneq(m1,m2); }
/*	INLINE bool operator <(const cmatrix &m1,const cmatrix &m2) throw() { return _mmless(m1,m2); }
	INLINE bool operator <=(const cmatrix &m1,const cmatrix &m2) throw() { return _mmleq(m1,m2); }
	INLINE bool operator >(const cmatrix &m1,const cmatrix &m2) throw() { return _mmless(m2,m1); }
	INLINE bool operator >=(const cmatrix &m1,const cmatrix &m2) throw() { return _mmleq(m2,m1); }
*/
	INLINE bool operator ==(const cmatrix &m1,const cmatrix_slice &ms) throw() { return _mmseq(m1,ms); }
	INLINE bool operator !=(const cmatrix &m1,const cmatrix_slice &ms) throw() { return _mmsneq(m1,ms); }
/*	INLINE bool operator <(const cmatrix &m1,const cmatrix_slice &ms) throw() { return _mmsless(m1,ms); }
	INLINE bool operator <=(const cmatrix &m1,const cmatrix_slice &ms) throw() { return _mmsleq(m1,ms); }
	INLINE bool operator >(const cmatrix &m1,const cmatrix_slice &ms) throw() { return _msmless(ms,m1); }
	INLINE bool operator >=(const cmatrix &m1,const cmatrix_slice &ms) throw() { return _msmleq(ms,m1); }
*/
	INLINE bool operator ==(const cmatrix_slice &m1,const cmatrix_slice &m2) throw() { return _msmseq(m1,m2); }
	INLINE bool operator !=(const cmatrix_slice &m1,const cmatrix_slice &m2) throw() { return _msmsneq(m1,m2); }
/*	INLINE bool operator <(const cmatrix_slice &m1,const cmatrix_slice &m2) throw() { return _msmsless(m1,m2); }
	INLINE bool operator <=(const cmatrix_slice &m1,const cmatrix_slice &m2) throw() { return _msmsleq(m1,m2); }
	INLINE bool operator >(const cmatrix_slice &m1,const cmatrix_slice &m2) throw() { return _msmsless(m2,m1); }
	INLINE bool operator >=(const cmatrix_slice &m1,const cmatrix_slice &m2) throw() { return _msmsleq(m2,m1); }
*/
	INLINE bool operator !(const cmatrix &ms) throw() { return _mnot(ms); }
	INLINE bool operator !(const cmatrix_slice &ms) throw() { return _msnot(ms); }
	INLINE std::ostream &operator <<(std::ostream &s,const cmatrix &r) throw() { return _mout(s,r); }
	INLINE std::ostream &operator <<(std::ostream &s,const cmatrix_slice &r) throw() { return _msout(s,r); }
	INLINE std::istream &operator >>(std::istream &s,cmatrix &r) throw() { return _min(s,r); }
	INLINE std::istream &operator >>(std::istream &s,cmatrix_slice &r) throw() { return _msin(s,r); }

        //! Computes permutation of matrix according to permutation vectors, C=PAQ
        INLINE cmatrix cmatrix::operator()(const intvector& p, const intvector& q) {
          cmatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
            for(int j=0 ; j<RowLen(A) ; j++)
              A[i+Lb(A,1)][j+Lb(A,2)] = (*this)[p[i+Lb(p)]+Lb(A,1)][q[j+Lb(q)]+Lb(A,2)];
          return A;
        }

        //! Computes permutation of matrix according to permutation vector, C=PA
        INLINE cmatrix cmatrix::operator()(const intvector& p) {
          cmatrix A(*this);
          for(int i=0 ; i<ColLen(A) ; i++)
              A[i+Lb(A,1)] = (*this)[p[i+Lb(p)]+Lb(A,1)];
          return A;
        }

        //! Computes permutation of matrix according to permutation matrix, C=PA
	INLINE cmatrix cmatrix::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

        //! Computes permutation of matrix according to permutation matrices, C=PAQ
        INLINE cmatrix cmatrix::operator()(const intmatrix& P, const intmatrix& Q) {
          intvector p = permvec(P);
          intvector q = perminv(permvec(Q));
          return (*this)(p,q);
        }

        //! Computes permutation of vector according to permutation matrix, C=Px
        INLINE cvector cvector::operator()(const intmatrix& P) {
          intvector p = permvec(P);
          return (*this)(p);
        }

} // namespace cxsc

#endif
