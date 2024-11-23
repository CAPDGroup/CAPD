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

/* CVS $Id: cxscmatr.hpp,v 1.23 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CXSCMATRIX_HPP_INCLUDED
#define _CXSC_CXSCMATRIX_HPP_INCLUDED

#include "cxscvect.hpp"
#include "except.hpp"
#include <iostream>

namespace cxsc {

enum ASPECT { ROW=1,  COL=2 };


class cxscmatrix_column
{
	private:
	int c;
	public:
	inline cxscmatrix_column(const int &i) throw():c(i) { }
	inline int col() const throw() { return c; }
};

inline int Row(const int &i) throw() { return i; }
inline cxscmatrix_column Col(const int &i) throw() { return cxscmatrix_column(i); }

template <class T>
class cxscmatrix_subv
{
	friend class cxscvector<T>;
	private:
	T *dat;
	int lb,ub;
	int size,start,offset; // start=first element index 0..n-1
	
	public:
	//----------------- Konstruktoren ----------------------------------

	cxscmatrix_subv<T> (T *d, const int &l, const int &u, const int &s, const int &st, const int &o) throw():dat(d),lb(l),ub(u),size(s),start(st),offset(o) { }

	//---------------------- Standardfunktionen ------------------------

	inline T &operator [](const int &i) throw() { return dat[start+((i-lb)*offset)]; }

};

template <class T>
class cxscmatrix_slice;

template <class T>
class cxscmatrix
{
	friend class cxscmatrix_slice<T>;
	friend class cxscmatrix_subv<T>;
	private:
	T *dat;
	int lb1,ub1,lb2,ub2,xsize,ysize;

	public:
	//--------------------------  Konstruktoren ----------------------------

	inline cxscmatrix<T>(const cxscmatrix &rm) throw()
	{
		lb1=rm.lb1; ub1=rm.ub1;
		lb2=rm.lb2; ub2=rm.ub2;
		xsize=rm.xsize;
		ysize=rm.ysize;
		dat=new T[xsize*ysize];
	}
	
	inline cxscmatrix<T>() throw()
	{
		dat=NULL;
		lb1=lb2=1;
		ub1=ub2=xsize=ysize=0;
	}

	explicit inline cxscmatrix<T>(const int &m, const int &n) throw()
	{
		lb1=lb2=1;
		ub1=ysize=m;
		ub2=xsize=n;
		dat=new T[m*n];
	}

	explicit inline cxscmatrix<T>(const int &m1, const int &n1, const int &m2, const int &n2) throw()
	{	
		lb1=m1; ub1=m2;
		lb2=n1; ub2=n2;
		xsize=n2-n1+1;
		ysize=m2-m2+1;
		dat=new T[xsize*ysize];
	}

	inline cxscmatrix<T>(const cxscmatrix_slice<T> &sl) throw()
	{
		int i,j;
		
		lb1=sl.ystart; ub1=sl.yend;
		lb2=sl.xstart; ub2=sl.xend;
		ysize=ub1-lb1+1;
		xsize=ub2-lb2+1;
		dat=new T[xsize*ysize];
		for (i=0;i<ysize;i++)
		{
			for(j=0;j<xsize;j++)
			{
				dat[i*xsize+j]=sl.dat[(lb1-sl.lb1+i)*xsize+lb2-sl.lb2+j];
			}
		}
	}
	
	inline cxscmatrix<T> &operator =(const T &r) throw()
	{
		for(int i=0;i<xsize*ysize;i++)
			dat[i]=r;
		return *this;
	}

	inline cxscmatrix<T> &operator =(const cxscmatrix<T> &m) throw()
	{
		T *ndat=new T[m.xsize*m.ysize];
		for(int i=0;i<m.xsize*m.ysize;i++)
			ndat[i]=m.dat[i];

		delete [] dat;
		dat=ndat;
		return *this;
	}

	cxscmatrix<T> &operator =(const cxscvector<T> &v) throw();

	//--------------------------- Destruktoren -----------------------------

	inline ~cxscmatrix() throw()
	{
		delete [] dat;
	}

	//------------------------- Standardfunktionen -------------------------

	friend inline int Lb(const cxscmatrix<T> &rm, const ASPECT &i) throw() { return (i==ROW?rm.lb1:rm.lb2); }
	friend inline int Ub(const cxscmatrix<T> &rm, const ASPECT &i) throw() { return (i==ROW?rm.ub1:rm.ub2); }
	friend inline void SetLb(cxscmatrix<T> &m, const ASPECT &i,const int &j) throw(WRONG_ROW_OR_COL)
	{
		switch(i)
		{
			case ROW:
				m.lb2=j;
				m.ub2=j+m.xsize-1;
				break;
			case COLUMN:
				m.lb1=j;
				m.ub1=j+m.ysize-1;
				break;
			default:
				throw WRONG_ROW_OR_COL();
		}
	}
	
	friend inline void SetUb(cxscmatrix<T> &m, const ASPECT &i,const int &j) throw(WRONG_ROW_OR_COL)
	{
		switch(i)
		{
			case ROW:
				m.ub2=j;
				m.lb2=j-m.xsize+1;
				break;
			case COLUMN:
				m.ub1=j;
				m.lb1=j-m.ysize+1;
				break;
			default:
				throw WRONG_ROW_OR_COL();
		}
	}
	
	friend inline void Resize(cxscmatrix<T> &A) throw()
	{
		A.xsize=A.ysize=A.ub1=A.ub2=0;
		A.lb1=A.lb2=1;
		delete [] A.dat;
		A.dat=NULL;
	}

	friend inline void Resize(cxscmatrix<T> &A,const int &m, const int &n) throw()
	{
		T *ndat=new T[m*n];
		int l1,u1,l2,u2;
		int i,j;
		l1=(A.lb1>1)?A.l1:1;
		u1=(A.ub1<m)?A.u1:m;
		l2=(A.lb2>1)?A.l2:1;
		u2=(A.ub2<n)?A.ub2:n;
		for(i=0;i<l1-1;i++)
		{
			for(j=0;j<n;j++)
				ndat[i*n+j]=(T)0;
		}
		for(;i<u1;i++)
		{
			for(j=0;j<l2-1;j++)
				ndat[i*n+j]=(T)0;
			for(;j<u2;j++)
				ndat[i*n+j]=A.dat[(i-A.lb1+1)*A.xsize+j];
			for(;j<n;j++)
				ndat[i*n+j]=(T)0;
		}
		for(;i<m;i++)
		{
			for(j=0;j<n;j++)
				ndat[i*n+j]=(T)0;
		}
		delete [] A.dat;
		A.dat=ndat;
		A.xsize=A.ub2=n;
		A.ysize=A.ub1=m;
		A.lb1=A.lb2=1;
	}

	friend inline void Resize(cxscmatrix<T> &A,const int &m1, const int &m2,const int &n1,const int &n2) throw()
	{
		int m,n;
		n=A.xsize=n2-n1+1;
		m=A.ysize=m2-m1+1;
		T *ndat=new T[A.xsize*A.ysize];
		int l1,u1,l2,u2;
		l1=(A.lb1>m1)?A.l1:m1;
		u1=(A.ub1<m2)?A.u1:m2;
		l2=(A.lb2>n1)?A.l2:n1;
		u2=(A.ub2<n2)?A.ub2:n2;
		for(i=0;i<l1-1;i++)
		{
			for(j=0;j<n;j++)
				ndat[i*n+j]=(T)0;
		}
		for(;i<u1;i++)
		{
			for(j=0;j<l2-1;j++)
				ndat[i*n+j]=(T)0;
			for(;j<u2;j++)
				ndat[i*n+j]=A.dat[(i-A.lb1+1)*A.xsize+j];
			for(;j<n;j++)
				ndat[i*n+j]=(T)0;
		}
		for(;i<m;i++)
		{
			for(j=0;j<n;j++)
				ndat[i*n+j]=(T)0;
		}
		delete [] A.dat;
		A.dat=ndat;
		A.lb1=m1;
		A.ub1=m2;
		A.lb2=n1;
		A.ub2=n2;
	}

	inline cxscmatrix_subv<T> operator [](const int &i) throw() { return cxscmatrix_subv<T>(dat, lb2, ub2, xsize, xsize*(i-lb1),1); }
	inline cxscmatrix_subv<T> operator [](const cxscmatrix_column &i) throw() { return cxscmatrix_subv<T>(dat, lb1, ub1, ysize, i.col()-lb2, xsize); }
	inline cxscmatrix_slice<T> operator ()() throw() { return cxscmatrix_slice<T>(*this,lb1,ub1,lb2,ub2); }
	inline cxscmatrix_slice<T> operator ()(const int &m, const int &n) throw() { return cxscmatrix_slice<T>(*this,1,ub1,1,ub2); }
	inline cxscmatrix_slice<T> operator ()(const int &m1, const int &m2, const int &n1, const int &n2) throw() { return cxscmatrix_slice<T>(*this,m1,m2,n1,n2); }
	friend std::ostream &operator <<(std::ostream &s,const cxscmatrix<T> &r) throw()
	{
		int i,j;
		for (i=0;i<r.ysize;i++)
		{
			for (j=0;j<r.xsize;j++)
			{
				s << r.dat[i*r.xsize+j]<<" ";
			}
			s<<std::endl;
		}
		return s;
	}

	friend inline cxscmatrix<T> operator +(const cxscmatrix<T> &m1,const cxscmatrix<T> &m2) throw()
	{
		int x=1,y=1;

		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }
		
		cxscmatrix<T> r(y,x,y+m1.ysize,x+m1.xsize);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]+m2.dat[i];

		return r;
	}
	
	friend inline cxscmatrix<T> operator -(const cxscmatrix<T> &m1,const cxscmatrix<T> &m2) throw()
	{
		int x=1,y=1;

		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }
		
		cxscmatrix<T> r(y,x,y+m1.ysize,x+m1.xsize);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]-m2.dat[i];

		return r;
	}
	
	friend inline cxscmatrix<T> operator |(const cxscmatrix<T> &m1,const cxscmatrix<T> &m2) throw()
	{
		int x=1,y=1;

		if((m1.lb1==m2.lb1)&&(m1.lb2==m2.lb2)) { y=m1.lb1; x=m1.lb2; }
		
		cxscmatrix<T> r(y,x,y+m1.ysize,x+m1.xsize);

		// assume that the matrices are of the same size, no bounds checking here
		for(int i=0;i<m1.xsize*m1.ysize;i++)
			r.dat[i]=m1.dat[i]|m2.dat[i];

		return r;
	}

	friend inline cxscmatrix<T> operator *(const cxscmatrix<T> &m1, const cxscmatrix<T> &m2) throw()
	{
		cxscmatrix<T> r(m1.lb1,m1.ub1,m2.lb2,m2.ub2);

		for(int i=0;i<m1.ysize;i++)
		{
			for(int j=0;j<m2.xsize;j++)
			{
				for(k=0;k<m1.xsize;k++)
					accumulate(dotaccu[0],m1.dat[i*m1.xsize+k],m2.dat[k*m2.xsize+j]);
				r.dat[i*m2.xsize+j]=rnd(dotaccu[0]);
			}
		}

		return r;
	}

	friend inline cxscmatrix<T> operator *(const cxscmatrix<T> &m,const T &c) throw() { return c*m; }
//	template <class S>//,class U>
//	friend inline cxscmatrix<T> _ms_mult(const cxscmatrix<T> &m,const S &c) throw()//,const U &u) throw()
//	{
//		cxscmatrix<T> r(m.lb1,m.ub1,m.lb2,m.ub2);
//
//		for(int i=0;i<m.xsize*m.ysize;m++)
//			r.dat[i]=c*m.dat[i];
//		return r;
//	}

	friend inline cxscmatrix<T> operator /(const cxscmatrix<T> &m,const T &c) throw()
	{
		cxscmatrix<T> r(m.lb1,m.ub1,m.lb2,m.ub2);

		for(int i=0;i<m.xsize*m.ysize;m++)
			r.dat[i]=m.dat[i]/c;
		return r;
	}

	friend cxscvector<T> &cxscvector<T>::operator =(const cxscmatrix<T> &m) throw();
	
};

template <class T>
class cxscmatrix_slice
{
	friend class cxscmatrix<T>;
	private:
	T *dat;
	int lb1,lb2,ub1,ub2,xsize,ysize; // matrix size
	int xstart,xend,ystart,yend,sxsize,sysize;     // slice size

	public:
	//--------------- Konstruktoren ----------------------------------------

	explicit inline cxscmatrix_slice<T>(cxscmatrix<T> &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),lb1(a.lb1),lb2(a.lb2),ub1(a.ub1),ub2(a.ub2),xsize(u2-l2+1),ysize(u1-l1+1),xstart(l2),xend(u2),ystart(l1),yend(u1),sxsize(u2-l2+1),sysize(u1-l1+1) { }
	explicit inline cxscmatrix_slice<T>(cxscmatrix_slice<T> &a,const int &l1,const int &u1,const int &l2, const int &u2) throw():dat(a.dat),lb1(a.lb1),lb2(a.lb2),ub1(a.ub1),ub2(a.ub2),xsize(u2-l2+1),ysize(u1-l1+1),xstart(l2),xend(u2),ystart(l1),yend(u1) { }

	//---------------- Standardfunktionen -----------------------------------

	friend inline int Lb(const cxscmatrix_slice<T> &rm, const int &i) throw() { return (i==ROW?rm.ystart:rm.xstart); }
	friend inline int Ub(const cxscmatrix_slice<T> &rm, const int &i) throw() { return (i==ROW?rm.yend:rm.xend); }
	inline cxscmatrix_subv<T> operator [](const int &i) throw() { return cxscmatrix_subv<T>(dat, lb2, ub2, xsize, xsize*(i-lb1),1); }
	inline cxscmatrix_subv<T> operator [](const cxscmatrix_column &i) throw() { return cxscmatrix_subv<T>(dat, lb1, ub1, ysize, i.col()-lb2, xsize); }
	inline cxscmatrix_slice<T> operator ()() throw() { return cxscmatrix_slice<T>(*this,lb1,ub1,lb2,ub2); }
	inline cxscmatrix_slice<T> operator ()(const int &m, const int &n) throw() { return cxscmatrix_slice<T>(*this,1,ub1,1,ub2); }
	inline cxscmatrix_slice<T> operator ()(const int &m1, const int &m2, const int &n1, const int &n2) throw() { return cxscmatrix_slice<T>(*this,m1,m2,n1,n2); }
	friend std::ostream &operator <<(std::ostream &s,const cxscmatrix_slice<T> &r) throw()
	{
		int i,j;
		for (i=r.ystart;i<=r.yend;i++)
		{
			for (j=r.xstart;j<=r.xend;j++)
			{
				s << r.dat[(i-r.lb1)*r.xsize+j-r.lb2]<<" ";
			}
			s<<std::endl;
		}
		return s;
	}

	friend inline cxscmatrix<T> operator +(const cxscmatrix_slice<T> &m1,const cxscmatrix_slice<T> &m2) throw()
	{
		int x=1,y=1;

		if((m1.xstart==m2.xstart)&&(m1.ystart==m2.ystart)) { y=m1.ystart; x=m1.xstart; }
		
		cxscmatrix<T> r(y,x,y+m1.sysize-1,x+m1.sxsize-1);

		// no bounds checking here
		for(int i=0;i<m1.sysize;i++)
		{
			for(int j=0;j<m1.sxsize;j++)
				r.dat[i*m1.sxsize+j]=m1.dat[(m1.ystart-m1.lb1+i)*m1.xsize+m1.xstart-m1.lb2+j]+m2.dat[(m2.ystart-m2.lb1+i)*m2.xsize+m2.xstart-m2.lb2+j];
		}

		return r;
	}

};

// uses components of class cxscmatrix_subv, therefore cannot be in cxscvect.hpp
template <class T>
cxscvector<T>::cxscvector<T>(const cxscmatrix_subv<T> &v) throw()
{
	l=v.l;
	u=v.u;
	size=v.size;
	dat=new T[size];
	for (int i=0, j=start;i<size;i++,j+=v.offset)
		dat[i]=v.dat[j];
}

template <class T>
cxscvector<T> &cxscvector<T>::operator =(const cxscmatrix<T> &m) throw()
{
	if(m.xsize!=1&&m.ysize!=1)
		cerr<<"Falsche Indexgrenzen bei cxscvect=cxscmatr\n"<<std::endl;
	else
		if(m.ysize==1)
		{
			T *ndat=new T[m.xsize];
			for(int i=0;i<m.xsize;i++)
				ndat[i]=m.dat[i];
			l=m.lb2;
			u=m.ub2;
			size=m.xsize;
			delete [] dat;
			dat=ndat;
		}
		else
		{
			T *ndat=new T[m.ysize];
			for(int i=0;i<m.ysize;i++)
				ndat[i]=m.dat[i];
			l=m.lb1;
			u=m.ub1;
			size=m.ysize;
			delete [] dat;
			dat=ndat;
		}
	return *this;
}

template <class T>
cxscmatrix<T> &cxscmatrix<T>::operator =(const cxscvector<T> &v) throw()
{
	T *ndat=new T[v.size];
	for(int i=0;i<v.size;i++)
		ndat[i]=v.dat[i];
	delete [] dat;
	dat=ndat;
	if(xsize==1&&ysize!=1)
	{
		lb1=ub1=ysize=1;
		lb2=v.l;
		ub2=v.u;
		xsize=v.size;
	}
	else
	{
		lb1=v.l;
		ub1=v.u;
		ysize=v.size;
		lb2=ub2=xsize=1;
	}
	return *this;
}

template <class T,class S,class U>
cxscmatrix<U> _ms_mult(const cxscmatrix<S> &m,const T &s,const U &u=0) throw()
{
	cxscmatrix<U> r(m.lb1,m.ub1,m.lb2,m.ub2);

	for(int i=0;i<m.ysize*m.xsize;i++)
		r.dat[i]=m.dat[i]*s;

	return r;
}

template <class T,class S,class R>
cxscvector<R> _mv_mult(const cxscmatrix<S> &m,const cxscvector<T> &v,R) throw()
{
	cxscvector<R> r(m.lb1,m.ub1);

	for(int i=0;i<m.ysize;i++)
	{
		for(j=0;j<v.size;j++)
			accumulate(dotaccu[0],m.dat[i*m.xsize+j],v.dat[j]);
		r.dat[i]=rnd(dotaccu[0]);
	}

	return r;
}

template <class T,class S, class R>
cxscvector<T> _vm_mult(const cxscvector<T> &v,const cxscmatrix<S> &m,R) throw()
{
	cxscvector<R> r(m.lb2,m.ub2);

	for(int i=0;i<m.xsize;i++)
	{
		for(j=0;j<v.size;j++)
			accumulate(dotaccu[0],m.dat[j*m.xsize+i],v.dat[j]);
		r.dat[i]=rnd(dotaccu[0]);
	}

	return r;
}


} // namespace cxsc 


#endif

