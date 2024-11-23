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

/* CVS $Id: cxscvect.hpp,v 1.23 2014/01/30 17:23:44 cxsc Exp $ */

#ifndef _CXSC_CXSCVECT_HPP_INCLUDED
#define _CXSC_CXSCVECT_HPP_INCLUDED

#include "dot.hpp"
#include "real.hpp" // used for declaration of Inf, Sup,...
//#include "cxscmatr.hpp"
#include <iostream>

namespace cxsc {

template <class T>
class cxscmatrix_subv;

template <class T>
class cxscmatrix;

template <class T>
class cxscvector_slice;

template <class T>
class cxscvector
{
	friend class cxscvector_slice<T>;
	private:
	T *dat;
	int l,u,size;

	public:
	//------ Konstruktoren ----------------------------------------------------
	inline cxscvector::cxscvector<T> () throw()
	{
		dat=NULL;
		size=u=0;
		l=1;
	}

	explicit inline cxscvector<T>(const int &i) throw()
	{
		dat=new T[i];
		l=1;
		size=u=i;
	}

	explicit inline cxscvector<T>(const int &i1,const int &i2) throw()
	{
		l=i1;
		u=i2;
		size=u-l+1;

		dat=new T[size];
	}

	inline cxscvector<T>(const cxscvector_slice<T> &rs) throw()
	{
		T *from=rs.dat;
		
		l=rs.start;
		u=rs.end;
		size=u-l+1;
		
		dat=new T[size];
		for(int i=0, to=l-rs.l;i<size;i++,to++)
			dat[i]=from[to];
	}

	inline cxscvector<T>(const cxscvector<T> &v) throw()
	{
		u=size=v.size;
		l=1;
		
		dat=new T[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}

	cxscvector<T>(const cxscmatrix_subv<T> &v) throw(); // defined in cxscmatr.hpp
	
	inline cxscvector<T> &cxscvector<T>::operator =(const cxscvector<T> &rv) throw()
	{
		T *ndat=new T[rv.size];
		
		for (int i=0;i<rv.size;i++)
			ndat[i]=rv.dat[i];

		delete [] dat;
		dat=ndat;

		return *this;
	}

	inline cxscvector<T> &operator =(const T &r) throw()
	{
		for (int i=0;i<size;i++)
			dat[i]=r;

	return *this;
	}

	cxscvector<T> &operator =(const cxscmatrix<T> &) throw();

	//--------- Destruktor ----------------------------------------------------
	inline ~cxscvector() { delete [] dat; }

	//------ Standardfunktionen -----------------------------------------------
	
//	friend cxscvector<real> Inf<T>(const cxscvector<T> &) throw();
	friend inline int Lb(const cxscvector<T> &rv) throw() { return rv.l; }
	friend inline int Ub(const cxscvector<T> &rv) throw() { return rv.u; }
	friend inline void SetLb(cxscvector<T> &rv, const int &l) throw() { rv.l=l; rv.u=l+rv.size-1; }
	friend inline void SetUb(cxscvector<T> &rv, const int &u) throw() { rv.u=u; rv.l=u-rv.size+1; }
	friend inline void Resize(cxscvector<T> &rv) throw()
	{
		rv.size=rv.u=0;
		rv.l=1;
		delete [] rv.dat;
		rv.dat=NULL;
	}

	friend inline void Resize(cxscvector<T> &rv, const int &len) throw()
	{
		if(rv.size==len)
			SetLb(rv,1);
		else
		{
			T *ndat=new T[len];
			int beg,end;
			beg=(rv.l>1)?rv.l:1;
			end=(rv.u<len)?rv.u:len;
			for(int i=beg-1;i<end;i++)
				ndat[i]=rv.dat[i-rv.l+1];
			delete [] rv.dat;
			rv.dat=ndat;
			rv.size=rv.u=len;
			rv.l=1;
		}
	}

	friend void Resize(cxscvector<T> &rv, const int &lb, const int &ub) throw()
	{
		if(rv.size==ub-lb+1)
			SetUb(rv,ub);
		else
		{		
			rv.size=ub-lb+1;
			T *ndat=new T[rv.size];
			int beg,end;
			beg=(rv.l>lb)?rv.l:lb;
			end=(rv.u<ub)?rv.u:ub;
			for(int i=beg;i<=end;i++)
				ndat[i-lb]=rv.dat[i-rv.l];
			delete [] rv.dat;
			rv.dat=ndat;
			rv.l=lb;
			rv.u=ub;
		}
	}

	inline T & operator [](const int &i) throw() { return dat[i-l]; }
	inline cxscvector<T> & operator ()() throw() { return *this; }
	inline cxscvector_slice<T> operator ()(const int &) throw() { return cxscvector_slice<T>(*this,1,i); }
	friend inline void accumulate(dotprecision &dp, const cxscvector<T> & rv1, const cxscvector<T> &rv2) throw()
	{
		for(int i=0;i<rv1.size;i++)
			accumulate(dp,rv1.dat[i],rv2.dat[i]);
	}

	friend inline void accumulate(dotprecision &dp, const cxscvector_slice<T> & sl, const cxscvector<T> &rv) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dp,sl.dat[j],rv.dat[i]);
	}

	friend inline void accumulate(dotprecision &dp, const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dp,sl.dat[j],rv.dat[i]);
	}

	friend inline cxscvector<T> operator /(const cxscvector<T> &rv, const T &s) throw()
	{
		cxscvector<T> p(rv.l,rv.u);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[i]/s;
		
		return p;
	}

	friend inline cxscvector<T> operator /(const cxscvector_slice<T> &sl, const T &s) throw()
	{
		cxscvector<T> p(sl.l,sl.u);

		for(int i=sl.start-sl.l;i<sl.size;i++)
			p.dat[i]=sl.dat[i]/s;
		
		return p;
	}

	friend inline cxscvector<T> operator *(const cxscvector<T> &rv, const T &s) throw()
	{
		cxscvector<T> p(rv.l,rv.u);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[i]*s;
		
		return p;
	}

	friend inline cxscvector<T> operator *(const T &s, const cxscvector<T> &rv) throw()
	{
		cxscvector<T> p(rv.l,rv.u);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[i]*s;
		
		return p;
	}

	friend inline cxscvector<T> operator *(const cxscvector_slice<T> &sl, const T &s) throw()
	{
		cxscvector<T> p(sl.l,sl.u);

		for(int i=sl.start-sl.l;i<sl.size;i++)
			p.dat[i]=sl.dat[i]*s;
		
		return p;
	}

	friend inline T operator *(const cxscvector<T> & rv1, const cxscvector<T> &rv2) throw()
	{
		dotakku[0]=0;
		for(int i=0;i<rv1.size;i++)
			accumulate(dotakku[0],rv1.dat[i],rv2.dat[i]);

		return rnd(dotakku[0]);
	}

	friend inline T operator *(const cxscvector_slice<T> & sl, const cxscvector<T> &rv) throw()
	{
		dotakku[0]=0;
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dotakku[0],sl.dat[j],rv.dat[i]);

		return rnd(dotakku[0]);
	}

	friend inline T operator *(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		dotakku[0]=0;
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dotakku[0],sl.dat[j],rv.dat[i]);

		return rnd(dotakku[0]);
	}

	friend inline cxscvector<T> operator +(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		cxscvector<T> sum(rv1.l,rv1.u);

		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[i]+rv2.dat[i];
		return sum;
	}

	friend inline cxscvector<T> operator +(const cxscvector<T> &rv,const cxscvector_slice<T> &sl) throw()
	{
		cxscvector<T> sum(rv.l,rv.u);

		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]+sl.dat[j];
		return sum;
	}

	friend inline cxscvector<T> operator +(const cxscvector_slice<T> &sl,const cxscvector<T> &rv) throw()
	{
		cxscvector<T> sum(rv.l,rv.u);

		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]+sl.dat[j];
		return sum;
	}

	friend inline cxscvector<T> operator +(const cxscvector<T> &rv) throw() { return rv; }
	friend inline cxscvector<T> operator +(const cxscvector_slice<T> &s1,const cxscvector_slice<T> &s2) throw()
	{
		cxscvector<T> sum(s1.l,s1.u);

		for(int i=s1.start-s1.l,j=s2.start-s2.l;i<s1.end;i++,j++)
			sum.dat[i]=s1.dat[i]+s2.dat[j];
		return sum;
	}

	friend inline cxscvector<T> & operator +=(cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]+=rv2.dat[i];
		return rv1;
	}

	friend inline cxscvector<T> & operator +=(cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]+=sl.dat[j];
		return rv;
	}

	friend inline cxscvector_slice<T> & operator +=(cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]+=rv.dat[i];
		return sl;
	}

	friend inline cxscvector<T> & operator -=(cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]-=rv2.dat[i];
		return rv1;
	}

	friend inline cxscvector<T> & operator -=(cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]-=sl.dat[j];
		return rv;
	}

	friend inline cxscvector_slice<T> & operator -=(cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]-=rv.dat[i];
		return sl;
	}

	friend inline cxscvector<T> operator -(const cxscvector<T> &rv) throw()
	{
		cxscvector<T> sum(rv.l,rv.u);

		for (int i=0;i<rv.size;i++)
			sum.dat[i]= -rv.dat[i];

		return sum;
	}

	friend inline cxscvector<T> operator -(const cxscvector_slice<T> &sl) throw()
	{
		cxscvector<T> sum(sl.l,sl.u);

		for (int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			sum.dat[i]= -sl.dat[j];

		return sum;
	}

	friend inline cxscvector<T> operator -(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		cxscvector<T> sum(rv1.l,rv1.u);

		for (int i=0;i<rv1.size;i++)
				sum.dat[i]=rv1.dat[i]-rv2.dat[i];
		return sum; 
	}

	friend inline cxscvector<T> operator -(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		cxscvector<T> sum(rv.l,rv.u);

		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]-sl.dat[j];

		return sum;
	}
	
	friend inline cxscvector<T> operator -(const cxscvector_slice<T> &sl,const cxscvector<T> &rv) throw()
	{
		cxscvector<T> sum(sl.l,sl.u);

		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=sl.dat[j]-rv.dat[i];
		return sum;
	}

	friend inline cxscvector<T> operator |(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		cxscvector<T> sum(rv1.l,rv1.u);

		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[i]|rv2.dat[i];
		return sum;
	}

	friend inline cxscvector<T> operator |(const cxscvector<T> &rv,const cxscvector_slice<T> &sl) throw()
	{
		cxscvector<T> sum(rv.l,rv.u);

		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]|sl.dat[j];
		return sum;
	}

	friend inline cxscvector<T> operator |(const cxscvector_slice<T> &sl,const cxscvector<T> &rv) throw()
	{
		cxscvector<T> sum(rv.l,rv.u);

		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]|sl.dat[j];
		return sum;
	}

	friend inline cxscvector<T> & operator |=(cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]|=rv2.dat[i];
		return rv1;
	}

	friend inline cxscvector<T> & operator |=(cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]|=sl.dat[j];
		return rv;
	}

	friend inline cxscvector_slice<T> & operator |=(cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]|=rv.dat[i];
		return sl;
	}

	friend inline bool operator ==(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);

		int i;
		for (i=0;i<rv1.size && rv1.dat[i]==rv2.dat[i];i++);
		
		return (i==rv1.size);
	}	

	friend inline bool operator ==(const cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);

		int i,j;
		for (i=0,j=sl.start-sl.l;i<rv.size && sl.dat[j]==rv.dat[i];i++,j++);
		
		return (i==rv.size);
	}	

	friend inline bool operator ==(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);

		int i,j;
		for (i=0,j=sl.start-sl.l;i<rv.size && sl.dat[j]==rv.dat[i];i++,j++);
		
		return (i==rv.size);
	}	

	friend inline bool operator !=(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(true);

		int i;
		for (i=0;i<rv1.size && rv1.dat[i]!=rv2.dat[i];i++);
		
		return (i!=rv1.size);
	}	

	friend inline bool operator !=(const cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		if(sl.size!=rv.size)
			return(true);

		int i,j;
		for (i=0,j=sl.start-sl.l;i<rv.size && sl.dat[j]!=rv.dat[i];i++,j++);
		
		return (i!=rv.size);
	}	

	friend inline bool operator !=(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		if(sl.size!=rv.size)
			return(true);

		int i,j;
		for (i=0,j=sl.start-sl.l;i<rv.size && sl.dat[j]!=rv.dat[i];i++,j++);
		
		return (i!=rv.size);
	}	

	friend inline bool operator <(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);
		
		int i;
		for(i=0;i<rv1.size&&(rv1.dat[i]<rv2.dat[i]);i++);

		return (i!=rv1.size);
	}

	friend inline bool operator <(const cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(sl.dat[i]<rv.dat[j]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator <(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(rv.dat[j]<sl.dat[i]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator <=(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);
		
		int i;
		for(i=0;i<rv1.size&&(rv1.dat[i]<=rv2.dat[i]);i++);

		return (i!=rv1.size);
	}

	friend inline bool operator <=(const cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(sl.dat[i]<=rv.dat[j]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator <=(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(rv.dat[j]<=sl.dat[i]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator >(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);
		
		int i;
		for(i=0;i<rv1.size&&(rv1.dat[i]>rv2.dat[i]);i++);

		return (i!=rv1.size);
	}

	friend inline bool operator >(const cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(sl.dat[i]>rv.dat[j]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator >(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(rv.dat[j]>sl.dat[i]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator >=(const cxscvector<T> &rv1, const cxscvector<T> &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);
		
		int i;
		for(i=0;i<rv1.size&&(rv1.dat[i]>=rv2.dat[i]);i++);

		return (i!=rv1.size);
	}

	friend inline bool operator >=(const cxscvector_slice<T> &sl, const cxscvector<T> &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(sl.dat[i]>=rv.dat[j]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator >=(const cxscvector<T> &rv, const cxscvector_slice<T> &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;i<sl.size&&(rv.dat[j]>=sl.dat[i]);i++,j++);

		return (i!=sl.size);
	}

	friend inline bool operator !(const cxscvector<T> &rv) throw()
	{
		int i;
		for(i=0;i<rv.size;i++)
		{
			if(rv.dat[i])
				return(false);
		}

		return true;
	}

	inline cxscvector<T>::operator void*() throw()
	{
		for(int i=0;i<size;i++)
		{
			if(dat[i]!=0.0)
				return (void *)1;
		}

		return (void *)0;
	}

	friend inline cxscvector<T> abs(const cxscvector<T> &rv) throw()
	{
		cxscvector<T> erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=abs(rv.dat[i]);

		return erg;
	}

	friend inline cxscvector<T> abs(const cxscvector_slice<T> &sl) throw()
	{
		cxscvector<T> erg(sl.l,sl.u);
		
		for(int i=0;i<sl.size;i++)
			erg.dat[i]=abs(sl.dat[i]);

		return erg;
	}

	friend std::ostream &operator <<(std::ostream &s, const cxscvector<T> &rv) throw()
	{
		for(int i=0;i<rv.size;i++)
			s<<rv.dat[i]<<std::endl;
		return s;
	}

	friend std::istream &operator >>(std::istream &s, cxscvector<T> &rv) throw()
	{
		for(int i=0;i<rv.size;i++)
			s>>rv.dat[i];
		return s;
	}
};


template <class T>
class cxscvector_slice
{
	friend class cxscvector<T>;
	private:
	T *dat;
	int l,u,size;
	int start,end;

	public:
	
	//--------------------- Konstruktoren -----------------------------------
	explicit inline cxscvector_slice<T>(cxscvector<T> &a, const int &l, const int &u) throw():dat(a.dat),l(a.l),u(a.u),size(u-l+1),start(l),end(u) { }
	inline cxscvector_slice<T> & operator =(const cxscvector_slice<T> &sl) throw()
	{
		for (int i=start-lb, j=sl.start-sl.lb;i<end;i++,j++)
			dat[i]=sl.dat[j];
		return *this;
	}
	
	inline cxscvector_slice<T> & operator =(const cxscvector<T> &rv) throw()
	{
		for (int i=start-lb, j=0;i<end;i++,j++)
			dat[i]=rv.dat[j];
		return *this;
	}

	inline cxscvector_slice<T> & operator =(const T &r) throw()
	{
		for (int i=start-lb;i<end;i++)
			dat[i]=r;
		return *this;
	}


	//--------------------- Standardfunktionen ------------------------------
	friend inline int Lb(const cxscvector_slice<T> &sl) throw() { return sl.start; }
	friend inline int Ub(const cxscvector_slice<T> &sl) throw() { return sl.end; }
	inline T & operator [](const int &i) throw() { return dat[i-l]; }
	friend inline void accumulate<T>(dotprecision &,const cxscvector_slice<T> &,const cxscvector<T> &) throw();
	friend inline void accumulate<T>(dotprecision &,const cxscvector<T> &,const cxscvector_slice<T> &) throw();
	friend inline void accumulate(dotprecision &dp, const cxscvector_slice<T> & sl1, const cxscvector_slice<T> &sl2) throw()
	{
		for(int i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size;i++,j++)
			accumulate(dp,sl1.dat[i],sl2.dat[j]);
	}

	friend inline cxscvector<T> operator /<T>(const cxscvector_slice<T> &sl, const T &s) throw();
	friend inline cxscvector<T> operator *<T>(const T &s, const cxscvector_slice<T> &sl) throw();
	friend inline cxscvector<T> operator *<T>(const cxscvector_slice<T> &sl, const T &s) throw();
//	{
//		cxscvector<T> p(sl.l,sl.u);
//
//		for(int i=sl.start-sl.l;i<sl.size;i++)
//			p.dat[i]=sl.dat[i]*s;
//		
//		return p;
//	}

	friend inline T operator *(const cxscvector_slice<T> & sl1, const cxscvector_slice<T> &sl2) throw()
	{
		dotakku[0]=0;
		for(int i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size;i++,j++)
			accumulate(dotakku[0],sl1.dat[i],sl2.dat[j]);

		return rnd(dotakku[0]);
	}

	friend inline T operator *<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline T operator *<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector<T> operator +<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector<T> operator +<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline cxscvector<T> operator +<T>(const cxscvector_slice<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector<T> operator +(const cxscvector_slice<T> &sl) throw() { return sl; }
	friend inline cxscvector<T> &operator +=<T>(cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector_slice<T> &operator +=<T>(cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline cxscvector_slice<T> & operator +=(cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		for(int i=sl1.start-sl1.l,j=sl2.start-sl1.l;i<sl1.end;i++,j++)
			sl1.dat[i]+=sl2.dat[j];
		return sl1;
	}

	friend inline cxscvector<T> &operator -=<T>(cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector_slice<T> &operator -=<T>(cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline cxscvector_slice<T> & operator -=(cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		for(int i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.end;i++,j++)
			sl1.dat[i]-=sl2.dat[j];
		return sl1;
	}

	friend inline cxscvector<T> operator -<T>(const cxscvector_slice<T> &) throw();
	friend inline cxscvector<T> operator -<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector<T> operator -<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline cxscvector<T> operator |<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector<T> operator |<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline cxscvector<T> &operator |=<T>(cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline cxscvector_slice<T> &operator |=<T>(cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline cxscvector_slice<T> & operator |=(cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		for(int i=sl1.start-sl1.l,j=sl2.start-sl1.l;i<sl1.end;i++,j++)
			sl1.dat[i]|=sl2.dat[j];
		return sl1;
	}

	friend inline bool operator ==(const cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);

		int i,j;
		for (i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size && sl1.dat[i]==sl2.dat[i];i++);
		
		return (i==sl1.size);
	}	

	friend inline bool operator ==<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline bool operator ==<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline bool operator !=(const cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(true);

		int i,j;
		for (i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size && sl1.dat[i]!=sl2.dat[i];i++);
		
		return (i!=sl1.size);
	}	

	friend inline bool operator !=<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline bool operator !=<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline bool operator <(const cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);
		
		int i,j;
		for(i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size&&(sl1.dat[i]<sl2.dat[j]);i++,j++);

		return (i!=sl1.size);
	}

	friend inline bool operator < <T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline bool operator < <T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline bool operator <=(const cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);
		
		int i,j;
		for(i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size&&(sl1.dat[i]<=sl2.dat[j]);i++,j++);

		return (i!=sl1.size);
	}

	friend inline bool operator <=<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline bool operator <=<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline bool operator >(const cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);
		
		int i,j;
		for(i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size&&(sl1.dat[i]>sl2.dat[j]);i++,j++);

		return (i!=sl1.size);
	}

	friend inline bool operator ><T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline bool operator ><T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline bool operator >=(const cxscvector_slice<T> &sl1, const cxscvector_slice<T> &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);
		
		int i,j;
		for(i=sl1.start-sl1.l,j=sl2.start-sl2.l;i<sl1.size&&(sl1.dat[i]>=sl2.dat[j]);i++,j++);

		return (i!=sl1.size);
	}

	friend inline bool operator >=<T>(const cxscvector_slice<T> &, const cxscvector<T> &) throw();
	friend inline bool operator >=<T>(const cxscvector<T> &, const cxscvector_slice<T> &) throw();
	friend inline bool operator !(const cxscvector_slice<T> &sl) throw()
	{
		int i;
		for(i=sl.start-sl.l;i<sl.size;i++)
		{
			if(sl.dat[i])
				return(false);
		}

		return true;
	}

	inline cxscvector_slice<T>::operator void*() throw()
	{
		for(int i=start-l;i<size;i++)
		{
			if(dat[i]!=0.0)
				return (void *)1;
		}

		return (void *)0;
	}

	friend inline cxscvector<T> abs<T>(const cxscvector_slice<T> &) throw();
};

} // namespace cxsc 

#endif

