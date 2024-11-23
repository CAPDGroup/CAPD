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

namespace cxsc {
	
	/*!
	Creation of a variable of type lx_civector with length \f$ n = 1 \f$ and index bounds \f$ lb = ub = 1 \f$. The value of the element is undefined.
	 */
	inline lx_civector::lx_civector () throw():dat(NULL),l(1),u(0),size(0)
	{ }
	
	/*!
	\param i Dimension of vector

	Creation of a variable of type lx_civector with length \f$ n = i \f$ and index bounds \f$ lb = 1 \f$, and \f$ ub = i \f$. The values of the elements are undefined.
	*/
	inline lx_civector::lx_civector(int i) throw():l(1),u(i),size(i)
	{
		dat=new lx_cinterval[i];
	}
	
		/*!
	\param i1 Starting dimension of vector
	\param i2 Ending dimension of vector

	Creation of a variable of type lx_civector with length \f$ n = i2 - i1 + 1 \f$ and index bounds \f$ lb = i1 \f$, and \f$ ub = i2 \f$. The values of the elements are undefined.
		 */
	inline lx_civector::lx_civector(int i1, int i2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_WRONG_BOUNDARIES,ERROR_IVECTOR_NO_MORE_MEMORY):
			                                            l(i1),u(i2),size(i2-i1+1)
#else
		throw():l(i1),u(i2),size(i2-i1+1)
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i1>i2) cxscthrow(ERROR_IVECTOR_WRONG_BOUNDARIES(
					"lx_civector(const int &i1,const int &i2)"));
#endif
		dat=new lx_cinterval[size];
	}
			
	inline lx_civector::lx_civector(const lx_cinterval &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}

	inline lx_civector::lx_civector(const l_cinterval &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
	
	inline lx_civector::lx_civector(const cinterval &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
	
	inline lx_civector::lx_civector(const lx_complex &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}

	inline lx_civector::lx_civector(const l_complex &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
	
	inline lx_civector::lx_civector(const complex &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
	
	inline lx_civector::lx_civector(const lx_interval &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
	
	inline lx_civector::lx_civector(const l_interval &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}	
	
	inline lx_civector::lx_civector(const interval &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
		
	inline lx_civector::lx_civector(const lx_real &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
		
	inline lx_civector::lx_civector(const l_real &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}
		
	inline lx_civector::lx_civector(const real &r) throw():l(1),u(1),size(1)
	{
		dat=new lx_cinterval[1];
		*dat=r;
	}	

	inline lx_civector::lx_civector(const lx_civector &v)
			                                   throw():l(v.l),u(v.u),size(v.size)
	{
		dat=new lx_cinterval[size];
		for (int i=0;i<size;i++)
			dat[i]=v.dat[i];
	}
	
	inline lx_civector &lx_civector::operator =(const lx_civector &rv) throw()
	{ 
		l = rv.l; u = rv.u; size = rv.size;
		dat=new lx_cinterval[size];
		for (int i=0;i<size;i++)
			dat[i]=rv.dat[i];
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const lx_cinterval &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const l_cinterval &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const cinterval &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const lx_complex &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const l_complex &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const complex &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}

	inline lx_civector &lx_civector::operator =(const lx_interval &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const l_interval &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const interval &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const lx_real &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const l_real &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_civector &lx_civector::operator =(const real &r) throw()
	{
		lx_cinterval *newdat = new lx_cinterval[size];
		for (int i=0;i<size;i++)
			newdat[i] = r;
		delete [] dat;
		dat = newdat;
		return *this;
	}
	
	inline lx_cinterval & lx_civector::operator [](const int &i)
#if(CXSC_INDEX_CHECK)
		throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
		throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l||i>u) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC(
					"lx_cinterval & lx_civector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}
		
	inline const lx_cinterval & lx_civector::operator [](const int &i) const
#if(CXSC_INDEX_CHECK)
			throw(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC)
#else
			throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(i<l || i>u) cxscthrow(ERROR_IVECTOR_ELEMENT_NOT_IN_VEC(
			"lx_cinterval & lx_civector::operator [](const int &i)"));
#endif
		return dat[i-l];
	}	

inline void Resize(lx_civector &rv, int len)
#if(CXSC_INDEX_CHECK)
			throw(ERROR__WRONG_BOUNDARIES<lx_civector>);
#else
	throw()
#endif
	{
		if (rv.size == len)
			SetLb(rv,1);
		else
		{
#if(CXSC_INDEX_CHECK)	
			if (len<0) 
				cxscthrow(ERROR__WRONG_BOUNDARIES(
							 "Resize(lx_civector &rv, int len)"));
#endif
			lx_cinterval *ndat = new lx_cinterval[len];
			int beg, end;
			beg = (rv.l>1)? rv.l : 1;
			end = (rv.u<len)? rv.u : len;
			for(int i=beg-1;i<end;i++)
				ndat[i]=rv.dat[i-rv.l+1];
			delete [] rv.dat;
			rv.dat=ndat;
			rv.size=rv.u=len;
			rv.l=1;
		}
	}
	
inline void Resize(lx_civector &rv, int lb, int ub)
#if(CXSC_INDEX_CHECK)
			throw(ERROR__WRONG_BOUNDARIES<lx_civector>)
#else
			throw()
#endif
{ 
	if (rv.size == ub-lb+1)
		SetUb(rv,ub);
	else
	{
		rv.size = ub-lb+1;
#if(CXSC_INDEX_CHECK)
	if (rv.size<0)
		cxscthrow(ERROR__WRONG_BOUNDARIES(
					 "Resize(lx_civector &rv, int lb, int ub)"));
#endif
		lx_cinterval *ndat = new lx_cinterval[rv.size];
		int beg, end;
		beg = (rv.l>lb)? rv.l : lb;
		end = (rv.u<ub)? rv.u : ub;
		for(int i=0;i<=rv.size-1;i++)
			ndat[i]=interval(0);
		for(int i=beg;i<=end;i++)
			ndat[i-lb]=rv.dat[i-rv.l];
		delete [] rv.dat;
		rv.dat=ndat;
		rv.l=lb;
		rv.u=ub;
	}
}

inline void DoubleSize(lx_civector& x) throw()
{
	int n = Lb(x);
	Resize(x,n,2*Ub(x)-n+1);
}

} // namespace cxsc
