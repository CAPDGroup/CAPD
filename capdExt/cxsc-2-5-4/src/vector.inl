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

/* CVS $Id: vector.inl,v 1.28 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_VECTOR_INL_INCLUDED
#define _CXSC_VECTOR_INL_INCLUDED

#include "except.hpp"

#ifdef CXSC_USE_BLAS
#include "cxsc_blas.hpp"
#endif

namespace cxsc {
int abs(int a);

	template <class V>
	TINLINE void _vresize(V &rv) throw()
	{
		rv.size=rv.u=0;
		rv.l=1;
		delete [] rv.dat;
		rv.dat=NULL;
	}

	template <class V,class S>
	TINLINE void _vresize(V &rv, const int &len)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>)
#else
	throw()
#endif
	{
		if(rv.size==len)
			SetLb(rv,1);
		else
		{
#if(CXSC_INDEX_CHECK)
			if(len<0) cxscthrow(ERROR__WRONG_BOUNDARIES<V>(" Resize("+nameof(rv)+" &, const int &)"));
#endif
			S *ndat=new S[len];
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

	template <class V,class S>
	TINLINE void _vresize(V &rv, const int &lb, const int &ub)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__WRONG_BOUNDARIES<V>)
#else
	throw()
#endif
	{
		if(rv.size==ub-lb+1)
			SetUb(rv,ub);
		else
		{		
			rv.size=ub-lb+1;
#if(CXSC_INDEX_CHECK)
			if(rv.size<0) cxscthrow(ERROR__WRONG_BOUNDARIES<V>("void Resize("+nameof(rv)+" &, const int &, const int &)"));
#endif
			S *ndat=new S[rv.size];
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

	template <class V1,class V2,class S>
	TINLINE V1 &_vvassign(V1 &rv1,const V2 &rv2) throw()
	{
		S *ndat=new S[rv2.size];
		rv1.l=rv2.l;
		rv1.u=rv2.u;
		rv1.size=rv2.size;
		
		for (int i=0;i<rv2.size;i++)
			ndat[i]=rv2.dat[i];
		delete [] rv1.dat;
		rv1.dat=ndat;

		return rv1;
	}

	template <class V,class S>
	TINLINE V & _vsassign(V &rv,const S &r) throw()
	{
		for (int i=0;i<rv.size;i++)
			rv.dat[i]=r;

		return rv;
	}

	template <class VS1,class VS2>
	TINLINE VS1 & _vsvsassign(VS1 &sl1,const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" "+nameof(sl1)+"::operator =(const "+nameof(sl2)+" &)"));
#endif
		for (int i=sl1.start-sl1.l, j=sl2.start-sl2.l;i<=sl1.end-sl1.l;i++,j++)
			sl1.dat[i]=sl2.dat[j];
		return sl1;
	}
	
	template <class V,class VS,class S>
	TINLINE V & _vvsassign(V &rv,const VS &sl) throw()
	{
		S *ndat=new S[sl.size];
		rv.l=sl.start;
		rv.u=sl.end;
		rv.size=sl.size;
		for (int i=sl.start-sl.l, j=0;j<rv.size;i++,j++)
			ndat[j]=sl.dat[i];
		delete [] rv.dat;
		rv.dat=ndat;
		return rv;
	}

	template <class VS,class V>
	TINLINE VS & _vsvassign(VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" "+nameof(sl)+"::operator =(const "+nameof(rv)+" &)"));
#endif
		for (int i=sl.start-sl.l, j=0;j<sl.size;i++,j++)
			sl.dat[i]=rv.dat[j];
		return sl;
	}

	template <class VS,class S>
	TINLINE VS & _vssassign(VS &sl,const S &r) throw()
	{
		for (int i=sl.start-sl.l;i<=sl.end-sl.l;i++)
			sl.dat[i]=r;
		return sl;
	}

	template <class DP,class V1,class V2>
	TINLINE void _vvaccu(DP &dp, const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(OP_WITH_WRONG_DIM("void accumulate("+nameof(dp)+" &, const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			accumulate(dp,rv1.dat[i],rv2.dat[i]);
	}

	template <class DP,class VS,class V>
	TINLINE void _vsvaccu(DP &dp, const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(OP_WITH_WRONG_DIM("void accumulate("+nameof(dp)+" &, const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dp,sl.dat[j],rv.dat[i]);
	}

	template <class V,class S,class E>
	TINLINE E _vsdiv(const V &rv, const S &s) throw()
	{
		E p(rv.l,rv.u);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[i]/s;
		
		return p;
	}

	template <class V,class S>
	TINLINE V &_vsdivassign(V &rv,const S &r) throw()
	{
		for(int i=0;i<rv.size;i++)
			rv.dat[i]/=r;
		return rv;
	}

	template <class VS,class S,class E>
	TINLINE E _vssdiv(const VS &sl, const S &s) throw()
	{
		E p(sl.start,sl.end);

		for(int i=sl.start-sl.l;i<sl.size;i++)
			p.dat[i]=sl.dat[i]/s;
		
		return p;
	}

	template <class V,class S,class E>
	TINLINE E _vsmult(const V &rv, const S &s) throw()
	{
		E p(rv.l,rv.u);

		for(int i=0;i<rv.size;i++)
			p.dat[i]=rv.dat[i]*s;
		
		return p;
	}

	template <class VS,class S,class E>
	TINLINE E _vssmult(const VS &sl, const S &s) throw()
	{
		E p(sl.start,sl.end);

		for(int i=sl.start-sl.l,k=0;k<sl.size;i++,k++)
			p.dat[k]=sl.dat[i]*s;
		
		return p;
	}

	template <class V1,class V2,class E>
	TINLINE E _vvlmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(E())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		dotprecision dot(0.0);
                accumulate(dot,rv1,rv2);
/*		  for(int i=0;i<rv1.size;i++)
			accumulate(dotakku[0],rv1.dat[i],rv2.dat[i]);*/

		return E(dot);
	}

	template <class VS,class V,class E>
	TINLINE E _vsvlmult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(E())+" operator *(const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		dotprecision dot(0.0);
                accumulate(dot,sl,rv);
/*		   for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dotakku[0],sl.dat[j],rv.dat[i]);*/

		return E(dot);
	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvslmult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(E())+" operator *(const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		dotprecision dot(0.0);
                accumulate(dot,sl1,sl2);
		  /*for(int i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size;i++,j++,k++)
			accumulate(dotakku[0],sl1.dat[i],sl2.dat[j]);*/

		return E(dot);
	}

	template <class V1,class V2,class E>
	TINLINE E _vvlimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(E())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		idotprecision idot(0.0);
		accumulate(idot,rv1,rv2);

		return E(idot);
	}

	template <class VS,class V,class E>
	TINLINE E _vsvlimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(E())+" operator *(const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		idotprecision idot(0.0);
		accumulate(idot,sl,rv);

		return E(idot);
	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvslimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(E())+" operator *(const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		idotprecision idot(0.0);
		accumulate(idot,sl1,sl2);

		return E(idot);
	}

	template <class V1,class V2,class E>
	TINLINE E _vvmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(E())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif

#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(rv1,rv2,ret);
			return ret;
		} 
                else
#endif
		{
			dotprecision dot(0.0);
                	dot.set_k(opdotprec);
                	accumulate_approx(dot,rv1,rv2);
			return rnd(dot);
		}
/*		   for(int i=0;i<rv1.size;i++)
			accumulate(dotakku[0],rv1.dat[i],rv2.dat[i]);*/

	}

	template <class VS,class V,class E>
	TINLINE E _vsvmult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(E())+" operator *(const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl,rv,ret);
			return ret;
		} 
                else
#endif
		{
			dotprecision dot(0.0);
        	        dot.set_k(opdotprec);
                	accumulate_approx(dot,sl,rv);
			return rnd(dot);
		}
		  /*for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			accumulate(dotakku[0],sl.dat[j],rv.dat[i]);*/

	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsmult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(E())+" operator *(const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl1,sl2,ret);
			return ret;
		} 
                else
#endif
		{
			dotprecision dot(0.0);
                	dot.set_k(opdotprec);
                	accumulate_approx(dot,sl1,sl2);
			return rnd(dot);
		}
/*		  for(int i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size;i++,j++,k++)
			accumulate(dotakku[0],sl1.dat[i],sl2.dat[j]);*/

	}

	template <class V1,class V2,class E>
	TINLINE E _vvimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(E())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(rv1,rv2,ret);
			return ret;
		} 
                else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
			accumulate(idot,rv1,rv2);
			return rnd(idot);
		}

	}

	template <class VS,class V,class E>
	TINLINE E _vsvimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(E())+" operator *(const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl,rv,ret);
			return ret;
		} 
                else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
			accumulate(idot,sl,rv);
			return rnd(idot);
		}

	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(E())+" operator *(const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl1,sl2,ret);
			return ret;
		} 
                else
#endif
		{
			idotprecision idot(0.0);
			idot.set_k(opdotprec);
			accumulate(idot,sl1,sl2);
			return rnd(idot);
		}

	}

	template <class V1,class V2,class E>
	TINLINE E _vvcmult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(E())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(rv1,rv2,ret);
			return ret;
		} 
                else
#endif
		{
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
			accumulate_approx(cdot,rv1,rv2);
			return rnd(cdot);
		}

	}

	template <class VS,class V,class E>
	TINLINE E _vsvcmult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(E())+" operator *(const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl,rv,ret);
			return ret;
		} 
                else
#endif
		{
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
			accumulate_approx(cdot,sl,rv);
			return rnd(cdot);
		}

	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvscmult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(E())+" operator *(const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl1,sl2,ret);
			return ret;
		} 
                else
#endif
		{
			cdotprecision cdot(0.0);
			cdot.set_k(opdotprec);
			accumulate_approx(cdot,sl1,sl2);
			return rnd(cdot);
		}

	}

	template <class V1,class V2,class E>
	TINLINE E _vvcimult(const V1 & rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(E())+" operator *(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(rv1,rv2,ret);
			return ret;
		} 
                else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
			accumulate(cidot,rv1,rv2);
			return rnd(cidot);
		}

	}

	template <class VS,class V,class E>
	TINLINE E _vsvcimult(const VS & sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl.size!=rv.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(E())+" operator *(const "+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl,rv,ret);
			return ret;
		} 
                else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
			accumulate(cidot,sl,rv);
			return rnd(cidot);
		}

	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvscimult(const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(E())+" operator *(const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
#ifdef CXSC_USE_BLAS
                if(opdotprec == 1) 
		{
			E ret;
			blasdot(sl1,sl2,ret);
			return ret;
		} 
                else
#endif
		{
			cidotprecision cidot(0.0);
			cidot.set_k(opdotprec);
			accumulate(cidot,sl1,sl2);
			return rnd(cidot);
		}

	}

	template <class V,class S>
	TINLINE V &_vsmultassign(V &rv,const S &r) throw()
	{
		for(int i=0;i<rv.size;i++)
			rv.dat[i]*=r;
		return rv;
	}

	template <class VS,class S>
	TINLINE VS &_vssmultassign(VS &rv,const S &r) throw()
	{
		for(int i=rv.start-rv.l;i<=rv.end-rv.l;i++)
			rv.dat[i]*=r;
		return rv;
	}

	template <class VS,class S>
	TINLINE VS &_vssdivassign(VS &rv,const S &r) throw()
	{
		for(int i=rv.start-rv.l;i<=rv.end-rv.l;i++)
			rv.dat[i]/=r;
		return rv;
	}

	template <class V1,class V2,class E>
	TINLINE E _vvplus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
		E sum(rv1.l,rv1.u);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" operator +(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[i]+rv2.dat[i];
		return sum;
	}

	template <class V,class VS,class E>
	TINLINE E _vvsplus(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
		E sum(rv.l,rv.u);

#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(sum)+" operator +(const "+nameof(rv)+" &,const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]+sl.dat[j];
		return sum;
	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsplus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
		E sum(s1.start,s1.end);

#if(CXSC_INDEX_CHECK)
		if(s1.size!=s2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sum)+" operator +(const "+nameof(s1)+" &,const "+nameof(s2)+" &)"));
#endif
		for(int i=s1.start-s1.l,j=s2.start-s2.l,k=0;k<s1.size;i++,j++,k++)
			sum.dat[k]=s1.dat[i]+s2.dat[j];
		return sum;
	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsminus(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
		E sum(s1.start,s1.end);

#if(CXSC_INDEX_CHECK)
		if(s1.size!=s2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sum)+" operator -(const "+nameof(s1)+" &,const "+nameof(s2)+" &)"));
#endif
		for(int i=s1.start-s1.l,j=s2.start-s2.l,k=0;k<s1.size;i++,j++,k++)
			sum.dat[k]=s1.dat[i]-s2.dat[j];
		return sum;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvplusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" & operator +=("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]+=rv2.dat[i];
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvsplusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" & operator +=("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]+=sl.dat[j];
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvplusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" & operator +=("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]+=rv.dat[i];
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsplusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" & operator +=("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			sl1.dat[j]+=sl2.dat[k];
		return sl1;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsminusassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" & operator +=("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			sl1.dat[j]-=sl2.dat[k];
		return sl1;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvminusassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" & operator -=("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]-=rv2.dat[i];
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvsminusassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" & operator -=("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]-=sl.dat[j];
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvminusassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" & operator -=("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]-=rv.dat[i];
		return sl;
	}

	template <class V>
	TINLINE V _vminus(const V &rv) throw()
	{
		V sum(rv.l,rv.u);

		for (int i=0;i<rv.size;i++)
			sum.dat[i]= -rv.dat[i];

		return sum;
	}

	template <class VS,class V>
	TINLINE V _vsminus(const VS &sl) throw()
	{
		V sum(sl.start,sl.end);

		for (int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			sum.dat[i]= -sl.dat[j];

		return sum;
	}

	template <class V1,class V2,class E>
	TINLINE E _vvminus(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.l,rv1.u);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator -(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
				sum.dat[i]=rv1.dat[i]-rv2.dat[i];
		return sum; 
	}

	template <class V,class VS,class E>
	TINLINE E _vvsminus(const V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv.l,rv.u);

#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator -(const "+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]-sl.dat[j];

		return sum;
	}
	
	template <class VS,class V,class E>
	TINLINE E _vsvminus(const VS &sl,const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(sl.start,sl.end);

#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator -(const "+nameof(sl)+" &,const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=sl.dat[j]-rv.dat[i];
		return sum;
	}

	template <class V1,class V2,class E>
	TINLINE E _vvconv(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv1.l,rv1.u);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[i]|rv2.dat[i];
		return sum;
	}

	template <class V,class VS,class E>
	TINLINE E _vvsconv(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv.l,rv.u);

#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(rv)+" &,const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]|sl.dat[j];
		return sum;
	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvsconv(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(s1.start,s1.end);

#if(CXSC_INDEX_CHECK)
		if(s1.size!=s2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(s1)+" &,const "+nameof(s2)+" &)"));
#endif
		for(int i=s1.start-s1.l,j=s2.start-s2.l,k=0;k<s1.size;i++,j++,k++)
			sum.dat[k]=s1.dat[i]|s2.dat[j];
		return sum;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvconvassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" & operator +=("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]|=rv2.dat[i];
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvsconvassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" & operator +=("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]|=sl.dat[j];
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvconvassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" & operator +=("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]|=rv.dat[i];
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsconvassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" & operator +=("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			sl1.dat[j]|=sl2.dat[k];
		return sl1;
	}

	template <class V1,class V2,class E>
	TINLINE E _vvsect(const V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
		E sum(rv1.l,rv1.u);

#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(sum)+" operator +(const "+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for (int i=0;i<rv1.size;i++)
			sum.dat[i]=rv1.dat[i]&rv2.dat[i];
		return sum;
	}

	template <class V,class VS,class E>
	TINLINE E _vvssect(const V &rv,const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(rv.l,rv.u);

#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(rv)+" &,const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sum.dat[i]=rv.dat[i]&sl.dat[j];
		return sum;
	}

	template <class VS1,class VS2,class E>
	TINLINE E _vsvssect(const VS1 &s1,const VS2 &s2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<E>)
#else
	throw()
#endif
	{
		E sum(s1.start,s1.end);

#if(CXSC_INDEX_CHECK)
		if(s1.size!=s2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<E>(nameof(sum)+" operator +(const "+nameof(s1)+" &,const "+nameof(s2)+" &)"));
#endif
		for(int i=s1.start-s1.l,j=s2.start-s2.l,k=0;k<s1.size;i++,j++,k++)
			sum.dat[k]=s1.dat[i]&s2.dat[j];
		return sum;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvsectassign(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" & operator +=("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			rv1.dat[i]&=rv2.dat[i];
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvssectassign(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" & operator +=("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			rv.dat[i]&=sl.dat[j];
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvsectassign(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" & operator +=("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			sl.dat[j]&=rv.dat[i];
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssectassign(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" & operator +=("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			sl1.dat[j]&=sl2.dat[k];
		return sl1;
	}

	template <class V1,class V2>
	TINLINE bool _vveq(const V1 &rv1, const V2 &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);

		int i;
		for (i=0;i<rv1.size && rv1.dat[i]==rv2.dat[i];i++);
		
		return (i==rv1.size);
	}	

	template <class VS,class V>
	TINLINE bool _vsveq(const VS &sl, const V &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);

		int i,j;
		for (i=0,j=sl.start-sl.l;i<rv.size && sl.dat[j]==rv.dat[i];i++,j++);
		
		return (i==rv.size);
	}	

	template <class V1,class V2>
	TINLINE bool _vvneq(const V1 &rv1, const V2 &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(true);

		int i;
		for (i=0;i<rv1.size && rv1.dat[i]==rv2.dat[i];i++);
		
		return (i!=rv1.size);
	}	

	template <class VS,class V>
	TINLINE bool _vsvneq(const VS &sl, const V &rv) throw()
	{
		if(sl.size!=rv.size)
			return(true);

		int i,j;
		for (i=0,j=sl.start-sl.l;i<rv.size && sl.dat[j]==rv.dat[i];i++,j++);
		
		return (i!=rv.size);
	}	

	template <class V1,class V2>
	TINLINE bool _vvless(const V1 &rv1, const V2 &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);
		
		int i;
		for(i=0;i<rv1.size&&(rv1.dat[i]<rv2.dat[i]);i++);

		return (i==rv1.size);
	}

	template <class VS,class V>
	TINLINE bool _vsvless(const VS &sl, const V &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;j<sl.size&&(sl.dat[i]<rv.dat[j]);i++,j++);

		return (j==sl.size);
	}

	template <class V1,class V2>
	TINLINE bool _vvleq(const V1 &rv1, const V2 &rv2) throw()
	{
		if(rv1.size!=rv2.size)
			return(false);
		
		int i;
		for(i=0;i<rv1.size&&(rv1.dat[i]<=rv2.dat[i]);i++);

		return (i==rv1.size);
	}

	template <class VS,class V>
	TINLINE bool _vsvleq(const VS &sl, const V &rv) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;j<sl.size&&(sl.dat[i]<=rv.dat[j]);i++,j++);

		return (j==sl.size);
	}

	template <class V,class VS>
	TINLINE bool _vvsless(const V &rv, const VS &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;j<sl.size&&(sl.dat[i]>rv.dat[j]);i++,j++);

		return (j==sl.size);
	}

	template <class V,class VS>
	TINLINE bool _vvsleq(const V &rv, const VS &sl) throw()
	{
		if(sl.size!=rv.size)
			return(false);
		
		int i,j;
		for(i=sl.start-sl.l,j=0;j<sl.size&&(sl.dat[i]>=rv.dat[j]);i++,j++);

		return (j==sl.size);
	}

	template <class V>
	TINLINE bool _vnot(const V &rv) throw()
	{
		int i;
		for(i=0;i<rv.size;i++)
		{
			if(!!rv.dat[i])
				return(false);
		}

		return true;
	}

	template <class V>
	TINLINE void *_vvoid(const V &rv) throw()
	{
		for(int i=0;i<rv.size;i++)
		{
			if(!!rv.dat[i])
				return (void *)1;
		}

		return (void *)0;
	}

	template <class V>
	TINLINE V _vconj(const V &rv) throw()
	{
		V erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=conj(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsconj(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=conj(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vabs(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=abs(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsabs(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=abs(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vdiam(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=diam(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsdiam(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=diam(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vmid(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=mid(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsmid(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=mid(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vinf(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=Inf(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsinf(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=Inf(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vsup(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=Sup(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vssup(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=Sup(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vre(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=Re(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsre(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=Re(sl.dat[j]);

		return erg;
	}

	template <class V,class E>
	TINLINE E _vim(const V &rv) throw()
	{
		E erg(rv.l,rv.u);
		
		for(int i=0;i<rv.size;i++)
			erg.dat[i]=Im(rv.dat[i]);

		return erg;
	}

	template <class VS,class E>
	TINLINE E _vsim(const VS &sl) throw()
	{
		E erg(sl.start,sl.end);
		
		for(int i=0,j=sl.start-sl.l;i<sl.size;i++,j++)
			erg.dat[i]=Im(sl.dat[j]);

		return erg;
	}

	template <class V,class S>
	TINLINE V &_vsusetsup(V &v, const S &s) throw()
	{
		for(int i=0;i<v.size;i++)
			UncheckedSetInf(v.dat[i],s);
		return v;
	}

	template <class V,class S>
	TINLINE V &_vsusetinf(V &v, const S &s) throw()
	{
		for(int i=0;i<v.size;i++)
			UncheckedSetInf(v.dat[i],s);
		return v;
	}

	template <class V,class S>
	TINLINE V &_vssetinf(V &v, const S &s) throw()
	{
		for(int i=0;i<v.size;i++)
			SetInf(v.dat[i],s);
		return v;
	}

	template <class V,class S>
	TINLINE V &_vssetsup(V &v, const S &s) throw()
	{
		for(int i=0;i<v.size;i++)
			SetSup(v.dat[i],s);
		return v;
	}

	template <class V,class S>
	TINLINE V &_vssetre(V &v, const S &s) throw()
	{
		for(int i=0;i<v.size;i++)
			SetRe(v.dat[i],s);
		return v;
	}

	template <class V,class S>
	TINLINE V &_vssetim(V &v, const S &s) throw()
	{
		for(int i=0;i<v.size;i++)
			SetIm(v.dat[i],s);
		return v;
	}

	template <class VS,class S>
	TINLINE VS &_vssusetsup(VS &vs, const S &s) throw()
	{
		for(int i=vs.start-vs.l;i<=vs.end-vs.l;i++)
			UncheckedSetInf(vs.dat[i],s);
		return vs;
	}

	template <class VS,class S>
	TINLINE VS &_vssusetinf(VS &vs, const S &s) throw()
	{
		for(int i=vs.start-vs.l;i<=vs.end-vs.l;i++)
			UncheckedSetInf(vs.dat[i],s);
		return vs;
	}

	template <class VS,class S>
	TINLINE VS &_vsssetinf(VS &vs, const S &s) throw()
	{
		for(int i=vs.start-vs.l;i<=vs.end-vs.l;i++)
			SetInf(vs.dat[i],s);
		return vs;
	}

	template <class VS,class S>
	TINLINE VS &_vsssetsup(VS &vs, const S &s) throw()
	{
		for(int i=vs.start-vs.l;i<=vs.end-vs.l;i++)
			SetSup(vs.dat[i],s);
		return vs;
	}

	template <class VS,class S>
	TINLINE VS &_vsssetre(VS &vs, const S &s) throw()
	{
		for(int i=vs.start-vs.l;i<=vs.end-vs.l;i++)
			SetRe(vs.dat[i],s);
		return vs;
	}

	template <class VS,class S>
	TINLINE VS &_vsssetim(VS &vs, const S &s) throw()
	{
		for(int i=vs.start-vs.l;i<=vs.end-vs.l;i++)
			SetIm(vs.dat[i],s);
		return vs;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvsetinf(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" & SetInf("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			SetInf(rv1.dat[i],rv2.dat[i]);
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvssetinf(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" & SetInf("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetInf(rv.dat[i],sl.dat[j]);
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvsetinf(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" & SetInf("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetInf(sl.dat[j],rv.dat[i]);
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetinf(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" &SetInf("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			SetInf(sl1.dat[j],sl2.dat[k]);
		return sl1;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvsetsup(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" &SetSup("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			SetSup(rv1.dat[i],rv2.dat[i]);
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvssetsup(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetSup("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetSup(rv.dat[i],sl.dat[j]);
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvsetsup(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" &SetSup("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetSup(sl.dat[j],rv.dat[i]);
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetsup(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" &SetSup("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			SetSup(sl1.dat[j],sl2.dat[k]);
		return sl1;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvusetinf(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" &UncheckedSetInf("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			UncheckedSetInf(rv1.dat[i],rv2.dat[i]);
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvsusetinf(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &UncheckedSetInf("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			UncheckedSetInf(rv.dat[i],sl.dat[j]);
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvusetinf(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" &UncheckedSetInf("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			UncheckedSetInf(sl.dat[j],rv.dat[i]);
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsusetinf(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" &UncheckedSetInf("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			UncheckedSetInf(sl1.dat[j],sl2.dat[k]);
		return sl1;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvusetsup(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" &UncheckedSetSup("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			UncheckedSetSup(rv1.dat[i],rv2.dat[i]);
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvsusetsup(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &UncheckedSetSup("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			UncheckedSetSup(rv.dat[i],sl.dat[j]);
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvusetsup(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" &UncheckedSetSup("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			UncheckedSetSup(sl.dat[j],rv.dat[i]);
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvsusetsup(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" &UncheckedSetSup("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			UncheckedSetSup(sl1.dat[j],sl2.dat[k]);
		return sl1;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvsetim(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" &SetIm("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			SetIm(rv1.dat[i],rv2.dat[i]);
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvssetim(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetIm("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetIm(rv.dat[i],sl.dat[j]);
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvsetim(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>("SetIm("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetIm(sl.dat[j],rv.dat[i]);
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetim(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" &SetIm("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			SetIm(sl1.dat[j],sl2.dat[k]);
		return sl1;
	}

	template <class V1,class V2>
	TINLINE V1 &_vvsetre(V1 &rv1, const V2 &rv2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv1.size!=rv2.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V1>(nameof(rv1)+" &SetRe("+nameof(rv1)+" &, const "+nameof(rv2)+" &)"));
#endif
		for(int i=0;i<rv1.size;i++)
			SetRe(rv1.dat[i],rv2.dat[i]);
		return rv1;
	}

	template <class V,class VS>
	TINLINE V &_vvssetre(V &rv, const VS &sl)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<V>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<V>(nameof(rv)+" &SetRe("+nameof(rv)+" &, const "+nameof(sl)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetRe(rv.dat[i],sl.dat[j]);
		return rv;
	}

	template <class VS,class V>
	TINLINE VS &_vsvsetre(VS &sl, const V &rv)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(rv.size!=sl.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS>(nameof(sl)+" &SetRe("+nameof(sl)+" &, const "+nameof(rv)+" &)"));
#endif
		for(int i=0,j=sl.start-sl.l;i<rv.size;i++,j++)
			SetRe(sl.dat[j],rv.dat[i]);
		return sl;
	}

	template <class VS1,class VS2>
	TINLINE VS1 &_vsvssetre(VS1 &sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(ERROR__OP_WITH_WRONG_DIM<VS1>)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl2.size!=sl1.size) cxscthrow(ERROR__OP_WITH_WRONG_DIM<VS1>(nameof(sl1)+" &SetRe("+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=0,j=sl1.start-sl1.l,k=sl2.start-sl2.l;i<sl2.size;i++,j++,k++)
			SetRe(sl1.dat[j],sl2.dat[k]);
		return sl1;
	}

	template <class DP,class VS1,class VS2>
	TINLINE void _vsvsaccu(DP &dp, const VS1 & sl1, const VS2 &sl2)
#if(CXSC_INDEX_CHECK)
		throw(OP_WITH_WRONG_DIM)
#else
	throw()
#endif
	{
#if(CXSC_INDEX_CHECK)
		if(sl1.size!=sl2.size) cxscthrow(OP_WITH_WRONG_DIM("void accumulate("+nameof(dp)+" &, const "+nameof(sl1)+" &, const "+nameof(sl2)+" &)"));
#endif
		for(int i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size;k++,i++,j++)
			accumulate(dp,sl1.dat[i],sl2.dat[j]);
	}

	template <class VS1,class VS2>
	TINLINE bool _vsvseq(const VS1 &sl1, const VS2 &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);

		int i,j,k;
		for (i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size && sl1.dat[i]==sl2.dat[j];k++,j++,i++);
		
		return (k==sl1.size);
	}	

	template <class VS1,class VS2>
	TINLINE bool _vsvsneq(const VS1 &sl1, const VS2 &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(true);

		int i,j,k;
		for (i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size && sl1.dat[i]==sl2.dat[j];i++,j++,k++);
		
		return (k!=sl1.size);
	}	

	template <class VS1,class VS2>
	TINLINE bool _vsvsless(const VS1 &sl1, const VS2 &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(false);

		int i,j,k;
		for (i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size && sl1.dat[i]<sl2.dat[j];k++,j++,i++);
		
		return (k==sl1.size);
	}	

	template <class VS1,class VS2>
	TINLINE bool _vsvsleq(const VS1 &sl1, const VS2 &sl2) throw()
	{
		if(sl1.size!=sl2.size)
			return(true);

		int i,j,k;
		for (i=sl1.start-sl1.l,j=sl2.start-sl2.l,k=0;k<sl1.size && sl1.dat[i]<=sl2.dat[j];i++,j++,k++);
		
		return (k==sl1.size);
	}	

	template <class VS>
	TINLINE bool _vsnot(const VS &sl) throw()
	{
		for(int i=sl.start-sl.l,k=0;k<sl.size;i++,k++)
		{
			if(!!sl.dat[i])
				return(false);
		}

		return true;
	}

	template <class VS>
	TINLINE void *_vsvoid(const VS &sl) throw()
	{
		for(int i=sl.start-sl.l,k=0;i<sl.size;i++,k++)
		{
			if(!!sl.dat[i])
				return (void *)1;
		}

		return (void *)0;
	}

	template <class V>
	std::ostream &_vout(std::ostream &s, const V &rv) throw()
	{
		for(int j=0;j<rv.size;j++)
			s<<rv.dat[j]<<std::endl;
		return s;
	}

	template <class V>
	std::istream &_vin(std::istream &s, V &rv) throw()
	{
		for(int j=0;j<rv.size;j++)
			s>>rv.dat[j];
		return s;
	}


	template <class V>
	std::ostream &_vsout(std::ostream &s, const V &rv) throw()
	{
		for(int j=rv.start;j<=rv.end;j++)
			s<<rv[j]<<std::endl;
		return s;
	}

	template <class V>
	std::istream &_vsin(std::istream &s, V &rv) throw()
	{
		for(int j=rv.start;j<=rv.end;j++)
			s>>rv[j];
		return s;
	}

} // namespace cxsc

#endif
