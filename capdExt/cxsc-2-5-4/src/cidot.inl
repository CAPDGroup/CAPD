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

/* CVS $Id: cidot.inl,v 1.28 2014/01/30 17:23:43 cxsc Exp $ */

namespace cxsc {
inline cidotprecision _cidotprecision(const complex &,const complex &) throw();
inline cidotprecision _cidotprecision(const complex &,const real &) throw();
inline cidotprecision _cidotprecision(const real &,const complex &) throw();
inline cidotprecision _cidotprecision(const interval &,const interval &) throw();
inline cidotprecision _cidotprecision(const interval &,const real &) throw();
inline cidotprecision _cidotprecision(const real &,const interval &) throw();
inline cidotprecision _cidotprecision(const real &) throw();
inline cidotprecision _cidotprecision(const complex &) throw();
inline cidotprecision _cidotprecision(const interval &) throw();
inline cidotprecision _cidotprecision(const cinterval &) throw();
      
inline cidotprecision _cidotprecision(const idotprecision &,const idotprecision &) throw();
inline cidotprecision _cidotprecision(const cdotprecision &,const cdotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
inline cidotprecision _cidotprecision(const idotprecision &,const dotprecision &) throw();
inline cidotprecision _cidotprecision(const cdotprecision &,const dotprecision &) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
inline cidotprecision _cidotprecision(const dotprecision &,const idotprecision &) throw();
inline cidotprecision _cidotprecision(const dotprecision &,const cdotprecision&) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL);
inline cidotprecision _cidotprecision(const cdotprecision &) throw();
inline cidotprecision _cidotprecision(const idotprecision &) throw();
inline cidotprecision _cidotprecision(const dotprecision &) throw();

inline cidotprecision::cidotprecision(const real & a) throw()
   : reinf(a), resup(reinf),
     iminf(0), imsup(iminf), k(0)
{
}

inline cidotprecision::cidotprecision(const dotprecision & a) throw()
   : reinf(a), resup(reinf),
     iminf(0), imsup(iminf)
{
	set_k(0);
}

inline cidotprecision::cidotprecision(const complex & a) throw()
   : reinf(Re(a)), resup(reinf),
     iminf(Im(a)), imsup(iminf), k(0)
{
}

inline cidotprecision::cidotprecision(const cdotprecision & a) throw()
   : reinf(Re(a)), resup(reinf),
     iminf(Im(a)), imsup(iminf)
{
	set_k(0);
}

inline cidotprecision::cidotprecision(const interval & a) throw()
   : reinf(Inf(a)), resup(Sup(a)),
     iminf(0),      imsup(iminf), k(0)
{
}

inline cidotprecision::cidotprecision(const idotprecision & a) throw()
   : reinf(Inf(a)), resup(Sup(a)),
     iminf(0),      imsup(iminf)
{
	set_k(0);
}

inline cidotprecision::cidotprecision(const cinterval & a) throw()
   : reinf(InfRe(a)), resup(SupRe(a)),
     iminf(InfIm(a)), imsup(SupIm(a)), k(0)
{
}

inline cidotprecision::cidotprecision(const cidotprecision & a) throw()
   : reinf(a.reinf), resup(a.resup),
     iminf(a.iminf), imsup(a.imsup)
{
	set_k(a.get_k());
}

inline cidotprecision::cidotprecision(const idotprecision & a, const idotprecision & b) throw()
   : reinf(Inf(a)), resup(Sup(a)),
     iminf(Inf(b)), imsup(Sup(b))
{
	set_k(0);
}

// ---- Typwandlungen ----
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const complex & a,const complex & b) throw()
{
   return _cidotprecision(_cinterval(a,b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const complex & a,const real & b) throw()
{
   return _cidotprecision(_cinterval(a,b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const real & a,const complex & b) throw()
{
   return _cidotprecision(_cinterval(a,b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const interval & a,const interval & b) throw()
{
   return _cidotprecision(_cinterval(a,b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const interval & a,const real &b) throw()
{
   return _cidotprecision(_cinterval(a,b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const real & a,const interval & b) throw()
{
   return _cidotprecision(_cinterval(a,b));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const real &)
*/
inline cidotprecision _cidotprecision(const real & a) throw()
{
   return _cidotprecision(_cinterval(a));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const complex &)
*/
inline cidotprecision _cidotprecision(const complex & a) throw()
{
   return _cidotprecision(_cinterval(a));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const interval &)
*/
inline cidotprecision _cidotprecision(const interval & a) throw()
{
   return _cidotprecision(_cinterval(a));
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const cinterval &)
*/
inline cidotprecision _cidotprecision(const cinterval & a) throw()
{
   cidotprecision tmp;
   tmp.reinf=Inf(Re(a));
   tmp.resup=Sup(Re(a));
   tmp.iminf=Inf(Im(a));
   tmp.imsup=Sup(Im(a));
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const idotprecision & a,const idotprecision & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=Inf(a);
   tmp.resup=Sup(a);
   tmp.iminf=Inf(b);
   tmp.imsup=Sup(b);
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const cdotprecision & a,const cdotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp;
   tmp.reinf=Re(a);
   tmp.resup=Re(b);
   tmp.iminf=Im(a);
   tmp.imsup=Im(b);
   if(tmp.reinf>tmp.resup || tmp.iminf>tmp.imsup)
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision _cidotprecision(const cdotprecision & a,const cdotprecision & b)"));
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const idotprecision & a,const dotprecision & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=Inf(a);
   tmp.resup=Sup(a);
   tmp.iminf=tmp.imsup=b;
   return tmp;   
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const cdotprecision & a,const dotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp;
   tmp.reinf=Re(a);
   tmp.iminf=Re(a);
   tmp.resup=b;
   tmp.imsup=0.0;
   if(tmp.reinf>tmp.resup || tmp.iminf>tmp.imsup)
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision _cidotprecision(const cdotprecision & a,const dotprecision & b)"));
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const dotprecision & a,const idotprecision & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=tmp.resup=a;
   tmp.iminf=Inf(b);
   tmp.imsup=Sup(b);
   return tmp;   
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const dotprecision & a,const cdotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp;
   tmp.reinf=a;
   tmp.iminf=0.0;
   tmp.resup=Re(b);
   tmp.imsup=Im(b);
   if(tmp.reinf>tmp.resup || tmp.iminf>tmp.imsup)
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision _cidotprecision(const dotprecision & a,const cdotprecision & b)"));
   return tmp;  
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _cidotprecision(const cdotprecision & a) throw()
{
   cidotprecision tmp;
   tmp.reinf=tmp.resup=Re(a);
   tmp.iminf=tmp.imsup=Im(a);
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const idotprecision &)
*/
inline cidotprecision _cidotprecision(const idotprecision & a) throw()
{
   cidotprecision tmp;
   tmp.reinf=Inf(a);
   tmp.resup=Sup(a);
   tmp.iminf=tmp.imsup=0.0;
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::cidotprecision::cidotprecision(const dotprecision &)
*/
inline cidotprecision _cidotprecision(const dotprecision & a) throw()
{
   cidotprecision tmp;
   tmp.reinf=tmp.resup=a;
   tmp.iminf=tmp.imsup=0.0;
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _unchecked_cidotprecision(const complex & a,const complex & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=Re(a);
   tmp.resup=Re(b);
   tmp.iminf=Im(a);
   tmp.imsup=Im(b);
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _unchecked_cidotprecision(const complex & a,const real & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=Re(a);
   tmp.iminf=Re(a);
   tmp.resup=b;
   tmp.imsup=0.0;
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _unchecked_cidotprecision(const real & a,const complex & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=a;
   tmp.iminf=0.0;
   tmp.resup=Re(b);
   tmp.imsup=Im(b);
   return tmp;  
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _unchecked_cidotprecision(const cdotprecision & a,const cdotprecision & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=Re(a);
   tmp.resup=Re(b);
   tmp.iminf=Im(a);
   tmp.imsup=Im(b);
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _unchecked_cidotprecision(const cdotprecision & a,const dotprecision & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=Re(a);
   tmp.iminf=Re(a);
   tmp.resup=b;
   tmp.imsup=0.0;
   return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline cidotprecision _unchecked_cidotprecision(const dotprecision & a,const cdotprecision & b) throw()
{
   cidotprecision tmp;
   tmp.reinf=a;
   tmp.iminf=0.0;
   tmp.resup=Re(b);
   tmp.imsup=Im(b);
   return tmp;  
}



// ---- Standardfunkt ---- (arithmetische Operatoren)
inline cidotprecision operator -(cidotprecision a) throw()
{
   dotprecision save=-a.reinf;
   a.reinf=-a.resup;
   a.resup=save;
   
   save=-a.iminf;
   a.iminf=-a.imsup;
   a.imsup=save;
   
   return a;
}

inline cidotprecision operator +(const cidotprecision &a) throw()
{
   return a;
}

inline cidotprecision operator +(const cidotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(a);
   tmp.reinf+=b.reinf;
   tmp.resup+=b.resup;
   tmp.iminf+=b.iminf;
   tmp.imsup+=b.imsup;
   return tmp;
}

inline cidotprecision operator -(const cidotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(a);
   tmp.reinf-=b.resup;
   tmp.resup-=b.reinf;
   tmp.iminf-=b.imsup;
   tmp.imsup-=b.iminf;
   return tmp;
}

inline cidotprecision operator |(const cidotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(a);
   if(tmp.reinf>b.reinf)
      tmp.reinf=b.reinf;
   if(tmp.iminf>b.iminf)
      tmp.iminf=b.iminf;
   if(tmp.resup<b.resup)
      tmp.resup=b.resup;
   if(tmp.imsup<b.imsup)
      tmp.imsup=b.imsup;
   return tmp;   
}

inline cidotprecision operator &(const cidotprecision & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   if(tmp.reinf<b.reinf)
      tmp.reinf=b.reinf;
   if(tmp.iminf<b.iminf)
      tmp.iminf=b.iminf;
   if(tmp.resup>b.resup)
      tmp.resup=b.resup;
   if(tmp.imsup>b.imsup)
      tmp.imsup=b.imsup;
      
   if (tmp.reinf >tmp.resup || tmp.iminf > tmp.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision operator &(const cidotprecision & a,const cidotprecision & b)"));

   return tmp;   
}

inline cidotprecision & operator +=(cidotprecision & a, const cidotprecision & b) throw()
{
   a.reinf+=b.reinf;
   a.resup+=b.resup;
   a.iminf+=b.iminf;
   a.imsup+=b.imsup;
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const cidotprecision & b) throw()
{
   a.reinf-=b.resup;
   a.resup-=b.reinf;
   a.iminf-=b.imsup;
   a.imsup-=b.iminf;
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const cidotprecision & b) throw()
{
   if(a.reinf>b.reinf)
      a.reinf=b.reinf;
   if(a.resup<b.resup)
      a.resup=b.resup;
   if(a.iminf>b.iminf)
      a.iminf=b.iminf;
   if(a.imsup<b.imsup)
      a.imsup=b.imsup;   
   return a;
}

inline cidotprecision & operator &=(cidotprecision & a, const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   if(a.reinf<b.reinf)
      a.reinf=b.reinf;
   if(a.resup>b.resup)
      a.resup=b.resup;
   if(a.iminf<b.iminf)
      a.iminf=b.iminf;
   if(a.imsup>b.imsup)
      a.imsup=b.imsup;   
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision operator &=(cidotprecision & a,const cidotprecision & b)"));

   return a;
}


// CID-R

inline cidotprecision operator +(const cidotprecision & a,const real & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const real & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const real & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const real & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const real & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=_cinterval(b);
}

inline cidotprecision operator |(const real & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=_cinterval(a);   
}

inline cidotprecision operator &(const cidotprecision & a,const real & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=_cinterval(b);   
}

inline cidotprecision operator &(const real & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=_cinterval(a);   
}

inline cidotprecision & operator +=(cidotprecision & a, const real & b) throw()
{
   a.reinf+=b;
   a.resup+=b;
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const real & b) throw()
{
   a.reinf-=b;
   a.resup-=b;
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const real & b) throw()
{
   return a|=_cinterval(b);
}

inline cidotprecision & operator &=(cidotprecision & a, const real & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   return a&=_cinterval(b);
}

// CID-C

inline cidotprecision operator +(const cidotprecision & a,const complex & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const complex & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const complex & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const complex & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const complex & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=_cinterval(b);
}

inline cidotprecision operator |(const complex & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=_cinterval(a);   
}

inline cidotprecision operator &(const cidotprecision & a,const complex & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=_cinterval(b);   
}

inline cidotprecision operator &(const complex & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=_cinterval(a);   
}

inline cidotprecision & operator +=(cidotprecision & a, const complex & b) throw()
{
   a.reinf+=Re(b);
   a.resup+=Re(b);
   a.iminf+=Im(b);
   a.imsup+=Im(b);
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const complex & b) throw()
{
   a.reinf-=Re(b);
   a.resup-=Re(b);
   a.iminf-=Im(b);
   a.iminf-=Im(b);
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const complex & b) throw()
{
   return a|=_cinterval(b);
}

inline cidotprecision & operator &=(cidotprecision & a, const complex & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   return a&=_cinterval(b);
}

// CID-I

inline cidotprecision operator +(const cidotprecision & a,const interval & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const interval & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const interval & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const interval & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const interval & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=_cinterval(b);
}

inline cidotprecision operator |(const interval & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=_cinterval(a);   
}

inline cidotprecision operator &(const cidotprecision & a,const interval & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=_cinterval(b);   
}

inline cidotprecision operator &(const interval & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=_cinterval(a);   
}

inline cidotprecision & operator +=(cidotprecision & a, const interval & b) throw()
{
   a.reinf+=Inf(b);
   a.resup+=Sup(b);
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const interval & b) throw()
{
   a.reinf-=Sup(b);
   a.resup-=Inf(b);
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const interval & b) throw()
{
   return a|=_cinterval(b);
}

inline cidotprecision & operator &=(cidotprecision & a, const interval & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   return a&=_cinterval(b);
}

// CID-CI

inline cidotprecision operator +(const cidotprecision & a,const cinterval & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const cinterval & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const cinterval & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const cinterval & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const cinterval & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=b;
}

inline cidotprecision operator |(const cinterval & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=a;   
}

inline cidotprecision operator &(const cidotprecision & a,const cinterval & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=b;   
}

inline cidotprecision operator &(const cinterval & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=a;   
}

inline cidotprecision & operator +=(cidotprecision & a, const cinterval & b) throw()
{
   a.reinf+=Inf(Re(b));
   a.resup+=Sup(Re(b));
   a.iminf+=Inf(Im(b));
   a.imsup+=Sup(Im(b));
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const cinterval & b) throw()
{
   a.reinf-=Sup(Re(b));
   a.resup-=Inf(Re(b));
   a.iminf-=Sup(Im(b));
   a.imsup-=Inf(Im(b));
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const cinterval & b) throw()
{
   if(Sup(Re(b))>a.resup)
      a.resup=Sup(Re(b));
   if(Inf(Re(b))<a.reinf)
      a.reinf=Inf(Re(b));
   if(Sup(Im(b))>a.imsup)
      a.imsup=Sup(Im(b));
   if(Inf(Im(b))<a.iminf)
      a.iminf=Inf(Im(b));
   return a;
}

inline cidotprecision & operator &=(cidotprecision & a, const cinterval & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   if(Sup(Re(b))<a.resup)
      a.resup=Sup(Re(b));
   if(Inf(Re(b))>a.reinf)
      a.reinf=Inf(Re(b));
   if(Sup(Im(b))<a.imsup)
      a.imsup=Sup(Im(b));
   if(Inf(Im(b))>a.iminf)
      a.iminf=Inf(Im(b));

   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision operator &=(cidotprecision & a,const cinterval & b)"));

   return a;
}

// CID-D

inline cidotprecision operator +(const cidotprecision & a,const dotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const dotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const dotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const dotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const dotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=b;
}

inline cidotprecision operator |(const dotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=a;
}

inline cidotprecision operator &(const cidotprecision & a,const dotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=b;
}

inline cidotprecision operator &(const dotprecision & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=a;
}

inline cidotprecision & operator +=(cidotprecision & a, const dotprecision & b) throw()
{
   a.reinf+=b;
   a.resup+=b;
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const dotprecision & b) throw()
{
   a.reinf-=b;
   a.resup-=b;
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const dotprecision & b) throw()
{
   if(b<a.reinf)
      a.reinf=b;
   if(b>a.resup)
      a.resup=b;
   if(0.<a.iminf)
      a.reinf=0.;
   if(0.>a.imsup)
      a.imsup=0.;
   return a;
}

inline cidotprecision & operator &=(cidotprecision & a, const dotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   if(b>a.reinf)
      a.reinf=b;
   if(b<a.resup)
      a.resup=b;
   if(0.>a.iminf)
      a.reinf=0.;
   if(0.<a.imsup)
      a.imsup=0.;
      
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision operator &=(cidotprecision & a,const dotprecision & b)"));

   return a;   
}

// CID-CD

inline cidotprecision operator +(const cidotprecision & a,const cdotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const cdotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const cdotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const cdotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const cdotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=b;
}

inline cidotprecision operator |(const cdotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=a;
}

inline cidotprecision operator &(const cidotprecision & a,const cdotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=b;
}

inline cidotprecision operator &(const cdotprecision & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=a;
}

inline cidotprecision & operator +=(cidotprecision & a, const cdotprecision & b) throw()
{
   a.reinf+=Re(b);
   a.resup+=Re(b);
   a.iminf+=Im(b);
   a.imsup+=Im(b);
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const cdotprecision & b) throw()
{
   a.reinf-=Re(b);
   a.resup-=Re(b);
   a.iminf-=Im(b);
   a.imsup-=Im(b);
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const cdotprecision & b) throw()
{
   if(Re(b)<a.reinf)
      a.reinf=Re(b);
   if(Re(b)>a.resup)
      a.resup=Re(b);
   if(Im(b)<a.iminf)
      a.reinf=Im(b);
   if(Im(b)>a.imsup)
      a.imsup=Im(b);
   return a;
}

inline cidotprecision & operator &=(cidotprecision & a, const cdotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   if(Re(b)>a.reinf)
      a.reinf=Re(b);
   if(Re(b)<a.resup)
      a.resup=Re(b);
   if(Im(b)>a.iminf)
      a.reinf=Im(b);
   if(Im(b)<a.imsup)
      a.imsup=Im(b);
      
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision operator &=(cidotprecision & a,const cdotprecision & b)"));

   return a;   
}


// CID-ID

inline cidotprecision operator +(const cidotprecision & a,const idotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp+=b;   
}

inline cidotprecision operator +(const idotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp+=a;
}

inline cidotprecision operator -(const cidotprecision & a,const idotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp-=b;
}

inline cidotprecision operator -(const idotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(-b);
   return tmp+=a;
}

inline cidotprecision operator |(const cidotprecision & a,const idotprecision & b) throw()
{
   cidotprecision tmp(a);
   return tmp|=b;
}

inline cidotprecision operator |(const idotprecision & a,const cidotprecision & b) throw()
{
   cidotprecision tmp(b);
   return tmp|=a;
}

inline cidotprecision operator &(const cidotprecision & a,const idotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(a);
   return tmp&=b;
}

inline cidotprecision operator &(const idotprecision & a,const cidotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   cidotprecision tmp(b);
   return tmp&=a;
}

inline cidotprecision & operator +=(cidotprecision & a, const idotprecision & b) throw()
{
   a.reinf+=Inf(b);
   a.resup+=Sup(b);
   return a;
}

inline cidotprecision & operator -=(cidotprecision & a, const idotprecision & b) throw()
{
   a.reinf-=Sup(b);
   a.resup-=Inf(b);
   return a;
}

inline cidotprecision & operator |=(cidotprecision & a, const idotprecision & b) throw()
{
   if(Inf(b)<a.reinf)
      a.reinf=Inf(b);
   if(Sup(b)>a.resup)
      a.resup=Sup(b);
   if(0.<a.iminf)
      a.reinf=0.;
   if(0.>a.imsup)
      a.imsup=0.;
   return a;
}

inline cidotprecision & operator &=(cidotprecision & a, const idotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   if(Inf(b)>a.reinf)
      a.reinf=Inf(b);
   if(Sup(b)<a.resup)
      a.resup=Sup(b);
   if(0.>a.iminf)
      a.reinf=0.;
   if(0.<a.imsup)
      a.imsup=0.;
      
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision operator &=(cidotprecision & a,const idotprecision & b)"));

   return a;   
}


// ---- Vergleichsop. ----
inline bool operator !(const cidotprecision & a) throw()
{
   return a.reinf<=0. && a.resup>=0. && a.iminf<=0. && a.imsup>=0.;
}

/*inline      cidotprecision::operator void *() const throw()
{
   if(reinf>0. || resup<0. || iminf>0. || imsup<0.)
      return (void *)1;

   return (void *)0;   
}*/

inline bool operator ==(const cidotprecision & a,const cidotprecision & b) throw()
{
   return a.reinf==b.reinf && a.resup==b.resup && a.iminf==b.iminf && a.imsup==b.imsup;
}

inline bool operator !=(const cidotprecision & a,const cidotprecision & b) throw()
{
   return a.reinf!=b.reinf || a.resup!=b.resup || a.iminf!=b.iminf || a.imsup!=b.imsup;
}

// CID-R

inline bool operator== (const cidotprecision & a, const real & b)    throw()
{
   return a.reinf==b && a.resup==b && a.iminf==0. && a.imsup==0.;
}

inline bool operator== (const real & a, const cidotprecision & b)    throw()
{
   return b.reinf==a && b.resup==a && b.iminf==0. && b.imsup==0.;
}

inline bool operator!= (const cidotprecision & a, const real & b)    throw()
{
   return a.reinf!=b || a.resup!=b || a.iminf!=0. || a.imsup!=0.;
}

inline bool operator!= (const real & a, const cidotprecision & b)    throw()
{
   return b.reinf!=a || b.resup!=a || b.iminf!=0. || b.imsup!=0.;
}

// CID-C

inline bool operator== (const cidotprecision & a, const complex & b)    throw()
{
   return a.reinf==Re(b) && a.resup==Re(b) && a.iminf==Im(b) && a.imsup==Im(b);
}

inline bool operator== (const complex & a, const cidotprecision & b)    throw()
{
   return b.reinf==Re(a) && b.resup==Re(a) && b.iminf==Im(a) && b.imsup==Im(a);
}

inline bool operator!= (const cidotprecision & a, const complex & b)    throw()
{
   return a.reinf!=Re(b) || a.resup!=Re(b) || a.iminf!=Im(a) || a.imsup!=Im(a);
}

inline bool operator!= (const complex & a, const cidotprecision & b)    throw()
{
   return b.reinf!=Re(a) || b.resup!=Re(a) || b.iminf!=Im(a) || b.imsup!=Im(a);
}

// CID-I

inline bool operator== (const cidotprecision & a, const interval & b)    throw()
{
   return a.reinf==Inf(b) && a.resup==Sup(b) && a.iminf==0. && a.imsup==0.;
}

inline bool operator== (const interval & a, const cidotprecision & b)    throw()
{
   return b.reinf==Inf(a) && b.resup==Sup(a) && b.iminf==0. && b.imsup==0.;
}

inline bool operator!= (const cidotprecision & a, const interval & b)    throw()
{
   return a.reinf!=Inf(b) || a.resup!=Sup(b) || a.iminf!=0. || a.imsup!=0.;
}

inline bool operator!= (const interval & a, const cidotprecision & b)    throw()
{
   return b.reinf!=Inf(a) || b.resup!=Sup(b) || b.iminf!=0. || b.imsup!=0.;
}

// CID-CI

inline bool operator== (const cidotprecision & a, const cinterval & b)    throw()
{
   return a.reinf==Inf(Re(b)) && a.resup==Sup(Re(b)) && a.iminf==Inf(Im(b)) && a.imsup==Sup(Im(b));
}

inline bool operator== (const cinterval & a, const cidotprecision & b)    throw()
{
   return b.reinf==Inf(Re(a)) && b.resup==Sup(Re(a)) && b.iminf==Inf(Im(a)) && b.imsup==Sup(Im(a));
}

inline bool operator!= (const cidotprecision & a, const cinterval & b)    throw()
{
   return a.reinf!=Inf(Re(b)) || a.resup!=Sup(Re(b)) || a.iminf!=Inf(Im(a)) || a.imsup!=Sup(Im(a));
}

inline bool operator!= (const cinterval & a, const cidotprecision & b)    throw()
{
   return b.reinf!=Inf(Re(a)) || b.resup!=Sup(Re(b)) || b.iminf!=Inf(Im(a)) || b.imsup!=Sup(Im(a));
}

// CID-D

inline bool operator== (const cidotprecision & a, const dotprecision & b)    throw()
{
   return a.reinf==b && a.resup==b && a.iminf==0. && a.imsup==0.;
}

inline bool operator== (const dotprecision & a, const cidotprecision & b)    throw()
{
   return b.reinf==a && b.resup==a && b.iminf==0. && b.imsup==0.;
}

inline bool operator!= (const cidotprecision & a, const dotprecision & b)    throw()
{
   return a.reinf!=b || a.resup!=b || a.iminf!=0. || a.imsup!=0.;
}

inline bool operator!= (const dotprecision & a, const cidotprecision & b)    throw()
{
   return b.reinf!=a || b.resup!=b || b.iminf!=0. || b.imsup!=0.;
}


// CID-CD

inline bool operator== (const cidotprecision & a, const cdotprecision & b)    throw()
{
   return a.reinf==Re(b) && a.resup==Re(b) && a.iminf==Im(b) && a.imsup==Im(b);
}

inline bool operator== (const cdotprecision & a, const cidotprecision & b)    throw()
{
   return b.reinf==Re(a) && b.resup==Re(a) && b.iminf==Im(a) && b.imsup==Im(a);
}

inline bool operator!= (const cidotprecision & a, const cdotprecision & b)    throw()
{
   return a.reinf!=Re(b) || a.resup!=Re(b) || a.iminf!=Im(a) || a.imsup!=Im(a);
}

inline bool operator!= (const cdotprecision & a, const cidotprecision & b)    throw()
{
   return b.reinf!=Re(a) || b.resup!=Re(b) || b.iminf!=Im(a) || b.imsup!=Im(a);
}

// CID-ID

inline bool operator== (const cidotprecision & a, const idotprecision & b)    throw()
{
   return a.reinf==Inf(b) && a.resup==Sup(b) && a.iminf==0. && a.imsup==0.;
}

inline bool operator== (const idotprecision & a, const cidotprecision & b)    throw()
{
   return b.reinf==Inf(a) && b.resup==Sup(a) && b.iminf==0. && b.imsup==0.;
}

inline bool operator!= (const cidotprecision & a, const idotprecision & b)    throw()
{
   return a.reinf!=Inf(b) || a.resup!=Sup(b) || a.iminf!=0. || a.imsup!=0.;
}

inline bool operator!= (const idotprecision & a, const cidotprecision & b)    throw()
{
   return b.reinf!=Inf(a) || b.resup!=Sup(b) || b.iminf!=0. || b.imsup!=0.;
}

// ---- Set Operators ----
inline bool operator  <(const cidotprecision & a,const cidotprecision & b) throw()
{
   return a.reinf>b.reinf && a.resup<b.resup && a.iminf>b.iminf && a.imsup<b.imsup;
}

inline bool operator  >(const cidotprecision & a,const cidotprecision & b) throw()
{
   return a.reinf<b.reinf && a.resup>b.resup && a.iminf<b.iminf && a.imsup>b.imsup;
}

inline bool operator <=(const cidotprecision & a,const cidotprecision & b) throw()
{
   return a.reinf>=b.reinf && a.resup<=b.resup && a.iminf>=b.iminf && a.imsup<=b.imsup;
}

inline bool operator >=(const cidotprecision & a,const cidotprecision & b) throw()
{
   return a.reinf<=b.reinf && a.resup>=b.resup && a.iminf<=b.iminf && a.imsup>=b.imsup;
}

// CID-R

inline bool operator  <(const real & a,const cidotprecision & b) throw() { return b>_cinterval(a); }
inline bool operator  >(const real & a,const cidotprecision & b) throw() { return false; }
inline bool operator <=(const real & a,const cidotprecision & b) throw() { return b>=_cinterval(a); }
inline bool operator >=(const real & a,const cidotprecision & b) throw() { return b<=_cinterval(a); }

inline bool operator  <(const cidotprecision & a,const real & b) throw() { return false; }
inline bool operator  >(const cidotprecision & a,const real & b) throw() { return a>_cinterval(b); }
inline bool operator <=(const cidotprecision & a,const real & b) throw() { return a<=_cinterval(b); }
inline bool operator >=(const cidotprecision & a,const real & b) throw() { return a>=_cinterval(b); }

// CID-C

inline bool operator  <(const complex & a,const cidotprecision & b) throw() { return b>_cinterval(a); }
inline bool operator  >(const complex & a,const cidotprecision & b) throw() { return false; }
inline bool operator <=(const complex & a,const cidotprecision & b) throw() { return b>=_cinterval(a); }
inline bool operator >=(const complex & a,const cidotprecision & b) throw() { return b<=_cinterval(a); }

inline bool operator  <(const cidotprecision & a,const complex & b) throw() { return false; }
inline bool operator  >(const cidotprecision & a,const complex & b) throw() { return a>_cinterval(b); }
inline bool operator <=(const cidotprecision & a,const complex & b) throw() { return a<=_cinterval(b); }
inline bool operator >=(const cidotprecision & a,const complex & b) throw() { return a>=_cinterval(b); }

// CID-I

inline bool operator  <(const interval & a,const cidotprecision & b) throw() { return b>_cinterval(a); }
inline bool operator  >(const interval & a,const cidotprecision & b) throw() { return false; }
inline bool operator <=(const interval & a,const cidotprecision & b) throw() { return b>=_cinterval(a); }
inline bool operator >=(const interval & a,const cidotprecision & b) throw() { return b<=_cinterval(a); }

inline bool operator  <(const cidotprecision & a,const interval & b) throw() { return false; }
inline bool operator  >(const cidotprecision & a,const interval & b) throw() { return a>_cinterval(b); }
inline bool operator <=(const cidotprecision & a,const interval & b) throw() { return a<=_cinterval(b); }
inline bool operator >=(const cidotprecision & a,const interval & b) throw() { return a>=_cinterval(b); }

// CID-CI

inline bool operator  <(const cinterval & a,const cidotprecision & b) throw() 
{ 
   return Inf(Re(a))>b.reinf && Sup(Re(a))<b.resup && Inf(Im(a))>b.iminf && Sup(Im(a))<b.imsup; 
}

inline bool operator  >(const cinterval & a,const cidotprecision & b) throw() 
{ 
   return Inf(Re(a))<b.reinf && Sup(Re(a))>b.resup && Inf(Im(a))<b.iminf && Sup(Im(a))>b.imsup; 
}

inline bool operator <=(const cinterval & a,const cidotprecision & b) throw() 
{ 
   return Inf(Re(a))>=b.reinf && Sup(Re(a))<=b.resup && Inf(Im(a))>=b.iminf && Sup(Im(a))<=b.imsup; 
}

inline bool operator >=(const cinterval & a,const cidotprecision & b) throw() 
{ 
   return Inf(Re(a))<=b.reinf && Sup(Re(a))>=b.resup && Inf(Im(a))<=b.iminf && Sup(Im(a))>=b.imsup; 
}

inline bool operator  <(const cidotprecision & a,const cinterval & b) throw() 
{ 
   return Inf(Re(b))<a.reinf && Sup(Re(b))>a.resup && Inf(Im(b))<a.iminf && Sup(Im(b))>a.imsup; 
}

inline bool operator  >(const cidotprecision & a,const cinterval & b) throw() 
{ 
   return Inf(Re(b))>a.reinf && Sup(Re(b))<a.resup && Inf(Im(b))>a.iminf && Sup(Im(b))<a.imsup; 
}

inline bool operator <=(const cidotprecision & a,const cinterval & b) throw() 
{
   return Inf(Re(b))<=a.reinf && Sup(Re(b))>=a.resup && Inf(Im(b))<=a.iminf && Sup(Im(b))>=a.imsup; 
}

inline bool operator >=(const cidotprecision & a,const cinterval & b) throw() 
{ 
   return Inf(Re(b))>=a.reinf && Sup(Re(b))<=a.resup && Inf(Im(b))>=a.iminf && Sup(Im(b))<=a.imsup; 
}

// CID-D

inline bool operator  <(const dotprecision & a,const cidotprecision & b) throw()
{
   return a>b.reinf && a<b.resup && 0.>b.iminf && 0.<b.imsup;
}

inline bool operator  >(const dotprecision & a,const cidotprecision & b) throw()
{
   return false;
}

inline bool operator <=(const dotprecision & a,const cidotprecision & b) throw()
{
   return a>=b.reinf && a<=b.resup && 0.>=b.iminf && 0.<=b.imsup;
}

inline bool operator >=(const dotprecision & a,const cidotprecision & b) throw()
{
   return a==b;
}

inline bool operator  <(const cidotprecision & a,const dotprecision & b) throw()
{
   return false;
}

inline bool operator  >(const cidotprecision & a,const dotprecision & b) throw()
{
   return b>a.reinf && b<a.resup && 0.>a.iminf && 0.<a.imsup;
}

inline bool operator <=(const cidotprecision & a,const dotprecision & b) throw()
{
   return a==b;
}

inline bool operator >=(const cidotprecision & a,const dotprecision & b) throw()
{
   return b>=a.reinf && b<=a.resup && 0.>=a.iminf && 0.<=a.imsup;
}


// CID-CD

inline bool operator  <(const cdotprecision & a,const cidotprecision & b) throw()
{
   return Re(a)>b.reinf && Re(a)<b.resup && Im(a)>b.iminf && Im(a)<b.imsup;
}

inline bool operator  >(const cdotprecision & a,const cidotprecision & b) throw()
{
   return false;
}

inline bool operator <=(const cdotprecision & a,const cidotprecision & b) throw()
{
   return Re(a)>=b.reinf && Re(a)<=b.resup && Im(a)>=b.iminf && Im(a)<=b.imsup;
}

inline bool operator >=(const cdotprecision & a,const cidotprecision & b) throw()
{
   return a==b;
}

inline bool operator  <(const cidotprecision & a,const cdotprecision & b) throw()
{
   return false;
}

inline bool operator  >(const cidotprecision & a,const cdotprecision & b) throw()
{
   return Re(b)>a.reinf && Re(b)<a.resup && Im(b)>a.iminf && Im(b)<a.imsup;
}

inline bool operator <=(const cidotprecision & a,const cdotprecision & b) throw()
{
   return a==b;
}

inline bool operator >=(const cidotprecision & a,const cdotprecision & b) throw()
{
   return Re(b)>=a.reinf && Re(b)<=a.resup && Im(a)>=a.iminf && Im(a)<=a.imsup;
}

// CID-ID

inline bool operator  <(const idotprecision & a,const cidotprecision & b) throw()
{
   return Inf(a)>b.reinf && Sup(a)<b.resup && 0.>b.iminf && 0.<b.imsup;
}

inline bool operator  >(const idotprecision & a,const cidotprecision & b) throw()
{
   return false;
}

inline bool operator <=(const idotprecision & a,const cidotprecision & b) throw()
{
   return Inf(a)>=b.reinf && Sup(a)<=b.resup && 0.>=b.iminf && 0.<=b.imsup;
}

inline bool operator >=(const idotprecision & a,const cidotprecision & b) throw()
{
   return Inf(a)<=b.reinf && Sup(a)>=b.resup && b.iminf==0. && b.imsup==0.;
}

inline bool operator  <(const cidotprecision & a,const idotprecision & b) throw()
{
   return false;
}

inline bool operator  >(const cidotprecision & a,const idotprecision & b) throw()
{
   return Inf(b)>a.reinf && Sup(b)<a.resup && 0.>a.iminf && 0.<a.imsup;
}

inline bool operator <=(const cidotprecision & a,const idotprecision & b) throw()
{
   return Inf(b)<=a.reinf && Sup(b)>=a.resup && a.iminf==0. && a.imsup==0.;
}

inline bool operator >=(const cidotprecision & a,const idotprecision & b) throw()
{
   return Inf(b)>=a.reinf && Sup(b)<=a.resup && 0.>=a.iminf && 0.<=a.imsup;
}

// ---- Funktionen    ----

cdotprecision   Inf(const cidotprecision & a) throw() { return cdotprecision(a.reinf,a.iminf); }
cdotprecision   Sup(const cidotprecision & a) throw() { return cdotprecision(a.resup,a.imsup); }

cidotprecision & SetInf(cidotprecision & a, const cdotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.reinf=Re(b);
   a.iminf=Im(b);
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetInf(cidotprecision & a,const cdotprecision & b)"));

   return a;
}

cidotprecision & SetSup(cidotprecision & a, const cdotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.resup=Re(b);
   a.imsup=Im(b);
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetSup(cidotprecision & a,const cdotprecision & b)"));

   return a;
}

cidotprecision & SetInf(cidotprecision & a, const dotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.reinf=b;
   a.iminf=b;
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetInf(cidotprecision & a,const dotprecision & b)"));

   return a;
}

cidotprecision & SetSup(cidotprecision & a, const dotprecision & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.resup=b;
   a.imsup=b;
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetInf(cidotprecision & a,const dotprecision & b)"));

   return a;
}

cidotprecision & SetInf(cidotprecision & a, const complex & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.reinf=Re(b);
   a.iminf=Im(b);
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetInf(cidotprecision & a,const complex & b)"));

   return a;
}

cidotprecision & SetSup(cidotprecision & a, const complex & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.resup=Re(b);
   a.imsup=Im(b);
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetSup(cidotprecision & a,const complex & b)"));

   return a;
}

cidotprecision & SetInf(cidotprecision & a, const real & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.reinf=b;
   a.iminf=b;
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetInf(cidotprecision & a,const real & b)"));

   return a;
}

cidotprecision & SetSup(cidotprecision & a, const real & b) throw(ERROR_CIDOTPRECISION_EMPTY_INTERVAL)
{
   a.resup=b;
   a.imsup=b;
   
   if (a.reinf >a.resup || a.iminf > a.imsup) 
      cxscthrow(ERROR_CIDOTPRECISION_EMPTY_INTERVAL("inline cidotprecision & SetInf(cidotprecision & a,const real & b)"));

   return a;
}

cidotprecision & UncheckedSetInf(cidotprecision & a, const cdotprecision & b) throw()
{
   a.reinf=Re(b);
   a.iminf=Im(b);
   return a;
}

cidotprecision & UncheckedSetSup(cidotprecision & a, const cdotprecision & b) throw()
{
   a.resup=Re(b);
   a.imsup=Im(b);
   return a;
}

cidotprecision & UncheckedSetInf(cidotprecision & a, const dotprecision & b) throw()
{
   a.reinf=b;
   a.iminf=b;
   return a;
}

cidotprecision & UncheckedSetSup(cidotprecision & a, const dotprecision & b) throw()
{
   a.resup=b;
   a.imsup=b;
   return a;
}

cidotprecision & UncheckedSetInf(cidotprecision & a, const complex & b) throw()
{
   a.reinf=Re(b);
   a.iminf=Im(b);
   return a;
}

cidotprecision & UncheckedSetSup(cidotprecision & a, const complex & b) throw()
{
   a.resup=Re(b);
   a.imsup=Im(b);
   return a;
}

cidotprecision & UncheckedSetInf(cidotprecision & a, const real & b) throw()
{
   a.reinf=b;
   a.iminf=b;
   return a;
}

cidotprecision & UncheckedSetSup(cidotprecision & a, const real & b) throw()
{
   a.resup=b;
   a.imsup=b;
   return a;
}

idotprecision    Re(const cidotprecision & a) throw() { return idotprecision(a.reinf,a.resup); }
idotprecision    Im(const cidotprecision & a) throw() { return idotprecision(a.iminf,a.imsup); }

inline const dotprecision & InfRe(const cidotprecision & a) throw() { return a.reinf; }
inline const dotprecision & InfIm(const cidotprecision & a) throw() { return a.iminf; }
inline const dotprecision & SupRe(const cidotprecision & a) throw() { return a.resup; }
inline const dotprecision & SupIm(const cidotprecision & a) throw() { return a.imsup; } 
      
inline       dotprecision & InfRe(cidotprecision & a) throw() { return a.reinf; }
inline       dotprecision & InfIm(cidotprecision & a) throw() { return a.iminf; } 
inline       dotprecision & SupRe(cidotprecision & a) throw() { return a.resup; } 
inline       dotprecision & SupIm(cidotprecision & a) throw() { return a.imsup; }

cidotprecision & SetRe(cidotprecision & a, const idotprecision & b) throw()
{
   a.reinf=Inf(b);
   a.resup=Sup(b);
   return a;
}

cidotprecision & SetIm(cidotprecision & a, const idotprecision & b) throw()
{
   a.iminf=Inf(b);
   a.imsup=Sup(b);
   return a;
}

cidotprecision & SetRe(cidotprecision & a, const dotprecision & b) throw()
{
   a.reinf=b;
   a.resup=b;
   return a;
}

cidotprecision & SetIm(cidotprecision & a, const dotprecision & b) throw()
{
   a.iminf=b;
   a.imsup=b;
   return a;
}

cidotprecision & SetRe(cidotprecision & a, const interval & b) throw()
{
   a.reinf=Inf(b);
   a.resup=Sup(b);
   return a;
}

cidotprecision & SetIm(cidotprecision & a, const interval & b) throw()
{
   a.iminf=Inf(b);
   a.imsup=Sup(b);
   return a;
}

cidotprecision & SetRe(cidotprecision & a, const real & b) throw()
{
   a.reinf=b;
   a.resup=b;
   return a;
}

cidotprecision & SetIm(cidotprecision & a, const real & b) throw()   
{
   a.iminf=b;
   a.imsup=b;
   return a;
}

void rnd(const cidotprecision & a,cinterval & b) throw()
{
   complex c;
   SetRe(c,rnd(a.reinf,RND_DOWN));
   SetIm(c,rnd(a.iminf,RND_DOWN));
   UncheckedSetInf(b,c);
   SetRe(c,rnd(a.resup,RND_UP));
   SetIm(c,rnd(a.imsup,RND_UP));
   UncheckedSetSup(b,c);
}

cinterval rnd(const cidotprecision & a) throw()
{
   cinterval tmp;
   rnd(a,tmp);
   return tmp;
}

inline void accumulate  (cidotprecision & a, const cinterval & b, const interval & c) throw()
{
   accumulate(a,b,cinterval(c));
}

inline void accumulate  (cidotprecision & a, const cinterval & b, const complex & c) throw()
{
   accumulate(a,b,cinterval(c));
}

inline void accumulate  (cidotprecision & a, const cinterval & b, const real & c) throw()
{
   accumulate(a,b,cinterval(c));
}

inline void accumulate  (cidotprecision & a, const interval & b,const cinterval & c) throw()
{
   accumulate(a,cinterval(b),c);
}

inline void accumulate  (cidotprecision & a, const complex & b,const cinterval & c) throw()
{
   accumulate(a,cinterval(b),c);
}

inline void accumulate  (cidotprecision & a, const real & b,const cinterval & c) throw()
{
   accumulate(a,cinterval(b),c);
}

inline void accumulate  (cidotprecision & a, const complex & b,const interval & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const interval & b,const complex & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const interval & b,const interval & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const interval & b,const real & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const real & b,const interval & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const complex & b,const complex & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const real & b,const complex & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const complex & b,const real & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

inline void accumulate  (cidotprecision & a, const real & b,const real & c) throw()
{
   accumulate(a,cinterval(b),cinterval(c));
}

} // namespace cxsc

