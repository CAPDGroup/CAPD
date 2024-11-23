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

/* CVS $Id: idot.inl,v 1.25 2014/01/30 17:23:45 cxsc Exp $ */

namespace cxsc {

// ---- Konstruktoren ----

inline idotprecision::idotprecision(const dotprecision &a,const dotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL)
         : inf(a), sup(b), k(0)
{
   if(a>b)
      cxscthrow(ERROR_IDOTPRECISION_EMPTY_INTERVAL("inline idotprecision::idotprecision(const dotprecision &a,const dotprecision &b)"));
   inf.set_k(0);
   sup.set_k(0);
}

inline idotprecision::idotprecision(const idotprecision &a)
         : inf(a.inf), sup(a.sup), k(a.k)
{
   // if (a > b) errmon (ERR_INTERVAL(EMPTY)); throw!
}

// ---- Typwandlungen ----
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::idotprecision::idotprecision(const real & a)
*/
inline idotprecision _idotprecision(const real & a) 
{
  return idotprecision (a,a);
}
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::idotprecision::idotprecision(const real & a,const real & b)
*/
inline idotprecision _idotprecision(const real & a, const real & b) 
{
  return idotprecision (a,b);
}
/*!
\deprecated use standard contructors for typecasting

\sa cxsc::idotprecision::idotprecision(const real & a,const real & b)
*/
inline idotprecision _unchecked_idotprecision(const real & a, const real & b) {
  idotprecision tmp;
  tmp.inf = a;
  tmp.sup = b;
  return tmp;
}

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::idotprecision::idotprecision(const dotprecision &r)
*/
inline idotprecision _idotprecision(const dotprecision & a) 
{
  return idotprecision (a,a);
}
/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline idotprecision _idotprecision(const dotprecision & a, const dotprecision & b) 
{
  return idotprecision (a,b);
}
/*!
\deprecated use standard contructors for typecasting

\sa ???
*/
inline idotprecision _unchecked_idotprecision(const dotprecision & a, const dotprecision & b) {
  idotprecision tmp;
  tmp.inf = a;
  tmp.sup = b;
  return tmp;
}


// ---- Standardfunkt ---- (arithmetische Operatoren)

inline idotprecision operator-(const idotprecision &a) throw() { return idotprecision (-a.sup, -a.inf); }
inline idotprecision operator+(const idotprecision &a) throw() { return a; }

inline idotprecision operator+(const idotprecision &a,const idotprecision &b) throw() { return idotprecision(a.inf+b.inf,a.sup+b.sup); }
inline idotprecision operator-(const idotprecision &a,const idotprecision &b) throw() { return idotprecision(a.inf-b.sup,a.sup-b.inf); }
inline idotprecision operator |(const idotprecision &a,const idotprecision &b) throw() 
{
   return idotprecision((a.inf<b.inf)?a.inf:b.inf,(a.sup>b.sup)?a.sup:b.sup);
}
// ID-D

inline idotprecision operator +(const idotprecision &a,const dotprecision &b) throw() { return idotprecision(a.inf+b,a.sup+b); }
inline idotprecision operator +(const dotprecision &b,const idotprecision &a) throw() { return idotprecision(a.inf+b,a.sup+b); } 
inline idotprecision operator -(const idotprecision &a,const dotprecision &b) throw() { return idotprecision(a.inf-b,a.sup-b); } 
inline idotprecision operator -(const dotprecision &a,const idotprecision &b) throw() { return idotprecision(a-b.sup,a-b.inf); }

inline idotprecision operator |(const dotprecision &a,const idotprecision &b) throw() 
{
   return idotprecision((a<b.inf)?a:b.inf,(a>b.sup)?a:b.sup);
}
inline idotprecision operator |(const idotprecision &a,const dotprecision &b) throw() 
{
   return idotprecision((a.inf<b)?a.inf:b,(a.sup>b)?a.sup:b);
}
inline idotprecision operator |(const dotprecision &a,const dotprecision &b) throw()
{
   if(a>b) return idotprecision(b,a);
   else    return idotprecision(a,b);
}
inline idotprecision operator &(const dotprecision &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a>b.inf)?a:b.inf,(a<b.sup)?a:b.sup);
}
inline idotprecision operator &(const idotprecision &a,const dotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b)?a.inf:b,(a.sup<b)?a.sup:b);
}

// ID-L

inline idotprecision operator +(const idotprecision &a,const long &b) throw() { return idotprecision(a.inf+b,a.sup+b); }
inline idotprecision operator +(const long &b,const idotprecision &a) throw() { return idotprecision(a.inf+b,a.sup+b); } 
inline idotprecision operator -(const idotprecision &a,const long &b) throw() { return idotprecision(a.inf-b,a.sup-b); } 
inline idotprecision operator -(const long &a,const idotprecision &b) throw() { return idotprecision(a-b.sup,a-b.inf); }

inline idotprecision operator |(const long &a,const idotprecision &b) throw() 
{
   return idotprecision((a<b.inf)?dotprecision(a):b.inf,(a>b.sup)?dotprecision(a):b.sup);
}
inline idotprecision operator |(const idotprecision &a,const long &b) throw() 
{
   return idotprecision((a.inf<b)?a.inf:dotprecision(b),(a.sup>b)?a.sup:dotprecision(b));
}
inline idotprecision operator &(const idotprecision &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b.inf)?a.inf:b.inf,(a.sup<b.sup)?a.sup:b.sup);
}
inline idotprecision operator &(const long &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a>b.inf)?dotprecision(a):b.inf,(a<b.sup)?dotprecision(a):b.sup);
}
inline idotprecision operator &(const idotprecision &a,const long &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b)?a.inf:dotprecision(b),(a.sup<b)?a.sup:dotprecision(b));
}

// ID-R

inline idotprecision operator +(const idotprecision &a,const real &b) throw() { return idotprecision(a.inf+b,a.sup+b); }
inline idotprecision operator +(const real &b,const idotprecision &a) throw() { return idotprecision(a.inf+b,a.sup+b); } 
inline idotprecision operator -(const idotprecision &a,const real &b) throw() { return idotprecision(a.inf-b,a.sup-b); } 
inline idotprecision operator -(const real &a,const idotprecision &b) throw() { return idotprecision(a-b.sup,a-b.inf); }

inline idotprecision operator |(const real &a,const idotprecision &b) throw() 
{
   return idotprecision((a<b.inf)?dotprecision(a):b.inf,(a>b.sup)?dotprecision(a):b.sup);
}
inline idotprecision operator |(const idotprecision &a,const real &b) throw() 
{
   return idotprecision((a.inf<b)?a.inf:dotprecision(b),(a.sup>b)?a.sup:dotprecision(b));
}
inline idotprecision operator &(const real &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a>b.inf)?dotprecision(a):b.inf,(a<b.sup)?dotprecision(a):b.sup);
}
inline idotprecision operator &(const idotprecision &a,const real &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b)?a.inf:dotprecision(b),(a.sup<b)?a.sup:dotprecision(b));
}


// ID-I

inline idotprecision operator +(const idotprecision &a,const interval &b) throw() { return idotprecision(a.inf+b.inf,a.sup+b.sup); }
inline idotprecision operator +(const interval &b,const idotprecision &a) throw() { return idotprecision(a.inf+b.inf,a.sup+b.sup); } 
inline idotprecision operator -(const idotprecision &a,const interval &b) throw() { return idotprecision(a.inf-b.sup,a.sup-b.inf); } 
inline idotprecision operator -(const interval &a,const idotprecision &b) throw() { return idotprecision(a.inf-b.sup,a.sup-b.inf); }

inline idotprecision operator |(const interval &a,const idotprecision &b) throw() 
{
   return idotprecision((a.inf<b.inf)?dotprecision(a.inf):b.inf,(a.sup>b.sup)?dotprecision(a.sup):b.sup);
}
inline idotprecision operator |(const idotprecision &a,const interval &b) throw() 
{
   return idotprecision((a.inf<b.inf)?a.inf:dotprecision(b.inf),(a.sup>b.sup)?a.sup:dotprecision(b.sup));
}
inline idotprecision operator &(const interval &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b.inf)?dotprecision(a.inf):b.inf,(a.sup<b.sup)?dotprecision(a.sup):b.sup);
}
inline idotprecision operator &(const idotprecision &a,const interval &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b.inf)?a.inf:dotprecision(b.inf),(a.sup<b.sup)?a.sup:dotprecision(b.sup));
}

// D-R

inline idotprecision operator |(const dotprecision &a,const real &b) throw() 
{
   if(a<b)
      return idotprecision(a,dotprecision(b));
   return idotprecision(dotprecision(b),a);
}
inline idotprecision operator |(const real &b,const dotprecision &a) throw() 
{
   if(a<b)
      return idotprecision(a,dotprecision(b));
   return idotprecision(dotprecision(b),a);
}



inline idotprecision & operator +=(idotprecision &a,const idotprecision &b) throw() { a.inf+=b.inf,a.sup+=b.sup; return a;}
inline idotprecision & operator +=(idotprecision &a,const dotprecision &b) throw() { a.inf+=b,a.sup+=b; return a;}
inline idotprecision & operator -=(idotprecision &a,const idotprecision &b) throw() { a.inf-=b.inf,a.sup-=b.sup; return a;}
inline idotprecision & operator -=(idotprecision &a,const dotprecision &b) throw() { a.inf-=b,a.sup-=b; return a;}
inline idotprecision & operator +=(idotprecision &a,const interval &b) throw() 
{
   a.inf+=Inf(b);
   a.sup+=Sup(b);
   return a;
}
inline idotprecision & operator -=(idotprecision &a,const interval &b) throw()
{
   a.inf-=Sup(b);
   a.sup-=Inf(b);
   return a;
}
inline idotprecision & operator +=(idotprecision &a,const real &b) throw()
{
   a.inf+=b;
   a.sup+=b;
   return a;
}
inline idotprecision & operator -=(idotprecision &a,const real &b) throw()
{
   a.inf-=b;
   a.sup-=b;
   return a;
}
      

inline idotprecision & operator |=(idotprecision &a,const idotprecision &b) throw() 
{
   if(b.inf<a.inf)
      a.inf=b.inf;
   if(b.sup>a.sup)
      a.sup=b.sup;
   return a;
}
inline idotprecision & operator &=(idotprecision &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   if(b.inf>a.inf)
      a.inf=b.inf;
   if(b.sup<a.sup)
      a.sup=b.sup;
      
   if(a.inf>a.sup)
      cxscthrow(ERROR_IDOTPRECISION_EMPTY_INTERVAL("inline idotprecision & operator &=(idotprecision &a,const idotprecision &b)"));
   return a;
}
inline idotprecision & operator |=(idotprecision &a,const dotprecision &b) throw() 
{
   if(b<a.inf)
      a.inf=b;
   if(b>a.sup)
      a.sup=b;
   return a;
}
inline idotprecision & operator &=(idotprecision &a,const dotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   if(b>a.inf)
      a.inf=b;
   if(b<a.sup)
      a.sup=b;
   if(a.inf>a.sup)
      cxscthrow(ERROR_IDOTPRECISION_EMPTY_INTERVAL("inline idotprecision & operator &=(idotprecision &a,const dotprecision &b)"));
   return a;
}

// --- Vergleichsoperationen ----
inline bool operator ==(const idotprecision &a,const idotprecision &b) throw() {   return(a.inf==b.inf && a.sup==b.sup); }
inline bool operator !=(const idotprecision &a,const idotprecision &b) throw() {   return(a.inf!=b.inf || a.sup!=b.sup); }
inline bool operator ==(const dotprecision &r,const idotprecision &a)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const dotprecision &r,const idotprecision &a)     throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const idotprecision &a,const dotprecision &r)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const idotprecision &a,const dotprecision &r)     throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const real &r,const idotprecision &a)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const real &r,const idotprecision &a)     throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const idotprecision &a,const real &r)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const idotprecision &a,const real &r)     throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const interval &a,const idotprecision &b) throw() {   return(Inf(a)==b.inf && Sup(a)==b.sup); }
inline bool operator !=(const interval &a,const idotprecision &b) throw() {   return(Inf(a)!=b.inf || Sup(a)!=b.sup); }
inline bool operator ==(const idotprecision &a,const interval &b) throw() {   return(a.inf==Inf(b) && a.sup==Sup(b)); }
inline bool operator !=(const idotprecision &a,const interval &b) throw() {   return(a.inf!=Inf(b) || a.sup!=Sup(b)); }

// --- Mengenvergleiche ---
// <,>,...
inline bool operator <=(const idotprecision &a,const idotprecision &b) throw()
{
   return(a.inf>=b.inf && a.sup<=b.sup);   
}
inline bool operator >=(const idotprecision &a,const idotprecision &b) throw()
{
   return(a.inf<=b.inf && a.sup>=b.sup);   
}
inline bool operator <(const idotprecision &a,const idotprecision &b) throw()
{
   return(a.inf>b.inf && a.sup<b.sup);   
}
inline bool operator >(const idotprecision &a,const idotprecision &b) throw()
{
   return(a.inf<b.inf && a.sup>b.sup);   
}

inline bool operator <=(const dotprecision &a,const idotprecision &b) throw()
{
   return(a>=b.inf && a<=b.sup);   
}
inline bool operator >=(const dotprecision &a,const idotprecision &b) throw()
{
   return(a<=b.inf && a>=b.sup);   
}
inline bool operator <(const dotprecision &a,const idotprecision &b) throw()
{
   return(a>b.inf && a<b.sup);   
}

inline bool operator <=(const idotprecision &a,const dotprecision &b) throw()
{
   return(a.inf>=b && a.sup<=b);   
}
inline bool operator >=(const idotprecision &a,const dotprecision &b) throw()
{
   return(a.inf<=b && a.sup>=b);   
}
inline bool operator >(const idotprecision &a,const dotprecision &b) throw()
{
   return(a.inf<b && a.sup>b);   
}

inline bool operator <=(const real &a,const idotprecision &b) throw()
{
   return(a>=b.inf && a<=b.sup);   
}
inline bool operator >=(const real &a,const idotprecision &b) throw()
{
   return(a<=b.inf && a>=b.sup);   
}
inline bool operator <(const real &a,const idotprecision &b) throw()
{
   return(a>b.inf && a<b.sup);   
}

inline bool operator <=(const idotprecision &a,const real &b) throw()
{
   return(a.inf>=b && a.sup<=b);   
}
inline bool operator >=(const idotprecision &a,const real &b) throw()
{
   return(a.inf<=b && a.sup>=b);   
}
inline bool operator >(const idotprecision &a,const real &b) throw()
{
   return(a.inf<b && a.sup>b);   
}

inline bool operator <=(const interval &a,const idotprecision &b) throw()
{
   return(Inf(a)>=b.inf && Sup(a)<=b.sup);   
}
inline bool operator >=(const interval &a,const idotprecision &b) throw()
{
   return(Inf(a)<=b.inf && Sup(a)>=b.sup);   
}
inline bool operator <(const interval &a,const idotprecision &b) throw()
{
   return(Inf(a)>b.inf && Sup(a)<b.sup);   
}
inline bool operator >(const interval &a,const idotprecision &b) throw()
{
   return(Inf(a)<b.inf && Sup(a)>b.sup);   
}

inline bool operator <=(const idotprecision &a,const interval &b) throw()
{
   return(a.inf>=Inf(b) && a.sup<=Sup(b));   
}
inline bool operator >=(const idotprecision &a,const interval &b) throw()
{
   return(a.inf<=Inf(b) && a.sup>=Sup(b));   
}
inline bool operator <(const idotprecision &a,const interval &b) throw()
{
   return(a.inf>Inf(b) && a.sup<Sup(b));   
}
inline bool operator >(const idotprecision &a,const interval &b) throw()
{
   return(a.inf<Inf(b) && a.sup>Sup(b));   
}

// ----- Funktionen -----

inline idotprecision & SetInf (idotprecision & a, const dotprecision & b)  throw()
{ // ggf. exception
   a.inf=b; 
   return a;
}
inline idotprecision & SetSup (idotprecision & a, const dotprecision & b) throw()
{
   a.sup=b;
   return a;
}
inline idotprecision & SetInf (idotprecision & a, const real & b)  throw()
{ // ggf. exception
   a.inf=b; 
   return a;
}
inline idotprecision & SetSup (idotprecision & a, const real & b) throw()
{
   a.sup=b;
   return a;
}
inline idotprecision & UncheckedSetInf (idotprecision & a, const dotprecision & b) throw()
{
   a.inf=b;
   return a;
}
inline idotprecision & UncheckedSetSup (idotprecision & a, const dotprecision & b) throw()
{
   a.sup=b;
   return a;
}
inline idotprecision & UncheckedSetInf (idotprecision & a, const real & b) throw()
{
   a.inf=b;
   return a;
}
inline idotprecision & UncheckedSetSup (idotprecision & a, const real & b) throw()
{
   a.sup=b;
   return a;
}

inline bool operator !(const idotprecision &a) throw() { return (a.inf <= dotprecision(0.0) && a.sup >= dotprecision(0.0)); }  

inline bool IsEmpty(const idotprecision &a) throw() { return (a.inf>a.sup); }

inline idotprecision abs(const idotprecision &a) throw()
{
   dotprecision h1  = abs(a.inf);
   dotprecision h2  = abs(a.sup);

   if (IsEmpty(a)) return a;
   if (!a)         
      return idotprecision(dotprecision(0.0), (h1 > h2) ? h1 : h2);
   if (h1 > h2)    
      return idotprecision(h2, h1);

   return idotprecision(h1, h2); 
}

inline void accumulate  (idotprecision & a, const interval & b, const real & c) throw() { accumulate(a,b,_interval(c)); }
inline void accumulate  (idotprecision & a, const real & b, const interval & c) throw() { accumulate(a,_interval(b),c); }
inline void accumulate  (idotprecision & a, const real & b, const real & c) throw() { accumulate(a,_interval(b),_interval(c)); }

} // namespace cxsc

