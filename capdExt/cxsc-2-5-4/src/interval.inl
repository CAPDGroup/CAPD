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

/* CVS $Id: interval.inl,v 1.31 2014/01/30 17:23:45 cxsc Exp $ */

namespace cxsc {

// ---- Konstruktoren ----

inline interval::interval(const real &a,const real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
                                                       : inf(a), sup(b)
{
   if (a > b) 
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval::interval(const real &a,const real &b)"));
}

/*inline interval::interval(int a,int b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
                                      : inf(a), sup(b)
{
   if (a > b) 
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval::interval(int a,int b)"));
}

inline interval::interval(const double & a,const double & b) throw(ERROR_INTERVAL_EMPTY_INTERVAL)
                                                             : inf(a), sup(b)
{
   if (a > b) 
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval::interval(const double & a,const double & b)"));
}
*/

inline interval& interval::operator= (const real& a)
{
   inf = sup = a;
   return *this;
}


// ---- Typwandlungen ----

/*!
\deprecated use standard contructors for typecasting

\sa cxsc::interval::interval(const real&, const real&)
*/
inline interval _unchecked_interval(const real& a, const real& b) 
{
  interval tmp;
  tmp.inf = a;
  tmp.sup = b;
  return tmp;
}


// ---- Standardfunkt ---- (arithmetische Operatoren)

inline interval operator-(const interval &a) throw() { return interval (-a.sup, -a.inf); }
inline interval operator+(const interval &a) throw() { return a; }

inline interval & operator +=(interval &a,const interval &b) throw() { return a=a+b; }
inline interval & operator +=(interval &a,const real &b)     throw() { return a=a+b; }
inline interval & operator -=(interval &a,const interval &b) throw() { return a=a-b; }
inline interval & operator -=(interval &a,const real &b)     throw() { return a=a-b; }
inline interval & operator *=(interval &a,const interval &b) throw() { return a=a*b; }
inline interval & operator *=(interval &a,const real &b)     throw() { return a=a*b; }
inline interval & operator /=(interval &a,const interval &b) throw() { return a=a/b; }
inline interval & operator /=(interval &a,const real &b)     throw() { return a=a/b; }


inline interval operator |(const interval &a,const interval &b) throw() 
{
   return interval((a.inf<b.inf)?a.inf:b.inf,(a.sup>b.sup)?a.sup:b.sup);
}
inline interval operator &(const interval &a,const interval &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   return interval((a.inf>b.inf)?a.inf:b.inf,(a.sup<b.sup)?a.sup:b.sup);
}
inline interval operator |(const real &a,const interval &b) throw() 
{
   return interval((a<b.inf)?a:b.inf,(a>b.sup)?a:b.sup);
}
inline interval operator &(const real &a,const interval &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   return interval((a>b.inf)?a:b.inf,(a<b.sup)?a:b.sup);
}
inline interval operator |(const interval &a,const real &b) throw() 
{
   return interval((a.inf<b)?a.inf:b,(a.sup>b)?a.sup:b);
}
inline interval operator |(const real &a,const real &b) throw()
{
   if(a>b) return interval(b,a);
   else    return interval(a,b);
}
inline interval operator &(const interval &a,const real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   return interval((a.inf>b)?a.inf:b,(a.sup<b)?a.sup:b);
}
inline interval & operator |=(interval &a,const interval &b) throw() 
{
   a.inf=(a.inf<b.inf)?a.inf:b.inf,a.sup=(a.sup>b.sup)?a.sup:b.sup;
   return a;
}
inline interval & operator &=(interval &a,const interval &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   a.inf=(a.inf>b.inf)?a.inf:b.inf,a.sup=(a.sup<b.sup)?a.sup:b.sup;
   if(a.inf>a.sup)
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval & operator &=(interval &a,const interval &b)"));
   return a;
}
inline interval & operator |=(interval &a,const real &b) throw() 
{
   a.inf=(a.inf<b)?a.inf:b,a.sup=(a.sup>b)?a.sup:b;
   return a;
}
inline interval & operator &=(interval &a,const real &b) throw(ERROR_INTERVAL_EMPTY_INTERVAL) 
{
   a.inf=(a.inf>b)?a.inf:b,a.sup=(a.sup<b)?a.sup:b;
   if(a.inf>a.sup)
      cxscthrow(ERROR_INTERVAL_EMPTY_INTERVAL("inline interval & operator &=(interval &a,const real &b)"));
   return a;
}

// --- Vergleichsoperationen ----
inline bool operator ==(const interval &a,const interval &b) throw() {   return(a.inf==b.inf && a.sup==b.sup); }
inline bool operator !=(const interval &a,const interval &b) throw() {   return(a.inf!=b.inf || a.sup!=b.sup); }
inline bool operator ==(const real &r,const interval &a)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const real &r,const interval &a)     throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const real &r)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const real &r)     throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const int &r,const interval &a)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const int &r,const interval &a)     throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const int &r)     throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const int &r)     throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const long &r,const interval &a)    throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const long &r,const interval &a)    throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const long &r)    throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const long &r)    throw() {   return(r!=a.inf || r!=a.sup); }

inline bool operator ==(const double &r,const interval &a)  throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const double &r,const interval &a)  throw() {   return(r!=a.inf || r!=a.sup); }
inline bool operator ==(const interval &a,const double &r)  throw() {   return(r==a.inf && r==a.sup); }
inline bool operator !=(const interval &a,const double &r)  throw() {   return(r!=a.inf || r!=a.sup); }

// --- Mengenvergleiche ---
// <,>,...
inline bool operator <=(const interval &a,const interval &b) throw()
{
   return(a.inf>=b.inf && a.sup<=b.sup);   
}
inline bool operator >=(const interval &a,const interval &b) throw()
{
   return(a.inf<=b.inf && a.sup>=b.sup);   
}
inline bool operator <(const interval &a,const interval &b) throw()
{
   return(a.inf>b.inf && a.sup<b.sup);   
}
inline bool operator >(const interval &a,const interval &b) throw()
{
   return(a.inf<b.inf && a.sup>b.sup);   
}

inline bool operator <=(const real &a,const interval &b) throw()
{
   return(a>=b.inf && a<=b.sup);   
}
inline bool operator >=(const real &a,const interval &b) throw()
{
   return(a<=b.inf && a>=b.sup);   
}
inline bool operator <(const real &a,const interval &b) throw()
{
   return(a>b.inf && a<b.sup);   
}

inline bool operator <=(const interval &a,const real &b) throw()
{
   return(a.inf>=b && a.sup<=b);   
}
inline bool operator >=(const interval &a,const real &b) throw()
{
   return(a.inf<=b && a.sup>=b);   
}
inline bool operator >(const interval &a,const real &b) throw()
{
   return(a.inf<b && a.sup>b);   
}

inline bool operator !(const interval &a) throw() { return (a.inf <= 0.0 && a.sup >= 0.0); }  

inline       real & Inf (interval& a)       throw() { return a.inf; }
inline const real & Inf (const interval &a) throw() { return a.inf; }
inline       real & Sup (interval& a)       throw() { return a.sup; }
inline const real & Sup (const interval &a) throw() { return a.sup; }

inline interval& SetInf (interval& a, const real& b)  throw() {a.inf=b; return a;}
inline interval& SetSup (interval& a, const real& b) throw()  {a.sup=b; return a;}
inline interval& UncheckedSetInf (interval& a, const real& b) throw() { a.inf=b; return a;}
inline interval& UncheckedSetSup (interval& a, const real& b) throw() { a.sup=b; return a;}

inline bool IsEmpty(const interval &a) throw() { return (a.inf>a.sup); }

inline interval abs(const interval &a) throw()
{
   real h1  = abs(a.inf);
   real h2  = abs(a.sup);

   if (IsEmpty(a)) return a;
   if (!a)         
      return interval(0.0, (h1 > h2) ? h1 : h2);
   if (h1 > h2)    
      return interval(h2, h1);

   return interval(h1, h2); 
}

inline real diam(const interval & a) throw()
{
   return subup(a.sup,a.inf); 
}

inline real Mid(const interval & a) throw()
{
   return addd( a.inf, subd(0.5*a.sup,0.5*a.inf) );
}

inline void times2pown(interval& x, const int& n) throw()
{
    real r1,r2;
    int j;
    r1 = x.inf;  r2 = x.sup;
    j = expo(r1) + n;
    if (j >= -1021) r1 = comp(mant(r1),j);
    else 
    {
	j += 1021;
	r1 = comp(mant(r1), -1021);
	if (j<-53)
	{
	    if (sign(r1)>=0) r1 = 0;
	    else r1 = -minreal;
	} else 
	    r1 = muld(r1,comp(0.5,j+1));
    }
    j = expo(r2) + n;
    if (j >= -1021) r2 = comp(mant(r2),j);
    else 
    {
	j += 1021;
	r2 = comp(mant(r2), -1021);
	if (j<-53)
	{
	    if (sign(r2)>0) r2 = minreal;
	    else r2 = 0;
	} else r2 = mulu(r2,comp(0.5,j+1));
    }
    x = _interval(r1,r2);
} // times2pown(...)

interval operator+ (const interval& a, const interval& b) throw()
  { return interval (adddown(a.inf,b.inf),addup(a.sup,b.sup)); }
interval operator+ (const interval& a, const real& b) throw() 
  { return interval (adddown(a.inf,b),addup(a.sup,b)); }
interval operator+ (const real& a, const interval& b) throw() 
  { return interval (adddown(a,b.inf),addup(a,b.sup)); }

interval operator-  (const interval& a, const interval& b) throw() 
  { return interval ( subdown(a.inf,b.sup), subup(a.sup,b.inf)); }
interval operator-  (const interval& a, const real& b) throw() 
  { return interval ( subdown(a.inf,b), subup(a.sup,b)); }
interval operator-  (const real& a, const interval& b) throw() 
  { return interval ( subdown(a,b.sup), subup(a,b.inf)); }

  //------------------------------------------------------------------------
  //  a * b  = ueber Entscheidungstabelle :
  //
  //                      bi,bs >= 0       bi < 0, bs >=0       bi,bs < 0
  // ----------------+------------------+------------------+----------------
  // ai,as >= 0      I   ai*bi, as*bs   I   as*bi, as*bs   I   as*bi, ai*bs
  // ----------------+------------------+------------------+----------------
  // ai < 0, as >= 0 I   ai*bs, as*bs   I      ....        I   as*bi, ai*bi
  // ----------------+------------------+------------------+----------------
  // ai,as < 0       I   ai*bs, as*bi   I   ai*bs, ai*bi   I   as*bs, ai*bi
  // ----------------+------------------+------------------+----------------
  //
  //  .... :  min(ai*bs, as*bi), max(ai*ai, as*as)
  //
interval operator *(const interval& a, const interval& b) throw()
{
   interval tmp;

   if (sign(a.inf) >= 0) 
   {                                   // 1. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0)
      {                                  //     1. Spalte: [ai*bi, as*bs]
         tmp.inf = multdown(a.inf,b.inf);
         tmp.sup = multup(a.sup,b.sup);
      } else if (sign(b.sup) >= 0) 
      {                                  //     2. Spalte: [as*bi, as*bs]
         tmp.inf = multdown(a.sup,b.inf);
         tmp.sup = multup(a.sup,b.sup);
      } else 
      {                                  //     3. Spalte: [as*bi, ai*bs]
         tmp.inf = multdown(a.sup,b.inf);
         tmp.sup = multup(a.inf,b.sup);
      }
   } else if (sign(a.sup) >= 0) 
   {                                   // 2. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0) 
      {                                  //     1. Spalte: [ai*bs, as*bs]
         tmp.inf = multdown(a.inf,b.sup);
         tmp.sup = multup(a.sup,b.sup);
      } else if (sign(b.sup) >= 0) 
      {                                  //   2. Spalte: [min(ai*bs, as*bi),
         real hlp;                       //                 max(ai*ai, as*as)]

         tmp.inf = multdown(a.inf,b.sup);
         hlp     = multdown(a.sup,b.inf);

         if (hlp < tmp.inf) 
            tmp.inf = hlp;

         tmp.sup = multup(a.inf,b.inf);
         hlp     = multup(a.sup,b.sup);
         
         if (hlp > tmp.sup) 
            tmp.sup = hlp;               
      } else 
      {                                  //     3. Spalte: [as*bi, ai*bi]
         tmp.inf = multdown(a.sup,b.inf);
         tmp.sup = multup(a.inf,b.inf);
      }
   } else 
   {                                   // 3. Zeile der Entscheidungstabelle
      if (sign(b.inf) >= 0) 
      {                                  //     1. Spalte: [ai*bs, as*bi]
         tmp.inf = multdown(a.inf,b.sup);
         tmp.sup = multup(a.sup,b.inf);
      } else if (sign(b.sup) >= 0) 
      {                                  //   2. Spalte: [ai*bs, ai*bi]
         tmp.inf = multdown(a.inf,b.sup);
         tmp.sup = multup(a.inf,b.inf);
      } else 
      {                                  //     3. Spalte: [as*bs, ai*bi]
         tmp.inf = multdown(a.sup,b.sup);
         tmp.sup = multup(a.inf,b.inf);
      }
   }
   return tmp;
}

  //------------------------------------------------------------------------
        //  a / b  = ueber Entscheidungstabelle :
  //
  //                      bi,bs > 0          bi,bs < 0
  // ----------------+------------------+------------------
  // ai,as >= 0      I   ai/bs, as/bi   I   as/bs, ai/bi
  // ----------------+------------------+------------------
  // ai < 0, as >= 0 I   ai/bi, as/bi   I   as/bs, ai/bs
  // ----------------+------------------+------------------
  // ai,as < 0       I   ai/bi, as/bs   I   as/bi, ai/bs
  // ----------------+------------------+------------------
  //
interval operator/  (const interval& a, const interval& b) throw(DIV_BY_ZERO)
{
   interval tmp;

   if ((sign(b.inf) <= 0) && (sign(b.sup) >= 0))
      cxscthrow(DIV_BY_ZERO("interval::interval operator/(const interval&,const interval&)"));   

   if (sign(a.inf) >= 0) 
   {                                    // 1. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) 
      {                                 //     1. Spalte: [ai/bs, as/bi]
         tmp.inf = divdown(a.inf,b.sup);
         tmp.sup = divup(a.sup,b.inf);
      } else 
      {                                 //     2. Spalte: [as/bs, ai/bi]
         tmp.inf = divdown(a.sup,b.sup);
         tmp.sup = divup(a.inf,b.inf);
      }
   } else if (sign(a.sup) >= 0) 
   {                                    // 2. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) 
      {                                 //     1. Spalte: [ai/bi, as/bi]
         tmp.inf = divdown(a.inf,b.inf);
         tmp.sup = divup(a.sup,b.inf);
      } else 
      {                                 //     2. Spalte: [as/bs, ai/bs]
         tmp.inf = divdown(a.sup,b.sup);
         tmp.sup = divup(a.inf,b.sup);
      }
   } else 
   {                                    // 3. Zeile der Entscheidungstabelle
      if (sign(b.inf) > 0) 
      {                                 //     1. Spalte: [ai/bi, as/bs]
         tmp.inf = divdown(a.inf,b.inf);
         tmp.sup = divup(a.sup,b.sup);
      } else 
      {                                 //     2. Spalte: [as/bi, ai/bs]
         tmp.inf = divdown(a.sup,b.inf);
         tmp.sup = divup(a.inf,b.sup);
      }
   }

   return tmp;
}

interval operator*  (const real& a, const interval& b) throw()
{
   interval tmp;
   if (sign(a) == 0) 
   {
      tmp.inf = 0.0;
      tmp.sup = 0.0;      
   } else if (sign(a) > 0) 
   {
      tmp.inf = multdown(a,b.inf);
      tmp.sup = multup(a,b.sup);
   } else // if (sign(a) < 0) 
   {
      tmp.inf = multdown(a,b.sup);
      tmp.sup = multup(a,b.inf);
   }
   return tmp;
}

interval operator*  (const interval& a, const real& b) throw()
{
   interval tmp;
   if (sign(b) == 0) 
   {
      tmp.inf = 0.0;
      tmp.sup = 0.0;      
   } else if (sign(b) > 0) 
   {
      tmp.inf = multdown(a.inf,b);
      tmp.sup = multup(a.sup,b);
   } else // if (sign(b) < 0) 
   {
      tmp.inf = multdown(a.sup,b);
      tmp.sup = multup(a.inf,b);
   }
   return tmp;
}

interval operator/  (const real& a, const interval& b)  throw() { return (interval(a) / b); }
interval operator/  (const interval& a, const real& b)  throw() { return (a / interval(b)); }


} // namespace cxsc

