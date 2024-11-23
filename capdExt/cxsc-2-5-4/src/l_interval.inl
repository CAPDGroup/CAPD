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

/* CVS $Id: l_interval.inl,v 1.22 2014/01/30 17:23:46 cxsc Exp $ */

namespace cxsc {

inline void l_interval::_allo(int nprec)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   prec=nprec;
#if (CXSC_INDEX_CHECK)
   if(nprec<=0)
      cxscthrow(ERROR_LINTERVAL_WRONG_STAGPREC("void l_interval::_allo(int nprec)"));
#endif
   data=new real[(prec)+1];
}

// ---- Konstruktoren ----

inline l_interval::l_interval() 
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   _allo(stagprec);
}

inline l_interval::l_interval(const l_interval & a) 
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   int i;
   _allo(a.prec);
   for(i=0;i<=prec;i++)
      data[i]=a.data[i];
}
      
inline l_interval::l_interval(const real &a, const real &b) 
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(1);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("inline l_interval::l_interval(const real &a, const real &b)"));
   elem(1)=a,elem(2)=b;
   
}

l_interval::l_interval(const real &a) throw()
{
   try { _allo(1); }
   catch(const ERROR_LINTERVAL_WRONG_STAGPREC &) {} // won't occur!
   elem(1)=a,elem(2)=a;
}

inline l_interval & l_interval::operator= (const real & a)        throw()
{
   return (*this)=_l_interval(a);
}

inline l_interval & l_interval::operator= (const l_real & a)      throw()
{
   return (*this)=_l_interval(a);
}

inline l_interval & l_interval::operator= (const interval & a)    throw()
{
   return (*this)=_l_interval(a);
}

inline l_interval::~l_interval() throw()
{
   delete [] data;
}
// ---- Typwandlungen ----

l_interval::l_interval(const interval & a) throw()
{
   try { _allo(1); }
   catch(const ERROR_LINTERVAL_WRONG_STAGPREC &) {} // won't occur!
   elem(1)=Inf(a),elem(2)=Sup(a);
}

l_interval::l_interval(const l_real & a) throw()
{
   try { _allo(StagPrec(a)); }
   catch(const ERROR_LINTERVAL_WRONG_STAGPREC &) {} // shouldn't occur!
   int i;
   for(i=1;i<=prec;i++)
      elem(i)=a[i];
   elem(prec+1)=a[prec];
}

inline interval _interval(const real & a, const l_real & b) throw() // Sollte in l_real!!!
{
   return _interval(_l_interval(_l_real(a),b));
}

inline interval _interval(const l_real & a, const real & b) throw() // Sollte in l_real!!!
{
   return _interval(_l_interval(a,_l_real(b)));
}

inline interval _interval(const l_real & a) throw() // l_real!
{
   return _interval(_l_interval(a,a));
}

inline interval _unchecked_interval(const l_real & a, const l_real & b) throw() // l_real!
{
   return _unchecked_interval(_real(a),_real(b));
}

// ---- Standardfunkt ---- (arithmetische Operatoren)

inline l_interval operator+(const l_interval &a) throw() { return a; }

// LI-LI
inline l_interval operator|(const l_interval & li1, const l_interval & li2) throw(ERROR_LINTERVAL_IN_EXACT_CH_OR_IS)
{  
   // liefert ConvexHull zweier Intervalle; Einschluss von aussen!
   l_interval li3, li4;    // innen, aussen

   ConvexHull(li1, li2, li3, li4);

   if (li3!=li4)
   {
      cxscthrow(ERROR_LINTERVAL_IN_EXACT_CH_OR_IS("inline l_interval operator|(const l_interval & li1, const l_interval & li2)"));
   }
   return li4;
}  

inline l_interval operator&(const l_interval & li1, const l_interval & li2) throw(ERROR_LINTERVAL_EMPTY_INTERVAL, ERROR_LINTERVAL_IN_EXACT_CH_OR_IS)
{  
   // liefert Intersection zweier Intervalle; Einschluss von aussen!
   l_interval li3, li4;    // innen, aussen
   
   try 
   {
      Intersection(li1, li2, li3, li4);
   }
   catch(const EMPTY_INTERVAL &)
   {
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("inline l_interval operator&(const l_interval & li1, const l_interval & li2)"));
   }
      
   if (li3!=li4) 
   {
      cxscthrow(ERROR_LINTERVAL_IN_EXACT_CH_OR_IS("inline l_interval operator&(const l_interval & li1, const l_interval & li2)")); 
   }
   return li4;
}   

inline l_interval & operator +=(l_interval &a,const l_interval &b) throw() { return a=a+b; }
inline l_interval & operator -=(l_interval &a,const l_interval &b) throw() { return a=a-b; }
inline l_interval & operator *=(l_interval &a,const l_interval &b) throw() { return a=a*b; }
inline l_interval & operator /=(l_interval &a,const l_interval &b) throw() { return a=a/b; }
inline l_interval & operator |=(l_interval &a,const l_interval &b) throw() { return a=a|b; }
inline l_interval & operator &=(l_interval &a,const l_interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a=a&b; }

// LI-LR
inline l_interval operator +(const l_interval &a,const l_real &b) throw() { return a+_l_interval(b); }
inline l_interval operator +(const l_real &a,const l_interval &b) throw() { return _l_interval(a)+b; }
inline l_interval operator -(const l_interval &a,const l_real &b) throw() { return a-_l_interval(b); }
inline l_interval operator -(const l_real &a,const l_interval &b) throw() { return _l_interval(a)-b; }
inline l_interval operator *(const l_interval &a,const l_real &b) throw() { return a*_l_interval(b); }
inline l_interval operator *(const l_real &a,const l_interval &b) throw() { return _l_interval(a)*b; }
inline l_interval operator /(const l_interval &a,const l_real &b) throw() { return a/_l_interval(b); }
inline l_interval operator /(const l_real &a,const l_interval &b) throw() { return _l_interval(a)/b; } 
inline l_interval operator |(const l_real &a,const l_interval &b) throw() { return _l_interval(a)|b; }
inline l_interval operator |(const l_interval &a,const l_real &b) throw() { return a|_l_interval(b); }
inline l_interval operator |(const l_real &a,const l_real &b)     throw() { return _l_interval(a)|_l_interval(b); }  // WARNING! For exact upper and lower bounds use the constructor!
inline l_interval operator &(const l_real &a,const l_interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return _l_interval(a)&b; }
inline l_interval operator &(const l_interval &a,const l_real &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a&_l_interval(b); }

inline l_interval & operator +=(l_interval &a,const l_real &b) throw() { return a=a+_l_interval(b); }      
inline l_interval & operator -=(l_interval &a,const l_real &b) throw() { return a=a-_l_interval(b); }
inline l_interval & operator *=(l_interval &a,const l_real &b) throw() { return a=a*_l_interval(b); }              
inline l_interval & operator /=(l_interval &a,const l_real &b) throw() { return a=a/_l_interval(b); } 
inline l_interval & operator |=(l_interval &a,const l_real &b) throw() { return a=a|_l_interval(b); }
inline l_interval & operator &=(l_interval &a,const l_real &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a=a&_l_interval(b); }
 
// LI-I
inline l_interval operator +(const l_interval &a,const interval &b) throw() { return a+_l_interval(b); }
inline l_interval operator +(const interval &a,const l_interval &b) throw() { return _l_interval(a)+b; }
inline l_interval operator -(const l_interval &a,const interval &b) throw() { return a-_l_interval(b); }
inline l_interval operator -(const interval &a,const l_interval &b) throw() { return _l_interval(a)-b; }
inline l_interval operator *(const l_interval &a,const interval &b) throw() { return a*_l_interval(b); }
inline l_interval operator *(const interval &a,const l_interval &b) throw() { return _l_interval(a)*b; }
inline l_interval operator /(const l_interval &a,const interval &b) throw() { return a/_l_interval(b); }
inline l_interval operator /(const interval &a,const l_interval &b) throw() { return _l_interval(a)/b; } 
inline l_interval operator |(const interval &a,const l_interval &b) throw() { return _l_interval(a)|b; }
inline l_interval operator |(const l_interval &a,const interval &b) throw() { return a|_l_interval(b); }
inline l_interval operator &(const interval &a,const l_interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return _l_interval(a)&b; }
inline l_interval operator &(const l_interval &a,const interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a&_l_interval(b); }

inline l_interval & operator +=(l_interval &a,const interval &b) throw() { return a=a+_l_interval(b); }      
inline l_interval & operator -=(l_interval &a,const interval &b) throw() { return a=a-_l_interval(b); }
inline l_interval & operator *=(l_interval &a,const interval &b) throw() { return a=a*_l_interval(b); }              
inline l_interval & operator /=(l_interval &a,const interval &b) throw() { return a=a/_l_interval(b); } 
inline l_interval & operator |=(l_interval &a,const interval &b) throw() { return a=a|_l_interval(b); }
inline l_interval & operator &=(l_interval &a,const interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a=a&_l_interval(b); }
 
// LI-R
inline l_interval operator +(const l_interval &a,const real &b) throw() { return a+_l_interval(b); }
inline l_interval operator +(const real &a,const l_interval &b) throw() { return _l_interval(a)+b; }
inline l_interval operator -(const l_interval &a,const real &b) throw() { return a-_l_interval(b); }
inline l_interval operator -(const real &a,const l_interval &b) throw() { return _l_interval(a)-b; }
inline l_interval operator *(const l_interval &a,const real &b) throw() { return a*_l_interval(b); }
inline l_interval operator *(const real &a,const l_interval &b) throw() { return _l_interval(a)*b; }
inline l_interval operator /(const l_interval &a,const real &b) throw() { return a/_l_interval(b); }
inline l_interval operator /(const real &a,const l_interval &b) throw() { return _l_interval(a)/b; } 
inline l_interval operator |(const real &a,const l_interval &b) throw() { return _l_interval(a)|b; }
inline l_interval operator |(const l_interval &a,const real &b) throw() { return a|_l_interval(b); }
inline l_interval operator &(const real &a,const l_interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return _l_interval(a)&b; }
inline l_interval operator &(const l_interval &a,const real &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a&_l_interval(b); }

inline l_interval & operator +=(l_interval &a,const real &b) throw() { return a=a+_l_interval(b); }      
inline l_interval & operator -=(l_interval &a,const real &b) throw() { return a=a-_l_interval(b); }
inline l_interval & operator *=(l_interval &a,const real &b) throw() { return a=a*_l_interval(b); }              
inline l_interval & operator /=(l_interval &a,const real &b) throw() { return a=a/_l_interval(b); } 
inline l_interval & operator |=(l_interval &a,const real &b) throw() { return a=a|_l_interval(b); }
inline l_interval & operator &=(l_interval &a,const real &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a=a&_l_interval(b); }

// LI-ID
inline idotprecision operator +(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)+b; }
inline idotprecision operator +(const idotprecision &a,const l_interval &b) throw() { return a+_idotprecision(b); }
inline idotprecision operator -(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)-b; }
inline idotprecision operator -(const idotprecision &a,const l_interval &b) throw() { return a-_idotprecision(b); }
inline idotprecision operator |(const idotprecision &a,const l_interval &b) throw() { return a|_idotprecision(b); }
inline idotprecision operator |(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)|b; }
inline idotprecision operator &(const idotprecision &a,const l_interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a&_idotprecision(b); }
inline idotprecision operator &(const l_interval &a,const idotprecision &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return _idotprecision(a)&b; }

inline idotprecision & operator +=(idotprecision &a,const l_interval &b) throw() { return a+=_idotprecision(b); }      
inline idotprecision & operator -=(idotprecision &a,const l_interval &b) throw() { return a-=_idotprecision(b); }
inline idotprecision & operator |=(idotprecision &a,const l_interval &b) throw() { return a|=_idotprecision(b); }
inline idotprecision & operator &=(idotprecision &a,const l_interval &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) { return a&=idotprecision(b); }

// LR-R
inline l_interval operator |(const real &a,const l_real &b) throw() { return l_real(a)|b; }
inline l_interval operator |(const l_real &a,const real &b) throw() { return a|l_real(b); }

// LR-D
inline idotprecision operator |(const dotprecision &a,const l_real &b) throw() { return a|dotprecision(b); }
inline idotprecision operator |(const l_real &a,const dotprecision &b) throw() { return dotprecision(a)|b; }

// LR-I
inline l_interval operator +(const l_real &a,const interval &b) throw() { return a+_l_interval(b); }
inline l_interval operator +(const interval &a,const l_real &b) throw() { return _l_interval(a)+b; }
inline l_interval operator -(const l_real &a,const interval &b) throw() { return a-_l_interval(b); }
inline l_interval operator -(const interval &a,const l_real &b) throw() { return _l_interval(a)-b; }
inline l_interval operator *(const l_real &a,const interval &b) throw() { return a*_l_interval(b); }
inline l_interval operator *(const interval &a,const l_real &b) throw() { return _l_interval(a)*b; }
inline l_interval operator /(const l_real &a,const interval &b) throw() { return a/_l_interval(b); }
inline l_interval operator /(const interval &a,const l_real &b) throw() { return _l_interval(a)/b; } 
inline l_interval operator |(const interval &a,const l_real &b) throw() { return _l_interval(a)|b; }
inline l_interval operator |(const l_real &a,const interval &b) throw() { return a|_l_interval(b); }
inline l_interval operator &(const interval &a,const l_real &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return _l_interval(a)&b; }
inline l_interval operator &(const l_real &a,const interval &b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a&_l_interval(b); }

// ---- Vergleichsop. ----
inline bool operator !=(const l_interval &a,const l_interval &b) throw() { return !(a==b); }

inline bool operator ==(const l_real &a,const l_interval &b) throw() { return _l_interval(a)==b; }
inline bool operator !=(const l_real &a,const l_interval &b) throw() { return _l_interval(a)!=b; }
inline bool operator ==(const l_interval &a,const l_real &b) throw() { return a==_l_interval(b); }
inline bool operator !=(const l_interval &a,const l_real &b) throw() { return a!=_l_interval(b); }

inline bool operator ==(const interval &a,const l_interval &b) throw() { return _l_interval(a)==b; }
inline bool operator !=(const interval &a,const l_interval &b) throw() { return _l_interval(a)!=b; }
inline bool operator ==(const l_interval &a,const interval &b) throw() { return a==_l_interval(b); }
inline bool operator !=(const l_interval &a,const interval &b) throw() { return a!=_l_interval(b); }

inline bool operator ==(const real &a,const l_interval &b) throw() { return _l_interval(a)==b; }
inline bool operator !=(const real &a,const l_interval &b) throw() { return _l_interval(a)!=b; }
inline bool operator ==(const l_interval &a,const real &b) throw() { return a==_l_interval(b); }
inline bool operator !=(const l_interval &a,const real &b) throw() { return a!=_l_interval(b); }

inline bool operator ==(const idotprecision &a,const l_interval &b) throw() { return a==_idotprecision(b); }
inline bool operator !=(const idotprecision &a,const l_interval &b) throw() { return a!=_idotprecision(b); }
inline bool operator ==(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)==b; }
inline bool operator !=(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)!=b; }

inline bool operator ==(const dotprecision &a,const l_interval &b) throw() { return a==_idotprecision(b); }
inline bool operator !=(const dotprecision &a,const l_interval &b) throw() { return a!=_idotprecision(b); }
inline bool operator ==(const l_interval &a,const dotprecision &b) throw() { return _idotprecision(a)==b; }
inline bool operator !=(const l_interval &a,const dotprecision &b) throw() { return _idotprecision(a)!=b; }

// ---- Mengenvergle. ----
inline bool operator  <(const l_real &a,const l_interval &b) throw() { return _l_interval(a)<b; }
inline bool operator  >(const l_real &a,const l_interval &b) throw() { return _l_interval(a)>b; }
inline bool operator <=(const l_real &a,const l_interval &b) throw() { return _l_interval(a)<=b; }
inline bool operator >=(const l_real &a,const l_interval &b) throw() { return _l_interval(a)>=b; }
inline bool operator  <(const l_interval &a,const l_real &b) throw() { return a<_l_interval(b); }
inline bool operator  >(const l_interval &a,const l_real &b) throw() { return a>_l_interval(b); }
inline bool operator <=(const l_interval &a,const l_real &b) throw() { return a<=_l_interval(b); }
inline bool operator >=(const l_interval &a,const l_real &b) throw() { return a>=_l_interval(b); }

inline bool operator  <(const interval &a,const l_interval &b) throw() { return _l_interval(a)<b; }
inline bool operator  >(const interval &a,const l_interval &b) throw() { return _l_interval(a)>b; }
inline bool operator <=(const interval &a,const l_interval &b) throw() { return _l_interval(a)<=b; }
inline bool operator >=(const interval &a,const l_interval &b) throw() { return _l_interval(a)>=b; }
inline bool operator  <(const l_interval &a,const interval &b) throw() { return a<_l_interval(b); }
inline bool operator  >(const l_interval &a,const interval &b) throw() { return a>_l_interval(b); }
inline bool operator <=(const l_interval &a,const interval &b) throw() { return a<=_l_interval(b); }
inline bool operator >=(const l_interval &a,const interval &b) throw() { return a>=_l_interval(b); }

inline bool operator  <(const real &a,const l_interval &b) throw() { return _l_interval(a)<b; }
inline bool operator  >(const real &a,const l_interval &b) throw() { return _l_interval(a)>b; }
inline bool operator <=(const real &a,const l_interval &b) throw() { return _l_interval(a)<=b; }
inline bool operator >=(const real &a,const l_interval &b) throw() { return _l_interval(a)>=b; }
inline bool operator  <(const l_interval &a,const real &b) throw() { return a<_l_interval(b); }
inline bool operator  >(const l_interval &a,const real &b) throw() { return a>_l_interval(b); }
inline bool operator <=(const l_interval &a,const real &b) throw() { return a<=_l_interval(b); }
inline bool operator >=(const l_interval &a,const real &b) throw() { return a>=_l_interval(b); }

inline bool operator  <(const idotprecision &a,const l_interval &b) throw() { return a<_idotprecision(b); }
inline bool operator  >(const idotprecision &a,const l_interval &b) throw() { return a>_idotprecision(b); }
inline bool operator <=(const idotprecision &a,const l_interval &b) throw() { return a<=_idotprecision(b); }
inline bool operator >=(const idotprecision &a,const l_interval &b) throw() { return a>=_idotprecision(b); }
inline bool operator  <(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)<b; }
inline bool operator  >(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)>b; }
inline bool operator <=(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)<=b; }
inline bool operator >=(const l_interval &a,const idotprecision &b) throw() { return _idotprecision(a)>=b; }

inline bool operator  <(const dotprecision &a,const l_interval &b) throw() { return a<_idotprecision(b); }
inline bool operator  >(const dotprecision &a,const l_interval &b) throw() { return a>_idotprecision(b); }
inline bool operator <=(const dotprecision &a,const l_interval &b) throw() { return a<=_idotprecision(b); }
inline bool operator >=(const dotprecision &a,const l_interval &b) throw() { return a>=_idotprecision(b); }
inline bool operator  <(const l_interval &a,const dotprecision &b) throw() { return _idotprecision(a)<b; }
inline bool operator  >(const l_interval &a,const dotprecision &b) throw() { return _idotprecision(a)>b; }
inline bool operator <=(const l_interval &a,const dotprecision &b) throw() { return _idotprecision(a)<=b; }
inline bool operator >=(const l_interval &a,const dotprecision &b) throw() { return _idotprecision(a)>=b; }

// ---- Funktionen    ----

static const int LI_Min_Exp_ = -1074,
                 LI_maxexpo1 = 1023;
      
//inline l_interval_Inf Inf (l_interval & a)  throw() { return l_interval_Inf(a); }
//inline l_interval_Sup Sup (l_interval & a)  throw() { return l_interval_Sup(a); }
inline l_real Inf (const l_interval &li) throw()
{
   // l_real in der Praezision des Intervals
   int save_stagprec=stagprec;
   stagprec=li.prec;
   l_real lr; // l_real in der Praez. des l_intervals angele

   for (int i=1; i<=stagprec; i++)
      lr.elem(i)=li.elem(i);   // vgl. Darstellung der l_intervals

   stagprec=save_stagprec;
        
   return lr;     
}
inline l_real Sup (const l_interval &li) throw()
{
   // l_real in der Praezision des Intervals
   int save_stagprec=stagprec;
   stagprec=li.prec;
   l_real lr; // l_real in der Praez. des l_intervals angele
   
   for (int i=1; i<stagprec; i++)
      lr.elem(i)=li.elem(i);   // vgl. Darstellung der l_intervals
        
   lr.elem(stagprec)=li.elem(stagprec+1);

   stagprec=save_stagprec;

   return lr;     
}

inline int StagPrec(const l_interval &a) throw() { return a.prec; }
      
inline l_interval & SetInf (l_interval & a, const l_real & b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a=_l_interval(b,Sup((const l_interval)a)); }
inline l_interval & SetSup (l_interval & a, const l_real & b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return a=_l_interval(Inf((const l_interval)a),b); }
inline l_interval & SetInf (l_interval & a, const real & b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return SetInf(a,_l_real(b)); }
inline l_interval & SetSup (l_interval & a, const real & b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL) { return SetSup(a,_l_real(b)); }

inline l_interval adjust(const l_interval & x) throw()
{  
   l_interval  y;

   if (x.prec == stagprec) 
      y = x; 
   else if (x.prec > stagprec)
      { 
           y = x + 0.0; 
      }
   else 
   { // x.prec < stagprec
      int i;
      for (i = 0; i <= stagprec-x.prec-1; i++)
         y.data[i] = 0.0;
      for (i = stagprec-x.prec; i <= stagprec; i++)
         y.data[i] = x.data[i-(stagprec-x.prec)];
   }
   return y;
}

inline l_interval & UncheckedSetInf (l_interval & a, const l_real & b) throw() { return a=_unchecked_l_interval(b,Sup(a)); }
inline l_interval & UncheckedSetSup (l_interval & a, const l_real & b) throw() { return a=_unchecked_l_interval(Inf(a),b); }
inline l_interval & UncheckedSetInf (l_interval & a, const real & b) throw() { return UncheckedSetInf(a,_l_real(b)); }
inline l_interval & UncheckedSetSup (l_interval & a, const real & b) throw() { return UncheckedSetSup(a,_l_real(b)); }

inline l_interval abs  (const l_interval & li1) throw()
{
   l_interval li2;   
   l_real i   = Inf(li1), s   = Sup(li1);
   l_real inf = abs(i),   sup = abs(s);

   if (i>s) 
      li2 = li1;
   else if (!li1)         
      li2 = l_interval(0.0, (inf > sup) ? inf : sup);
   else if (inf > sup)    
      li2 = l_interval(sup, inf);
   else 
      li2 = l_interval(inf, sup);

   return li2;     
}

inline l_real     diam (const l_interval & li) throw()
{
   return _l_real(diam(_interval(li.elem(li.prec),li.elem(li.prec+1))));
}

inline void accumulate(idotprecision & a, const real & b, const l_interval & c) throw() { accumulate(a,_l_interval(b),c); }
inline void accumulate(idotprecision & a, const l_interval & b, const real & c) throw() { accumulate(a,b,_l_interval(c)); }
inline void accumulate(idotprecision & a, const interval & b, const l_real & c) throw() { accumulate(a,_l_interval(b),_l_interval(c)); }
inline void accumulate(idotprecision & a, const l_real & b, const interval & c) throw() { accumulate(a,_l_interval(b),_l_interval(c)); }    
inline void accumulate(idotprecision & a, const l_interval & b, const l_real & c) throw() { accumulate(a,b,_l_interval(c)); }
inline void accumulate(idotprecision & a, const l_real & b, const l_interval & c) throw() { accumulate(a,_l_interval(b),c); }
inline void accumulate(idotprecision & a, const l_interval & b, const interval & c) throw() { accumulate(a,b,_l_interval(c)); }
inline void accumulate(idotprecision & a, const interval & b, const l_interval & c) throw() { accumulate(a,_l_interval(b),c); }


inline real & l_interval::operator[](int a)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG)
#else
   throw()
#endif
{
#if (CXSC_INDEX_CHECK)
   if(a<1 || a>prec+1)
      cxscthrow(ERROR_LINTERVAL_ELEMENT_NOT_IN_LONG("inline real & l_interval::operator[](int a)"));
#endif
   return data[a-1];
}

inline void l_interval::_clear(int p) throw()
{ 
   int i;
   // fuellt l_interval ab Stelle p bis zum Ende mit Null.
   // funktioniert immer, da for-Schleife abweisend!
   for (i=p; i<=prec+1; i++)
      this->elem(i)=0.0;
}

inline bool point_intv(const l_interval &a ) // bool delivers: a is a point interval;
{
   int k(a.prec);
   return a.elem(k+1) == a.elem(k);
}

inline bool zero_(const l_interval& li) throw()
{  // returns only true if all li.elem(i) == 0; Blomquist, 27.11.02; 
    int i=1, p=StagPrec(li)+1;
    bool tmp = true;
    do
    {
	if (sign(li.elem(i))!=0) tmp = false;
	i++;
    }  while(tmp && i <= p );
    return tmp;
}



    /*!
    \param x The value for which to calculate
    \return The result of the calculation

    \sa expo_gr(const l_real&)
    */
    inline int expo_gr(const l_interval& x)
    // Calculating expo(x[k]) of the greatest |x[k]|.
      	{
	    int k(1),p(x.prec),ex1,ex2;
	    l_interval y(x);
	    while (y.elem(k)==0 && k<p) k++;
            ex1 = expo(y.elem(k));
            ex2 = expo(y.elem(k+1));
            if (ex2>ex1) ex1 = ex2; // ex1: Maximum;
	    return ex1;
	}
    /*!
    \param x The value for which to calculate
    \return The result of the calculation

    \sa expo_sm(const l_real&)
    */
    inline int expo_sm(const l_interval& x)
    // Calculating expo(x[k]) of the smallest |x[k]|<>0.
	{
	    int k(x.prec+1),ex1,ex2;
	    l_interval y(x);

	    while (y.elem(k)==0 && k>2) k--;
            ex1 = expo(y.elem(k));
            ex2 = expo(y.elem(k-1));
            if (ex2<ex1) 
	    {
	       k = ex1;  ex1 = ex2; ex2 = k;    
	    } // ex1: Minimum;
            if (ex1<-100000) ex1 = ex2;			
	    return ex1;
	}


} // namespace cxsc

