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

/* CVS $Id: l_interval.cpp,v 1.21 2014/01/30 17:23:46 cxsc Exp $ */

#include <math.h>
#include "l_interval.hpp"
#include "idot.hpp"
#include "RtsFunc.h"

namespace cxsc {

#define CXSC_Zero 0.

l_interval & l_interval::operator= (const l_interval & a)    throw()
{
   real *newdata=new real[a.prec+1];
   int i;
   for(i=0;i<=a.prec;i++)
      newdata[i]=a.data[i];
   delete [] data;
   data=newdata;
   prec=a.prec;
   return *this;
}   

// ---- Typwandlungen ----

interval::interval(const l_interval & a) throw()
{
   idotprecision idot(0.0);
   a._akku_add(idot);
   *this=rnd(idot);
}

interval & interval::operator =(const l_interval & a) throw()
{
   idotprecision idot(0.0);
   a._akku_add(idot);
   return *this=rnd(idot);
}


interval _interval(const l_interval & a) throw()
{
   idotprecision idot(0.0);
   a._akku_add(idot);
   return rnd(idot);
}

l_interval::l_interval(const dotprecision & a)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   _allo(stagprec);
   idotprecision idot(a);
   _akku_out(idot);   
}

l_interval::l_interval(const dotprecision & a,const dotprecision & b)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const dotprecision&,const dotprecision&)"));
   _allo(stagprec);
   idotprecision idot;
   UncheckedSetInf(idot,a);
   UncheckedSetSup(idot,b);
   _akku_out(idot);   
}

l_interval & l_interval::operator =(const dotprecision & a)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   if(stagprec!=prec)
   {
      delete [] data;
      _allo(stagprec);
   }
   idotprecision idot(a);
   _akku_out(idot);   
   return *this;
}

l_interval::l_interval(const idotprecision & a)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   _allo(stagprec);
   idotprecision idot(a);
   _akku_out(idot);   
}

l_interval & l_interval::operator =(const idotprecision & a)
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC)
#else
   throw()
#endif
{
   if(stagprec!=prec)
   {
      delete [] data;
      _allo(stagprec);
   }
   idotprecision idot(a);
   _akku_out(idot);   
   return *this;
}

l_interval::l_interval(const l_real &a, const l_real &b) 
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(stagprec);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const l_real &a, const l_real &b)"));
   dotprecision dot1, dot2;
   dot1 = a;
   dot2 = b;
   idotprecision idot(dot1,dot2);
   _akku_out(idot);
}

l_interval::l_interval(const real &a, const l_real &b) 
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(stagprec);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const l_real &a, const l_real &b)"));
   dotprecision dot1, dot2;
   dot1 = a;
   dot2 = b;
   idotprecision idot(dot1,dot2);
   _akku_out(idot);
}

void l_realz2l_interval(const l_real& lr, const interval& z,
			                        l_interval& li) throw()
{  // converts lr+z to li of type l_interval; Blomquist, 12.10.06;
   // lr+z is included by li.
   // prec(lr) <= prec(li)  must be realized!
    int p = StagPrec(lr);  int q = StagPrec(li);
    if(p>q) 
    { 
	std::cerr << "l_realz2l_interval(const l_real& lr, const interval& z, l_interval& li): incorrect precisions of lr,li !" 
		  << std:: endl; 
	exit (1); 
    }

    for (int i=1; i<=q-1; i++) li[i] = 0; 
    li[q]   = Inf(z);
    li[q+1] = Sup(z);
    if(p<q) for (int i=1; i<=p; i++) li[i] = lr[i];
    else // p=q
    {
	p--;
	for (int i=1; i<=p; i++) li[i] = lr[i];
	li[q]   = addd(lr[p+1],Inf(z)); 
	li[q+1] = addu(lr[p+1],Sup(z));
    }
} // l_realz2l_interval

l_interval::l_interval(const l_real &a, const real &b) 
#if (CXSC_INDEX_CHECK)
   throw(ERROR_LINTERVAL_WRONG_STAGPREC,ERROR_LINTERVAL_EMPTY_INTERVAL)
#else
   throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
#endif
{
   _allo(stagprec);
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval::l_interval(const l_real &a, const l_real &b)"));
   dotprecision dot1, dot2;
   dot1 = a;
   dot2 = b;
   idotprecision idot(dot1,dot2);
   _akku_out(idot);
}


/*
l_interval _l_interval(const l_real & a, const l_real & b) throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
   if(a>b)
      cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("l_interval _l_interval(const l_real & a, const l_real & b)"));

   dotakku[0]=0.0;
   dotakku[1]=0.0;
   a._akku_add(dotakku[0]);
   b._akku_add(dotakku[1]);
   idotakku[0]=_idotprecision(dotakku[0],dotakku[1]);
   l_interval tmp;
   tmp._akku_out();
   return tmp;
}   
*/

/*
l_interval _unchecked_l_interval(const l_real & lr1, const l_real & lr2) throw()
{
   real inf, sup, tmp;          // fuer Naeherungen, entspricht Interval z
   int i=1;
   l_interval li;
   dotprecision dot1(0.0);
   dotprecision dot2(0.0);
   lr1._akku_add(dot1);
   lr2._akku_add(dot2);
   inf = rnd(dot1,RND_DOWN);
   sup = rnd(dot2,RND_UP);
   while (i<li.prec && !(inf<=0. && sup>=0.)) 
   {
      tmp = inf + (sup-inf)/2.0; // middle(interval)
      li.elem(i) = tmp;
      dot1 -= tmp;
      dot2 -= tmp;
      inf = rnd(dot1,RND_DOWN);
      sup = rnd(dot2,RND_UP);
      i++;
   }
   li.elem(li.prec)=inf;           // Intervall in die letzten Stellen
   li.elem(li.prec+1)=sup;         // schreiben

   return li;
}
*/

l_interval _unchecked_l_interval(const l_real & lr1, const l_real & lr2) throw()
{
   dotprecision dot1, dot2;
   dot1 = lr1;
   dot2 = lr2;
   idotprecision idot(dot1,dot2);
   l_interval li;
   li._akku_out(idot);  
   return li;
}


idotprecision _idotprecision(const l_interval & a) throw()
{
   return idotprecision(a);
}

idotprecision::idotprecision(const l_interval & a) throw() : inf(0),
                                                             sup(0)
{ 
   a._akku_add(*this);
}

idotprecision & idotprecision::operator =(const l_interval & a) throw()
{
   *this=0;
   a._akku_add(*this);
   return *this;
}

// ---- Standardfunkt ---- (arithmetische Operatoren)

l_interval operator-(const l_interval & a) throw()
{
   int precsave=stagprec;
   stagprec=a.prec;

   l_interval tmp;

   stagprec=precsave;

   int i;
   for(i=0;i<a.prec-1;i++)
      tmp.data[i]=-a.data[i];

   tmp.data[a.prec-1]=-a.data[a.prec];
   tmp.data[a.prec]=-a.data[a.prec-1];

   return tmp;
} 

// LI-LI

l_interval operator+(const l_interval & li1, const l_interval & li2) throw()
{
   l_interval li3;
   idotprecision idot(0.0);
   li1._akku_add(idot);
   li2._akku_add(idot);
   li3._akku_out(idot);
   return li3;
}

l_interval operator-(const l_interval & li1, const l_interval & li2) throw()
{
   l_interval li3;
   idotprecision idot(0.0);
   li1._akku_add(idot);
   li2._akku_sub(idot);
   li3._akku_out(idot);
   return li3;
} 
 
/* l_interval operator*(const l_interval & li1, const l_interval & li2) throw()
{
   l_interval li3;
   interval stdmul;
   
   stdmul = _interval(li1)*_interval(li2);

   if (abs(Inf(stdmul)) < MinReal) 
      li3 = _l_interval(stdmul) + 0.0; // Was das +0.0 soll ist mir ein Raetsel SR
   else 
   {               // if eingef�gt am 30.07.92 von Frederic Toussaint
      idotakku[0]=0.0;
      accumulate(idotakku[0], li1, li2);
      li3._akku_out();
      li3 = li3 & _l_interval(stdmul);  // Nachkorrektur
      // Wie bitte? (SR)
   }
   return li3;
}  */

l_interval operator * (const l_interval& li1, const l_interval& li2) throw()
{  // Blomquist, Neue Version, 21.11.02;
    l_interval li3;
    interval stdmul;
    stdmul = _interval(li1) * _interval(li2); // grobe Einschlie�ung
    idotprecision idot(0.0);
    accumulate(idot,li1,li2);
    li3._akku_out(idot);       // li3: meist feinere Einschlie�ung! Aber bei zu
             // breiten li1,li2 ist dieses Ergebnis li3 breiter als stdmul!
    li3 = li3 & _l_interval(stdmul);   // Das optimale Ergebnis liegt daher 
              // im Durchschnitt der groben und der feineren Einschlie�ung!
    return li3;
}    

l_interval operator/(const l_interval & li1, const l_interval & li2) throw(ERROR_LINTERVAL_DIV_BY_ZERO)
{ 
// ge�ndert am 29.07.92 von Frederic Toussaint
// 13.5.93 AW: neuer Algorithmus nach der Beschreibung von W. Kraemer

   int         i, j, k, Z_sign;
   real        dzmitte;
   interval    dn = _interval(li2),
               dz = _interval(li1),
               stddiv = dz/dn;
   idotprecision idot;
   dotprecision dot1,dot2;
   l_interval  li3;

   li3._clear(1);

   if(!li2) 
   {
      cxscthrow(ERROR_LINTERVAL_DIV_BY_ZERO("l_interval operator/(const l_interval & li1, const l_interval & li2)"));
   } else 
   {
      dz = dz/dn;
      if (stagprec > 1) 
      {         // bei == 1 passiert nix!
         idot = 0.0;          // X, also Zaehler
         li1._akku_add(idot);
         k = 1;
         Z_sign = (Inf(dz) > 0.0);
         do {
            dzmitte = (Inf(dz) + Sup(dz)) / 2.0;  // mid(interval); grob!
            if (!!dz)                           
               if (Sup(abs(dz)) > MinReal)
                  if (diam(dz) < 1e-14 * abs(dzmitte) ) 
                  {
                     li3.elem(k) = dzmitte;
                     // Bestimme Rest in dotakku[2], [3] (Inf, Sup):
                     //
                     dot1 = 0.0;
                     for (i=1; i<=k; i++)                   // Rumpf
                        for (j=1; j<=li2.prec-1; j++)
                           accumulate(dot1, -li2.elem(j), li3.elem(i));

                     dot2 = dot1;        // Rumpf ist fuer Inf, Sup gleich
                     dot1 += Inf(idot);  // Inf += Inf(X)
                     dot2 += Sup(idot);  // Sup += Sup(X)
                     if (Z_sign)                      // INCL(Z) > 0.0
                     for (i=1; i<=k; i++) 
                     {
                        accumulate(dot1, -li2.elem(li2.prec+1), li3.elem(i));
                        accumulate(dot2, -li2.elem(li2.prec),   li3.elem(i));
                     } else                             // INCL(Z) < 0.0
                        for (i=1; i<=k; i++) 
                        {
                           accumulate(dot1, -li2.elem(li2.prec),   li3.elem(i));
                           accumulate(dot2, -li2.elem(li2.prec+1), li3.elem(i));
                        }
                     //
                     // Rausrunden in dz
                     dz = _interval(rnd(dot1, RND_DOWN),
                                    rnd(dot2, RND_UP));
                     dz = dz / dn;
                  }
            k++;
         } while (k < stagprec);
      } // if (stagprec > 1)
      li3.elem(stagprec) = Inf(dz);     // Intervall in die letzten Stelle
      li3.elem(stagprec+1) = Sup(dz);
      li3 = li3 & _l_interval(stddiv);  // Nachkorrektur
   } // if keine Null im Nenner
   return li3;
}

// ---- Vergleichsop. ----
bool operator!(const l_interval & li) throw()
{
   idotprecision idot(0.0);
   li._akku_add(idot);
   return (!idot);
}

/*l_interval::operator void *(void) throw()
{
   idotakku[0]=0.0;
   _akku_add(idotakku[0]);
   return (idotakku[0]);
}*/

bool operator==(const l_interval & li1, const l_interval & li2) throw()
{
   idotprecision idot1(0.0);
   idotprecision idot2(0.0);
   li1._akku_add(idot1);
   li2._akku_add(idot2);
   return (idot1==idot2);
}   

// ---- Mengenvergleiche ----

bool operator<(const l_interval & li1, const l_interval & li2) throw()
{
   idotprecision idot1(0.0);
   idotprecision idot2(0.0);
   li1._akku_add(idot1);
   li2._akku_add(idot2);
   return (idot1<idot2);
} 

bool operator>(const l_interval & li1, const l_interval & li2) throw()
{
   idotprecision idot1(0.0);
   idotprecision idot2(0.0);
   li1._akku_add(idot1);
   li2._akku_add(idot2);
   return (idot1>idot2);
} 

bool operator<=(const l_interval & li1, const l_interval & li2) throw()
{
   idotprecision idot1(0.0);
   idotprecision idot2(0.0);
   li1._akku_add(idot1);
   li2._akku_add(idot2);
   return (idot1<=idot2);
} 

bool operator>=(const l_interval & li1, const l_interval & li2) throw()
{ 
   idotprecision idot1(0.0);
   idotprecision idot2(0.0);
   li1._akku_add(idot1);
   li2._akku_add(idot2);
   return (idot1>=idot2);
}

void ConvexHull(const l_interval & li1, const l_interval & li2, l_interval & li3, l_interval & li4) throw()
{
   if(li1<=li2) 
   {                                      // Trivialfall 1
      li3=li2;
      li4=li2;
   } else if(li2<=li1) 
   {                                 // Trivialfall 2
      li3=li1;
      li4=li1;
   } else 
   {                                              // rechne
      idotprecision idot1(0.0);
      idotprecision idot2(0.0);
      li1._akku_add(idot1);
      li2._akku_add(idot2);

      idot1 |= idot2;
      // nun steht das richtige Ergebnis in idotakku[0]
      idot2 = idot1;                              // sichern
      li3._akku_out_inn(idot1);                              // Rundung nach innen
      idot1 = idot2;                              // und wiederherstellen
      li4._akku_out(idot1);                                  // ... nach aussen
   }
}

void Intersection(const l_interval & li1, const l_interval & li2, l_interval & li3, l_interval & li4) throw(ERROR_LINTERVAL_EMPTY_INTERVAL)
{
   if(li1<=li2) 
   {                                      // Trivialfall 1
      li3=li1;
      li4=li1;
   } else if(li2<=li1) 
   {                                 // Trivialfall 2
      li3=li2;
      li4=li2;
   } else 
   {                                              // rechne
      idotprecision idot1(0.0);
      idotprecision idot2(0.0);
      idotprecision idot;
      li1._akku_add(idot1);
      li2._akku_add(idot2);
      try
      {
         idot=(idot1 &= idot2);
      }
      catch(const EMPTY_INTERVAL &)
      {
         cxscthrow(ERROR_LINTERVAL_EMPTY_INTERVAL("void Intersection(const l_interval & li1, const l_interval & li2, l_interval & li3, l_interval & li4)"));
      }
      // nun steht das richtige Ergebnis in idotakku[0]
      li3._akku_out_inn(idot);                              // Rundung nach innen
      idot = idot1;
      li4._akku_out(idot);                                  // ... nach aussen
   }
}

l_real     mid  (const l_interval & li) throw()
{
   l_real lr;

   // --------------------------------------------------------
   //   dotakku[0] = Inf(li) + Sup(li)

   dotprecision dot(0.0);
   for (int j=1; j<=li.prec-1; j++)
      dot += li.elem(j);

   dot += dot; // mal 2
   dot += li.elem(li.prec);
   dot += li.elem(li.prec+1);

   //  ------------------------------------------------------
   //  Division nur bei ungleich 0

   if (dot != 0.0) 
   {
      Dotprecision ptr = *dot.ptr();

      // --------------------------------------------------------
      //  Dividiere dotakku[0] durch 2, mittels 1 Rechtsshift

      ptr[(a_intg)++ptr[A_END]] = ZERO;
          
      b_shr1 (&ptr[(a_intg)ptr[A_BEGIN]], (a_intg)(ptr[A_END]-ptr[A_BEGIN]+1));
          
      if (ptr[(a_intg)ptr[A_END]]   == ZERO) 
         --ptr[A_END];
      if (ptr[(a_intg)ptr[A_BEGIN]] == ZERO) 
         ++ptr[A_BEGIN];

      // --------------------------------------------------------
   }

   lr._akku_out(dot);

   return lr;
}
 
/* void accumulate(idotprecision & d, const l_interval & li1, const l_interval & li2) throw()
{
   // �nderungen am 24.9.92 von F. Toussaint wegen Underflow-Fehlern

   int  i, j, n = 0;
   real tmp;

   for (i=1; i<=li1.prec-1; i++)
      for (j=1; j<=li2.prec-1; j++) 
      {
         tmp = abs(li1.elem(i)*li2.elem(j));
         if (tmp < MinReal && tmp != CXSC_Zero) 
            n++;
         else 
            accumulate(d, _interval(li1.elem(i)), _interval(li2.elem(j)));
      }
      for (i=1; i<=li2.prec-1; i++) 
      {
         tmp = abs(Inf(_interval(li1.elem(li1.prec), li1.elem(li1.prec+1)) * _interval(li2.elem(i))));
         if (tmp < MinReal && tmp != CXSC_Zero) 
            n++;
         else 
            accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                          _interval(li2.elem(i)));
      }
      for (i=1; i<=li1.prec-1; i++) 
      {
         tmp = abs(Inf(_interval(li2.elem(li2.prec), li2.elem(li2.prec+1)) * _interval(li1.elem(i))));
         if (tmp < MinReal && tmp != CXSC_Zero) 
            n++;
         else 
            accumulate(d, _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)),
                          _interval(li1.elem(i)));
      }
      tmp = abs(Inf(_interval(li1.elem(li1.prec), li1.elem(li1.prec+1)) *
                    _interval(li2.elem(li2.prec), li2.elem(li2.prec+1))));
      if (tmp < MinReal && tmp != CXSC_Zero) 
         n++;
      else 
         accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                       _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)));
      if (n > 0) 
         accumulate(d, _interval(_real(n)), 
                       _interval(-MinReal, MinReal));
} */

void accumulate(idotprecision & d, const l_interval & li1, 
                                   const l_interval & li2) throw()
{   // Blomquist, Neue Version vom 21.11.02;
    // Die �nderungen von Toussaint wurden r�ckg�ngig gemacht.
   int  i,j;

   for (i=1; i<=li1.prec-1; i++)
      for (j=1; j<=li2.prec-1; j++) 
      {
            accumulate(d, _interval(li1.elem(i)), _interval(li2.elem(j)));
      }
      for (i=1; i<=li2.prec-1; i++) 
      {
            accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                          _interval(li2.elem(i)));
      }
      for (i=1; i<=li1.prec-1; i++) 
      {
            accumulate(d, _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)),
                          _interval(li1.elem(i)));
      }
      accumulate(d, _interval(li1.elem(li1.prec), li1.elem(li1.prec+1)),
                    _interval(li2.elem(li2.prec), li2.elem(li2.prec+1)));
}


/* void l_interval::_akku_out() throw()
{  
   // ein im idotakku[0] liegendes Zwischenergebnis
   // wird in der entsprechenden Praezision in das aufrufende l_interval
   // gerundet. Rundung, also Einschluss von aussen!!
   // Alle MinReal wurden ge�ndert in minreal, Blomquist 19.11.02;
   _clear(1);
   interval z;
   real tmp, r, s;
   int i = 1,
   
   n = 0;
   z = rnd(idotakku[0]);

   while (i<prec && !(CXSC_Zero <= z)) 
   {       
      // z ist Intervall, deshalb <=
      if (succ(succ(Inf(z))) < Sup(z))
         break;                            // bei zu grobem Einschluss:
                                           // Intervall direkt in
                                           // Einschlusskomponente
                                           // schreiben
      tmp = Inf(z) + (Sup(z)-Inf(z))/2.0; // middle(interval)
      if (abs(Inf(z)) < minreal) 
      {
         if (abs(Sup(z)) < minreal) 
         {
            if (tmp != 0.0) 
            {
               tmp = 0.0;
               n++;
            }
            i = prec;
         }
      }
          
      this->elem(i) = tmp;
      idotakku[0] -= tmp;
      z = rnd(idotakku[0]);
      i++;
   } // while
   r = Inf(z);
   if (abs(r) < minreal) 
   {
      if (r < 0.0) 
         r = -minreal;
      else         
         r = 0.0;
   }
   s = Sup(z);
   if (abs(s) < minreal) 
   {
      if (s <= 0.0) 
         s = 0.0;
      else          
         s = minreal;
   }
   if (n > 0) 
   {  
      this->elem(prec) = r-_real(n+1)*minreal;   // Intervall in die letzten
      this->elem(prec+1) = s+_real(n+1)*minreal; // Stellen schreiben
      // (n+1), da Rundungsrichtung nicht bestimmt werden kann!
   } else 
   {
      this->elem(prec) = r;           // Intervall in die letzten Stellen
      this->elem(prec+1) = s;         // schreiben
   }
}  */

void l_interval::_akku_out(idotprecision& idot) throw()
{  
   // ein im idotakku[0] liegendes Zwischenergebnis
   // wird in der entsprechenden Praezision in das aufrufende l_interval
   // so gerundet, dass idotakku[0] eingeschlossen wird.
   // Neueste Version: Blomquist 21.07.03;
   _clear(1);
   interval z;
   real tmp;
   int i = 1;
   
   z = rnd(idot);
   while (i<prec && !(CXSC_Zero <= z)) 
   {       
      // z ist Intervall, deshalb <=
      if (succ(succ(Inf(z))) < Sup(z))
         break;  // bei zu grobem Einschluss: Intervall direkt in
                 // Einschlusskomponente schreiben.
      tmp = mid(z);  // middle(interval)
      this->elem(i) = tmp;
      idot -= tmp;
      z = rnd(idot);
      i++;
   } // while
   this->elem(prec) = Inf(z);      // Intervall in die letzten Stellen
   this->elem(prec+1) = Sup(z);    // schreiben
} // _akku_out()

void l_interval::_akku_out_inn(idotprecision& idot) throw()
{ 
   // ein im idotakku[0] liegendes Zwischenergebnis
   // wird in der entsprechenden Praezision in die aufrufende l_interval Zahl
   // li so gerundet, dass gilt:  li enthalten in idotakku[0].
   _clear(1);
   real inf, sup, tmp;          // fuer Naeherungen, entspricht Interval z
   int i=1;
   inf = rnd(Inf(idot),RND_UP);
   sup = rnd(Sup(idot),RND_DOWN);
   
   if (inf>sup) 
      inf = sup;      // Vorsichtsmassnahme
        
   while (i<prec && !(inf<=CXSC_Zero&&sup>=CXSC_Zero)) 
   {
      tmp = inf + (sup-inf)/2.0; // middle(interval)
      this->elem(i) = tmp;
      idot -= tmp;
      inf = rnd(Inf(idot),RND_UP);
      sup = rnd(Sup(idot),RND_DOWN);
      if (inf>sup) 
         inf = sup;    // Vorsichtsmassnahme
      i++;
   }
   this->elem(prec)=inf;           // Intervall in die letzten Stellen
   this->elem(prec+1)=sup;         // schreiben
}

/* void l_interval::_akku_add(idotprecision& d) const throw()
{ 
   // addiert aufrufenden l_interval auf iakku d.
   // �nderung am 23.9.92 von F. Toussaint, da Fehler im Unterlaufbereich
   // Ausgeblendet von Blomquist, 20.11.02, da Fehler im Unterlaufbereich
   // nicht erkennbar. Meine neue Version direkt anschliessend!
   int        n = 0;
   real       r, s;

   for (int i=1; i<=prec-1; i++) 
   {
      r = this->elem(i);
      if (abs(r) < MinReal) 
      {
         if (r != 0.0) 
            n++;
      } else 
         d += _interval(r);
   }
   r = this->elem(prec);
   if (abs(r) < MinReal) 
   {
      if (r < 0.0) 
         r = -MinReal;
      else         
         r = 0.0;
   }
   s = this->elem(prec+1);
   if (abs(s) < MinReal) 
   {
      if (s <= 0.0) 
         s = 0.0;
      else          
         s = MinReal;
   }
   if (r != CXSC_Zero || s != CXSC_Zero) 
      d += _interval(r, s);
   if (n > 0) 
   {
      d += _interval(-_real(n+1)*MinReal, _real(n+1)*MinReal);
      // (n+1), da Rundungsrichtung nicht bestimmt werden kann!
   }
} */

void l_interval::_akku_add(idotprecision& d) const throw()
{ 
    // Addiert aufrufendes Intervall vom Typ l_interval auf d.
    // Meine neue Version; Blomquist, 20.11.02;
   real r,s;
   for (int i=1; i<=prec-1; i++)
   {   
       r = this->elem(i);
       if (sign(r) != CXSC_Zero) // Addition nur, falls n�tig!
	   d += _interval(r); 
   }
   r = this->elem(prec);
   s = this->elem(prec+1);
   if (r != CXSC_Zero || s != CXSC_Zero) // Addition nur, falls n�tig!
      d += _interval(r, s);
}

/* void l_interval::_akku_sub(idotprecision& d) const throw()
{ 
   // Subtrahiert aufrufendes Intervall vom Typ l_interval von d.
   // Intervallsubtraktion!!
   // �nderung am 23.9.92 von F. Toussaint, da Fehler im Unterlaufbereich
   // Blomquist: Mir sind keine solchen Fehler bekannt, daher wurde diese
   // Version von mir ausgeklammert, 20.11.02; Neue Version anschliessend!

   int        n = 0;
   real       r, s;

   for (int i=1; i<=prec-1; i++) 
   {
      r = this->elem(i);
      if (abs(r) < MinReal) 
      {
         if (r != 0.0) 
            n++;
      } else 
         d -= _interval(r);
   }
   r = this->elem(prec);
   if (abs(r) < MinReal) 
   {
      if (r < 0.0) 
         r = -MinReal;
      else         
         r = 0.0;
   }
   s = this->elem(prec+1);
   if (abs(s) < MinReal) 
   {
      if (s <= 0.0) 
         s = 0.0;
      else          
         s = MinReal;
   }
   d -= _interval(r, s);
   if (n > 0) 
   {
      d -= _interval(-_real(n+1)*MinReal, _real(n+1)*MinReal);
      // (n+1), da Rundungsrichtung nicht bestimmt werden kann!
   }
} */

void l_interval::_akku_sub(idotprecision& d) const throw()
{ 
    // Subtrahiert aufrufendes Intervall vom Typ l_interval von d.
    // Meine neue Version; Blomquist, 20.11.02;
   real r,s;

   for (int i=1; i<=prec-1; i++)
   { 
       r = this->elem(i);
       if (sign(r) != CXSC_Zero) // Subtraktion nur, falls n�tig!
	   d -= _interval(r); 
   }
   r = this->elem(prec);
   s = this->elem(prec+1);
   if (r != CXSC_Zero || s != CXSC_Zero) // Subtraktion nur, falls n�tig!
      d -= _interval(r, s);
}

// ---- Ausgabefunkt. ---------------------------------------

std::ostream & operator << (std::ostream &s, const l_interval & a) throw()
{
   idotprecision idot(0.0);
   a._akku_add(idot);
   s << idot;
   return s;
}
std::string & operator << (std::string &s, const l_interval & a) throw()
{
   idotprecision idot(0.0);
   a._akku_add(idot);
   s << idot;
   return s;
}

std::istream & operator >> (std::istream &s, l_interval & a) throw()
{
   idotprecision idot;
   s >> idot;
   a._akku_out(idot);
   return s;
}

std::string & operator >> (std::string &s, l_interval & a) throw()
{
   idotprecision idot;
   s >> idot;
   a._akku_out(idot);
   return s;

}

void operator >>(const std::string &s,l_interval & a) throw()
{
   std::string r(s);
   idotprecision idot;
   r >> idot;
   a._akku_out(idot);
}
void operator >>(const char *s,l_interval & a) throw()
{
   std::string r(s);
   idotprecision idot;
   r>>idot;
   a._akku_out(idot);
}

int in ( const real& x, const l_interval& y )         // Contained-in relation
  { return ( (Inf(y) <= x) && (x <= Sup(y)) ); }      //----------------------

int in ( const l_real& x, const l_interval& y )       // Contained-in relation
  { return ( (Inf(y) <= x) && (x <= Sup(y)) ); }      //----------------------

int in ( const interval& x, const l_interval& y )     // Contained-in relation
{                                                     //----------------------
  return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) );
}

int in ( const l_interval& x, const l_interval& y )   // Contained-in relation
{                                                     //----------------------
  return ( (Inf(y) < Inf(x)) && (Sup(x) < Sup(y)) );
}

l_interval Blow (const l_interval& x, const real& eps )
{  
    int n;
    l_interval y;
    l_real lr1,lr2;
    y = x + interval(-eps,eps) * diam(x);
    lr1 = Inf(y);
    n = StagPrec(lr1);
    lr1[n] = pred(lr1[n]);
    lr2 = Sup(y);
    lr2[n] = succ(lr2[n]);
    return l_interval(lr1,lr2);
}

int Disjoint (const l_interval& a, const l_interval& b ) 
                                                 // Test for disjointedness
{                                                //------------------------
    return (Inf(a) > Sup(b)) || (Inf(b) > Sup(a));
}

l_real AbsMin ( const l_interval& x )              // Absolute minimum of
{                                                  // an interval
   if ( in(0.0,x) ) return l_real(0.0);
   else 
   { 
       l_real y(Inf(x));
       if (y > 0.0) return y;
       else return -Sup(x);
   }

}

l_real AbsMax (const l_interval& x )            // Absolute maximum of
{                                               // an interval
    l_real a, b;                                //--------------------

    a = abs(Inf(x));
    b = abs(Sup(x));

    if (a > b)
	return a;
    else
	return b;
}

l_real RelDiam ( const l_interval x )                     // Relative diameter
{                                                         // of an interval
  if ( in(0.0,x) )                                        //------------------
    return diam(x);
  else
    return( Sup( l_interval(diam(x)) / l_interval(AbsMin(x)) ) );
}

inline void times2pown(l_interval& x, int n) throw()
{ // x = x * 2^n;  Blomquist, 28.11.02;
    real mt,t;
    interval z;
    int p = StagPrec(x);
    if ( n<LI_Min_Exp_ || n>LI_maxexpo1 ) 
    { std::cerr << "Error in:  " 
           << "times2pown(l_interval& x, const int& n): " << std::endl
           << " -1074 <= n <= +1023 not fulfilled" << std::endl; exit(0); 
    }
    real F = comp(0.5,n+1);
    z = _interval(x[p],x[p+1]);
    times2pown(z,n);  // Scaling the original interval;

    for (int i=1; i<=p-1; i++)
    {
	mt = mant(x[i]);
	t = x[i];
	times2pown(x[i],n);
	if ( mt != mant(x[i]) ) 
	{
	    x[i] = 0;
	    z += _interval(t) * F;
        }
    }
    x[p]   = Inf(z);
    x[p+1] = Sup(z);
} // times2pown(...)


void Times2pown(l_interval& a, const real& p) throw()
// The first parameter delivers an inclusion of a * 2^p;
// For p in [-2100,+2100] p must be an integer value.
// This condition is NOT tested in this function!
// For p outside [-2100,+2100] an inclusion of a * 2^p is
// calculated for any p of type real, unless an overflow occurs.
// If the function is activated with the second parameter of type int, 
// then the first parameter delivers correct inclusions of a * 2^p, 
// unless an overflow occurs.
// Blomquist, 28.01.2008;
{
    const int c2 = 2100;
    int ex(expo_gr(a)),fac,rest,n;
    double dbl;

    if (ex > -1000000)
    {
	if (p>=0)
	    if (p>c2)
		times2pown(a,c2); // Produces an error
	    else // 0 <= p <= 2100
	    {
		dbl = _double(p);
		n = (int) dbl;
		fac = n/LI_maxexpo1;
		rest = n%LI_maxexpo1;
		for (int k=1; k<=fac; k++)
		    times2pown(a,LI_maxexpo1);
		times2pown(a,rest);
	    }
	else // p<0; No overflow or underflow!
	    if (p<-c2)
	    {
		if (0<a)
		    a = l_interval(-minreal,minreal);
		else
		    if (Inf(a)>=0)
			a = l_interval(0,minreal);
		    else a = l_interval(-minreal,0);
	    }
	    else // -2100 <= p < 0
	    {
		dbl = _double(p);
		n = (int) dbl;
		fac = n/LI_Min_Exp_;
		rest = n%LI_Min_Exp_;
		for (int k=1; k<=fac; k++)
		    times2pown(a,LI_Min_Exp_);
		times2pown(a,rest);
	    }
    }
} // Times2pown(...)

} // namespace cxsc

