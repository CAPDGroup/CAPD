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

/* CVS $Id: dot.cpp,v 1.28 2014/01/30 17:23:44 cxsc Exp $ */

#include "dot.hpp"
#include "real.hpp"
#include "interval.hpp"

#include "RtsTyp.h" 
#include "RtsFunc.h"
#include "ioflags.hpp"

// Fuer memcpy... Kann auch <string.h> sein.
#include <memory.h> 
// #include <sstream>

namespace cxsc {

#ifdef CXSC_USE_TLS_PREC

#ifdef _WIN32
__declspec(thread) unsigned int opdotprec = 0;
#else
__thread unsigned int opdotprec = 0;
#endif

#else

unsigned int opdotprec = 0;

#endif

//----------------------------------------------------------------------
// global verf�gbare Dotprecision Variablen
//
//  dotakku[0..3] stehen fuer Matrix, Langzahl u.a. Pakete zur
//                Verfuegung
//  dotakku[4]    wird in den skalaren Paketen intern verwendet
//dotprecision dotakku[MAXDOTAKKU];

dotprecision::dotprecision(void) throw() : akku(new a_btyp[A_LENGTH]), err(0.0), k(0)
{
   memset(akku,0,BUFFERSIZE); // d_clr(&akku);
   // d_clr(&akku); // Da d_init calloc verwendet, nicht noetig!?!?
   // d_clr  ist definiert in p88rts.h, also RtsTyp.h
}

dotprecision::dotprecision(const dotprecision &from) throw()
                                 : akku(new a_btyp[A_LENGTH]), err(from.err), k(from.k)
{
   memcpy((void *)akku,(void *)from.akku,BUFFERSIZE);
}

dotprecision::dotprecision(const real &r) throw()
                                          : akku(new a_btyp[A_LENGTH]), err(0.0), k(0)
{
   real dummy(r);
   memset(akku,0,BUFFERSIZE); // d_clr(&akku);

   d_radd(&akku,*((a_real *)&dummy));   
}

dotprecision & dotprecision::operator =(const dotprecision &from) throw()
{
   if(&from==this) // Handle self-assignment
      return *this; // Hier kein new, bzw delete, deswegen diese Methode

   memcpy((void *)akku,(void *)from.akku,BUFFERSIZE);
   err = from.err;
   //k = from.k;

   return *this;
}

dotprecision & dotprecision::operator =(const real &r) throw()
{
   real dummy(r);
   memset(akku,0,BUFFERSIZE); 
   d_radd(&akku,*((a_real *)&dummy));   // KRANK!
   err = 0.0;
   return *this;
}

dotprecision::~dotprecision()
{
   delete [] akku;
}

dotprecision operator +(const dotprecision &a,const dotprecision &b) throw()
{
   dotprecision tmp(a); 
   d_dadd(&(tmp.akku),(Dotprecision)b.akku);
   tmp.err = addu(a.err,b.err);
   //tmp.k = 0;
   return tmp;
}

dotprecision operator -(const dotprecision &a,const dotprecision &b) throw()
{
   dotprecision tmp(b);
   tmp.negdot();
   d_dadd(&(tmp.akku),(Dotprecision)a.akku);
   tmp.err = addu(a.err,b.err);
   //tmp.k = 0;
   return tmp;
}

dotprecision operator +(const dotprecision &d,const real &r) throw()
{
   dotprecision erg(d);
   d_radd(&(erg.akku),*((a_real *)&r));
   return erg;
}

dotprecision operator +(const real &r,const dotprecision &d) throw()
{
   dotprecision erg(d);
   d_radd(&(erg.akku),*((a_real *)&r));
   return erg;
}

dotprecision operator -(const dotprecision &d,const real &r) throw()
{
   dotprecision erg(d);
   d_radd(&(erg.akku),-*((a_real *)&r));
   return erg;
}

dotprecision operator -(const real &r,const dotprecision &d) throw()
{
   dotprecision erg(d);
   erg.negdot();
   d_radd(&(erg.akku),*((a_real *)&r));
   return erg;
}

dotprecision & operator +=(dotprecision &d,const real &r) throw()
{
   d_radd(&d.akku,*((a_real *)&r));
   return d;
}

dotprecision & operator +=(dotprecision &d,const dotprecision &d2) throw()
{
   d_dadd(&(d.akku),(Dotprecision)(d2.akku));
   d.err = addu(d.err,d2.err);
   return d;
}

dotprecision & operator -=(dotprecision &d,const real &r) throw()
{
   d_radd(&d.akku,-*((a_real *)&r));
   return d;
}

dotprecision & operator -=(dotprecision &d,const dotprecision &d2) throw()
{
   dotprecision tmp(d2);
   tmp.negdot();
   d_dadd(&(d.akku),(Dotprecision)(tmp.akku));
   d.err = addu(d.err,d2.err);
   return d;
}

dotprecision operator -(const dotprecision &a) throw()
{
   dotprecision tmp(a);
   tmp.negdot();
   return tmp;
}

dotprecision operator +(const dotprecision &a) throw()
{
   return a;
}

// ---------------------------------------------------------------------------
// Alle Vergleichsoperatoren werden auf die folgenden zwei (==, <=) zurueckgefuehrt!

bool operator ==(const dotprecision &a,const dotprecision &b)  throw()  // Uebernommen aus dot.cpp "testequal"
{
   int res = true;

   int leftsign = sign(a);
   if (leftsign != sign(b) )
      res = false;
   else 
   {
      int leftstart  = (a_intg)a.akku[A_BEGIN];
      int leftend    = (a_intg)a.akku[A_END];
      int rightstart = (a_intg)b.akku[A_BEGIN];
      int rightend   = (a_intg)b.akku[A_END];

      if (leftend < rightstart || rightend < leftstart) 
         res = false;
      else 
      {
         res = true;
         while (res && leftstart < rightstart)
            res = (a.akku[leftstart++] == ZERO);

         while (res && rightstart < leftstart)
            res = (b.akku[rightstart++] == ZERO);

         while (res && leftstart <= leftend && leftstart <= rightend) 
         {
            res = (a.akku[leftstart++] == b.akku[rightstart++]);
         }

         while (res && leftstart <= leftend)
            res = (a.akku[leftstart++] == ZERO);
         while (res && rightstart <= rightend)
            res = (b.akku[rightstart++] == ZERO);
      }
   } 

   return res && (a.err == b.err);
}

bool operator <=(const dotprecision &x,const dotprecision &y) throw()   // Uebernommen aus dot.cpp "testlessequal"
{
   int res = true, cont;
   // Dotprecision dotakku = ((dotprecision*) &dot)->akku; b.akku

   dotprecision a = x + x.err;
   dotprecision b = y - y.err;

   int leftsign = sign(a);
   int rightsign = sign(b);
   if (leftsign != rightsign) 
      res = (leftsign < rightsign);
   else if (leftsign == 0) 
      res = true;
   else 
   {
      int leftstart  = (a_intg)a.akku[A_BEGIN];
      int leftend    = (a_intg)a.akku[A_END];
      int rightstart = (a_intg)b.akku[A_BEGIN];
      int rightend   = (a_intg)b.akku[A_END];

      if (leftend < rightstart) 
         res = (leftsign == -1);
      else if (rightend < leftstart) 
         res = (leftsign != -1);
      else 
      {
         cont = true;

         // ---------------------------------------------------------
         //  hat a Mantissenelemente vor b ==> false

         while (cont && leftstart < rightstart) 
         {
            cont = (a.akku[leftstart++] == ZERO);
            if (!cont) 
               res = false;
         }
         // ---------------------------------------------------------
         //  hat b Mantissenelemente vor a ==> true

         while (cont && rightstart < leftstart) 
         {
            cont = (b.akku[rightstart++] == ZERO);
            if (!cont) 
               res = true;
         }

         // ---------------------------------------------------------
         //  ein gemeinsames Mantissenelement a <= b  ?

         while (cont && leftstart <= leftend && leftstart <= rightend) 
         {
            cont = (a.akku[leftstart] == b.akku[rightstart]);
            if (!cont) 
            {
               res = (a.akku[leftstart] <= b.akku[rightstart]);
            }
            leftstart++,rightstart++;
         }

         // ---------------------------------------------------------
         //  hat a weitere Mantissenelemente ==> false

         while (cont && leftstart <= leftend) 
         {
            cont = (a.akku[leftstart++] == ZERO);
            if (!cont) 
               res = false;
         }
         // ---------------------------------------------------------
         //  hat b weitere Mantissenelemente ==> true

         while (cont && rightstart <= rightend) 
         {
            cont = (b.akku[rightstart++] == ZERO);
            if (!cont) 
               res = true;
         }

         if (cont)            // --------------------------------
            res = true;       // Mantissen waren gleich
         else                 // --------------------------------
            if (leftsign == -1)// Mantissen waren unterschiedlich
               res = !res;    // negiere Vergleichsresultat falls
                              // Akku's negativ waren
      }
   }

   return res;
}

bool operator !=(const dotprecision &a,const dotprecision &b) throw() { return !(a==b); }
bool operator  <(const dotprecision &a,const dotprecision &b) throw() { return !(b<=a); }
bool operator  >(const dotprecision &a,const dotprecision &b) throw() { return !(a<=b); }
bool operator >=(const dotprecision &a,const dotprecision &b) throw() { return (b<=a);  }

bool operator ==(const real &r,const dotprecision &d) throw() { return dotprecision(r)==d; }
bool operator !=(const real &r,const dotprecision &d) throw() { return dotprecision(r)!=d; }
bool operator  <(const real &r,const dotprecision &d) throw() { return dotprecision(r)<d;  }  
bool operator  >(const real &r,const dotprecision &d) throw() { return dotprecision(r)>d;  }
bool operator <=(const real &r,const dotprecision &d) throw() { return dotprecision(r)<=d; }
bool operator >=(const real &r,const dotprecision &d) throw() { return dotprecision(r)>=d; }

bool operator ==(const dotprecision &d,const real &r) throw() { return d==dotprecision(r); }
bool operator !=(const dotprecision &d,const real &r) throw() { return d!=dotprecision(r); }
bool operator  <(const dotprecision &d,const real &r) throw() { return d<dotprecision(r);  }
bool operator  >(const dotprecision &d,const real &r) throw() { return d>dotprecision(r);  }
bool operator <=(const dotprecision &d,const real &r) throw() { return d<=dotprecision(r); }
bool operator >=(const dotprecision &d,const real &r) throw() { return d>=dotprecision(r); }

bool operator  !(const dotprecision &a) throw()               { return sign(a)==0; }


void rnd (const dotprecision& d, real& r, rndtype type) throw() 
{
   if(type==RND_NEXT) {
      *((a_real *)(&r))=d_stan((Dotprecision)d.akku);
   } else if(type==RND_UP) {
      *((a_real *)(&r))=d_stau((Dotprecision)d.akku);
      r = addu(r,d.err);
   } else {
      *((a_real *)(&r))=d_stad((Dotprecision)d.akku); 
      r = subd(r,d.err);
   }
}

void rnd (const dotprecision& d, real& l, real& u) throw() 
{
   *((a_real *)(&l))=d_stad(d.akku);
   *((a_real *)(&u))=d_stau(d.akku);

   l = subd(l,d.err);
   u = addu(u,d.err);
}

void rnd (const dotprecision& d, interval& x) throw() 
{
    real a,b;
    rnd(d,a,b);
    x = interval(a,b);
}

real rnd (const dotprecision& d, rndtype type) throw() 
{
    real r;
   if(type==RND_NEXT) {
      *((a_real *)(&r))=d_stan((Dotprecision)d.akku);
   } else if(type==RND_UP) {
      *((a_real *)(&r))=d_stau((Dotprecision)d.akku);
      r = addu(r,d.err);
   } else {
      *((a_real *)(&r))=d_stad((Dotprecision)d.akku); 
      r = subd(r,d.err);
   }
   return r;
}

/*std::string  & operator << (std::string  & s, const dotprecision & d) throw() 
{
   //std::ostringstream o(s);
   
   if (ioflags.isset(IOFlags::realformat))
   {   
      
      real rl,ru; 
      rnd (d, rl,ru);    
      s += "dot(";
      //s << SaveOpt << RndDown;
      // o << rl;
      //sprintf (&sh[strlen(sh)] << rl, ", "); sh << RndUp;
      //
      s+=",";
      //      sprintf (&sh[strlen(sh)] << ru, ")");  sh << RestoreOpt;
      // o << ru;
      
      s+=")";
      //          strcpy (s, sh);
   } else
   {
      //real r=rnd(d);

      //o << r;
      s+="<zahl>";
             
      / *
      rndtype rnd;
      
      unsigned long length, formatflag, addblanks;
      unsigned long FracDigits = dotdigits;
      
      if (ioflags.isset(ioflags::rndup)) rnd = RND_UP;
      else if (ioflags.isset(ioflags::rnddown)) rnd = RND_DOWN;
      else rnd = RND_NEXT;
      
      if (ioflags.isset(IOFlags::variable))
         formatflag = dotwidth; 
      else if (ioflags.isset(IOFlags::varfixwidth))
         formatflag = dotwidth, FracDigits = -FracDigits;
      else
         formatflag = (ioflags.isset(IOFlags::fixed)) ? 0 : -1;
         
      // d_outp (str = dm, this->akku, formatflag, digits, rnd, &length);
      {
         long   dexpo,bdp,len;
         long   expo,i,digits,IntDigits,DecPlaces,vz;
         long   HoldWidth = (FracDigits < 0);
         
         if (HoldWidth) FracDigits = -FracDigits;
             
         // Kehre noetigenfalls Rundungsrichtung um
         if (d.sign() < 0) 
         {
              rnd = -rnd;
         }
         
         bdp = 2+A_I_DIGITS;
         len = bdp+1;
         if (c[A_END] > A_D_P) len += B_LENGTH * ((a_intg)c[A_END] - A_D_P);

      
      }
      * /
          
   }
   //s=o.str();                             
   return s;
}


std::ostream & operator << (std::ostream & o, const dotprecision & d) throw() 
{
   std::string buffer;
   buffer << d;
   o << d;
   return o;
   
/ *
   // nur real-Ausgabe   
   real rl,ru;
   *((a_real *)(&rl))=d_stad(d.akku);
   *((a_real *)(&ru))=d_stau(d.akku);

   o << "dot(" << rl << "," << ru << ") " << ru-rl;
* /
}

std::istream & operator >> (std::istream & i,dotprecision & d) throw() 
{
   // Dies ist noch _schlecht_!
   double b;
   i >> b;
   d=b;
   return i;
}
*/
dotprecision & dotprecision::negdot(void) throw() 
{
   this->akku[A_SIGN] = 1 - this->akku[A_SIGN];
   return *this;
}

int sign(const dotprecision &a) throw() 
{
   if(a.akku[A_BEGIN]==ZERO)
      return 0;
   return (a.akku[A_SIGN] ? -1 : 1);
}

dotprecision abs(const dotprecision &a) throw() 
{
   if(sign(a)>=0)
      return a;
   dotprecision tmp(a);
   tmp.negdot();
   return tmp;
}

dotprecision & accumulate (dotprecision& dot, const real& a, const real& b) throw() 
{
   // nur dann accumulieren, wenn noetig
   if (!!a && !!b)
      d_padd (&dot.akku,*(a_real*)&a,*(a_real*)&b);
   return dot;  
}

} // namespace cxsc

