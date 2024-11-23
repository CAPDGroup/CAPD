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

/* CVS $Id: l_real.cpp,v 1.34 2014/01/30 17:23:46 cxsc Exp $ */

#include <memory.h>

#include "l_real.hpp"
#include "dot.hpp"
#include "idot.hpp"
#include "interval.hpp"
#include "l_interval.hpp"

namespace cxsc {

#ifdef CXSC_USE_TLS_PREC

#ifdef _WIN32
__declspec(thread) int stagprec = 2;
#else
__thread int stagprec = 2;
#endif

#else

int stagprec = 2;

#endif


interval::interval(const l_real & a) throw()
{
   dotprecision dot(a);
   rnd(dot,inf,sup);
}

interval & interval::operator =(const l_real & a) throw()
{
   dotprecision dot(a);
   rnd(dot,inf,sup);
   return *this;
}

idotprecision::idotprecision(const l_real & a) throw() : inf(0),
                                                         sup(0)
{ 
   *this+=a;
}

idotprecision::idotprecision(const l_real & a,const l_real & b) : 
                                                         inf(0),
                                                         sup(0)
{ 
   inf+=a;
   sup+=b;
}



void l_real::_clear(int p) throw()
{  
   // filling a l_real number from element p up to the end with zero.
   for (int i=p; i<=prec; i++)
      this->elem(i)=0.0;
}


void l_real::_akku_out(const dotprecision& d) throw()
{   
    // The dotprecision value is rounded to the activated l_real number
    // in its own precision.
    dotprecision dot(d);
    _clear(1);
    this->elem(1) = rnd(dot);
    if (prec>1) 
    {
       int i=2, weiter;
       do {
           dot -= this->elem(i-1);
           this->elem(i) = rnd(dot);
           i++;
           weiter = this->elem(i-1) != 0;  // Blomquist 02.10.02.
       } while (weiter && (i <= prec));
    }
}

void l_real::_akku_out_up(const dotprecision& d) throw()
{   
    // The dotprecision value is rounded up to the activated l_real
    // number in its own precision.
    // Blomquist, 20.11.2006;
    bool weiter;
    dotprecision dot(d);
    _clear(1);
    this->elem(1) = (prec==1)? rnd(dot,RND_UP) : rnd(dot);
    if (prec>1) 
    {
       int i=2;
       do {
           dot -= this->elem(i-1);
	   weiter = (sign(dot) != 0);
	   if (weiter)
	       this->elem(i) = (i==prec)? 
		   rnd(dot,RND_UP) : rnd(dot); 
           i++;
       } while (weiter && (i <= prec));
    }
}

void l_real::_akku_out_down(const dotprecision& d) throw()
{   
    // The dotprecision value is rounded up to the activated l_real
    // number in its own precision.
    // Blomquist, 20.11.2006;
    bool weiter;
    dotprecision dot(d);
    _clear(1);
    this->elem(1) = (prec==1)? rnd(dot,RND_DOWN) : rnd(dot);
    if (prec>1) 
    {
       int i=2;
       do {
           dot -= this->elem(i-1);
	   weiter = (sign(dot) != 0);
	   if (weiter)
	       this->elem(i) = (i==prec)? 
		   rnd(dot,RND_DOWN) : rnd(dot); 
           i++;
       } while (weiter && (i <= prec));
    }
}

void l_real::_akku_add(dotprecision& d) const throw()
{ 
   // adding the activated l_real number to the accumulator d.
   for (int i=1; i<=prec; i++) 
   {
       if (this->elem(i) != 0) d += this->elem(i); 
   }
}

void l_real::_akku_sub(dotprecision& d) const throw()
{ 
   // subtracting the activated l_real number from the accumulator d.
   for (int i=1; i<=prec; i++) 
   {
      if (this->elem(i) != 0) d -= this->elem(i);
   }
}

//---------------------------------------------------------------------------
// Constructors and destructors
//
l_real::l_real() throw()
{
   data=new real[stagprec];      
   if(data)
      prec=stagprec;
   else
      prec=0;
}

l_real::~l_real() throw()
{
   delete [] data;
}

l_real::l_real(const l_real& lr) throw() 
                                 : data(new real[lr.prec])
{
   if(data)
   {
      prec=lr.prec;
      memcpy(data,lr.data,sizeof(real)*prec);
   } else
      prec=0;
}

l_real::l_real(const real & r) throw()
                             : prec(1), 
                               data(new real[1])
{
    data[0]=r;
}

l_real::l_real(const double & db) throw()
                             : prec(1), data(new real[1])
{  
   data[0] = db;   // Blomquist 10.09.02. Deklaration in l_real.hpp
}

l_real::l_real(int i) throw()
                    : prec(1),
                      data(new real[1])
{
   data[0]=i;
}

l_real::l_real(long i) throw()
                    : prec(1),
                      data(new real[1])
{
   data[0]=i;
}
//---------------------------------------------------------------------------
// Type transfers

real::real(const l_real& lr) throw()
{
   dotprecision dot(0.0);
   lr._akku_add(dot);
   w=rnd(dot).w;
}

dotprecision::dotprecision(const l_real& lr) throw() : akku(new a_btyp[A_LENGTH])
{
   memset(akku,0,BUFFERSIZE);
   lr._akku_add(*this);
}

dotprecision & dotprecision::operator =(const l_real & lr) throw()
{
   memset(akku,0,BUFFERSIZE);
   lr._akku_add(*this);
   return *this;
}

l_real::l_real(const dotprecision & d) throw() : prec(stagprec),
                                                 data(new real[prec])
{
   _akku_out(d);
}

//---------------------------------------------------------------------------
// assignments 
//
l_real& l_real::operator=(const real& r) throw()
{ 
   // Siehe _l_real(real..)

   if(prec!=1)
   {
      delete [] data;
      data=new real[1];
      if(data)
         prec=1;
      else
         prec=0;
   }
   elem(1)=r; 
   return *this;
}

l_real& l_real::operator=(const l_real& lr) throw()
{ 
   // Siehe c++-FAQ assignment-operators.html#[12.1]

   // This code gracefully (albeit implicitly) handles self assignment
   real* tmp = new real[lr.prec];
   // It would be OK if an exception got thrown here

   if((prec=lr.prec)>0)
      memcpy(tmp,lr.data,sizeof(real)*prec);

   delete [] data;
   data = tmp;
   return *this;
}

l_real& l_real::operator=(const dotprecision & d) throw()
{ 
   if(prec!=stagprec)
   {
      delete [] data;
      data=new real[prec=stagprec];
   }
   _akku_out(d);   
   return *this;
}


real& l_real::operator[](int i) const throw() 
{
        if (i<1 || i>prec) 
        {
                // throw!
                i=1;
        }
        return data[i-1];
}

int StagPrec(const l_real& lr) throw()
{
        return lr.prec;
}

std::istream& operator>>(std::istream& s, l_real& lr) throw()
{
        dotprecision dot;
        s >> dot;
        lr._akku_out(dot);
        return s;
}

std::ostream& operator<<(std::ostream& s, const l_real& lr) throw()
{
	dotprecision dot(0.0);
        lr._akku_add(dot);
        s << dot;
        return s;
}
std::string& operator>>(std::string& s, l_real& lr) throw()
{
	dotprecision dot;
        s >> dot;
        lr._akku_out(dot);
        return s;
}

std::string& operator<<(std::string& s, const l_real& lr) throw()
{
        dotprecision dot(0.0);
        lr._akku_add(dot);
        s << dot;
        return s;
}
void operator>>(const std::string& s, l_real& lr) throw()
{
   std::string r(s);
   r>>lr;
}
void operator>>(const char *s, l_real& lr) throw()
{
   std::string r(s);
   r>>lr;
}


void accumulate(dotprecision& d, const real& r, const l_real& lr) throw()
{
// accumulate(d, l_real(r), lr);       Blomquist: Old version, not from me!
    for (int i=1; i<=lr.prec; i++)  // Blomquist: My new version, 24.09.02
	accumulate(d, lr.elem(i), r);
}

void accumulate(dotprecision& d, const l_real& lr, const real& r) throw()
{
//  accumulate(d, lr, l_real(r));      Blomquist: Old version, not from me!
    for (int i=1; i<=lr.prec; i++)  // Blomquist: My new version, 24.09.02
	accumulate(d, lr.elem(i), r);
}

void accumulate(dotprecision& d, const l_real& lr1, const l_real&lr2) throw()
{
   int i, j;

   for (i=1; i<=lr1.prec; i++)
      for (j=1; j<=lr2.prec; j++)
//  if (abs(lr1.elem(i) * lr2.elem(j)) > MinReal) // alte Zeile von Toussaint
            accumulate(d, lr1.elem(i), lr2.elem(j));
}

void accumulate(idotprecision & a, const real & b, const l_real & c) throw() { accumulate(a,l_interval(b),l_interval(c)); }
void accumulate(idotprecision & a, const l_real & b, const real & c) throw() { accumulate(a,l_interval(b),l_interval(c)); }
void accumulate(idotprecision & a, const l_real & b, const l_real & c) throw() { accumulate(a,l_interval(b),l_interval(c)); }

l_real rnd_up(const dotprecision& a)
{ // Blomquist, 20.11.2006;
    l_real lr;
    lr._akku_out_up(a);
    return lr;
}

l_real rnd_down(const dotprecision& a)
{ // Blomquist, 20.11.2006;
    l_real lr;
    lr._akku_out_down(a);
    return lr;
}

l_real  operator-(const l_real& lr, const real& r) throw() { return lr-_l_real(r); }
l_real  operator-(const real& r, const l_real& lr) throw() { return _l_real(r)-lr; }
l_real  operator+(const l_real& lr, const real& r) throw() { return lr+_l_real(r); }
l_real  operator+(const real& r, const l_real& lr) throw() { return _l_real(r)+lr; }
l_real  operator*(const l_real& lr, const real& r) throw() { return lr*_l_real(r); }
l_real  operator*(const real& r, const l_real& lr) throw() { return _l_real(r)*lr; }
l_real  operator/(const l_real& lr, const real& r) throw() { return lr/_l_real(r); }
l_real  operator/(const real& r, const l_real& lr) throw() { return _l_real(r)/lr; }

dotprecision operator-(const l_real& lr, const dotprecision& r) throw() 
{ 
   return _dotprecision(lr)-r; 
}
dotprecision operator-(const dotprecision& r, const l_real& lr) throw() 
{ 
   return r-_dotprecision(lr); 
}
dotprecision operator+(const l_real& lr, const dotprecision& r) throw() 
{ 
   return _dotprecision(lr)+r; 
}
dotprecision operator+(const dotprecision& r, const l_real& lr) throw() 
{ 
   return r+_dotprecision(lr); 
}


l_real& operator-=(l_real& lr, const real& r) throw()
{ lr = lr-_l_real(r); return lr; }
    
l_real& operator+=(l_real& lr, const real& r) throw()
{ lr = lr+_l_real(r); return lr; }

l_real& operator*=(l_real& lr, const real& r) throw()
{ lr = lr*_l_real(r); return lr; }

l_real& operator/=(l_real& lr, const real& r) throw()
{ lr = lr/_l_real(r); return lr; }

real& operator-=(real& r, const l_real& lr) throw() { r = r-_real(lr); return r; }
real& operator+=(real& r, const l_real& lr) throw() { r = r+_real(lr); return r; }
real& operator*=(real& r, const l_real& lr) throw() { r = r*_real(lr); return r; }
real& operator/=(real& r, const l_real& lr) throw() { r = r/_real(lr); return r; }

bool operator==(const l_real& lr, const real& r) throw() { return lr==_l_real(r); }
bool operator==(const real& r, const l_real& lr) throw() { return _l_real(r)==lr; }
bool operator!=(const l_real& lr, const real& r) throw() { return lr!=_l_real(r); }
bool operator!=(const real& r, const l_real& lr) throw() { return _l_real(r)!=lr; }

bool operator<=(const l_real& lr, const real& r) throw() { return lr<=_l_real(r); }
bool operator<=(const real& r, const l_real& lr) throw() { return _l_real(r)<=lr; }
bool operator>=(const l_real& lr, const real& r) throw() { return lr>=_l_real(r); }
bool operator>=(const real& r, const l_real& lr) throw() { return _l_real(r)>=lr; }

bool operator<(const l_real& lr, const real& r) throw() { return lr<_l_real(r); }
bool operator<(const real& r, const l_real& lr) throw() { return _l_real(r)<lr; }
bool operator>(const l_real& lr, const real& r) throw() { return lr>_l_real(r); }
bool operator>(const real& r, const l_real& lr) throw() { return _l_real(r)>lr; }

bool operator==(const l_real& lr, const dotprecision& d) throw() { return _dotprecision(lr)==d; }
bool operator==(const dotprecision& d, const l_real& lr) throw() { return d==_dotprecision(lr); }
bool operator!=(const l_real& lr, const dotprecision& d) throw() { return _dotprecision(lr)!=d; }
bool operator!=(const dotprecision& d, const l_real& lr) throw() { return d!=_dotprecision(lr); }

bool operator<=(const l_real& lr, const dotprecision& d) throw() { return _dotprecision(lr)<=d; }
bool operator<=(const dotprecision& d, const l_real& lr) throw() { return d<=_dotprecision(lr); }
bool operator>=(const l_real& lr, const dotprecision& d) throw() { return _dotprecision(lr)>=d; }
bool operator>=(const dotprecision& d, const l_real& lr) throw() { return d>=_dotprecision(lr); }

bool operator<(const l_real& lr, const dotprecision& d) throw() { return _dotprecision(lr)<d; }
bool operator<(const dotprecision& d, const l_real& lr) throw() { return d<_dotprecision(lr); }
bool operator>(const l_real& lr, const dotprecision& d) throw() { return _dotprecision(lr)>d; }
bool operator>(const dotprecision& d, const l_real& lr) throw() { return d>_dotprecision(lr); }

l_real operator-(const l_real& lr1) throw()
{
   l_real lr2(lr1);
   for (int i=1; i<=lr1.prec; i++)
      lr2.elem(i) =  -(lr1.elem(i));
   return lr2;
}
l_real operator+(const l_real& lr1) throw()
{
   return lr1;
}

l_real operator-(const l_real& lr1, const l_real& lr2) throw()
{
   l_real lr3;
   dotprecision dot(0.0);
   lr1._akku_add(dot);
   lr2._akku_sub(dot);
   lr3._akku_out(dot);
   return lr3;
}

l_real operator+(const l_real& lr1, const l_real& lr2) throw()
{
   l_real lr3;
   dotprecision dot(0.0);
   lr1._akku_add(dot);
   lr2._akku_add(dot);
   lr3._akku_out(dot);

   return lr3;
}

l_real operator*(const l_real& lr1, const l_real& lr2) throw()
{
   l_real lr3;
   dotprecision dot(0.0);
   accumulate(dot, lr1, lr2);
   lr3._akku_out(dot);
   return lr3;
}

l_real operator/(const l_real& lr1, const l_real& lr2) throw(DIV_BY_ZERO)
{  // Blomquist, 09.12.02; throw() ---> throw(DIV_BY_ZERO)
   real a, b;
   l_real lr3;
   lr3._clear(1);
   dotprecision dot1(0.0);
   dotprecision dot2(0.0);
   lr1._akku_add(dot1);
   lr2._akku_add(dot2);
   a = rnd(dot1, RND_DOWN);
   b = rnd(dot2, RND_UP);
   if (!b) 
   { // Blomquist: cxscthrow(DIV_BY_ ... ) 
     cxscthrow(DIV_BY_ZERO("l_real operator/(const l_real&, const l_real&)"));
   } else 
   {
      lr3.elem(1) = a/b;
      for (int i=2; i<=stagprec; i++) 
      {
         if (!a) 
            break;
         for (int j=1; j<=lr2.prec; j++)
            if (!!lr3.elem(i-1) && !!lr2.elem(j) )
               accumulate(dot1, lr3.elem(i-1), -lr2.elem(j));
         a = rnd(dot1, RND_DOWN);
         lr3.elem(i) = a/b;
      }
   } // no division by zero
   return lr3;
}


l_real & operator-=(l_real& lr1, const l_real& lr2) throw() { return lr1 = lr1-lr2; }
l_real & operator+=(l_real& lr1, const l_real& lr2) throw() { return lr1 = lr1+lr2; }
l_real & operator*=(l_real& lr1, const l_real& lr2) throw() { return lr1 = lr1*lr2; }
l_real & operator/=(l_real& lr1, const l_real& lr2) throw() { return lr1 = lr1/lr2; }

bool operator==(const l_real& lr1, const l_real& lr2) throw()
{
   dotprecision dot1(0.0);
   dotprecision dot2(0.0);
   lr1._akku_add(dot1);
   lr2._akku_add(dot2);
   return (dot1==dot2);
}

bool operator<=(const l_real& lr1, const l_real& lr2) throw()
{
   dotprecision dot1(0.0);
   dotprecision dot2(0.0);
   lr1._akku_add(dot1);
   lr2._akku_add(dot2);
   return (dot1<=dot2);
}

bool operator!(const l_real& lr) throw()
{
   dotprecision dot(0.0);
   lr._akku_add(dot);
   return (!dot);
}

bool operator!=(const l_real& lr1, const l_real& lr2) throw() { return !(lr1==lr2); }
bool operator< (const l_real& lr1, const l_real& lr2) throw() { return !(lr2<=lr1); }
bool operator> (const l_real& lr1, const l_real& lr2) throw() { return !(lr1<=lr2); }
bool operator>=(const l_real& lr1, const l_real& lr2) throw() { return (lr2<=lr1);  }

// ID-LR

 idotprecision operator +(const idotprecision &a,const l_real &b) throw() 
{ 
   return idotprecision(a.inf+b,a.sup+b); 
}

 idotprecision operator +(const l_real &b,const idotprecision &a) throw() 
{ 
   return idotprecision(a.inf+b,a.sup+b); 
} 

 idotprecision operator -(const idotprecision &a,const l_real &b) throw() 
{ 
   return idotprecision(a.inf-b,a.sup-b); 
} 

 idotprecision operator -(const l_real &a,const idotprecision &b) throw() 
{ 
   return idotprecision(a-b.sup,a-b.inf); 
}

idotprecision operator |(const l_real &a,const idotprecision &b) throw() 
{
   return idotprecision((a<b.inf)?_dotprecision(a):b.inf,(a>b.sup)?_dotprecision(a):b.sup);
}

 idotprecision operator |(const idotprecision &a,const l_real &b) throw() 
{
   return idotprecision((a.inf<b)?a.inf:_dotprecision(b),(a.sup>b)?a.sup:_dotprecision(b));
}

 idotprecision operator &(const l_real &a,const idotprecision &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a>b.inf)?_dotprecision(a):b.inf,(a<b.sup)?_dotprecision(a):b.sup);
}

 idotprecision operator &(const idotprecision &a,const l_real &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) 
{
   return idotprecision((a.inf>b)?a.inf:_dotprecision(b),(a.sup<b)?a.sup:_dotprecision(b));
}

// LR-ID
idotprecision & operator +=(idotprecision &a,const l_real &b) throw() { return a+=dotprecision(b); }      
idotprecision & operator -=(idotprecision &a,const l_real &b) throw() { return a-=dotprecision(b); }
idotprecision & operator |=(idotprecision &a,const l_real &b) throw() { return a|=dotprecision(b); }
idotprecision & operator &=(idotprecision &a,const l_real &b) throw(ERROR_IDOTPRECISION_EMPTY_INTERVAL) { return a&=dotprecision(b); }


bool operator ==(const interval &i,const l_real &r) throw() { return Inf(i)==r && Sup(i)==r; }
bool operator !=(const interval &i,const l_real &r) throw() { return Inf(i)!=r || Sup(i)!=r; }
bool operator  <(const interval &i,const l_real &r) throw() { return false; }
bool operator  >(const interval &i,const l_real &r) throw() { return Inf(i)<r && Sup(i)>r; }
bool operator <=(const interval &i,const l_real &r) throw() { return Inf(i)==r && Sup(i)==r; }
bool operator >=(const interval &i,const l_real &r) throw() { return Inf(i)<=r && Sup(i)>=r; }

bool operator ==(const l_real &r,const interval &i) throw() { return Inf(i)==r && Sup(i)==r; }
bool operator !=(const l_real &r,const interval &i) throw()   { return Inf(i)!=r || Sup(i)!=r; }
bool operator  <(const l_real &r,const interval &i) throw()   { return Inf(i)<r && Sup(i)>r; }
bool operator  >(const l_real &r,const interval &i) throw()   { return false; }
bool operator <=(const l_real &r,const interval &i) throw()   { return Inf(i)<=r && Sup(i)>=r; }
bool operator >=(const l_real &r,const interval &i) throw()   { return Inf(i)==r && Sup(i)==r; }

bool operator ==(const idotprecision &i,const l_real &r) throw() { return Inf(i)==r && Sup(i)==r; }
bool operator !=(const idotprecision &i,const l_real &r) throw() { return Inf(i)!=r || Sup(i)!=r; }
bool operator  <(const idotprecision &i,const l_real &r) throw() { return false; }
bool operator  >(const idotprecision &i,const l_real &r) throw() { return Inf(i)<r && Sup(i)>r; }
bool operator <=(const idotprecision &i,const l_real &r) throw() { return Inf(i)==r && Sup(i)==r; }
bool operator >=(const idotprecision &i,const l_real &r) throw() { return Inf(i)<=r && Sup(i)>=r; }

bool operator ==(const l_real &r,const idotprecision &i) throw() { return Inf(i)==r && Sup(i)==r; }
bool operator !=(const l_real &r,const idotprecision &i) throw()   { return Inf(i)!=r || Sup(i)!=r; }
bool operator  <(const l_real &r,const idotprecision &i) throw()   { return Inf(i)<r && Sup(i)>r; }
bool operator  >(const l_real &r,const idotprecision &i) throw()   { return false; }
bool operator <=(const l_real &r,const idotprecision &i) throw()   { return Inf(i)<=r && Sup(i)>=r; }
bool operator >=(const l_real &r,const idotprecision &i) throw()   { return Inf(i)==r && Sup(i)==r; }

real & real::operator = (const l_real& a) throw()
{  // Blomquist, 12.11.2008;
	real x(a);
	return *this = x;
}

//---------------------------------------------------------------------------
// Other functions:
//
l_real  abs(const l_real& lr1) throw()
{
   l_real lr2;
   dotprecision dot(0.0);
   lr1._akku_add(dot);

   if (dot < 0.0) 
      dot = -dot;

   lr2._akku_out(dot);

   return lr2;
}

int sign(const l_real& lr) throw()
{
   dotprecision dot(0.0);
   lr._akku_add(dot);
   return sign(dot);
}

l_real adjust(const l_real & x) throw()
{
   l_real  y;
  
   if (x.prec == stagprec) 
      y = x;
   else if (x.prec > stagprec) 
   {
      dotprecision dot(0.0);
      x._akku_add(dot);
      y._akku_out(dot);
   } else 
   {
      int i;
      for (i = 0; i <= stagprec-x.prec-1; i++)
         y.data[i] = 0;
      for (i = stagprec-x.prec; i <= stagprec-1; i++)
         y.data[i] = x.data[i-(stagprec-x.prec)];
   }
        
   return y;
}

/*!
\param x The value for which to calculate
\return The result of the calculation

Result for a multiple-precisionnumber \f$ x = \sum \limits_{i=1}^n x_i \f$ .

\f[
\mbox{expo}_{ \mbox{sm} }(x) = \mbox{expo}( \mbox{ min } \{ | x_i | x_i \not= 0 \;,\; i = 1,...,n\})
\f]

*/
int expo_sm(const l_real& x)
// Calculating expo(x[k]) of the smallest |x[k]|<>0.
{
    int k(x.prec);
    l_real y(x);

    while (y.elem(k)==0 && k>1) k--;
    return expo(y.elem(k));
}

/*!
\param x The value for which to calculate
\return The result of the calculation

Result for a multiple-precisionnumber \f$ x = \sum \limits_{i=1}^n x_i \f$.
\f[
\mbox{expo}_{ \mbox{gr} } (x) = \mbox{expo}( \mbox{ max } \{ | x_i | x_i \not= 0 \;,\; i = 1,...,n\})
\f]
*/
int expo_gr(const l_real& x)
// Calculating expo(x[k]) of the greatest |x[k]|.
{
    int k(1),p(x.prec);
    l_real y(x);

    while (y.elem(k)==0 && k<p) k++;
    return expo(y.elem(k));
}


} // namespace cxsc

