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

/* CVS $Id: testsklr.hpp,v 1.25 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTSKLR_HPP_INCLUDED
#define _CXSC_TESTSKLR_HPP_INCLUDED

#include <fstream>

// ------------- Scalar - Tests ---------------------------

namespace cxsc {

template <class T,class S>
class scalarassign : public testclass
{
   public:
      scalarassign(void)
      {
         T x;
         S s;
         
         s=1;

         testing(nameof(x),"operator =("+nameof(x)+","+nameof(s)+")");
         if(!tested())
         {
            x=s;
            
            if(x==s && !(x!=s))
               ok();
            else
               error();
         }
      }
};

template <class T>
class allscalarassign : public testclass
{
   public:
      allscalarassign(void)
      {
         scalarassign<T,int>    assign_from_int;
         scalarassign<T,float>  assign_from_float;
         scalarassign<T,double> assign_from_double;
         scalarassign<T,real>   assign_from_real;

         if(!assign_from_int || !assign_from_float || !assign_from_double ||
            !assign_from_real)
            setfail();
      }
};

template <class A,class B>
class scalaraddsub : public testclass
{
   // Tests A+B, B+A, A-B, B-A, A+=B, A-=B
   //       -A, +A
   public:
      scalaraddsub(void)
      {
         A a,e;
         B b;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator +("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a= 123,b=1023,e=123+1023;
            if(a+b==e && !(a+b!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator +("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a= 123,b=1023,e=123+1023;
            if(b+a==e && !(b+a!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator +=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            a= 123,b=1023,e=123+1023;
            if((a+=b)==e && a==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a= 123,b=1023,e=123-1023;
            if(a-b==e && !(a-b!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a= 123,b=1023,e=1023-123;
            if(b-a==e && !(b-a!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            a= 123,b=1023,e=123-1023;
            if((a-=b)==e && a==e)
               ok();
            else
               error();
         }

         testing(nameof(e),"operator -("+nameof(a)+")");
         if(!tested())
         {
            a= 784,e= -784; 
            if((-a==e) && !(-a!=e))
               ok();
            else
               error();
         }

         testing(nameof(e),"operator +("+nameof(a)+")");
         if(!tested())
         {
            a= 784,e= 784; 
            if((+a==e) && !(+a!=e))
               ok();
            else
               error();
         }
      }
};

template <class A,class B>
class scalarmuldiv : public testclass
{
   // Tests A*B, B*A, A/B, B/A, A*=B, A/=B
   public:
      scalarmuldiv(void)
      {
         A a,e;
         B b;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator *("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a= 123,b=77,e=123*77L;
            if(a*b==e && !(a*b!=e) && (-a)*b==(-e) && a*(-b)==(-e) && (-a)*(-b)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator *("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a= 123,b=77,e=123*77L;
            if(b*a==e && !(b*a!=e) && (-b)*a==(-e) && b*(-a)==(-e) && (-b)*(-a)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator *=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            a= 123,b=77,e=123*77L;
            if((a*=b)==e && a==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=2048,b=512,e=4;
            if(a/b==e && !(a/b!=e) && (-a)/b==(-e) && a/(-b)==(-e) && (-a)/(-b)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a=512,b=2048,e=4;
            if(b/a==e && !(b/a!=e) && (-b)/a==(-e) && b/(-a)==(-e) && (-b)/(-a)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            a=2048,b=512,e=4;
            if((a/=b)==e && a==e)
               ok();
            else
               error();
         }
      }
};

template <class A,class B>
class scalarcompare : public testclass
{
   // Tests A<B, B<A, A>B, B>A, A<=B, A>=B, A==B, B==A, A!=B, B!=A
   //       !A, (A)
   public:
      scalarcompare(void)
      {
         A a1,a2;
         B b1,b2;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         a1=123,b1=123;
         a2=424,b2=424;
         
         testing(nameof(a1),"operator <("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( !(a1<b1) && !(a2<b2) && !(a2<b1) && (a1<b2) )
               ok();
            else
               error();
         }
         
         testing(nameof(b1),"operator <("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( !(b1<a1) && !(b2<a2) && !(b2<a1) && (b1<a2) )
               ok();
            else
               error();
         }
         
         testing(nameof(a1),"operator <=("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( (a1<=b1) && (a2<=b2) && !(a2<=b1) && (a1<=b2) )
               ok();
            else
               error();
         }
         
         testing(nameof(b1),"operator <=("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( (b1<=a1) && (b2<=a2) && !(b2<=a1) && (b1<=a2) )
               ok();
            else
               error();
         }
         
         testing(nameof(a1),"operator >=("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( (a1>=b1) && (a2>=b2) && (a2>=b1) && !(a1>=b2) )
               ok();
            else
               error();
         }
         
         testing(nameof(b1),"operator >=("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( (b1>=a1) && (b2>=a2) && (b2>=a1) && !(b1>=a2) )
               ok();
            else
               error();
         }
         
         testing(nameof(a1),"operator >("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( !(a1>b1) && !(a2>b2) && (a2>b1) && !(a1>b2) )
               ok();
            else
               error();
         }
         
         testing(nameof(b1),"operator >("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( !(b1>a1) && !(b2>a2) && (b2>a1) && !(b1>a2) )
               ok();
            else
               error();
         }
         
         testing(nameof(a1),"operator ==("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( (a1==b1) && (a2==b2) && !(a2==b1) && !(a1==b2) )
               ok();
            else
               error();
         }

         testing(nameof(a1),"operator !=("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( !(a1!=b1) && !(a2!=b2) && (a2!=b1) && (a1!=b2) )
               ok();
            else
               error();
         }
         
         testing(nameof(b1),"operator ==("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( (b1==a1) && (b2==a2) && !(b2==a1) && !(b1==a2) )
               ok();
            else
               error();
         }

         testing(nameof(b1),"operator !=("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( !(b1!=a1) && !(a2!=a2) && (b2!=a1) && (b1!=a2) )
               ok();
            else
               error();
         }
         
         a1=123,a2=0;
         
         testing(nameof(a1),"operator !("+nameof(a1)+")");
         if(!tested())
         {
            if( !!a1 && !a2 )
               ok();
            else
               error();
         }
         
/*         testing(nameof(a1),"operator (void *) ("+nameof(a1)+")");
         if(!tested())
         {
            if(a1)
               if(a2)
                  error();
               else
                  ok();
            else
               error();
         }*/
         
      }
};

template <class A,class B>
class scalar_eq_compare : public testclass
{
   // Tests A==B, B==A, A!=B, B!=A
   //       !A, (A)
   public:
      scalar_eq_compare(void)
      {
         A a1,a2;
         B b1,b2;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         a1=123,b1=123;
         a2=424,b2=424;
         
         testing(nameof(a1),"operator ==("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( (a1==b1) && (a2==b2) && !(a2==b1) && !(a1==b2) )
               ok();
            else
               error();
         }

         testing(nameof(a1),"operator !=("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( !(a1!=b1) && !(a2!=b2) && (a2!=b1) && (a1!=b2) )
               ok();
            else
               error();
         }
         
         testing(nameof(b1),"operator ==("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( (b1==a1) && (b2==a2) && !(b2==a1) && !(b1==a2) )
               ok();
            else
               error();
         }

         testing(nameof(b1),"operator !=("+nameof(b1)+","+nameof(a1)+")");
         if(!tested())
         {
            if( !(b1!=a1) && !(a2!=a2) && (b2!=a1) && (b1!=a2) )
               ok();
            else
               error();
         }
         
         a1=123,a2=0;
         
         testing(nameof(a1),"operator !("+nameof(a1)+")");
         if(!tested())
         {
            if( !!a1 && !a2 )
               ok();
            else
               error();
         }
         
/*         testing(nameof(a1),"operator (void *) ("+nameof(a1)+")");
         if(!tested())
         {
            if(a1)
               if(a2)
                  error();
               else
                  ok();
            else
               error();
         }*/
         
      }
};

template <class T>
class scalarabssgn : public testclass
{
   public:
      scalarabssgn(void)
      {
         T x1,x2,x3,y1,y2,y3;
         
         testing(nameof(x1),"abs("+nameof(x1)+")");
         if(!tested())
         {
            x1=3,y1=3;
            x2=-2,y2=2;
            x3=0,y3=0;

            if(abs(x1)==y1 && abs(x2)==y2 && abs(x3)==y3)
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sign("+nameof(x1)+")");
         if(!tested())
         {
            x1=3,y1=1;
            x2=-2,y2=-1;
            x3=0,y3=0;

            if(sign(x1)==y1 && sign(x2)==y2 && sign(x3)==y3)
               ok();
            else
               error();
         }
      }
};

template <class T>
class scalarstdfunc : public testclass
{
   public:
      scalarstdfunc(void)
      {
         T x1;
         
         T pi,pi_2,pi_4,e,e1d,ee;
         CSTRING_PI>>pi;
         CSTRING_PI_2>>pi_2;
         CSTRING_PI_4>>pi_4;
         CSTRING_E>>e;
         e1d=T(1.)/e;
         ee=e*e;  
         
         testing(nameof(x1),"sqr("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(11),sqr(T(11)),T(121))
            && compare(T(0),sqr(T(0)),T(0))
            && compare(T(-9),sqr(T(-9)),T(81)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sqrt("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(121),sqrt(T(121)),T(11))
            && compare(T(0),sqrt(T(0)),T(0))
            && compare(T(81),sqrt(T(81)),T(9)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sqrt("+nameof(x1)+",int)");
         if(!tested())
         {
            if(compare(T(27),T(3),sqrt(T(27),int(3)),T(3))
            && compare(T(0),T(4),sqrt(T(0),int(4)),T(0))
            && compare(T(1024),T(10),sqrt(T(1024),int(10)),T(2)) )
               ok();
            else
               error();
         }
         
/*         testing(nameof(x1),"sqrtm1("+nameof(x1)+") funktionioniert nicht?");
         if(!tested())
         {
            x1=121,e1=10;
            x2=4,  e2=1;
            x3=81, e3=8;
            
            if(sqrtm1(x1)==e1 && sqrtm1(x2)==e2 && sqrtm1(x3)==e3)
               ok();
            else
               error();
         } Spezielle real-Funktion */
         
         testing(nameof(x1),"sin("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),sin(T(0)),T(0))
            && compare(T(pi_2),sin(T(pi_2)),T(1))
            && compare(T(-pi),sin(T(-pi)),T(0)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"cos("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),cos(T(0)),T(1))
            && compare(T(pi_2),cos(T(pi_2)),T(0))
            && compare(T(-pi),cos(T(-pi)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"tan("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),tan(T(0)),T(0))
            && compare(T(pi_4),tan(T(pi_4)),T(1))
            && compare(T(-pi_4),tan(T(-pi_4)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"cot("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(pi_2),cot(T(pi_2)),T(0))
            && compare(T(pi_4),cot(T(pi_4)),T(1))
            && compare(T(-pi_4),cot(T(-pi_4)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"asin("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),asin(T(0)),T(0))
            && compare(T(1),asin(T(1)),T(pi_2))
            && compare(T(-1),asin(T(-1)),T(-pi_2)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acos("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),acos(T(0)),T(pi_2))
            && compare(T(1),acos(T(1)),T(0))
            && compare(T(-1),acos(T(-1)),T(pi)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"atan("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),atan(T(0)),T(0))
            && compare(T(1),atan(T(1)),T(pi_4))
            && compare(T(-1),atan(T(-1)),T(-pi_4)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acot("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),acot(T(0)),T(pi_2))
            && compare(T(1),acot(T(1)),T(pi_4))
            && compare(T(-1),acot(T(-1)),T(pi_4+pi_2)) )
               ok();
            else
               error();
         }
         /*testing(nameof(x1),"expm1("+nameof(x1)+")");
         if(!tested())
         {
            x1=0,e1=0;
            x2=1,e2=eM1;
            x3=-1, e3=e1d-1;
            if(expm1(x1)==e1 && expm1(x2)==e2 && expm1(x3)==e3)
               ok();
            else
               error();
         }
         testing(nameof(x1),"lnp1("+nameof(x1)+")");
         if(!tested())
         {
            x1=e,e1=2;
            x2=1,e2=1;
            x3=ee, e3=3;
            if(lnp1(x1)==e1 && lnp1(x2)==e2 && lnp1(x3)==e3)
               ok();
            else
               error();
         }beides spezielle real-Funktionen */
         
         testing(nameof(x1),"exp("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),exp(T(0)),T(1))
            && compare(T(1),exp(T(1)),T(e))
            && compare(T(-1),exp(T(-1)),T(e1d)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"ln("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(e),ln(T(e)),T(1))
            && compare(T(1),ln(T(1)),T(0))
            && compare(T(ee),ln(T(ee)),T(2)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sinh("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),sinh(T(0)),T(0))
            && compare(T(1),sinh(T(1)),T( (e-1/e)/2 ))
            && compare(T(-1),sinh(T(-1)),T( (1/e-e)/2 )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"cosh("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),cosh(T(0)),T(1))
            && compare(T(1),cosh(T(1)),T( (e+1/e)/2 ))
            && compare(T(-1),cosh(T(-1)),T( (e+1/e)/2 )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"tanh("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),tanh(T(0)),T(0))
            && compare(T(1),tanh(T(1)),T( (e-1/e)/(e+1/e) ))
            && compare(T(-1),tanh(T(-1)),T( (1/e-e)/(e+1/e) )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"coth("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(2),coth(T(2)),T( cosh(T(2))/sinh(T(2)) )) // =:}
            && compare(T(1),coth(T(1)),T( (e+1/e)/(e-1/e) ))
            && compare(T(-1),coth(T(-1)),T( (e+1/e)/(1/e-e) )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"asinh("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),asinh(T(0)),T(0))
            && compare(T((e-1/e)/2),asinh(T((e-1/e)/2)),T(1))
            && compare(T((1/e-e)/2),asinh(T((1/e-e)/2)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acosh("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(1),acosh(T(1)),T(0))
            && compare(T((e+1/e)/2),acosh(T((e+1/e)/2)),T(1))
            && compare(T(cosh(T(2))),acosh(T(cosh(T(2)))),T(2)) ) // =:}
               ok();
            else
               error();
         }
         testing(nameof(x1),"atanh("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(0),atanh(T(0)),T(0))
            && compare(T((e-1/e)/(e+1/e)),atanh(T((e-1/e)/(e+1/e))),T(1))
            && compare(T((1/e-e)/(e+1/e)),atanh(T((1/e-e)/(e+1/e))),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acoth("+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T((e+1/e)/(e-1/e)),acoth(T((e+1/e)/(e-1/e))),T(1))
            && compare(T((e+1/e)/(1/e-e)),acoth(T((e+1/e)/(1/e-e))),T(-1))
            && compare(T(coth(T(2))),acoth(T(coth(T(2)))),T(2)) )
               ok();
            else
               error();
         }

         
         testing(nameof(x1),"pow("+nameof(x1)+","+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(2),T(2),pow(T(2),T(2)),T(4))
            && compare(T(4),T(5),pow(T(4),T(5)),T(1024))
            && compare(T(+2),T(3),pow(T(+2),T(3)),T(+8)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"power("+nameof(x1)+",int)");
         if(!tested())
         {
            if(compare(T(2),T(2),power(T(2),int(2)),T(4))
            && compare(T(4),T(5),power(T(4),int(5)),T(1024))
            && compare(T(+2),T(3),power(T(+2),int(3)),T(+8)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"pow("+nameof(x1)+","+nameof(x1)+")");
         if(!tested())
         {
            if(compare(T(2),T(2),pow(T(2),T(2)),T(4))
            && compare(T(4),T(5),pow(T(4),T(5)),T(1024))
            && compare(T(+2),T(3),pow(T(+2),T(3)),T(+8)) )
               ok();
            else
               error();
         }
                                                   
         
      }
};


template <class T>
class allscalaraddsub : public testclass
{
   public:
      allscalaraddsub(void)
      {
         scalaraddsub<T,int>    add_int;
         scalaraddsub<T,long>   add_long;
         scalaraddsub<T,float>  add_float;
         scalaraddsub<T,double> add_double;
         scalaraddsub<T,real>   add_real;
         
         if(!add_int || !add_float || !add_double ||
            !add_real)
            setfail();
      }
};

template <class T>
class allscalarmuldiv : public testclass
{
   public:
      allscalarmuldiv(void)
      {
         scalarmuldiv<T,int>    mul_int;
         scalarmuldiv<T,long>   mul_long;
         scalarmuldiv<T,float>  mul_float;
         scalarmuldiv<T,double> mul_double;
         scalarmuldiv<T,real>   mul_real;
         
         if(!mul_int || !mul_long || !mul_float || !mul_double ||
            !mul_real)
            setfail();
      }
};

template <class T>
class allscalarcompare : public testclass
{
   public:
      allscalarcompare(void)
      {
         scalarcompare<T,int>    comp_int;
         scalarcompare<T,long>   comp_long;
         scalarcompare<T,float>  comp_float;
         scalarcompare<T,double> comp_double;
         scalarcompare<T,real>   comp_real;
         scalarcompare<T,dotprecision> comp_dotprecision;
         scalarcompare<T,l_real> comp_l_real; 
         
         if(!comp_int || !comp_long || !comp_float || !comp_double ||
            !comp_real || !comp_dotprecision || !comp_l_real)
            setfail();
      }
};

template <class T>
class allscalar_eq_compare : public testclass
{
   public: 
      allscalar_eq_compare(void)
      {
         scalar_eq_compare<T,int>    comp_int;
         scalar_eq_compare<T,long>   comp_long;
         scalar_eq_compare<T,float>  comp_float;
         scalar_eq_compare<T,double> comp_double;
         scalar_eq_compare<T,real>   comp_real;
         scalar_eq_compare<T,dotprecision> comp_dotprecision;
         
         if(!comp_int || !comp_long || !comp_float || !comp_double ||
            !comp_real || !comp_dotprecision)
            setfail();
      }
};

template <class D>
class cast_to_real : public testclass
{
   public:
      cast_to_real(void)
      {
         real a;
         D    d;

         testing(nameof(a),"real("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=real(d);
            if(a!=12 || real(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_real("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_real(d);
            if(a!=12 || _real(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=d;
            if(a!=12)
               error();
            else 
               ok();
         }
      }
};

template <class D>
class cast_to_l_real : public testclass
{
   public:
      cast_to_l_real(void)
      {
         l_real a;
         D    d;

         testing(nameof(a),"l_real("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=l_real(d);
            if(a!=12 || l_real(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_l_real("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_l_real(d);
            if(a!=12 || _l_real(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=d;
            if(a!=12)
               error();
            else 
               ok();
         }
      }
};

template <class D>
class cast_to_double : public testclass
{
   public:
      cast_to_double(void)
      {
         double a;
         D    d;

         testing(nameof(a),"_double("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_double(d);
            if(a!=12 || _double(d)!=12)
               error();
            else 
               ok();
         }
         
      }
};

class test_realio : public testclass
{
   public:
      test_realio(void)
      {
         real a;
         
         testing(nameof(a),"operator >>(string,"+nameof(a)+")");
         if(!tested())
         {
            int nok=0;
            a=0;
            
            string("1234.5678") >> a;
            if(a==1234.5678)
               nok++;
            
            string("-0.01234") >> a;
            if(a==-0.01234)
               nok++;
               
            string("2.34e+56") >> a;
            if(a==2.34e+56)
               nok++;
               
            string("-2.34e-5") >> a;
            if(a==-2.34e-5)
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator >>(char *,"+nameof(a)+")");
         if(!tested())
         {
            int nok=0;
            a=0;
            
            "1234.5678" >> a;
            
            if(a==1234.5678)
               nok++;
            
            (("-0.01234") >> a);
            
            if(a==-0.01234)
               nok++;
               
            (("2.34e+56") >> a);
            
            if(a==2.34e+56)
               nok++;
               
            (("-2.34e-5") >> a);
            if(a==-2.34e-5)
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator >>(std::istream,"+nameof(a)+")");
         if(!tested())
         {
            int nok=0;
            {
               std::ofstream out("test.tmp");
               out << "1234.5678" << std::endl 
                   << "-0.01234" << std::endl
                   << "2.34e+56" << std::endl
                   << "-2.34e-5" << std::endl;
            }
            std::ifstream in("test.tmp");
            
            a=0;
            
            in >> a;
            if(a==1234.5678)
               nok++;
            
            in >> a;
            if(a==-0.01234)
               nok++;
               
            in >> a;
            if(a==2.34e+56)
               nok++;
               
            in >> a;
            if(a==-2.34e-5)
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();
         }         
         
         testing(nameof(a),"operator <<(string,"+nameof(a)+")");
         if(!tested())
         {
            int nok=0;
            string s;
            
            cout << SetPrecision(19,11);
            s=">";
            
            CSTRING_PI >> a;
            s << a;

            if(s==">      3.14159265359")
               nok++;
            else 
               cout << ">"<<s<<"<"<< std::endl;
            
            a=-1024;
            s="<";
            s << a;
            if(s=="<  -1024.00000000000")
               nok++;
            else
               cout << ">" << s << "<" << std::endl;
               
            a=-123.456e+78;
            s="";
            s << a;
            if(s=="-1.23456000000E+080")
               nok++;
            else
               cout << ">" << s << "<" << std::endl;
            
            a=123.456e-78;
            s="";
            s<<a;
            if(s==" 1.23456000000E-076")
               nok++;
            else
               cout << ">" << s << "<" << std::endl;
            
            if(nok!=4)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator <<(std::ostream,"+nameof(a)+")");
         if(!tested())
         {
            int nok=0;
            string s;
            std::ofstream out("test.tmp");
            
            out << SetPrecision(19,11);
            
            CSTRING_PI >> a;
            out << a << std::endl;
            
            a=-1024;
            out << a << std::endl;
            
            a=-123.456e+78;
            out << a << std::endl;
            
            a=123.456e-78;
            out << a << std::endl;
            
            out.close();
            std::ifstream in("test.tmp");
            
            in >> s;
            
            if(s=="3.14159265359")
               nok++;
            
            in >> s;
            
            if(s=="-1024.00000000000")
               nok++;
               
            in >> s;
               
            if(s=="-1.23456000000E+080")
               nok++;
            
            in >> s;
            
            if(s=="1.23456000000E-076")
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();
         }
         
      }
};

template <class D>
class testpredsucc : public testclass
{
   public:
      testpredsucc(void)
      {
         D    d;

         testing(nameof(d),"succ("+nameof(d)+")");
         if(!tested())
         {
            int nok(0);
            D e;
            
            d=1e+20;
            e=(succ(d)-d);
            if(e>D(0))
               nok++;
            
            d=-1e-20;
            e=(succ(d)-d);
            if(e>D(0))
               nok++;
            
            
            d=-1;
            e=(succ(d)-d);
            if(e>D(0))
               nok++;
            
            
            d=0;
            e=(succ(d)-d);
            if(e>D(0))
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();  
         }

         testing(nameof(d),"pred("+nameof(d)+")");
         if(!tested())
         {
            int nok(0);
            D e;
            
            d=1e+20;
            e=(d-pred(d));
            
            if(e>D(0))
               nok++;
            
            d=-1e-20;
            e=(d-pred(d));
            
            if(e>D(0))
               nok++;
            
            
            d=-1;
            e=(d-pred(d));

            if(e>D(0))
               nok++;
            
            
            d=0;
            e=(d-pred(d));
            
            if(e>D(0))
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();  
         }
         
      }
};

template <class D>
class testaddsubmultdivupdown : public testclass
{
   public:
      testaddsubmultdivupdown(void)
      {
         D    d;

         testing(nameof(d),"addup("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=addup(a,b);
            
            if(e!=a+b && e!=succ(a+b))
               error();
            else 
               ok();  
         }
         testing(nameof(d),"adddown("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=adddown(a,b);
            
            if(e!=a+b && e!=pred(a+b))
               error();
            else 
               ok();  
         }
         testing(nameof(d),"subup("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(-10);
            e=subup(a,b);
            
            if(e!=(a-b) && e!=succ(a-b))
               error();
            else 
               ok();  
         }
         testing(nameof(d),"subdown("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=subdown(a,b);
            
            if(e!=a-b && e!=pred(a-b))
               error();
            else 
               ok();  
         }
         
         testing(nameof(d),"multup("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=multup(a,b);
            
            if(e!=a*b && e!=succ(a*b))
               error();
            else 
               ok();  
         }
         testing(nameof(d),"multdown("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=multdown(a,b);
            
            if(e!=a*b && e!=pred(a*b))
               error();
            else 
               ok();  
         }
         testing(nameof(d),"divup("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=divup(a,b);
            
            if(e!=(a/b) && e!=succ(a/b))
               error();
            else 
               ok();  
         }
         testing(nameof(d),"divdown("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            D e;
            D a,b;
            
            a=D(1e+30);
            b=D(10);
            e=divdown(a,b);
            
            if(e!=a/b && e!=pred(a/b))
               error();
            else 
               ok();  
         }
         
      }
};

template <class D>
class testminmax : public testclass
{
   public:
      testminmax(void)
      {
         D    a1,a2,a3,d;
         
         a1=-100;
         a2=0;
         a3=123;
         

         testing(nameof(d),"min("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            if(min(a1,a2)!=a1 || min(a1,a3)!=a1 || min(a2,a3)!=a2)
               error();
            else 
               ok();
         }
         testing(nameof(d),"max("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            if(max(a1,a2)!=a2 || max(a1,a3)!=a3 || max(a2,a3)!=a3)
               error();
            else 
               ok();
         }
         
      }
};

template <class D>
class testnullconstructor : public testclass
{
   public:
      testnullconstructor(void)
      {
         D d;
         testing(nameof(d),nameof(d)+"()");
         if(!tested())
         {
            ok();
         }
      }
};

} // namespace cxsc 

#endif // _CXSC_TESTSKLR_HPP_INCLUDED
