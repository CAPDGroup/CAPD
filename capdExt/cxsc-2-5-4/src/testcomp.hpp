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

/* CVS $Id: testcomp.hpp,v 1.25 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTCOMP_HPP_INCLUDED
#define _CXSC_TESTCOMP_HPP_INCLUDED

#include <fstream>

namespace cxsc {

// ------------- Complex - Scalar - Tests ---------------------------

template <class A,class B>
class complexaddsub : public testclass
{
   // Tests A+B, B+A, A-B, B-A, A+=B, A-=B
   //       -A, +A
   public:
      complexaddsub(void)
      {
         A a,e;
         B b;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator +("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=complex(123,12),b=complex(423,11),e=complex(546,23);
            if(a+b==e && !(a+b!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator +("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a=complex(123,12),b=complex(423,11),e=complex(546,23);
            if(b+a==e && !(b+a!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator +=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=complex(123,12),b=complex(423,11),e=complex(546,23);
            if((a+=b)==e && a==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=complex(123,12),b=complex(423,11),e=complex(-300,1);
            if(a-b==e && !(a-b!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a=complex(123,12),b=complex(423,11),e=complex(300,-1);

            if(b-a==e && !(b-a!=e))
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=complex(123,12),b=complex(423,11),e=complex(-300,1);
            if((a-=b)==e && a==e)
               ok();
            else
               error();
         }

         testing(nameof(e),"operator -("+nameof(a)+")");
         if(!tested())
         {
            a= complex(784,-17),e= complex(-784,17); 
            if((-a==e) && !(-a!=e))
               ok();
            else
               error();
         }

         testing(nameof(e),"operator +("+nameof(a)+")");
         if(!tested())
         {
            a= complex(784,34),e= complex(784,34); 
            if((+a==e) && !(+a!=e))
               ok();
            else
               error();
         }
      }
};

template <class A,class B>
class complexmuldiv : public testclass
{
   // Tests A*B, B*A, A/B, B/A, A*=B, A/=B
   public:
      complexmuldiv(void)
      {
         A a,e;
         B b;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator *("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a= complex(2,3),b=complex(4,5),e=complex(-7,22);
            if(a*b==e && !(a*b!=e) && (-a)*b==(-e) && a*(-b)==(-e) && (-a)*(-b)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator *("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a= complex(2,3),b=complex(4,5),e=complex(-7,22);
            if(b*a==e && !(b*a!=e) && (-b)*a==(-e) && b*(-a)==(-e) && (-b)*(-a)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator *=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a= complex(2,3),b=complex(4,5),e=complex(-7,22);
            if((a*=b)==e && a==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=complex(-7,22),b=complex(2,3),e=complex(4,5);
            if(a/b==e && !(a/b!=e) && (-a)/b==(-e) && a/(-b)==(-e) && (-a)/(-b)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            b=complex(-7,22),a=complex(4,5),e=complex(2,3);
            if(b/a==e && !(b/a!=e) && (-b)/a==(-e) && b/(-a)==(-e) && (-b)/(-a)==e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=complex(-7,22),b=complex(2,3),e=complex(4,5);
            if((a/=b)==e && a==e)
               ok();
            else
               error();
         }
      }
};

template <class A,class B>
class complexcompare : public testclass
{
   //       !A, (A)
   public:
      complexcompare(void)
      {
         A a1,a2,a3,a4;
         B b1,b2,b3,b4;

         a1=complex(123,17),b1=complex(123,17);
         a2=complex(424,23),b2=complex(424,23);
         a3=complex(123,23),b3=complex(123,23);
         a4=complex(424,17),b4=complex(424,17);
         
         testing(nameof(a1),"operator ==("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if(  (a1==b1) &&  (a2==b2) && !(a2==b1) && !(a1==b2) 
             &&  (a3==b3) &&  (a4==b4) && !(a4==b3) && !(a3==b4)
             && !(a1==b3) && !(a1==b4) && !(a2==b3) && !(a2==b4)
             && !(a3==b1) && !(a3==b2) && !(a4==b1) && !(a4==b2) )
               ok();
            else
               error();
         }

         testing(nameof(a1),"operator !=("+nameof(a1)+","+nameof(b1)+")");
         if(!tested())
         {
            if( !(a1!=b1) && !(a2!=b2) &&  (a2!=b1) &&  (a1!=b2) 
             && !(a3!=b3) && !(a4!=b4) &&  (a4!=b3) &&  (a3!=b4)
             &&  (a1!=b3) &&  (a1!=b4) &&  (a2!=b3) &&  (a2!=b4)
             &&  (a3!=b1) &&  (a3!=b2) &&  (a4!=b1) &&  (a4!=b2) )
               ok();
            else
               error();
         }
         
         a1=complex(0,0),a2=complex(123,0),a3=complex(0,232),a4=complex(12,42);
         
         testing(nameof(a1),"operator !("+nameof(a1)+")");
         if(!tested())
         {
            if( !a1 && !!a2 && !!a3 && !!a4)
               ok();
            else
               error();
         }
         
/*         testing(nameof(a1),"operator (void *) ("+nameof(a1)+")");
         if(!tested())
         {
            if(a1)
               error();
            else if(a2)
               if(a3)
                  if(a4)
                     ok();
                  else 
                     error();
               else 
                  error();
            else
               error();            
         }*/
         
      }
};

template <class T,class E>
class complexabs : public testclass
{
   public:
      complexabs(void)
      {
         T x1,x2,x3;
         E y1,y2,y3;
         
         testing(nameof(x1),"abs("+nameof(x1)+")");
         if(!tested())
         {
            x1=complex(3,4),y1=5;
            x2=complex(-4,-3),y2=5;
            x3=complex(0,0),y3=0;

            if(abs(x1)==y1 && abs(x2)==y2 && abs(x3)==y3)
               ok();
            else
               error();
         }
      }
};


/*
template <class D,class A,class B>
class scalaraccumulate : public testclass
{
   // A sollte der genauere Typ sein
   public:
      scalar_accumulate(void)
      {
         D d;
         A a;
         B b;

         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=434;
            b=134;
            d=0;
            accumulate(d,a,b);
            if(d!=434*134 || d!=a*b)
               error();
            else 
               ok();
         }
         
         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a=432;
            b=135;
            d=0;
            accumulate(d,b,a);
            if(d!=432*135 || d!=a*b)
               error();
            else 
               ok();
         }
         
      }
};
*/
template <class D>
class cast_c_to_complex : public testclass
{
   public:
      cast_c_to_complex(void)
      {
         complex a;
         D       d;

         testing(nameof(a),"complex("+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=complex(12,343);
            a=complex(d);
            if(a!=d)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_complex("+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=complex(12,343);
            a=_complex(d);
            if(a!=d)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=complex(12,343);
            a=d;
            if(a!=d)
               error();
            else 
               ok();         
         }
      }
};

template <class D>
class cast_r_to_complex : public testclass
{
   public:
      cast_r_to_complex(void)
      {
         complex a;
         D       d,e;

         testing(nameof(a),"complex("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=D(12),e=D(343);
            a=complex(d,e);
            if(Re(a)!=d || Im(a)!=e)
               error();
            else 
               ok();
         }

         testing(nameof(a),"_complex("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=D(12),e=D(343);
            a=_complex(d,e);
            if(Re(a)!=d || Im(a)!=e)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"=("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=D(12);
            a=d;
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=D(12);
            a=complex(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         testing(nameof(a),"_"+nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=complex(434,232);
            d=D(12);
            a=_complex(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }         
      }
};

template <class D>
class cast_r_to_cdotprecision : public testclass
{
   public:
      cast_r_to_cdotprecision(void)
      {
         cdotprecision a;
         D       d,e;

         testing(nameof(a),"cdotprecision("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12),e=D(343);
            a=cdotprecision(d,e);
            if(Re(a)!=d || Im(a)!=e)
               error();
            else 
               ok();
         }

         testing(nameof(a),"_cdotprecision("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12),e=D(343);
            a=cdotprecision(d,e);
            if(Re(a)!=d || Im(a)!=e)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"=("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12);
            a=d;
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12);
            a=cdotprecision(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         testing(nameof(a),"_"+nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12);
            a=_cdotprecision(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }         
      }
};

template <class D>
class cast_c_to_cdotprecision : public testclass
{
   public:
      cast_c_to_cdotprecision(void)
      {
         cdotprecision a;
         D       d;

         testing(nameof(a),"cdotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12,14);
            a=cdotprecision(d);
            if(Re(a)!=12 || Im(a)!=14)
               error();
            else 
               ok();
         }

         testing(nameof(a),"_cdotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12,343);
            a=_cdotprecision(d);
            if(Re(a)!=12 || Im(a)!=343)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"=("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=cdotprecision(434,232);
            d=D(12,14);
            a=d;
            if(Re(a)!=12 || Im(a)!=14)
               error();
            else 
               ok();
         }
      }
};

template <class D,class G>
class cast_scalar_to_complex : public testclass
{
   public:
      cast_scalar_to_complex(void)
      {
         complex a;
         D    d,e;
         
         testing(nameof(a),"complex("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=complex(d);
            if(a!=12 || complex(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_complex("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_complex(d);
            if(a!=12 || _complex(d)!=12)
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
         
         testing(nameof(a),"complex("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=complex(d,e);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_complex("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=_complex(d,e);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         G g;
         testing(nameof(a),"complex("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=complex(d,g);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_complex("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=_complex(d,g);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         testing(nameof(a),"complex("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=complex(g,d);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_complex("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=_complex(g,d);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
      }
};

template <class D,class G>
class cast_scalar_to_cdotprecision : public testclass
{
   public:
      cast_scalar_to_cdotprecision(void)
      {
         cdotprecision a;
         D    d,e;
         
         testing(nameof(a),"cdotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=cdotprecision(d);
            if(a!=12 || cdotprecision(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_cdotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_cdotprecision(d);
            if(a!=12 || _cdotprecision(d)!=12)
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
         
         testing(nameof(a),"cdotprecision("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=cdotprecision(d,e);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_cdotprecision("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=_cdotprecision(d,e);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         G g;
         testing(nameof(a),"cdotprecision("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=cdotprecision(d,g);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_cdotprecision("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=_cdotprecision(d,g);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         testing(nameof(a),"cdotprecision("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=cdotprecision(g,d);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_cdotprecision("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=_cdotprecision(g,d);
            if(Re(a)!=12 || Im(a)!=19)
               error();
            else 
               ok();
         }
         
      }
};



template <class D>
class cast_c_to_cinterval : public testclass
{
   public:
      cast_c_to_cinterval(void)
      {
         cinterval a;
         D       d;

         testing(nameof(a),"cinterval("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=complex(12,343);
            a=cinterval(d);
            if(a!=d)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_cinterval("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=complex(12,343);
            a=_cinterval(d);
            if(a!=d)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=complex(12,343);
            a=d;
            if(a!=d)
               error();
            else 
               ok();         
         }
      }
};

template <class D>
class cast_r_to_cinterval : public testclass
{
   public:
      cast_r_to_cinterval(void)
      {
         cinterval a;
         D       d,e;

         testing(nameof(a),"cinterval("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12),e=D(343);
            a=cinterval(d,e);
            if(Re(a)!=d || Im(a)!=e)
               error();
            else 
               ok();
         }

         testing(nameof(a),"_cinterval("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12),e=D(343);
            a=_cinterval(d,e);
            if(Re(a)!=d || Im(a)!=e)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"=("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12);
            a=d;
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12);
            a=cinterval(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         testing(nameof(a),"_"+nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12);
            a=_cinterval(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }         
      }
};


template <class D>
class cast_c_to_cidotprecision : public testclass
{
   public:
      cast_c_to_cidotprecision(void)
      {
         cidotprecision a;
         D       d;

         testing(nameof(a),"cidotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=complex(12,343);
            a=cidotprecision(d);
            if(a!=d)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_cidotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=complex(12,343);
            a=_cidotprecision(d);
            if(a!=d)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=complex(12,343);
            a=d;
            if(a!=d)
               error();
            else 
               ok();         
         }
      }
};

template <class D>
class cast_r_to_cidotprecision : public testclass
{
   public:
      cast_r_to_cidotprecision(void)
      {
         cidotprecision a;
         D       d,e;

         testing(nameof(a),"=("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12);
            a=d;
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12);
            a=cidotprecision(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }
         testing(nameof(a),"_"+nameof(a)+"("+nameof(d)+")");
         if(!tested())
         {
            a=cinterval(interval(1,2),interval(3,4));
            d=D(12);
            a=_cidotprecision(d);
            if(Re(a)!=d || Im(a)!=0)
               error();
            else 
               ok();
         }         
      }
};


/*
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
            
            a=-1024;
            s="<";
            s << a;
            if(s=="<  -1024.00000000000")
               nok++;
               
            a=123.456e+78;
            s="";
            s << a;
            if(s==" 1.23456000000E+080")
               nok++;
            
            a=123.456e-78;
            s="";
            s<<a;
            if(s=="-1.23456000000E-076");
               nok++;
            
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
            
            a=123.456e+78;
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
               
            if(s=="1.23456000000E+080")
               nok++;
            
            in >> s;
            
            if(s=="-1.23456000000E-076");
               nok++;
            
            if(nok!=4)
               error();
            else 
               ok();
         }
         
      }
};
*/

template <class T,class E>
class complexreimconj : public testclass
{
   public:
      complexreimconj(void)
      {
         T x1,x2,x3;
         E y1,y3;
         T y2;
         
         testing(nameof(x1),"Re("+nameof(x1)+")");
         if(!tested())
         {
            x1=complex(3,4);y1=3;
            
            x2=complex(-4,-3);
            Re(x2)=-7;
            y2=complex(-7,-3);
            
            x3=complex(0,0);
            Re(x3)=12;
            y3=12;

            if(Re(x1)==y1 && x2==y2 && Re(x3)==y3)
               ok();
            else
               error();
         }
         testing(nameof(x1),"Im("+nameof(x1)+")");
         if(!tested())
         {
            x1=complex(3,4);y1=4;
            x2=complex(-4,-3);
            Im(x2)=-7;
            y2=complex(-4,-7);
            x3=complex(0,0);
            Im(x3)=12;
            y3=12;

            if(Im(x1)==y1 && x2==y2 && Im(x3)==y3)
               ok();
            else
               error();
         }
         testing(nameof(x1),"SetRe("+nameof(x1)+","+nameof(y1)+")");
         if(!tested())
         {
            x1=complex(3,4);
            y1=4;
            SetRe(x1,y1);
            
            x2=complex(-4,-3);
            SetRe(x2,E(-7));
            y2=complex(-7,-3);
            
            x3=complex(0,0);
            SetRe(x3,E(12));
            y3=12;

            if(Re(x1)==y1 && x2==y2 && Re(x3)==y3)
               ok();
            else
               error();
         }
         testing(nameof(x1),"SetIm("+nameof(x1)+","+nameof(y1)+")");
         if(!tested())
         {
            x1=complex(3,4);
            y1=4;
            SetIm(x1,y1);
            
            x2=complex(-4,-3);
            SetIm(x2,E(-7));
            y2=complex(-4,-7);
            
            x3=complex(0,0);
            SetIm(x3,E(12));
            y3=12;

            if(Im(x1)==y1 && x2==y2 && Im(x3)==y3)
               ok();
            else
               error();
         }
         testing(nameof(x1),"conj("+nameof(x1)+")");
         if(!tested())
         {
            x1=complex(3,4);
            y1=-4;
            
            x2=complex(-4,-3);
            y2=complex(-4,3);
            
            x3=complex(0,0);
            SetIm(x3,E(12));
            y3=0;

            if(Im(conj(x1))==y1 && conj(x2)==y2 && Re(conj(x3))==y3)
               ok();
            else
               error();
         }
         
      }
};

template <class D,class A,class B>
class complexaccumulate : public testclass
{
   // A sollte der genauere Typ sein
   public:
      complexaccumulate(void)
      {
         D d;
         A a,e;
         B b;

         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            a=complex(1,2);
            b=complex(3,4);
            e=complex(-5,10);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
            else
               cout << a << "*" << b << "=" << d << " " << e << " " << a*b << endl;
           
            
            a=complex(-1,2);
            b=complex(3,4);
            e=complex(-11,2);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(2,-1);
            b=complex(3,4);
            e=complex(10,5);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=complex(-2,-1);
            b=complex(3,4);
            e=complex(-2,-11);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(1,2);
            b=complex(-3,4);
            e=complex(-11,-2);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            if(okz==5)
               ok();
            else
               error();
         }
      }
};


template <class D,class I,class R>
class complexmixaccumulate : public testclass
{
   // I Intervalltyp, R Realtyp
   public:
      complexmixaccumulate(void)
      {
         D d;
         I a,e;
         R b;

         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            a=complex(1,2);
            b=3;
            e=complex(3,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-1,2);
            b=3;
            e=complex(-3,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
            
            a=complex(-2,1);
            b=3;
            e=complex(-6,3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
               
            a=complex(-2,-1);
            b=3;
            e=complex(-6,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(1,2);
            b=-3;
            e=complex(-3,-6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-1,2);
            b=-3;
            e=complex(3,-6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-2,1);
            b=-3;
            e=complex(6,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-2,-1);
            b=-3;
            e=complex(6,3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            if(okz==8)
               ok();
            else
               error();
         }
         
         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            a=complex(1,2);
            b=3;
            e=complex(3,6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-1,2);
            b=3;
            e=complex(-3,6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
            
            a=complex(-2,1);
            b=3;
            e=complex(-6,3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
               
            a=complex(-2,-1);
            b=3;
            e=complex(-6,-3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(1,2);
            b=-3;
            e=complex(-3,-6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-1,2);
            b=-3;
            e=complex(3,-6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-2,1);
            b=-3;
            e=complex(6,-3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=complex(-2,-1);
            b=-3;
            e=complex(6,3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            if(okz==8)
               ok();
            else
               error();
         }
         
      }
};

} // namespace cxsc 

#endif // _CXSC_TESTCOMP_HPP_INCLUDED
