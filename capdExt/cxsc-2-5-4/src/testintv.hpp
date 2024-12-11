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

/* CVS $Id: testintv.hpp,v 1.24 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTINTV_HPP_INCLUDED
#define _CXSC_TESTINTV_HPP_INCLUDED

namespace cxsc {

template <class T,class S>
class intervalscalarassign : public testclass
{
   public:
      intervalscalarassign(void)
      {
         T x;
         S s,e;
         
         s=1;
         e=2;
         
         x=0;
         
         testing(nameof(x),"Inf("+nameof(x)+").operator =("+nameof(s)+")");
         if(!tested())
         {
            Sup(x)=s;
            Inf(x)=s;
            Sup(x)=e;
            
            if(Inf(x)==s)
               ok();
            else
               error();
         }
         
         x=0;
         
         s=10;
         e=5;
         
         testing(nameof(x),"Sup("+nameof(x)+").operator =("+nameof(s)+")");
         if(!tested())
         {
            Sup(x)=s;
            Inf(x)=e;
            
            if(Sup(x)==s)
               ok();
            else
               error();
         }
         x=0;
         
         s=10;
         e=5;
         
         testing(nameof(x),"SetInf("+nameof(x)+","+nameof(s)+")");
         if(!tested())
         {
            Sup(x)=e;
            SetInf(x,e);
            Sup(x)=s;
            
            if(Inf(x)==e)
               ok();
            else
               error();
         }
         
         x=0;
         
         s=10;
         e=5;
         
         testing(nameof(x),"SetSup("+nameof(x)+","+nameof(s)+")");
         if(!tested())
         {
            SetSup(x,s);
            Inf(x)=e;
            
            if(Sup(x)==s)
               ok();
            else
               error();
         }
         
         x=interval(1,2);
         
         /* testing(nameof(x),"exception EMPTY_INTERVAL");
         if(!tested())
         {
            int n=0;
            try                            { SetSup(x,-1); n--; }
            catch ( EMPTY_INTERVAL ) { n++; }
            try                            { SetInf(x,3); n--; }
            catch ( EMPTY_INTERVAL ) { n++; }
            try                            { SetSup(x,3); n++;  }
            catch ( EMPTY_INTERVAL ) { n--; }
            try                            { SetInf(x,-1); n++; }
            catch ( EMPTY_INTERVAL ) { n--; }
            x=interval(1,2);
            try                            { Sup(x)=S(-1); n--; }
            catch ( EMPTY_INTERVAL ) { n++; }
            try                            { Inf(x)=S(3); n--; }
            catch ( EMPTY_INTERVAL ) { n++; }
            try                            { Sup(x)=S(3); n++;  }
            catch ( EMPTY_INTERVAL ) { n--; }
            try                            { Inf(x)=S(-1); n++; }
            catch ( EMPTY_INTERVAL ) { n--; }
            
            if(n==8)
               ok();
            else
               error();
               
         }
         */         
      }
};

template <class T>
class intervalallscalarassign : public testclass
{
   public:
      intervalallscalarassign(void)
      {
         intervalscalarassign<T,int>    with_int;
         intervalscalarassign<T,long>   with_long;
         intervalscalarassign<T,float>  with_float;
         intervalscalarassign<T,double> with_double;
         intervalscalarassign<T,real>   with_real;

         if(!with_int || !with_float || !with_double ||
            !with_real || !with_long)
            setfail();
      }
};



template <class A,class B>
class intervaladdsub : public testclass
{
   // Tests A+B, B+A, A-B, B-A, A+=B, A-=B
   //       -A, +A
   public:
      intervaladdsub(void)
      {
         A a,e;
         B b;

         intervalscalarassign<A,int> ai;
         intervalscalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         
         testing(nameof(e),"operator +("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=interval(10,20);
            b=interval(13,17);
            e=interval(23,37);
            if(a+b>=e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator +("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a=interval(10,20);
            b=interval(13,17);
            e=interval(23,37);
            if(b+a>=e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator +=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            a=interval(10,20);
            b=interval(13,17);
            e=interval(23,37);

            if((a+=b)>=e && a>=e)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator -("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            a=interval(10,20);
            b=interval(13,16);
            e=interval(-3,4);
            
            if(a-b>=e)
               ok();
            else
               error(),cerr << a << "-" << b << "=" << a-b << "!=" << e << endl;
         }
         
         testing(nameof(e),"operator -("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            a=interval(10,20);
            b=interval(13,16);
            e=interval(-7,6);
            
            if(b-a>=e)
               ok();
            else
               error(),cerr << b << "-" << a << "=" << b-a << "!=" << e << endl;
         }
         
         testing(nameof(e),"operator -=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            a=interval(10,20);
            b=interval(13,16);
            e=interval(-3,4);
            
            if((a-=b)>=e && a>=e)
               ok();
            else
               error();
         }

         testing(nameof(e),"operator -("+nameof(a)+")");
         if(!tested())
         {
            a=interval(10,20);
            e=interval(-20,-10);
            
            if((-a)==e)
               ok();
            else
               error();
         }

         testing(nameof(e),"operator +("+nameof(a)+")");
         if(!tested())
         {
            a=interval(10,20);
            e=interval(10,20);
            
            if((+a==e))
               ok();
            else
               error();
         }
      }
};

template <class A,class B>
class intervalmuldiv : public testclass
{
   // Tests A*B, B*A, A/B, B/A, A*=B, A/=B
   public:
      intervalmuldiv(void)
      {
         A a,e;
         B b;

         intervalscalarassign<A,int> ai;
         intervalscalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator *("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            a=interval(1,2);
            b=interval(3,4);
            e=interval(3,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(3,4);
            e=interval(-4,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(3,4);
            e=interval(-8,4);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,-1);
            b=interval(3,4);
            e=interval(-8,-3);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,4);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-3,4);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-4,3);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,3);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(1,2);
            b=interval(-4,-3);
            e=interval(-8,-3);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-1,2);
            b=interval(-4,-3);
            e=interval(-8,4);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-4,-3);
            e=interval(-4,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,-3);
            e=interval(3,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            if(okz==16)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator *("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            A b;
            B a; // Just vice versa... ;)
            
            a=interval(1,2);
            b=interval(3,4);
            e=interval(3,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(3,4);
            e=interval(-4,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(3,4);
            e=interval(-8,4);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,-1);
            b=interval(3,4);
            e=interval(-8,-3);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,4);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-3,4);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-4,3);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,3);
            e=interval(-6,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(1,2);
            b=interval(-4,-3);
            e=interval(-8,-3);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-1,2);
            b=interval(-4,-3);
            e=interval(-8,4);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-4,-3);
            e=interval(-4,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,-3);
            e=interval(3,8);

            if(a*b==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            if(okz==16)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator *=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(1,2);
            b=interval(3,4);
            e=interval(3,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(3,4);
            e=interval(-4,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(3,4);
            e=interval(-8,4);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,-1);
            b=interval(3,4);
            e=interval(-8,-3);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,4);
            e=interval(-8,6);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-3,4);
            e=interval(-8,6);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-4,3);
            e=interval(-6,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,3);
            e=interval(-6,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(1,2);
            b=interval(-4,-3);
            e=interval(-8,-3);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-1,2);
            b=interval(-4,-3);
            e=interval(-8,4);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-4,-3);
            e=interval(-4,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,-3);
            e=interval(3,8);

            if((a*=b)==e && a==e)
               okz++;
            else cerr << a << "*" << b << "=" << a*b << "!=" << e << endl;
               
            if(okz==16)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            a=interval(1,2);
            b=interval(4,8);
            e=interval(1./8,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(4,8);
            e=interval(-1./4,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(4,8);
            e=interval(-1./2,1./4);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-2,-1);
            b=interval(4,8);
            e=interval(-1./2,-1./8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            /*a=interval(1,2);
            b=interval(-4,8);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-4,8);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-4,8);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,8);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-8,4);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-8,4);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-8,4);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-8,4);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            */   
            a=interval(1,2);
            b=interval(-8,-4);
            e=interval(-1./2,-1./8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-1,2);
            b=interval(-8,-4);
            e=interval(-1./2,1./4);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-8,-4);
            e=interval(-1./4,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-8,-4);
            e=interval(1./8,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            if(okz==8) // Noch ohne Exceptions fuer DivZero
               ok();
            else
               error();
         }
         
         
         testing(nameof(e),"operator /("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            A b;
            B a; // et vice versa
            
            a=interval(1,2);
            b=interval(4,8);
            e=interval(1./8,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(4,8);
            e=interval(-1./4,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(4,8);
            e=interval(-1./2,1./4);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-2,-1);
            b=interval(4,8);
            e=interval(-1./2,-1./8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            /*a=interval(1,2);
            b=interval(-4,8);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-4,8);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-4,8);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,8);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-8,4);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-8,4);
            e=interval(-8,6);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-8,4);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-8,4);
            e=interval(-6,8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            */   
            a=interval(1,2);
            b=interval(-8,-4);
            e=interval(-1./2,-1./8);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-1,2);
            b=interval(-8,-4);
            e=interval(-1./2,1./4);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-8,-4);
            e=interval(-1./4,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-8,-4);
            e=interval(1./8,1./2);

            if(a/b==e)
               okz++;
            else cerr << a << "/" << b << "=" << a/b << "!=" << e << endl;
               
            if(okz==8) // Noch ohne Exceptions fuer DivZero
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator /=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(1,2);
            b=interval(4,8);
            e=interval(1./8,1./2);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(4,8);
            e=interval(-1./4,1./2);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(4,8);
            e=interval(-1./2,1./4);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
               
            a=interval(-2,-1);
            b=interval(4,8);
            e=interval(-1./2,-1./8);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            /*a=interval(1,2);
            b=interval(-4,8);
            e=interval(-6,8);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-4,8);
            e=interval(-6,8);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-4,8);
            e=interval(-8,6);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-4,8);
            e=interval(-8,6);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(1,2);
            b=interval(-8,4);
            e=interval(-8,6);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-1,2);
            b=interval(-8,4);
            e=interval(-8,6);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-8,4);
            e=interval(-6,8);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-8,4);
            e=interval(-6,8);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            */   
            a=interval(1,2);
            b=interval(-8,-4);
            e=interval(-1./2,-1./8);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
               
            a=interval(-1,2);
            b=interval(-8,-4);
            e=interval(-1./2,1./4);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
               
            a=interval(-2,1);
            b=interval(-8,-4);
            e=interval(-1./4,1./2);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,-1);
            b=interval(-8,-4);
            e=interval(1./8,1./2);

            if((a/=b)==e && a==e)
               okz++;
            else cerr << "a/=" << b << "=" << a << "!=" << e << endl;
               
            if(okz==8) // Noch ohne Exceptions fuer DivZero
               ok();
            else
               error();
         }
      }
};

template <class A,class B>
class intervalsetops : public testclass
{
   // Tests A|B, B|A, A&B, B&A, A|=B, A&=B
   public:
      intervalsetops(void)
      {
         A a,e;
         B b;

         intervalscalarassign<A,int> ai;
         intervalscalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator |("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=interval(-4,-3);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,-1);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,4);
            e=interval(-4,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=interval(-2,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(1,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(3,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,-3);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,-1);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,4);
            e=interval(-4,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-1,1);
            e=interval(-2,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(1,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(3,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            if(okz==12)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator |("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            B a; 
            A b; // And once again...
            
            int okz=0;
            
            a=interval(-2,2);
            b=interval(-4,-3);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl; 
                      
            
            a=interval(-2,2);
            b=interval(-4,-1);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,4);
            e=interval(-4,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=interval(-2,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(1,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(3,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,-3);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,-1);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,4);
            e=interval(-4,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-1,1);
            e=interval(-2,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(1,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(3,4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            if(okz==12)
               ok();
            else
               error();
         }  
                
         testing(nameof(e),"operator |=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=interval(-4,-3);
            e=interval(-4,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,-1);
            e=interval(-4,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,4);
            e=interval(-4,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=interval(-2,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(1,4);
            e=interval(-2,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(3,4);
            e=interval(-2,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,-3);
            e=interval(-4,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,-1);
            e=interval(-4,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,4);
            e=interval(-4,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-1,1);
            e=interval(-2,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(1,4);
            e=interval(-2,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(3,4);
            e=interval(-2,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            if(okz==12)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator &("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            /* a=interval(-2,2);
            b=interval(-4,-3);
            e=interval(-4,2); // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
            a=interval(-2,2);
            b=interval(-4,-1);
            e=interval(-2,-1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,4);
            e=interval(-2,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=interval(-1,1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(1,4);
            e=interval(1,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            /*a=interval(-2,2);
            b=interval(3,4);
            e=interval(-2,4);  // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
            
            /* b=interval(-2,2);
            a=interval(-4,-3);
            e=interval(-4,2); // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
                        
            b=interval(-2,2);
            a=interval(-4,-1);
            e=interval(-2,-1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,4);
            e=interval(-2,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-1,1);
            e=interval(-1,1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(1,4);
            e=interval(1,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            /*b=interval(-2,2);
            a=interval(3,4);
            e=interval(-2,4); // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
            if(okz==8)
               ok();
            else
               error();
         }
                  
         testing(nameof(e),"operator &("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            B a; 
            A b;
            
            /* a=interval(-2,2);
            b=interval(-4,-3);
            e=interval(-4,2); // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
            a=interval(-2,2);
            b=interval(-4,-1);
            e=interval(-2,-1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,4);
            e=interval(-2,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=interval(-1,1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(1,4);
            e=interval(1,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            /*a=interval(-2,2);
            b=interval(3,4);
            e=interval(-2,4);  // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
            
            /* b=interval(-2,2);
            a=interval(-4,-3);
            e=interval(-4,2); // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
                        
            b=interval(-2,2);
            a=interval(-4,-1);
            e=interval(-2,-1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,4);
            e=interval(-2,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-1,1);
            e=interval(-1,1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(1,4);
            e=interval(1,2);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            /*b=interval(-2,2);
            a=interval(3,4);
            e=interval(-2,4); // empty

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            */
            if(okz==8)
               ok();
            else
               error();
         }         

         testing(nameof(e),"operator &=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            /* a=interval(-2,2);
            b=interval(-4,-3);
            e=interval(-4,2); // empty

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            */
            a=interval(-2,2);
            b=interval(-4,-1);
            e=interval(-2,-1);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-4,4);
            e=interval(-2,2);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=interval(-1,1);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(1,4);
            e=interval(1,2);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            /*a=interval(-2,2);
            b=interval(3,4);
            e=interval(-2,4);  // empty

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            */
            
            /* b=interval(-2,2);
            a=interval(-4,-3);
            e=interval(-4,2); // empty

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            */
                        
            b=interval(-2,2);
            a=interval(-4,-1);
            e=interval(-2,-1);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-4,4);
            e=interval(-2,2);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(-1,1);
            e=interval(-1,1);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            b=interval(-2,2);
            a=interval(1,4);
            e=interval(1,2);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            
            /*b=interval(-2,2);
            a=interval(3,4);
            e=interval(-2,4); // empty

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << a << "!=" << e << endl;
            */
            if(okz==8)
               ok();
            else
               error();
         }       
      }
};


template <class A,class B>
class intervalmixsetops : public testclass
{
   // Tests A|B, B|A, A&B, B&A, A|=B, A&=B
   // B is scalar-type
   public:
      intervalmixsetops(void)
      {
         A a,e;
         B b;

         intervalscalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator |("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=B(-4);
            e=interval(-4,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(1);
            e=interval(-2,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(4);
            e=interval(-2,4);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            if(okz==3)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator |("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=B(-4);
            e=interval(-4,2);

            if((b|a)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (b|a) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(1);
            e=interval(-2,2);

            if((b|a)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (b|a) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(4);
            e=interval(-2,4);

            if((b|a)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (b|a) << "!=" << e << endl;
            
            if(okz==3)
               ok();
            else
               error();
         }
                
         testing(nameof(e),"operator |=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=B(-4);
            e=interval(-4,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(1);
            e=interval(-2,2);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(4);
            e=interval(-2,4);

            if((a|=b)==e && a==e)
               okz++;
            else cerr << "a|=" << b << "=" << a << "!=" << e << endl;
            
            if(okz==3)
               ok();
            else
               error();
         }
         
         testing(nameof(e),"operator &("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=B(-4);
            e=interval(-4,2);// empty

            try 
            {
               a&b;
            }
            catch (EMPTY_INTERVAL)
            {
               okz++;
            }


            a=interval(-2,2);
            b=B(1);
            e=interval(1,1);

            if((a&b)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (a&b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(4);
            e=interval(-2,2); // empty

            try {
               a&b;
            }
            catch( EMPTY_INTERVAL )
            {
               okz++;
            }
            
            if(okz==3)
               ok();
            else
               error();
         }
                  
         testing(nameof(e),"operator &("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=B(-4);
            e=interval(-4,2);// empty

            try 
            {
               b&a;
            }
            catch ( EMPTY_INTERVAL )
            {
               okz++;
            }


            a=interval(-2,2);
            b=B(1);
            e=interval(1,1);

            if((b&a)==e)
               okz++;
            else cerr << a << "&" << b << "=" << (b&a) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(4);
            e=interval(-2,2); // empty

            try {
               b&a;
            }
            catch( EMPTY_INTERVAL )
            {
               okz++;
            }
            
            if(okz==3)
               ok();
            else
               error();
         }
         testing(nameof(e),"operator &=("+nameof(e)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-2,2);
            b=B(-4);
            e=interval(-4,2);// empty

            try 
            {
               a&=b;
            }
            catch ( EMPTY_INTERVAL )
            {
               okz++;
            }


            a=interval(-2,2);
            b=B(1);
            e=interval(1,1);

            if((a&=b)==e && a==e)
               okz++;
            else cerr << "a&=" << b << "=" << (a) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(4);
            e=interval(-2,2); // empty

            try {
               a&=b;
            }
            catch( EMPTY_INTERVAL )
            {
               okz++;
            }
            
            if(okz==3)
               ok();
            else
               error();         
         }       
      }
};

template <class A,class B>
class scalarmixsetops : public testclass
{
   // Tests A|B, B|A
   // A and B are scalar-type
   public:
      scalarmixsetops(void)
      {
         A a;
         B b;
         interval e;

         scalarassign<A,int> ai;
         scalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(e),"operator |("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=A(-2);
            b=B(-4);
            e=interval(-4,-2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=A(-2);
            b=B(-2);
            e=interval(-2,-2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            a=A(-2);
            b=B(2);
            e=interval(-2,2);

            if((a|b)==e)
               okz++;
            else cerr << a << "|" << b << "=" << (a|b) << "!=" << e << endl;
            
            if(okz==3)
               ok();
            else
               error();
         }

         testing(nameof(e),"operator |("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=A(-2);
            b=B(-4);
            e=interval(-4,-2);

            if((b|a)==e)
               okz++;
            else cerr << b << "|" << a << "=" << (b|a) << "!=" << e << endl;
            
            a=A(-2);
            b=B(-2);
            e=interval(-2,-2);

            if((b|a)==e)
               okz++;
            else cerr << b << "|" << a << "=" << (b|a) << "!=" << e << endl;
            
            a=A(-2);
            b=B(2);
            e=interval(-2,2);

            if((b|a)==e)
               okz++;
            else cerr << b << "|" << a << "=" << (b|a) << "!=" << e << endl;
            
            if(okz==3)
               ok();
            else
               error();
         }
      }
};

template <class T>
class intervalallscalarmixsetops : public testclass
{
   public:
      intervalallscalarmixsetops(void)
      {
         intervalmixsetops<T,int>    with_int;
         intervalmixsetops<T,long>   with_long;
         intervalmixsetops<T,float>  with_float;
         intervalmixsetops<T,double> with_double;
         intervalmixsetops<T,real>   with_real;

         if(!with_int || !with_float || !with_double ||
            !with_real || !with_long)
            setfail();
      }
};

template <class A,class B>
class intervalsetcompops : public testclass
{
   // Tests A<B, A>B, A<=B, A>=B, A==B, A!=B
   //       B<A, B>A, B<=A, B>=A, B==A, B!=A
   
   public:
      intervalsetcompops(void)
      {
         A a;
         B b;
         bool e;

         intervalscalarassign<A,int> ai;
         intervalscalarassign<B,int> bi;
         
         setfail(!ai || !bi);
         
         testing(nameof(a),"operator <("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=true;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator >("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=true;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator <=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator >=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator ==("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=true;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a==b)==e)
               okz++;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         
         testing(nameof(a),"operator !=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=false;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=true;

            if((a!=b)==e)
               okz++;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator <("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            B a;
            A b;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=true;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator >("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            B a;
            A b;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=true;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator <=("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            B a;
            A b;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator >=("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            B a;
            A b;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(b),"operator ==("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=true;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=false;

            if((b==a)==e)
               okz++;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         
         testing(nameof(b),"operator !=("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=interval(-1,2);
            e=false;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,1);
            b=interval(-3,2);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,1);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-1,2);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,1);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-2,3);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=interval(-3,2);
            e=true;

            if((b!=a)==e)
               okz++;
            
            if(okz==7)
               ok();
            else
               error();
         }
      }
};


template <class A,class B>
class intervalscalarsetcompops : public testclass
{
   // Tests A<B, A>B, A<=B, A>=B, A==B
   //       B<A, B>A, B<=A, B>=A, B==A
   // where B is skalar
   public:
      intervalscalarsetcompops(void)
      {
         A a;
         B b;
         bool e;

         testing(nameof(a),"operator <("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=false;

            if((a<b)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator >("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;
   
            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=true;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=false;

            if((a>b)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();

         }
         
         testing(nameof(a),"operator <=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=true;

            if((a<=b)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();

         }
         
         testing(nameof(a),"operator >=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;
   
            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=true;

            if((a>=b)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();
         }

         testing(nameof(a),"operator ==("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;
   
            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(0);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((a==b)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(-1);
            e=true;

            if((a==b)==e)
               okz++;
             
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator !=("+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=true;
   
            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(-2);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(0);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(2);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(3);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(1);
            e=true;

            if((a!=b)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(-1);
            e=false;

            if((a!=b)==e)
               okz++;
             
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator <("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=true;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=false;

            if((b<a)==e)
               okz++;
            else cerr << a << "<" << b << "=" << (a<b) << "!=" << e << endl;
            
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(a),"operator >("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;
   
            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=false;

            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=false;

            if((b>a)==e)
               okz++;
            else cerr << a << ">" << b << "=" << (a>b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();

         }
         
         testing(nameof(a),"operator <=("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=true;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=true;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=true;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=true;

            if((b<=a)==e)
               okz++;
            else cerr << a << "<=" << b << "=" << (a<=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();

         }
         
         testing(nameof(a),"operator >=("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;
   
            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(0);
            e=false;

            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            a=interval(-1,-1);
            b=B(-1);
            e=true;

            if((b>=a)==e)
               okz++;
            else cerr << a << ">=" << b << "=" << (a>=b) << "!=" << e << endl;
            
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(b),"operator ==("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=false;
   
            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(-2);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(0);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(2);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(3);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(1);
            e=false;

            if((b==a)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(-1);
            e=true;

            if((b==a)==e)
               okz++;
             
            if(okz==7)
               ok();
            else
               error();
         }
         
         testing(nameof(b),"operator !=("+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            
            a=interval(-1,2);
            b=B(-2);
            e=true;
   
            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(-2);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(0);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(2);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-2,2);
            b=B(3);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(1);
            e=true;

            if((b!=a)==e)
               okz++;
            
            a=interval(-1,-1);
            b=B(-1);
            e=false;

            if((b!=a)==e)
               okz++;
             
            if(okz==7)
               ok();
            else
               error();
         }
         
      }
};




template <class T>
class intervalstdfunc : public testclass
{
   public:
      intervalstdfunc(void)
      {
         T x1;
         testing(nameof(x1),"sqr("+nameof(x1)+")");
         
         if(!tested())
         {
            if(comparei(T(11),sqr(T(11)),T(121))
            && comparei(T(0),sqr(T(0)),T(0))
            && comparei(T(-9),sqr(T(-9)),T(81)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sqrt("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(121),sqrt(T(121)),T(11))
            && comparei(T(0),sqrt(T(0)),T(0))
            && comparei(T(81),sqrt(T(81)),T(9)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sqrt("+nameof(x1)+",int)");
         if(!tested())
         {
            if(comparei(T(27),T(3),sqrt(T(27),int(3)),T(3))
            && comparei(T(0),T(4),sqrt(T(0),int(4)),T(0))
            && comparei(T(1024),T(10),sqrt(T(1024),int(10)),T(2)) )
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
         
         T pi,pi_2,pi_4,e,e1d,ee;
         CSTRING_PI>>pi;
         CSTRING_PI_2>>pi_2;
         CSTRING_PI_4>>pi_4;
         CSTRING_E>>e;
         e1d=T(1.)/e;
         ee=e*e;

         testing(nameof(x1),"sin("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),sin(T(0)),T(0))
            && comparei(T(pi_2),sin(T(pi_2)),T(1))
            && comparei(T(-pi),sin(T(-pi)),T(0)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"cos("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),cos(T(0)),T(1))
            && comparei(T(pi_2),cos(T(pi_2)),T(0))
            && comparei(T(-pi),cos(T(-pi)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"tan("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),tan(T(0)),T(0))
            && comparei(T(pi_4),tan(T(pi_4)),T(1))
            && comparei(T(-pi_4),tan(T(-pi_4)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"cot("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(pi_2),cot(T(pi_2)),T(0))
            && comparei(T(pi_4),cot(T(pi_4)),T(1))
            && comparei(T(-pi_4),cot(T(-pi_4)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"asin("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),asin(T(0)),T(0))
            && comparei(T(1),asin(T(1)),T(pi_2))
            && comparei(T(-1),asin(T(-1)),T(-pi_2)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acos("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),acos(T(0)),T(pi_2))
            && comparei(T(1),acos(T(1)),T(0))
            && comparei(T(-1),acos(T(-1)),T(pi)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"atan("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),atan(T(0)),T(0))
            && comparei(T(1),atan(T(1)),T(pi_4))
            && comparei(T(-1),atan(T(-1)),T(-pi_4)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acot("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),acot(T(0)),T(pi_2))
            && comparei(T(1),acot(T(1)),T(pi_4))
            && comparei(T(-1),acot(T(-1)),T(pi_4+pi_2)) )
               ok();
            else
               error();
         }
         /*testing(nameof(x1),"expm1("+nameof(x1)+")");
         if(!tested())
         {
            x1=0,e1=0;
            x2=1,e2=EM1;
            x3=-1, e3=E1D-1;
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
            if(comparei(T(0),exp(T(0)),T(1))
            && comparei(T(1),exp(T(1)),T(e))
            && comparei(T(-1),exp(T(-1)),T(e1d)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"ln("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(e),ln(T(e)),T(1))
            && comparei(T(1),ln(T(1)),T(0))
            && comparei(T(ee),ln(T(ee)),T(2)) )
               ok();
            else
               error();
         }
         
         testing(nameof(x1),"sinh("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),sinh(T(0)),T(0))
            && comparei(T(1),sinh(T(1)),T( (e-1/e)/2 ))
            && comparei(T(-1),sinh(T(-1)),T( (1/e-e)/2 )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"cosh("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),cosh(T(0)),T(1))
            && comparei(T(1),cosh(T(1)),T( (e+1/e)/2 ))
            && comparei(T(-1),cosh(T(-1)),T( (e+1/e)/2 )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"tanh("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),tanh(T(0)),T(0))
            && comparei(T(1),tanh(T(1)),T( (e-1/e)/(e+1/e) ))
            && comparei(T(-1),tanh(T(-1)),T( (1/e-e)/(e+1/e) )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"coth("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(2),coth(T(2)),T( cosh(T(2))/sinh(T(2)) )) // =:}
            && comparei(T(1),coth(T(1)),T( (e+1/e)/(e-1/e) ))
            && comparei(T(-1),coth(T(-1)),T( (e+1/e)/(1/e-e) )) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"asinh("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),asinh(T(0)),T(0))
            && comparei(T((e-1/e)/2),asinh(T((e-1/e)/2)),T(1))
            && comparei(T((1/e-e)/2),asinh(T((1/e-e)/2)),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acosh("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(1),acosh(T(1)),T(0))
            && comparei(T((e+1/e)/2),acosh(T((e+1/e)/2)),T(1))
            && comparei(T(cosh(T(2))),acosh(T(cosh(T(2)))),T(2)) ) // =:}
               ok();
            else
               error();
         }
         testing(nameof(x1),"atanh("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(0),atanh(T(0)),T(0))
            && comparei(T((e-1/e)/(e+1/e)),atanh(T((e-1/e)/(e+1/e))),T(1))
            && comparei(T((1/e-e)/(e+1/e)),atanh(T((1/e-e)/(e+1/e))),T(-1)) )
               ok();
            else
               error();
         }
         testing(nameof(x1),"acoth("+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T((e+1/e)/(e-1/e)),acoth(T((e+1/e)/(e-1/e))),T(1))
            && comparei(T((e+1/e)/(1/e-e)),acoth(T((e+1/e)/(1/e-e))),T(-1))
            && comparei(T(coth(T(2))),acoth(T(coth(T(2)))),T(2)) )
               ok();
            else
               error();
         }
         


         testing(nameof(x1),"pow("+nameof(x1)+","+nameof(x1)+")");
         if(!tested())
         {
            if(comparei(T(2),T(2),pow(T(2),T(2)),T(4))
            && comparei(T(4),T(5),pow(T(4),T(5)),T(1024))
            && comparei(T(+2),T(3),pow(T(+2),T(3)),T(+8)) )// Negativ geht noch nicht
               ok();
            else
               error();
         }
         testing(nameof(x1),"power("+nameof(x1)+",int)");
         if(!tested())
         {
            if(comparei(T(2),T(2),power(T(2),int(2)),T(4))
            && comparei(T(4),T(5),power(T(4),int(5)),T(1024))
            && comparei(T(+2),T(3),power(T(+2),int(3)),T(+8)) ) // Negativ geht noch nicht
               ok();
            else
               error();
         }
                                                   
         
      }
};

template <class D,class A,class B>
class intervalaccumulate : public testclass
{
   // A sollte der genauere Typ sein
   public:
      intervalaccumulate(void)
      {
         D d;
         A a,e;
         B b;

         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            a=interval(1,2);
            b=interval(3,4);
            e=interval(3,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
            else
               cout << a << "*" << b << "=" << d << " " << e << " " << a*b << endl;
           
            
            a=interval(-1,2);
            b=interval(3,4);
            e=interval(-4,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,1);
            b=interval(3,4);
            e=interval(-8,4);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-2,-1);
            b=interval(3,4);
            e=interval(-8,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,1);
            b=interval(-3,4);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=interval(-3,4);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-2,1);
            b=interval(-4,3);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=interval(-4,3);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(1,2);
            b=interval(-4,-3);
            e=interval(-8,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-1,2);
            b=interval(-4,-3);
            e=interval(-8,4);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-2,1);
            b=interval(-4,-3);
            e=interval(-4,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=interval(-4,-3);
            e=interval(3,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            if(okz==16)
               ok();
            else
               error();
         }
         
         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(b)+","+nameof(a)+")");
         if(!tested())
         {
            int okz=0;
            B a;
            A b;
            
            a=interval(1,2);
            b=interval(3,4);
            e=interval(3,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=interval(3,4);
            e=interval(-4,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,1);
            b=interval(3,4);
            e=interval(-8,4);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-2,-1);
            b=interval(3,4);
            e=interval(-8,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=interval(-3,4);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,1);
            b=interval(-3,4);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=interval(-3,4);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=interval(-4,3);
            e=interval(-8,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-2,1);
            b=interval(-4,3);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=interval(-4,3);
            e=interval(-6,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(1,2);
            b=interval(-4,-3);
            e=interval(-8,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-1,2);
            b=interval(-4,-3);
            e=interval(-8,4);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
               
            a=interval(-2,1);
            b=interval(-4,-3);
            e=interval(-4,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=interval(-4,-3);
            e=interval(3,8);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            if(okz==16)
               ok();
            else
               error();
         }
         
      }
};


template <class D,class I,class R>
class intervalmixaccumulate : public testclass
{
   // I Intervalltyp, R Realtyp
   public:
      intervalmixaccumulate(void)
      {
         D d;
         I a,e;
         R b;

         testing(nameof(d),"accumulate("+nameof(d)+","+nameof(a)+","+nameof(b)+")");
         if(!tested())
         {
            int okz=0;
            a=interval(1,2);
            b=3;
            e=interval(3,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=3;
            e=interval(-3,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
            
            a=interval(-2,1);
            b=3;
            e=interval(-6,3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
               
            a=interval(-2,-1);
            b=3;
            e=interval(-6,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(1,2);
            b=-3;
            e=interval(-6,-3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=-3;
            e=interval(-6,3);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,1);
            b=-3;
            e=interval(-3,6);

            d=0;
            accumulate(d,a,b);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=-3;
            e=interval(3,6);

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
            a=interval(1,2);
            b=3;
            e=interval(3,6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=3;
            e=interval(-3,6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
            
            a=interval(-2,1);
            b=3;
            e=interval(-6,3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
               
            a=interval(-2,-1);
            b=3;
            e=interval(-6,-3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(1,2);
            b=-3;
            e=interval(-6,-3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-1,2);
            b=-3;
            e=interval(-6,3);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,1);
            b=-3;
            e=interval(-3,6);

            d=0;
            accumulate(d,b,a);
            if(d==e && d==a*b)
               okz++;
           
            
            a=interval(-2,-1);
            b=-3;
            e=interval(3,6);

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


template <class D,class G>
class cast_scalar_to_interval : public testclass
{
   public:
      cast_scalar_to_interval(void)
      {
         interval a;
         D    d,e;
         
         testing(nameof(a),"interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=interval(d);
            if(a!=12 || interval(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_interval(d);
            if(a!=12 || _interval(d)!=12)
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
         
         testing(nameof(a),"interval("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=interval(d,e);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_interval("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=_interval(d,e);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         G g;
         testing(nameof(a),"interval("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=interval(d,g);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_interval("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=_interval(d,g);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         testing(nameof(a),"interval("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=interval(g,d);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_interval("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=_interval(g,d);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
      }
};

template <class D>
class cast_scalar_to_l_interval : public testclass
{
   public:
      cast_scalar_to_l_interval(void)
      {
         l_interval a;
         D    d,e;
         
         testing(nameof(a),"l_interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=l_interval(d);
            if(a!=12 || l_interval(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_l_interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_l_interval(d);
            if(a!=12 || _l_interval(d)!=12)
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
         
         testing(nameof(a),"l_interval("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=l_interval(d,e);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_l_interval("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=_l_interval(d,e);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
      }
};

template <class D,class G>
class cast_scalar_to_idotprecision : public testclass
{
   public:
      cast_scalar_to_idotprecision(void)
      {
         idotprecision a;
         D    d,e;
         
         testing(nameof(a),"idotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=idotprecision(d);
            if(a!=12 || idotprecision(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_idotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_idotprecision(d);
            if(a!=12 || _idotprecision(d)!=12)
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
         
         testing(nameof(a),"idotprecision("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=idotprecision(d,e);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_idotprecision("+nameof(d)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=12;
            e=19;
            a=_idotprecision(d,e);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         G g;
         testing(nameof(a),"idotprecision("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=idotprecision(d,g);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_idotprecision("+nameof(d)+","+nameof(g)+")");
         if(!tested())
         {
            a=123;
            d=12;
            g=19;
            a=_idotprecision(d,g);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         testing(nameof(a),"idotprecision("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=idotprecision(g,d);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_idotprecision("+nameof(g)+","+nameof(d)+")");
         if(!tested())
         {
            a=123;
            d=19;
            g=12;
            a=_idotprecision(g,d);
            if(Inf(a)!=12 || Sup(a)!=19)
               error();
            else 
               ok();
         }
         
      }
};

template <class D>
class cast_to_interval : public testclass
{
   public:
      cast_to_interval(void)
      {
         interval a;
         D    d;

         testing(nameof(a),"interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=interval(d);
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=_interval(d);
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();         
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=d;
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();         
         }
      }
};


template <class D>
class cast_to_l_interval : public testclass
{
   public:
      cast_to_l_interval(void)
      {
         l_interval a;
         D    d;

         testing(nameof(a),"l_interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=l_interval(d);
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_l_interval("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=_l_interval(d);
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();         
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=d;
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();         
         }
      }
};


template <class D>
class cast_to_idotprecision : public testclass
{
   public:
      cast_to_idotprecision(void)
      {
         idotprecision a;
         D    d;

         testing(nameof(a),"idotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=idotprecision(d);
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_idotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=_idotprecision(d);
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();         
         }
         
         testing(nameof(a),"operator =("+nameof(a)+","+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=D(17,99);
            a=d;
            if(Inf(a)!=17 || Sup(a)!=99)
               error();
            else 
               ok();         
         }
      }
};

template <class D,class E>
class testabsmiddiam : public testclass
{
   public:
      testabsmiddiam(void)
      {
         D    d;
         E    e;

         testing(nameof(d),"abs("+nameof(d)+")");
         if(!tested())
         {
            D x1(-2,-1),y1(1,2);
            D x2(-1,2),y2(0,2);
            D x3(3,4),y3(3,4);
            D x4(-3,1),y4(0,3);

            if(abs(x1)!=y1 || abs(x2)!=y2 || abs(x3)!=y3 || abs(x4)!=y4)
               error();
            else 
               ok();
         }
         
         testing(nameof(e),"mid("+nameof(d)+")");
         if(!tested())
         {
            D x1(-4,-2),x2(-4,2),x3(-2,4),x4(2,4);
            E y1(-3)   ,y2(-1)  ,y3(1)   ,y4(3);

            if(mid(x1)!=y1 || mid(x2)!=y2 || mid(x3)!=y3 || mid(x4)!=y4)
               error();
            else 
               ok();
         }
         
         testing(nameof(e),"diam("+nameof(d)+")");
         if(!tested())
         {
            D x1(-4,-2),x2(-4,2),x3(-2,4),x4(2,4);
            E y1(2)   ,y2(6)  ,y3(6)   ,y4(2);

            if(diam(x1)!=y1 || diam(x2)!=y2 || diam(x3)!=y3 || diam(x4)!=y4)
               error();
            else 
               ok();
         }
         
      }
};

template <class D,class E>
class test_unchecked : public testclass
{
   public:
      test_unchecked(void)
      {
         D    d;
         E    e;

         testing(nameof(d),"UncheckedSetInf("+nameof(d)+","+nameof(e)+")");
         if(!tested())
         {  
            d=interval(12,75);
            e=345;
            
            if(Inf(UncheckedSetInf(d,e))!=e || Inf(d)!=e || Sup(d)!=75)
               error();
            else 
               ok();   
         }
         testing(nameof(d),"UncheckedSetSup("+nameof(d)+","+nameof(e)+")");
         if(!tested())
         {  
            d=interval(12,75);
            e=5;
            
            if(Sup(UncheckedSetSup(d,e))!=e || Sup(d)!=e || Inf(d)!=12)
               error();
            else 
               ok();   
         }
         testing(nameof(d),"IsEmpty("+nameof(d)+")");
         if(!tested())
         {  
            D i1(-5,6),i2(-5,6),i3(-5,6),i4(-5,6);

            UncheckedSetInf(i1,8);
            UncheckedSetInf(i2,2);
            UncheckedSetSup(i3,-10);
            UncheckedSetSup(i4,-5);
            
            if(!IsEmpty(i1) || IsEmpty(i2) || !IsEmpty(i3) || IsEmpty(i4))
               error();
            else 
               ok();   
         }
      }
};

} // namespace cxsc 

#endif // _CXSC_TESTINTV_HPP_INCLUDED
