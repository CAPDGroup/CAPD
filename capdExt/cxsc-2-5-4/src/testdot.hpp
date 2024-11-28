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

/* CVS $Id: testdot.hpp,v 1.25 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTDOT_HPP_INCLUDED
#define _CXSC_TESTDOT_HPP_INCLUDED

#include <fstream>

namespace cxsc {

// ------------- Dot - Tests ---------------------------


template <class D,class A,class B>
class scalaraccumulate : public testclass
{
   // A sollte der genauere Typ sein
   public:
      scalaraccumulate(void)
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


template <class D>
class cast_to_dotprecision : public testclass
{
   public:
      cast_to_dotprecision(void)
      {
         dotprecision a;
         D    d;

         testing(nameof(a),"dotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=dotprecision(d);
            if(a!=12 || dotprecision(d)!=12)
               error();
            else 
               ok();
         }
         
         testing(nameof(a),"_dotprecision("+nameof(d)+")");
         if(!tested())
         {
            a=434;
            d=12;
            a=_dotprecision(d);
            if(a!=12 || _dotprecision(d)!=12)
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
class test_dotio : public testclass
{
   public:                
      test_dotio(void)
      {
         D a;
         
         testing(nameof(a),"operator >>(string,"+nameof(a)+")");
         if(!tested())
         {
            int nok=0;
            
            string("1234.5678") >> a;
            if(abs(a-1234.5678)<1)
               nok++;
            
            string("-0.01234") >> a;
            if(abs(a+0.01234)<1e-5)
               nok++;
               
            string("2.34e+56") >> a;
            if(abs(a-2.34e+56)<1e+50)
               nok++;
               
            string("-2.34e-5") >> a;
            if(abs(a+2.34e-5)<1e-10)
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
            
            "1234.5678" >> a;
            
            if(abs(a-1234.5678)<1)
               nok++;
            
            (("-0.01234") >> a);
            
            if(abs(a+0.01234)<1e-5)
               nok++;
               
            (("2.34e+56") >> a);
            
            if(abs(a-2.34e+56)<1e+50)
               nok++;
               
            (("-2.34e-5") >> a);
            if(abs(a+2.34e-5)<1e-10)
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
            
            in >> a;
            if(abs(a-1234.5678)<1)
               nok++;
            
            in >> a;
            if(abs(a+0.01234)<1e-5)
               nok++;
               
            in >> a;
            if(abs(a-2.34e+56)<1e+40)
               nok++;
               
            in >> a;
            if(abs(a+2.34e-5)<1e-10)
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
            
            cout << SetDotPrecision(30,20);
            
            s=">";
            
            CSTRING_PI >> a;
            s << a;

            if(s==">        3.14159265358979323846")
               nok++;
            else
               cout << ">" << s << "<" << std::endl;
            
            "-1024" >> a;
            s="<";
            s << a;
            if(s=="<    -1024.00000000000000000000")
               nok++;
            else
               cout << ">" << s << "<" << std::endl;             
               
            "5.5e+78" >> a;
            s="";
            s << a;
            if(s=="  5.50000000000000000000E+0078")
               nok++;
            else
               cout << ">" << s << "<" << std::endl;               
               
            "5.5e-78" >> a;
            s="";
            s<<a;
            if(s=="  5.50000000000000000000E-0078")
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
            
            out << SetDotPrecision(30,20);
            
            CSTRING_PI >> a;
            out << a << std::endl;
            
            "-1024" >> a;
            out << a << std::endl;
            
            "5.5e+78" >> a;
            out << a << std::endl;
            
            "5.5e-78" >> a;
            out << a << std::endl;
            
            out.close();
            std::ifstream in("test.tmp");
            
            in >> s;
            
            if(s=="3.14159265358979323846")
               nok++;
            in >> s;
            
            if(s=="-1024.00000000000000000000")
               nok++;
            in >> s;
               
            if(s=="5.50000000000000000000E+0078")
               nok++;
            
            in >> s;
            
            if(s=="5.50000000000000000000E-0078")
               nok++;
               
            if(nok!=4)
               error();
            else 
               ok();
         }
         
      }
};

template <class D,class A>
class testrnd : public testclass
{
   // tests A rnd(D,RND_....)
   public:
      testrnd(void)
      {
         D d;
         A a;
         

         
         testing(nameof(a),"rnd("+nameof(d)+",RND_NEXT)");
         if(!tested())
         {
            a=10000;
            ("0.1")>>d;
            a=rnd(d,RND_NEXT);
            if(a==d || a-d>0.01 || d-a>0.01)
               error();
            else 
               ok();
         }
         testing(nameof(a),"rnd("+nameof(d)+",RND_UP)");
         if(!tested())
         {
            a=10000;
            ("0.1")>>d;
            a=rnd(d,RND_UP);
            if(a==d || a-d>0.01 || d-a>0.01 || a<d)
               error();
            else 
               ok();
         }
         testing(nameof(a),"rnd("+nameof(d)+",RND_DOWN)");
         if(!tested())
         {
            a=10000;
            ("0.1")>>d;
            a=rnd(d,RND_DOWN);
            if(a==d || a-d>0.01 || d-a>0.01 || a>d)
               error();
            else 
               ok();
         }
      }
};

} // namespace cxsc 

#endif // _CXSC_TESTDOT_HPP_INCLUDED
