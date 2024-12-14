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

/* CVS $Id: test.hpp,v 1.24 2014/01/30 17:23:49 cxsc Exp $ */

#define CXSC_INDEX_CHECK 1

#include "testclss.hpp"

#ifdef TEST_DOTPRECISION
   #include "real.hpp"
   #include "l_real.hpp"
   #include "dot.hpp"
   #include "testsklr.hpp"
   #include "testdot.hpp"
   #include "testintv.hpp"
#endif
#ifdef TEST_IDOTPRECISION
   #include "interval.hpp"
   #include "l_interv.hpp"
   #include "idot.hpp"
   #include "testsklr.hpp"
   #include "testdot.hpp"
   #include "testintv.hpp"
#endif
#ifdef TEST_REAL
   #include "real.hpp"
   #include "rmath.hpp"
   #include "ioflags.hpp"   
   
   #include "testsklr.hpp"
#endif
#ifdef TEST_L_REAL
   #include "l_real.hpp"
   #include "l_rmath.hpp"
   
   #include "testsklr.hpp"
   #include "testintv.hpp"
   #include "testdot.hpp"
#endif
#ifdef TEST_INTERVAL
   #include "interval.hpp"
   #include "imath.hpp"
   
   #include "testsklr.hpp"
   #include "testintv.hpp"
#endif
#ifdef TEST_L_INTERVAL
   #include "l_interv.hpp"
   #include "l_imath.hpp"
   
   #include "testsklr.hpp"
   #include "testintv.hpp"
#endif
#ifdef TEST_COMPLEX
   #include "complex.hpp"
   #include "cdot.hpp"

   #include "testsklr.hpp"
   #include "testcomp.hpp"
#endif
#ifdef TEST_CDOTPRECISION
   #include "complex.hpp"
   #include "cdot.hpp"

   #include "testsklr.hpp"
   #include "testcomp.hpp"
#endif
#ifdef TEST_CINTERVAL
   #include "cinterva.hpp"
   #include "cidot.hpp"

   #include "testsklr.hpp"
   #include "testcomp.hpp"
#endif
#ifdef TEST_CDOTPRECISION
   #include "complex.hpp"
   #include "cdot.hpp"

   #include "testsklr.hpp"
   #include "testcomp.hpp"
#endif
#ifdef TEST_RVECTOR
	#include "dot.hpp"
	#include "idot.hpp"
	#include "cdot.hpp"
	#include "cidot.hpp"
   #include "rvector.hpp"
   #include "rmatrix.hpp"
   
   #include "testvect.hpp"
#endif
#ifdef TEST_RMATRIX
   #include "rvector.hpp"
   #include "rmatrix.hpp"
   
   #include "testmatr.hpp"
#endif
#ifdef TEST_IVECTOR
	#include "idot.hpp"
	#include "cidot.hpp"
   #include "ivector.hpp"
   #include "imatrix.hpp"
   
   #include "testvect.hpp"
#endif
#ifdef TEST_IMATRIX
   #include "ivector.hpp"
   #include "imatrix.hpp"
   
   #include "testmatr.hpp"
#endif
#ifdef TEST_CVECTOR
	#include "cdot.hpp"
	#include "cidot.hpp"
   #include "cvector.hpp"
   #include "cmatrix.hpp"
   
   #include "testvect.hpp"
#endif
#ifdef TEST_CMATRIX
   #include "cvector.hpp"
   #include "cmatrix.hpp"
   
   #include "testmatr.hpp"
#endif
#ifdef TEST_CIVECTOR
	#include "cidot.hpp"
   #include "civector.hpp"
   #include "cimatrix.hpp"
   
   #include "testvect.hpp"
#endif
#ifdef TEST_CIMATRIX
   #include "civector.hpp"
   #include "cimatrix.hpp"
   
   #include "testmatr.hpp"
#endif
#ifdef TEST_LRVECTOR
	#include "lrvector.hpp"
	#include "testvect.hpp"
#endif
#ifdef TEST_LRMATRIX
	#include "lrmatrix.hpp"
	#include "testmatr.hpp"
#endif
#ifdef TEST_LIVECTOR
	#include "livector.hpp"
	#include "testvect.hpp"
#endif
#ifdef TEST_LIMATRIX
	#include "limatrix.hpp"
	#include "testmatr.hpp"
#endif

namespace cxsc {

template <class T>
class test : public testclass
{
   public:
      test(void)
      {
         T unused; 
         testing(nameof(unused),"Testfunctions not yet implemented");
         tested();
      }
};

#ifdef TEST_REAL
template <>
class test<real> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<real>(void)
      {
         cast_to_double<real>  a0;
         allscalarassign<real> a;
         cast_to_real<int>     b0;
         cast_to_real<long>    b1;
         cast_to_real<float>   b2;
         cast_to_real<double>  b3;
         allscalaraddsub<real> c;
         allscalarcompare<real> d;
         allscalarmuldiv<real> e;
         scalarabssgn<real>     f;
         scalarstdfunc<real>    g;
         test_realio            realio;
         testpredsucc<real>     predsucc;
         testminmax<real>       minmax;         
         testaddsubmultdivupdown<real> addsubmuldivupdown; 
         setfail(!a0 || !a || !b0 || !b1 || !b2 || !b3 
                     || !c || !d  || !e  || !f  || !g
                     || !realio || !predsucc || !addsubmuldivupdown
                     || !minmax);
      }
};
#endif
#ifdef TEST_L_REAL
template <>
class test<l_real> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<l_real>(void)
      {
         cast_to_l_real<int>           a0;
         cast_to_l_real<long>          a1;
         cast_to_l_real<float>         a2;
         cast_to_l_real<double>        a3;
         cast_to_l_real<real>          a4;
#ifdef TEST_DOTPRECISION
         cast_to_l_real<dotprecision>  a5;
         cast_to_dotprecision<l_real>  a6;
         scalaraddsub<dotprecision,l_real> a7;
#endif
         allscalarassign<l_real> a;
         allscalaraddsub<l_real> c;
         scalaraddsub<l_real,l_real> c2;
         allscalarcompare<l_real> d;
         allscalarmuldiv<l_real>  e;
         scalarmuldiv<l_real,l_real> c3;
         scalarabssgn<l_real>     f;
         scalarstdfunc<l_real>    g;
         scalarmixsetops<l_real,real>          k1;
         scalarmixsetops<l_real,int>          k3;
         scalarmixsetops<l_real,long>          k4;
         scalarmixsetops<l_real,float>          k5;
         scalarmixsetops<l_real,double>          k6;

         int stagprecsave=stagprec;
         
         stagprec=8;

         test_dotio<l_real>       k2;

         stagprec=stagprecsave;

         setfail(!a0 || !a1 || !a2 || !a3 || !a4 || !k1 || !k2
                     || !k3 || !k4 || !k5 || !k6 || !c2 || !c3
#ifdef TEST_DOTPRECISION
                     || !a5 || !a6 || !a7
#endif
                     || !a  || !c  || !d  || !e  || !f  || !g);
      }
};
#endif
#ifdef TEST_DOTPRECISION
template <>
class test<dotprecision> : public testclass
{
   // should test everything a dotprecision should be able to do
   public:
      test<dotprecision>(void)
      {
         cast_to_dotprecision<real>    a0;
         cast_to_dotprecision<int>     a1;
         cast_to_dotprecision<long>    a2;
         cast_to_dotprecision<float>   a3;
         cast_to_dotprecision<double>  a4;
         
         allscalarassign<dotprecision> a;
         allscalaraddsub<dotprecision> c;
         scalaraddsub<dotprecision,dotprecision> c2;
         allscalarcompare<dotprecision> d;
         
         scalarmixsetops<dotprecision,real>          k1;
         
         scalarabssgn<dotprecision>     e;
         scalaraccumulate<dotprecision,real,real> f;
         scalaraccumulate<dotprecision,real,l_real> g;
         scalaraccumulate<dotprecision,l_real,real> h;
         scalaraccumulate<dotprecision,l_real,l_real> i;
         testrnd<dotprecision,real>     j;
         test_dotio<dotprecision>       k;
         setfail(!a0 || !c2 || !a1 || !a2 || !a4 || !a3 || !k
              || !a || !c || !d || !e || !f || !g || !h || !i || !j);
      }
};
#endif
#ifdef TEST_INTERVAL
template <>
class test<interval> : public testclass
{
   // should test everything a interval should be able to do
   public:
      test<interval>(void)
      {
         interval rtstest(282429536481.);
         rtstest*=rtstest;
         if(Sup(rtstest)==Inf(rtstest))
            cerr << "The Run-Time-System has to be recompiled." << endl
                 << "The roundings do not work properly!" << endl;
      
         cast_scalar_to_interval<int,long> a0;
         cast_scalar_to_interval<int,float> a1;
         cast_scalar_to_interval<int,double> a2;
         cast_scalar_to_interval<long,float> a3;
         cast_scalar_to_interval<long,double> a4;
         cast_scalar_to_interval<float,double> a5;
         cast_scalar_to_interval<double,double> a6;
         cast_scalar_to_interval<real,real> a7;
         
#ifdef TEST_L_REAL
//         cast_scalar_to_interval<l_real,l_real> a8;
#endif
#ifdef TEST_DOTPRECISION
//         cast_scalar_to_interval<dotprecision,dotprecision> a9;
#endif
         
         intervalallscalarassign<interval> scalarassign;
         
         testnullconstructor<interval> nullconstructor;
      
         allscalarassign<interval> a; // Punktintervalle
         // intervalassign... gibts noch nicht
         intervaladdsub<interval,interval> b;
         allscalaraddsub<interval>         c;
         intervalmuldiv<interval,interval> d;
         allscalarmuldiv<interval>         e;
         intervalsetops<interval,interval> f;

         scalarmixsetops<real,int>            e0;
         scalarmixsetops<real,long>           e1;
         scalarmixsetops<real,float>          e2;
         scalarmixsetops<real,double>         e3;
         scalarmixsetops<int,real>            e4;
         scalarmixsetops<long,real>           e5;
         scalarmixsetops<float,real>          e6;
         scalarmixsetops<double,real>         e7;
         scalarmixsetops<real,real>           e8; 
         
         intervalsetcompops<interval,interval> f0;
         intervalscalarsetcompops<interval,real> f1;
         intervalscalarsetcompops<interval,int> f2;
         intervalscalarsetcompops<interval,long> f3;
         intervalscalarsetcompops<interval,float> f4;
         intervalscalarsetcompops<interval,double> f5;
#ifdef TEST_DOTPRECISION
         intervalscalarsetcompops<interval,dotprecision> f6;
#endif
#ifdef TEST_L_REAL
         intervalscalarsetcompops<interval,l_real> f7;
#endif
         intervalallscalarmixsetops<interval>  g;
         intervalstdfunc<interval>         h;
         testabsmiddiam<interval,real>     i;
         test_unchecked<interval,real>     j;
         setfail( !a0 || !a1 || !a2 || !a3 || !a4 || !a5 || !a6 || !a7 
               || !a  || !b  || !c  || !d  || !e  || !nullconstructor
               || !f  || !f0 || !f1 || !g  || !h  || !scalarassign
               || !i  || !f2 || !f3 || !f4 || !f5 || !j
               || !e0 || !e1 || !e2 || !e3 || !e4 || !e5 || !e6 || !e7 || !e8
#ifdef TEST_DOTPRECISION
               || !f6 
//               || !a9
#endif               
#ifdef TEST_L_REAL
//               || !a8 
               || !f7
#endif
               );
      }
};
#endif
#ifdef TEST_L_INTERVAL
template <>
class test<l_interval> : public testclass
{
   // should test everything a l_interval should be able to do
   public:
      test<l_interval>(void)
      {
         cast_scalar_to_l_interval<int> a0;
         cast_scalar_to_l_interval<long> a1;
         cast_scalar_to_l_interval<float> a2;
         cast_scalar_to_l_interval<double> a3;
         cast_scalar_to_l_interval<real> a4;

         cast_to_l_interval<interval> a5;
         cast_to_l_interval<idotprecision> a6;
         
         cast_to_interval<l_interval> a7;
#ifdef TEST_IDOTPRECISION
         cast_to_idotprecision<l_interval> a8;
         cast_to_l_interval<idotprecision> a9;
#endif      
         cast_scalar_to_l_interval<dotprecision>    a10;
         cast_scalar_to_l_interval<l_real>          a11;

         allscalarassign<l_interval> a; // Punktintervalle
         intervaladdsub<l_interval,l_interval> b;
         intervaladdsub<l_interval,interval>   c;
         allscalaraddsub<l_interval>           d;
         intervalmuldiv<l_interval,l_interval> e;
         intervalmuldiv<l_interval,interval>   f;
         allscalarmuldiv<l_interval>           g;
         intervalsetops<l_interval,l_interval> h;
         intervalsetops<l_interval,interval>   i;
         intervalallscalarmixsetops<l_interval>  j;
         intervalmixsetops<l_interval,l_real>  k;
         
         intervalsetcompops<l_interval,l_interval> f0;
         intervalsetcompops<l_interval,interval>   f1;
         intervalscalarsetcompops<l_interval,real> f2;
         
         intervalstdfunc<l_interval>           l;
         testabsmiddiam<l_interval,l_real>     m;
         setfail( !a0 || !a1 || !a2 || !a3 || !a4
               || !a5 || !a6 || !a10 || !a11
#ifdef TEST_IDOTPRECISION
               || !a8 || !a9
#endif
               || !f0 || !f1 || !f2
               || !a  || !b  || !c  || !d  || !e  || !f || !g || !h || !i || !j || !k || !l || !m);
      }
};
#endif
#ifdef TEST_IDOTPRECISION
template <>
class test<idotprecision> : public testclass
{
   // should test everything a idotprecision should be able to do
   public:
      test<idotprecision>(void)
      {
         testnullconstructor<idotprecision> nullconstructor;
      
         allscalarassign<idotprecision> a0;
         intervalscalarassign<idotprecision,dotprecision> a1;
         
         cast_scalar_to_idotprecision<int,long> b6;
         cast_scalar_to_idotprecision<int,float> b7;
         cast_scalar_to_idotprecision<int,double> b8;
         cast_scalar_to_idotprecision<long,float> b9;
         cast_scalar_to_idotprecision<long,double> b10;
         cast_scalar_to_idotprecision<float,double> b11;
         cast_scalar_to_idotprecision<double,double> b12;
         cast_scalar_to_idotprecision<real,real> b0;
         cast_scalar_to_idotprecision<dotprecision,dotprecision> b13;
         cast_scalar_to_idotprecision<l_real,l_real> b14;
         
         
         intervalsetops<idotprecision,idotprecision> a;
         intervalmixsetops<idotprecision,dotprecision> b;
         intervalmixsetops<idotprecision,int>          b1;
         intervalmixsetops<idotprecision,long>         b2;
//         intervalmixsetops<idotprecision,float>        b3;
//         intervalmixsetops<idotprecision,double>       b4;
         intervalmixsetops<idotprecision,real>         b5;
         
         
         intervalsetcompops<idotprecision,idotprecision> f0;
         
         intervalscalarsetcompops<idotprecision,int> f11;
         intervalscalarsetcompops<idotprecision,long> f12;
         intervalscalarsetcompops<idotprecision,float> f13;
         intervalscalarsetcompops<idotprecision,double> f14;
         
         intervalscalarsetcompops<idotprecision,real> f1;
         intervalscalarsetcompops<idotprecision,dotprecision> f2;
         intervalscalarsetcompops<idotprecision,l_real> f3;
         intervalsetcompops<idotprecision,interval>   f4;
         
         
         intervaladdsub<idotprecision,idotprecision> c;
         intervaladdsub<idotprecision,interval>      c1;
         intervaladdsub<idotprecision,l_interval>    c2;
         
         allscalaraddsub<idotprecision>              c3;       
         scalaraddsub<idotprecision,dotprecision>    c4;  
//         scalaraddsub<idotprecision,l_real>          c5;
         
         scalaraccumulate<idotprecision,real,real> l;
         scalaraccumulate<idotprecision,real,l_real> m;
         scalaraccumulate<idotprecision,l_real,real> n;
         scalaraccumulate<idotprecision,l_real,l_real> o;
         
         intervalaccumulate<idotprecision,interval,interval> d;
         intervalaccumulate<idotprecision,interval,l_interval> e;
         intervalaccumulate<idotprecision,l_interval,interval> f;
         intervalaccumulate<idotprecision,l_interval,l_interval> ff;
         intervalmixaccumulate<idotprecision,interval,real> g;
         intervalmixaccumulate<idotprecision,l_interval,real> h;
         intervalmixaccumulate<idotprecision,l_interval,l_real> i;
  
         test_unchecked<idotprecision,dotprecision> j;
//         testabsmiddiam<idotprecision,dotprecision> k;
         
         setfail(!a0 || !c1 || !c2 || !a1
              || !b0 || !b1 || !b2 || !b5
              || !b6 || !b7 || !b8 || !b9 || !b10 || !b11
              || !b12 || !b13 || !c3 || !c4 
//              || !c5
              || !nullconstructor
              || !f0 || !f1 || !f2 || !f3 || !ff || !f11 || !f12 || !f13 || !f14
              || !a  || !b || !c || !e || !f || !g || !h || !i || !j 
              || !l  || !m  || !n || !o
//              || !k
              );
      }
};
#endif

#ifdef TEST_COMPLEX
template <>
class test<complex> : public testclass
{
   // should test everything a complex should be able to do
   public:
      test<complex>(void)
      {
         allscalarassign<complex> a;
         scalarassign<complex,complex> b;
         cast_c_to_complex<cdotprecision> b2;
         cast_r_to_complex<int>     b3;
         cast_r_to_complex<long>    b4;
         cast_r_to_complex<float>   b5;
         cast_r_to_complex<double>  b6;
         cast_r_to_complex<real>    b7;
         
         cast_scalar_to_complex<int,long> b8;
         cast_scalar_to_complex<int,float> b9;
         cast_scalar_to_complex<int,double> b10;
         cast_scalar_to_complex<long,float> b11;
         cast_scalar_to_complex<long,double> b12;
         cast_scalar_to_complex<float,double> b13;
         cast_scalar_to_complex<double,double> b14;
         cast_scalar_to_complex<real,real> b15;
         
         testnullconstructor<complex> nullconstructor;
         
         allscalaraddsub<complex> c;
         complexaddsub<complex,complex> d;
         allscalar_eq_compare<complex> d2;
         scalar_eq_compare<complex,real> d3;
         complexcompare<complex,complex> e;
         complexcompare<complex,cdotprecision> e3;
         complexcompare<cdotprecision,complex> e4;
         allscalarmuldiv<complex> e2;
         complexmuldiv<complex,complex> f;
         complexabs<complex,real> g;
         complexreimconj<complex,real> h;
   //      test_realio            realio;
         setfail( !a || !b || !b2 || !b3 || !b4 || !b5 || !b6 || !b7
                     || !c || !d || !e || !e2 || !f || !g
                     || !d2  || !d3 || !e3 || !e4 || !h
                     || !b8 || !b9 || !b10 || !b11 || !b12 
                     || !b13 || !b14 || !b15 || !nullconstructor
                );
      }
};
#endif

#ifdef TEST_CDOTPRECISION
template <>
class test<cdotprecision> : public testclass
{
   // should test everything a cdotprecision should be able to do
   public:
      test<cdotprecision>(void)
      {
         allscalarassign<cdotprecision> a;
         scalarassign<cdotprecision,cdotprecision> b;
         // cast_c_to_cdotprecision<cdotprecision> b2;
         cast_r_to_cdotprecision<int>     b3;
         cast_r_to_cdotprecision<long>    b4;
         cast_r_to_cdotprecision<float>   b5;
         cast_r_to_cdotprecision<double>  b6;
         cast_r_to_cdotprecision<real>    b7;
         cast_c_to_cdotprecision<complex> b8;
         
         testnullconstructor<cdotprecision> nullconstructor;
         
         cast_scalar_to_cdotprecision<int,long> b17;
         cast_scalar_to_cdotprecision<int,float> b9;
         cast_scalar_to_cdotprecision<int,double> b10;
         cast_scalar_to_cdotprecision<long,float> b11;
         cast_scalar_to_cdotprecision<long,double> b12;
         cast_scalar_to_cdotprecision<float,double> b13;
         cast_scalar_to_cdotprecision<double,double> b14;
         cast_scalar_to_cdotprecision<real,real> b15;
         cast_scalar_to_cdotprecision<dotprecision,dotprecision> b16;
         
         
         allscalaraddsub<cdotprecision> c;
         scalaraddsub<cdotprecision,dotprecision> c2;
         
         complexaddsub<cdotprecision,cdotprecision> d;
         complexaddsub<cdotprecision,complex> d4;
         allscalar_eq_compare<cdotprecision> d2;
         scalar_eq_compare<cdotprecision,real> d3;
         complexcompare<cdotprecision,cdotprecision> e;
         complexcompare<complex,cdotprecision> e3;
         complexcompare<cdotprecision,complex> e4;
         complexreimconj<cdotprecision,dotprecision> h;
         
         complexaccumulate<cdotprecision,complex,complex> i;
         complexmixaccumulate<cdotprecision,complex,real> j;

         
   //      test_realio            realio;
         setfail( !a || !b || !b3 || !b4 || !b5 || !b6 || !b7
                     || !c || !c2 || !d || !e || !d4
                     || !d2  || !d3 || !e3 || !e4 || !h 
                     || !b8 || !b9 || !b10 || !b11 || !b12 || !b13
                     || !b14 || !b15 || !b16 || !nullconstructor
                     || !i || !j
                );
      }
};
#endif

#ifdef TEST_CINTERVAL
template <>
class test<cinterval> : public testclass
{
   // should test everything a cinterval should be able to do
   public:
      test<cinterval>(void)
      {
         allscalarassign<cinterval> a;
         scalarassign<cinterval,cinterval> b;
         cast_c_to_cinterval<complex> b2;
         cast_c_to_cinterval<cdotprecision> b3;        
         
/*         cast_scalar_to_cinterval<int,long> b8;
         cast_scalar_to_cinterval<int,float> b9;
         cast_scalar_to_cinterval<int,double> b10;
         cast_scalar_to_cinterval<long,float> b11;
         cast_scalar_to_cinterval<long,double> b12;
         cast_scalar_to_cinterval<float,double> b13;
         cast_scalar_to_cinterval<double,double> b14;
         cast_scalar_to_cinterval<real,real> b15;
*/         
         testnullconstructor<cinterval> nullconstructor;
         
         allscalaraddsub<cinterval> c;
//         cintervaladdsub<cinterval,cinterval> d;
         allscalar_eq_compare<cinterval> d2;
         scalar_eq_compare<cinterval,real> d3;
//         cintervalcompare<cinterval,cinterval> e;
//         cintervalcompare<cinterval,cdotprecision> e3;
//         cintervalcompare<cdotprecision,cinterval> e4;
         allscalarmuldiv<cinterval> e2;
//         cintervalmuldiv<cinterval,cinterval> f;
//         cintervalabs<cinterval,real> g;
//         cintervalreimconj<cinterval,real> h;
   //      test_realio            realio;
         setfail( !a || !b || !b2 || !b3 
                     || !c || !e2 
                     || !d2  || !d3
                     || !nullconstructor
                );
      }
};
#endif

#ifdef TEST_CIDOTPRECISION
template <>
class test<cidotprecision> : public testclass
{
   // should test everything a cidotprecision should be able to do
   public:
      test<cidotprecision>(void)
      {
         allscalarassign<cidotprecision> a;
         scalarassign<cidotprecision,cidotprecision> b;
         // cast_c_to_cidotprecision<cidotprecision> b2;
         cast_r_to_cidotprecision<int>     b3;
         cast_r_to_cidotprecision<long>    b4;
         cast_r_to_cidotprecision<float>   b5;
         cast_r_to_cidotprecision<double>  b6;
         cast_r_to_cidotprecision<real>    b7;
         cast_c_to_cidotprecision<complex> b8;
         
         testnullconstructor<cidotprecision> nullconstructor;
         
         allscalaraddsub<cidotprecision> c;
         scalaraddsub<cidotprecision,dotprecision> c2;
         
         complexaddsub<cidotprecision,cidotprecision> d;
         complexaddsub<cidotprecision,complex> d4;
         allscalar_eq_compare<cidotprecision> d2;
         scalar_eq_compare<cidotprecision,real> d3;
         complexcompare<cidotprecision,cidotprecision> e;
         complexcompare<complex,cidotprecision> e3;
         complexcompare<cidotprecision,complex> e4;
         complexreimconj<cidotprecision,dotprecision> h;
         
         complexaccumulate<cidotprecision,complex,complex> i;
         complexmixaccumulate<cidotprecision,complex,real> j;

         
   //      test_realio            realio;
         setfail( !a || !b || !b3 || !b4 || !b5 || !b6 || !b7
                     || !c || !c2 || !d || !e || !d4
                     || !d2  || !d3 || !e3 || !e4 || !h 
                     || !b8 
                     || !nullconstructor
                     || !i || !j
                );
      }
};
#endif



#ifdef TEST_RVECTOR
template <>
class test<rvector> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<rvector>(void)
      {
			rmatrix m(50,50),n(50,50);
			rvector v(50),w(50);
			
			vectorconstr<rvector,rvector,rmatrix,real,50> c(v,w,m);
			vectorconstr<rvector,rmatrix_subv,rmatrix,real,50> a26(v,m[9],n);
			
         vectorassign<rvector,rvector,rmatrix,50,real> a(v,w,m);
         vectorassign<rvector_slice,rvector,rmatrix,50,real> d(v(50),w,m);
         vectorassign<rmatrix_subv,rvector,rmatrix,50,real> b(m[4],w,n);
			
         vectoraddsub<rvector,rvector,50,rvector> e(v,w);
			vectoraddsub<rvector,rmatrix_subv,50,rvector> g(v,m[7]);
			vectoraddsub<rmatrix_subv,rvector,50,rvector> h(m[8],v);
			vectoraddsub<rmatrix_subv,rmatrix_subv,50,rvector> i(m[9],m[10]);
			vectoraddsub<rvector_slice,rvector,50,rvector> l(v(50),w);
			vectoraddsub<rvector,rvector_slice,50,rvector> a1(v,w(50));
			vectoraddsub<rvector_slice,rvector_slice,50,rvector> a2(v(50),w(50));
			vectoraddsub<rvector_slice,rmatrix_subv,50,rvector> a3(v(50),m[11]);
			vectoraddsub<rmatrix_subv,rvector_slice,50,rvector> a4(m[12],w(50));

         vectoraddsubassign<rvector,rvector,50,rvector> a27(v,w);
			vectoraddsubassign<rvector,rmatrix_subv,50,rvector> a28(v,m[7]);
			vectoraddsubassign<rmatrix_subv,rvector,50,rvector> a29(m[8],v);
			vectoraddsubassign<rmatrix_subv,rmatrix_subv,50,rvector> a30(m[9],m[10]);
			vectoraddsubassign<rvector_slice,rvector,50,rvector> a31(v(50),w);
			vectoraddsubassign<rvector,rvector_slice,50,rvector> a32(v,w(50));
			vectoraddsubassign<rvector_slice,rvector_slice,50,rvector> a33(v(50),w(50));
			vectoraddsubassign<rvector_slice,rmatrix_subv,50,rvector> a34(v(50),m[11]);
			vectoraddsubassign<rmatrix_subv,rvector_slice,50,rvector> a35(m[12],w(50));

			vectorscalar<rvector,real,50> a20(v);
			vectorscalar<rvector_slice,real,50> a21(v(50));
			vectorscalar<rmatrix_subv,real,50> a22(m[13]);

			vectorscalarassign<rvector,real,rvector,50> a23(v,w);
			vectorscalarassign<rvector_slice,real,rvector,50> a24(v(50),w);
			vectorscalarassign<rmatrix_subv,real,rvector,50> a25(m[14],w);

			vectormult<rvector,rmatrix_subv,50,real> j(v,m[11]);
			vectormult<rmatrix_subv,rmatrix_subv,50,real> k(m[2],m[3]);
			vectormult<rvector,rvector,50,real> f(v,w);
			vectormult<rvector,rvector_slice,50,real> a5(v,w(50));
			vectormult<rvector_slice,rvector_slice,50,real> a6(v(50),w(50));
			vectormult<rvector_slice,rmatrix_subv,50,real> a7(v(1,50),m[13]);
			
			vectoraccu<rvector,rvector,50,dotprecision> a8(v,w);
			vectoraccu<rvector,rmatrix_subv,50,dotprecision> a9(v,m[11]);
			vectoraccu<rmatrix_subv,rmatrix_subv,50,dotprecision> a10(m[2],m[3]);
			vectoraccu<rvector,rvector_slice,50,dotprecision> a11(v,w(50));
			vectoraccu<rvector_slice,rvector_slice,50,dotprecision> a12(v(50),w(50));
			vectoraccu<rvector_slice,rmatrix_subv,50,dotprecision> a13(v(1,50),m[13]);
			
			vectoraccu<rvector,rvector,50,idotprecision> a14(v,w);
			vectoraccu<rvector,rmatrix_subv,50,idotprecision> a15(v,m[11]);
			vectoraccu<rmatrix_subv,rmatrix_subv,50,idotprecision> a16(m[2],m[3]);
			vectoraccu<rvector,rvector_slice,50,idotprecision> a17(v,w(50));
			vectoraccu<rvector_slice,rvector_slice,50,idotprecision> a18(v(50),w(50));
			vectoraccu<rvector_slice,rmatrix_subv,50,idotprecision> a19(v(1,50),m[13]);
			
			vectoraccu<rvector,rvector,50,cdotprecision> a36(v,w);
			vectoraccu<rvector,rmatrix_subv,50,cdotprecision> a37(v,m[11]);
			vectoraccu<rmatrix_subv,rmatrix_subv,50,cdotprecision> a38(m[2],m[3]);
			vectoraccu<rvector,rvector_slice,50,cdotprecision> a47(v,w(50));
			vectoraccu<rvector_slice,rvector_slice,50,cdotprecision> a39(v(50),w(50));
			vectoraccu<rvector_slice,rmatrix_subv,50,cdotprecision> a40(v(1,50),m[13]);
			
			vectoraccu<rvector,rvector,50,cidotprecision> a41(v,w);
			vectoraccu<rvector,rmatrix_subv,50,cidotprecision> a42(v,m[11]);
			vectoraccu<rmatrix_subv,rmatrix_subv,50,cidotprecision> a43(m[2],m[3]);
			vectoraccu<rvector,rvector_slice,50,cidotprecision> a44(v,w(50));
			vectoraccu<rvector_slice,rvector_slice,50,cidotprecision> a45(v(50),w(50));
			vectoraccu<rvector_slice,rmatrix_subv,50,cidotprecision> a46(v(1,50),m[13]);
			
         setfail(!a || !b || !c || !d || !e || !f || !g || !h || !i||!j||!k || !l || !a1||!a2||!a3||!a4||!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a20||!a21||!a22||!a23||!a24||!a25||!a14||!a15||!a16||!a17||!a18||!a19||!a26||!a27||!a28||!a29||!a30||!a31||!a32||!a33||!a34||!a35||!a36||!a37||!a38||!a39||!a40||!a41||!a42||!a43||!a44||!a45||!a46||!a47);
      }
};
#endif
#ifdef TEST_RMATRIX
template <>
class test<rmatrix> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<rmatrix>(void)
      {
			rvector v(50);
			rmatrix m1(50,50),m2(50,50);

			matrixconstr<rmatrix,rvector,rmatrix,real,50> a1(m1,v,m2);
			matrixconstr<rmatrix,rmatrix_subv,rmatrix,real,50> a2(m1,m2[1],m2);
			
         matrixassign<rmatrix,rvector,rmatrix,real,50> a3(m1,v,m2);
         matrixassign<rmatrix_slice,rvector,rmatrix,real,50> a4(m1(50,50),v,m2);
         matrixassign<rmatrix,rvector_slice,rmatrix,real,50> a5(m1,v(50),m2);
         matrixassign<rmatrix_slice,rvector_slice,rmatrix,real,50> a6(m1(50,50),v(50),m2);
         matrixassign<rmatrix,rmatrix_subv,rmatrix,real,50> a7(m1,m2[Col(7)],m2);
         matrixassign<rmatrix_slice,rmatrix_subv,rmatrix,real,50> a8(m1(50,50),m2[Col(7)],m2);
         
			matrixaddsub<rmatrix,rmatrix,rmatrix,50> e(m1,m2);
			
			matrixscalarmult<rmatrix,real,rmatrix,50> f(m1);
			matrixscalarmult<rmatrix_slice,real,rmatrix,40> n(m1(1,40,1,40));
			
			matrixmult<rmatrix,rmatrix,rmatrix,dotprecision,50> o(m1,m2);
			matrixmult<rmatrix,rmatrix_slice,rmatrix,dotprecision,50> p(m1,m2(50,50));
			matrixmult<rmatrix_slice,rmatrix_slice,rmatrix,dotprecision,40> q(m1(1,40,1,40),m2(1,40,1,40));
         
			setfail( !e || !f || !n || !o|| !p || !q||!a1||!a2||!a3||!a4|!a5||!a6||!a7||!a8);
//         setfail(!e);
      }
};
#endif
#ifdef TEST_IVECTOR
template <>
class test<ivector> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<ivector>(void)
      {
			imatrix m(50,50),n(50,50);
			ivector v(50),w(50);
			rmatrix rm(50,50);
			rvector rv(50),rw(50);
			
			vectorconstr<ivector,ivector,imatrix,interval,50> c(v,w,m);
			vectorconstr<ivector,imatrix_subv,imatrix,interval,50> a26(v,m[9],n);
			
         vectorassign<ivector,ivector,imatrix,50,interval> a(v,w,m);
         vectorassign<ivector_slice,ivector,imatrix,50,interval> d(v(50),w,m);
         vectorassign<imatrix_subv,ivector,imatrix,50,interval> b(m[4],w,n);

         vectoraddsub<ivector,ivector,50,ivector> e(v,w);
			vectoraddsub<ivector,imatrix_subv,50,ivector> g(v,m[7]);
			vectoraddsub<imatrix_subv,ivector,50,ivector> h(m[8],v);
			vectoraddsub<imatrix_subv,imatrix_subv,50,ivector> i(m[9],m[10]);
			vectoraddsub<ivector_slice,ivector,50,ivector> l(v(50),w);
			vectoraddsub<ivector,ivector_slice,50,ivector> a1(v,w(50));
			vectoraddsub<ivector_slice,ivector_slice,50,ivector> a2(v(50),w(50));
			vectoraddsub<ivector_slice,imatrix_subv,50,ivector> a3(v(50),m[11]);
			vectoraddsub<imatrix_subv,ivector_slice,50,ivector> a4(m[12],w(50));

         vectoraddsubassign<ivector,ivector,50,ivector> c1(v,w);
			vectoraddsubassign<ivector,imatrix_subv,50,ivector> c2(v,m[7]);
			vectoraddsubassign<imatrix_subv,ivector,50,ivector> c3(m[8],v);
			vectoraddsubassign<imatrix_subv,imatrix_subv,50,ivector> c4(m[9],m[10]);
			vectoraddsubassign<ivector_slice,ivector,50,ivector> c5(v(50),w);
			vectoraddsubassign<ivector,ivector_slice,50,ivector> c6(v,w(50));
			vectoraddsubassign<ivector_slice,ivector_slice,50,ivector> c7(v(50),w(50));
			vectoraddsubassign<ivector_slice,imatrix_subv,50,ivector> c8(v(50),m[11]);
			vectoraddsubassign<imatrix_subv,ivector_slice,50,ivector> c9(m[12],w(50));

			vectorscalar<ivector,interval,50> a20(v);
			vectorscalar<ivector_slice,interval,50> a21(v(50));
			vectorscalar<imatrix_subv,interval,50> a22(m[13]);

			vectorscalarassign<ivector,interval,ivector,50> a23(v,w);
			vectorscalarassign<ivector_slice,interval,ivector,50> a24(v(50),w);
			vectorscalarassign<imatrix_subv,interval,ivector,50> a25(m[14],w);

			vectormult<ivector,imatrix_subv,50,interval> j(v,m[11]);
			vectormult<imatrix_subv,imatrix_subv,50,interval> k(m[2],m[3]);
			vectormult<ivector,ivector,50,interval> f(v,w);
			vectormult<ivector,ivector_slice,50,interval> a5(v,w(50));
			vectormult<ivector_slice,ivector_slice,50,interval> a6(v(50),w(50));
			vectormult<ivector_slice,imatrix_subv,50,interval> a7(v(1,50),m[13]);
			
			vectoraccu<ivector,ivector,50,idotprecision> a8(v,w);
			vectoraccu<ivector,imatrix_subv,50,idotprecision> a9(v,m[11]);
			vectoraccu<imatrix_subv,imatrix_subv,50,idotprecision> a10(m[2],m[3]);
			vectoraccu<ivector,ivector_slice,50,idotprecision> a11(v,w(50));
			vectoraccu<ivector_slice,ivector_slice,50,idotprecision> a12(v(50),w(50));
			vectoraccu<ivector_slice,imatrix_subv,50,idotprecision> a13(v(1,50),m[13]);
			
			vectoraccu<ivector,ivector,50,cidotprecision> a14(v,w);
			vectoraccu<ivector,imatrix_subv,50,cidotprecision> a15(v,m[11]);
			vectoraccu<imatrix_subv,imatrix_subv,50,cidotprecision> a16(m[2],m[3]);
			vectoraccu<ivector,ivector_slice,50,cidotprecision> a17(v,w(50));
			vectoraccu<ivector_slice,ivector_slice,50,cidotprecision> a18(v(50),w(50));
			vectoraccu<ivector_slice,imatrix_subv,50,cidotprecision> a19(v(1,50),m[13]);

         vectorconv<ivector,ivector,50,ivector> a63(v,w);
			vectorconv<imatrix_subv,ivector,50,ivector> a64(m[8],v);
			vectorconv<imatrix_subv,imatrix_subv,50,ivector> a65(m[9],m[10]);
			vectorconv<ivector_slice,ivector,50,ivector> a66(v(50),w);
			vectorconv<ivector_slice,ivector_slice,50,ivector> a67(v(50),w(50));
			vectorconv<ivector_slice,imatrix_subv,50,ivector> a68(v(50),m[11]);

         vectorconvassign<ivector,ivector,50,ivector> a69(v,w);
			vectorconvassign<imatrix_subv,ivector,50,ivector> a70(m[8],v);
			vectorconvassign<imatrix_subv,imatrix_subv,50,ivector> a71(m[9],m[10]);
			vectorconvassign<ivector_slice,ivector,50,ivector> a72(v(50),w);
			vectorconvassign<ivector_slice,ivector_slice,50,ivector> a73(v(50),w(50));
			vectorconvassign<ivector_slice,imatrix_subv,50,ivector> a74(v(50),m[11]);

         vectorsect<ivector,ivector,50,ivector,real> a75(v,w);
			vectorsect<imatrix_subv,ivector,50,ivector,real> a76(m[8],v);
			vectorsect<imatrix_subv,imatrix_subv,50,ivector,real> a77(m[9],m[10]);
			vectorsect<ivector_slice,ivector,50,ivector,real> a78(v(50),w);
			vectorsect<ivector_slice,ivector_slice,50,ivector,real> a79(v(50),w(50));
			vectorsect<ivector_slice,imatrix_subv,50,ivector,real> a80(v(50),m[11]);

         vectorsectassign<ivector,ivector,50,ivector,real> a81(v,w);
			vectorsectassign<imatrix_subv,ivector,50,ivector,real> a82(m[8],v);
			vectorsectassign<imatrix_subv,imatrix_subv,50,ivector,real> a83(m[9],m[10]);
			vectorsectassign<ivector_slice,ivector,50,ivector,real> a84(v(50),w);
			vectorsectassign<ivector_slice,ivector_slice,50,ivector,real> a85(v(50),w(50));
			vectorsectassign<ivector_slice,imatrix_subv,50,ivector,real> a86(v(50),m[11]);

		//-------- real ------------	
			vectorconstr<ivector,rvector,rmatrix,real,50> a27(v,rv,rm);
			vectorconstr<ivector,rmatrix_subv,rmatrix,real,50> a29(v,rm[3],rm);
			
			vectorassign<ivector,rvector,rmatrix,50,real> a28(v,rv,rm);
			vectorassign<ivector_slice,rvector,rmatrix,50,real> a30(v(50),rv,rm);
			vectorassign<imatrix_subv,rvector,rmatrix,50,real> a31(m[5],rv,rm);
			
         vectoraddsub<ivector,rvector,50,ivector> b1(v,rw);
         vectoraddsub<rvector,ivector,50,ivector> b2(rv,w);
			vectoraddsub<rvector,imatrix_subv,50,ivector> b3(rv,m[7]);
			vectoraddsub<ivector,rmatrix_subv,50,ivector> b4(v,rm[7]);
			vectoraddsub<rmatrix_subv,ivector,50,ivector> b5(rm[8],v);
			vectoraddsub<imatrix_subv,rvector,50,ivector> b6(m[8],rv);
			vectoraddsub<rmatrix_subv,imatrix_subv,50,ivector> b7(rm[9],m[10]);
			vectoraddsub<imatrix_subv,rmatrix_subv,50,ivector> b8(m[9],rm[10]);
			vectoraddsub<rvector_slice,ivector,50,ivector> b9(rv(50),w);
			vectoraddsub<ivector_slice,rvector,50,ivector> b10(v(50),rw);
			vectoraddsub<rvector,ivector_slice,50,ivector> b11(rv,w(50));
			vectoraddsub<ivector,rvector_slice,50,ivector> b12(v,rw(50));
			vectoraddsub<rvector_slice,ivector_slice,50,ivector> b13(rv(50),w(50));
			vectoraddsub<ivector_slice,rvector_slice,50,ivector> b14(v(50),rw(50));
			vectoraddsub<rvector_slice,imatrix_subv,50,ivector> b15(rv(50),m[11]);
			vectoraddsub<ivector_slice,rmatrix_subv,50,ivector> b16(v(50),rm[11]);
			vectoraddsub<rmatrix_subv,ivector_slice,50,ivector> b17(rm[12],w(50));
			vectoraddsub<imatrix_subv,rvector_slice,50,ivector> b18(m[12],rw(50));

         vectoraddsubassign<ivector,rvector,50,ivector> c10(v,rw);
			vectoraddsubassign<ivector,rmatrix_subv,50,ivector> c11(v,rm[7]);
			vectoraddsubassign<imatrix_subv,rvector,50,ivector> c12(m[8],rv);
			vectoraddsubassign<imatrix_subv,rmatrix_subv,50,ivector> c13(m[9],rm[10]);
			vectoraddsubassign<ivector_slice,rvector,50,ivector> c14(v(50),rw);
			vectoraddsubassign<ivector,rvector_slice,50,ivector> c15(v,rw(50));
			vectoraddsubassign<ivector_slice,rvector_slice,50,ivector> c16(v(50),rw(50));
			vectoraddsubassign<ivector_slice,rmatrix_subv,50,ivector> c17(v(50),rm[11]);
			vectoraddsubassign<imatrix_subv,rvector_slice,50,ivector> c18(m[12],rw(50));

			vectormult<rvector,imatrix_subv,50,interval> b19(rv,m[11]);
			vectormult<ivector,rmatrix_subv,50,interval> b20(v,rm[11]);
			vectormult<imatrix_subv,rvector,50,interval> c19(m[11],rv);
			vectormult<rmatrix_subv,ivector,50,interval> c20(rm[11],v);
			vectormult<rmatrix_subv,imatrix_subv,50,interval> b21(rm[2],m[3]);
			vectormult<imatrix_subv,rmatrix_subv,50,interval> b22(m[2],rm[3]);
			vectormult<rvector,ivector,50,interval> b23(rv,w);
			vectormult<ivector,rvector,50,interval> b24(v,rw);
			vectormult<rvector,ivector_slice,50,interval> b25(rv,w(50));
			vectormult<ivector,rvector_slice,50,interval> b26(v,rw(50));
			vectormult<rvector_slice,ivector,50,interval> c23(rv(50),w);
			vectormult<ivector_slice,rvector,50,interval> c24(v(50),rw);
			vectormult<rvector_slice,ivector_slice,50,interval> b27(rv(50),w(50));
			vectormult<ivector_slice,rvector_slice,50,interval> b28(v(50),rw(50));
			vectormult<rvector_slice,imatrix_subv,50,interval> b29(rv(1,50),m[13]);
			vectormult<ivector_slice,rmatrix_subv,50,interval> b30(v(1,50),rm[13]);
			vectormult<imatrix_subv,rvector_slice,50,interval> c22(m[11],rv(1,50));
			vectormult<rmatrix_subv,ivector_slice,50,interval> c21(rm[11],v(1,50));
			
			vectoraccu<ivector,rvector,50,idotprecision> a39(v,rw);
			vectoraccu<ivector,rmatrix_subv,50,idotprecision> a40(v,rm[11]);
			vectoraccu<imatrix_subv,rmatrix_subv,50,idotprecision> a41(m[2],rm[3]);
			vectoraccu<ivector,rvector_slice,50,idotprecision> a42(v,rw(50));
			vectoraccu<ivector_slice,rvector_slice,50,idotprecision> a43(v(50),rw(50));
			vectoraccu<ivector_slice,rmatrix_subv,50,idotprecision> a44(v(1,50),rm[13]);
			
			vectoraccu<rvector,ivector,50,idotprecision> a45(rv,w);
			vectoraccu<rvector,imatrix_subv,50,idotprecision> a46(rv,m[11]);
			vectoraccu<rmatrix_subv,imatrix_subv,50,idotprecision> a47(rm[2],m[3]);
			vectoraccu<rvector,ivector_slice,50,idotprecision> a48(rv,w(50));
			vectoraccu<rvector_slice,ivector_slice,50,idotprecision> a49(rv(50),w(50));
			vectoraccu<rvector_slice,imatrix_subv,50,idotprecision> a50(rv(1,50),m[13]);
			
			vectoraccu<ivector,rvector,50,cidotprecision> a51(v,rw);
			vectoraccu<ivector,rmatrix_subv,50,cidotprecision> a52(v,rm[11]);
			vectoraccu<imatrix_subv,rmatrix_subv,50,cidotprecision> a53(m[2],rm[3]);
			vectoraccu<ivector,rvector_slice,50,cidotprecision> a54(v,rw(50));
			vectoraccu<ivector_slice,rvector_slice,50,cidotprecision> a55(v(50),rw(50));
			vectoraccu<ivector_slice,rmatrix_subv,50,cidotprecision> a56(v(1,50),rm[13]);

			vectoraccu<rvector,ivector,50,cidotprecision> a57(rv,w);
			vectoraccu<rvector,imatrix_subv,50,cidotprecision> a58(rv,m[11]);
			vectoraccu<rmatrix_subv,imatrix_subv,50,cidotprecision> a59(rm[2],m[3]);
			vectoraccu<rvector,ivector_slice,50,cidotprecision> a60(rv,w(50));
			vectoraccu<rvector_slice,ivector_slice,50,cidotprecision> a61(rv(50),w(50));
			vectoraccu<rvector_slice,imatrix_subv,50,cidotprecision> a62(rv(1,50),m[13]);

         vectorconv<rvector,rvector,50,ivector> a111(rv,rw);
			vectorconv<rmatrix_subv,rvector,50,ivector> a112(rm[8],rv);
			vectorconv<rmatrix_subv,rmatrix_subv,50,ivector> a113(rm[9],rm[10]);
			vectorconv<rvector_slice,rvector,50,ivector> a114(rv(50),rw);
			vectorconv<rvector_slice,rvector_slice,50,ivector> a115(rv(50),rw(50));
			vectorconv<rvector_slice,rmatrix_subv,50,ivector> a116(rv(50),rm[11]);

         vectorconv<ivector,rvector,50,ivector> a87(v,rw);
			vectorconv<imatrix_subv,rvector,50,ivector> a88(m[8],rv);
			vectorconv<imatrix_subv,rmatrix_subv,50,ivector> a89(m[9],rm[10]);
			vectorconv<ivector_slice,rvector,50,ivector> a90(v(50),rw);
			vectorconv<ivector_slice,rvector_slice,50,ivector> a91(v(50),rw(50));
			vectorconv<ivector_slice,rmatrix_subv,50,ivector> a92(v(50),rm[11]);

         vectorconvassign<ivector,rvector,50,ivector> a93(v,rw);
			vectorconvassign<imatrix_subv,rvector,50,ivector> a94(m[8],rv);
			vectorconvassign<imatrix_subv,rmatrix_subv,50,ivector> a95(m[9],rm[10]);
			vectorconvassign<ivector_slice,rvector,50,ivector> a96(v(50),rw);
			vectorconvassign<ivector_slice,rvector_slice,50,ivector> a97(v(50),rw(50));
			vectorconvassign<ivector_slice,rmatrix_subv,50,ivector> a98(v(50),rm[11]);

         vectorsect<ivector,rvector,50,ivector,real> a99(v,rw);
			vectorsect<imatrix_subv,rvector,50,ivector,real> a100(m[8],rv);
			vectorsect<imatrix_subv,rmatrix_subv,50,ivector,real> a101(m[9],rm[10]);
			vectorsect<ivector_slice,rvector,50,ivector,real> a102(v(50),rw);
			vectorsect<ivector_slice,rvector_slice,50,ivector,real> a103(v(50),rw(50));
			vectorsect<ivector_slice,rmatrix_subv,50,ivector,real> a104(v(50),rm[11]);

         vectorsectassign<ivector,rvector,50,ivector,real> a105(v,rw);
			vectorsectassign<imatrix_subv,rvector,50,ivector,real> a106(m[8],rv);
			vectorsectassign<imatrix_subv,rmatrix_subv,50,ivector,real> a107(m[9],rm[10]);
			vectorsectassign<ivector_slice,rvector,50,ivector,real> a108(v(50),rw);
			vectorsectassign<ivector_slice,rvector_slice,50,ivector,real> a109(v(50),rw(50));
			vectorsectassign<ivector_slice,rmatrix_subv,50,ivector,real> a110(v(50),rm[11]);

         setfail(!a || !b || !c || !d || !e || !f || !g || !h || !i||!j||!k || !l || !a1||!a2||!a3||!a4||!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a20||!a21||!a22||!a23||!a24||!a25||!a14||!a15||!a16||!a17||!a18||!a19||!a26||!b1||!b2||!b3||!b4||!b5||!b6||!b7||!b8||!b9||!b10||!b11||!b12||!b13||!b14||!b15||!b16||!b17||!b18||!b19||!b20||!b21||!b22||!b23||!b24||!b25||!b26||!b27||!b28||!b29||!b30||!c1||!c2||!c3||!c4||!c5||!c6||!c7||!c8||!c9||!c10||!c11||!c12||!c13||!c14||!c15||!c16||!c17||!c18||!c19||!c20||!c21||!c22||!c23||!c24||!a27||!a28||!a29||!a30||!a31||!a39||!a40||!a41||!a42||!a43||!a44||!a45||!a46||!a47||!a48||!a49||!a50||!a51||!a52||!a53||!a54||!a55||!a56||!a57||!a58||!a59||!a60||!a61||!a62||!a63||!a64||!a65||!a66||!a67||!a68||!a69||!a70||!a71||!a72||!a73||!a74||!a75||!a76||!a77||!a78||!a79||!a80||!a81||!a82||!a83||!a84||!a85||!a86||!a87||!a88||!a89||!a90||!a91||!a92||!a93||!a94||!a95||!a96||!a97||!a98||!a99||!a100||!a101||!a102||!a103||!a104||!a105||!a106||!a107||!a108||!a109||!a110||!a111||!a112||!a113||!a114||!a115||!a116);
      }
};
#endif
#ifdef TEST_IMATRIX
template <>
class test<imatrix> : public testclass
{
   // should test everything a interval should be able to do
   public:
      test<imatrix>(void)
      {
			ivector v(50);
			imatrix m1(50,50),m2(50,50);
			rvector rv(50);
			rmatrix rm1(50,50),rm2(50,50);

			matrixconstr<imatrix,ivector,imatrix,interval,50> a1(m1,v,m2);
			matrixconstr<imatrix,imatrix_subv,imatrix,interval,50> a2(m1,m2[1],m2);
			
         matrixassign<imatrix,ivector,imatrix,interval,50> a3(m1,v,m2);
         matrixassign<imatrix_slice,ivector,imatrix,interval,50> a4(m1(50,50),v,m2);
         matrixassign<imatrix,ivector_slice,imatrix,interval,50> a5(m1,v(50),m2);
         matrixassign<imatrix_slice,ivector_slice,imatrix,interval,50> a6(m1(50,50),v(50),m2);
         matrixassign<imatrix,imatrix_subv,imatrix,interval,50> a7(m1,m2[Col(7)],m2);
         matrixassign<imatrix_slice,imatrix_subv,imatrix,interval,50> a8(m1(50,50),m2[Col(7)],m2);
         
			matrixaddsub<imatrix,imatrix,imatrix,50> e(m1,m2);
			
			matrixscalarmult<imatrix,interval,imatrix,50> f(m1);
			matrixscalarmult<imatrix_slice,interval,imatrix,40> n(m1(1,40,1,40));
			
			matrixmult<imatrix,imatrix,imatrix,idotprecision,50> o(m1,m2);
			matrixmult<imatrix,imatrix_slice,imatrix,idotprecision,50> p(m1,m2(50,50));
			matrixmult<imatrix_slice,imatrix_slice,imatrix,idotprecision,40> q(m1(1,40,1,40),m2(1,40,1,40));
			// Real ----------------------------------
			matrixconstr<imatrix,rvector,rmatrix,real,50> a9(m1,rv,rm2);
			matrixconstr<imatrix,rmatrix_subv,rmatrix,real,50> a10(m1,rm2[1],rm2);
         
         matrixassign<imatrix,rvector,rmatrix,real,50> a11(m1,rv,rm2);
         matrixassign<imatrix_slice,rvector,rmatrix,real,50> a12(m1(50,50),rv,rm2);
         matrixassign<imatrix,rvector_slice,rmatrix,real,50> a13(m1,rv(50),rm2);
         matrixassign<imatrix_slice,rvector_slice,rmatrix,real,50> a14(m1(50,50),rv(50),rm2);
         matrixassign<imatrix,rmatrix_subv,rmatrix,real,50> a15(m1,rm2[Col(7)],rm2);
         matrixassign<imatrix_slice,rmatrix_subv,rmatrix,real,50> a16(m1(50,50),rm2[Col(7)],rm2);
         
			setfail( !e || !f || !n || !o|| !p || !q||!a1||!a2||!a3||!a4|!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a14||!a15||!a16);
//         setfail(!e);
      }
};
#endif
#ifdef TEST_CVECTOR
template <>
class test<cvector> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<cvector>(void)
      {
			cmatrix m(50,50),n(50,50);
			rmatrix rm(50,50);
			rvector rv(50),rw(50);
			cvector v(50),w(50);
			
			vectorconstr<cvector,cvector,cmatrix,complex,50> c(v,w,m);
			vectorconstr<cvector,cmatrix_subv,cmatrix,complex,50> a26(v,m[9],n);
			
         vectorassign<cvector,cvector,cmatrix,50,complex> a(v,w,m);
         vectorassign<cvector_slice,cvector,cmatrix,50,complex> d(v(50),w,m);
         vectorassign<cmatrix_subv,cvector,cmatrix,50,complex> b(m[4],w,n);

         vectoraddsub<cvector,cvector,50,cvector> e(v,w);
			vectoraddsub<cvector,cmatrix_subv,50,cvector> g(v,m[7]);
			vectoraddsub<cmatrix_subv,cvector,50,cvector> h(m[8],v);
			vectoraddsub<cmatrix_subv,cmatrix_subv,50,cvector> i(m[9],m[10]);
			vectoraddsub<cvector_slice,cvector,50,cvector> l(v(50),w);
			vectoraddsub<cvector,cvector_slice,50,cvector> a1(v,w(50));
			vectoraddsub<cvector_slice,cvector_slice,50,cvector> a2(v(50),w(50));
			vectoraddsub<cvector_slice,cmatrix_subv,50,cvector> a3(v(50),m[11]);
			vectoraddsub<cmatrix_subv,cvector_slice,50,cvector> a4(m[12],w(50));

         vectoraddsubassign<cvector,cvector,50,cvector> c1(v,w);
			vectoraddsubassign<cvector,cmatrix_subv,50,cvector> c2(v,m[7]);
			vectoraddsubassign<cmatrix_subv,cvector,50,cvector> c3(m[8],v);
			vectoraddsubassign<cmatrix_subv,cmatrix_subv,50,cvector> c4(m[9],m[10]);
			vectoraddsubassign<cvector_slice,cvector,50,cvector> c5(v(50),w);
			vectoraddsubassign<cvector,cvector_slice,50,cvector> c6(v,w(50));
			vectoraddsubassign<cvector_slice,cvector_slice,50,cvector> c7(v(50),w(50));
			vectoraddsubassign<cvector_slice,cmatrix_subv,50,cvector> c8(v(50),m[11]);
			vectoraddsubassign<cmatrix_subv,cvector_slice,50,cvector> c9(m[12],w(50));

			vectorscalar<cvector,complex,50> a20(v);
			vectorscalar<cvector_slice,complex,50> a21(v(50));
			vectorscalar<cmatrix_subv,complex,50> a22(m[13]);

			vectorscalarassign<cvector,complex,cvector,50> a23(v,w);
			vectorscalarassign<cvector_slice,complex,cvector,50> a24(v(50),w);
			vectorscalarassign<cmatrix_subv,complex,cvector,50> a25(m[14],w);

			vectormult<cvector,cmatrix_subv,50,complex> j(v,m[11]);
			vectormult<cmatrix_subv,cmatrix_subv,50,complex> k(m[2],m[3]);
			vectormult<cvector,cvector,50,complex> f(v,w);
			vectormult<cvector,cvector_slice,50,complex> a5(v,w(50));
			vectormult<cvector_slice,cvector_slice,50,complex> a6(v(50),w(50));
			vectormult<cvector_slice,cmatrix_subv,50,complex> a7(v(1,50),m[13]);
			
			vectoraccu<cvector,cvector,50,cdotprecision> a8(v,w);
			vectoraccu<cvector,cmatrix_subv,50,cdotprecision> a9(v,m[11]);
			vectoraccu<cmatrix_subv,cmatrix_subv,50,cdotprecision> a10(m[2],m[3]);
			vectoraccu<cvector,cvector_slice,50,cdotprecision> a11(v,w(50));
			vectoraccu<cvector_slice,cvector_slice,50,cdotprecision> a12(v(50),w(50));
			vectoraccu<cvector_slice,cmatrix_subv,50,cdotprecision> a13(v(1,50),m[13]);
			
			vectoraccu<cvector,cvector,50,cidotprecision> a14(v,w);
			vectoraccu<cvector,cmatrix_subv,50,cidotprecision> a15(v,m[11]);
			vectoraccu<cmatrix_subv,cmatrix_subv,50,cidotprecision> a16(m[2],m[3]);
			vectoraccu<cvector,cvector_slice,50,cidotprecision> a17(v,w(50));
			vectoraccu<cvector_slice,cvector_slice,50,cidotprecision> a18(v(50),w(50));
			vectoraccu<cvector_slice,cmatrix_subv,50,cidotprecision> a19(v(1,50),m[13]);

		//-------- real ------------	
			vectorconstr<cvector,rvector,rmatrix,real,50> a27(v,rv,rm);
			vectorconstr<cvector,rmatrix_subv,rmatrix,real,50> a29(v,rm[3],rm);
			
			vectorassign<cvector,rvector,rmatrix,50,real> a28(v,rv,rm);
			vectorassign<cvector_slice,rvector,rmatrix,50,real> a30(v(50),rv,rm);
			vectorassign<cmatrix_subv,rvector,rmatrix,50,real> a31(m[5],rv,rm);
			
         vectoraddsub<cvector,rvector,50,cvector> b1(v,rw);
         vectoraddsub<rvector,cvector,50,cvector> b2(rv,w);
			vectoraddsub<rvector,cmatrix_subv,50,cvector> b3(rv,m[7]);
			vectoraddsub<cvector,rmatrix_subv,50,cvector> b4(v,rm[7]);
			vectoraddsub<rmatrix_subv,cvector,50,cvector> b5(rm[8],v);
			vectoraddsub<cmatrix_subv,rvector,50,cvector> b6(m[8],rv);
			vectoraddsub<rmatrix_subv,cmatrix_subv,50,cvector> b7(rm[9],m[10]);
			vectoraddsub<cmatrix_subv,rmatrix_subv,50,cvector> b8(m[9],rm[10]);
			vectoraddsub<rvector_slice,cvector,50,cvector> b9(rv(50),w);
			vectoraddsub<cvector_slice,rvector,50,cvector> b10(v(50),rw);
			vectoraddsub<rvector,cvector_slice,50,cvector> b11(rv,w(50));
			vectoraddsub<cvector,rvector_slice,50,cvector> b12(v,rw(50));
			vectoraddsub<rvector_slice,cvector_slice,50,cvector> b13(rv(50),w(50));
			vectoraddsub<cvector_slice,rvector_slice,50,cvector> b14(v(50),rw(50));
			vectoraddsub<rvector_slice,cmatrix_subv,50,cvector> b15(rv(50),m[11]);
			vectoraddsub<cvector_slice,rmatrix_subv,50,cvector> b16(v(50),rm[11]);
			vectoraddsub<rmatrix_subv,cvector_slice,50,cvector> b17(rm[12],w(50));
			vectoraddsub<cmatrix_subv,rvector_slice,50,cvector> b18(m[12],rw(50));

         vectoraddsubassign<cvector,rvector,50,cvector> c10(v,rw);
			vectoraddsubassign<cvector,rmatrix_subv,50,cvector> c11(v,rm[7]);
			vectoraddsubassign<cmatrix_subv,rvector,50,cvector> c12(m[8],rv);
			vectoraddsubassign<cmatrix_subv,rmatrix_subv,50,cvector> c13(m[9],rm[10]);
			vectoraddsubassign<cvector_slice,rvector,50,cvector> c14(v(50),rw);
			vectoraddsubassign<cvector,rvector_slice,50,cvector> c15(v,rw(50));
			vectoraddsubassign<cvector_slice,rvector_slice,50,cvector> c16(v(50),rw(50));
			vectoraddsubassign<cvector_slice,rmatrix_subv,50,cvector> c17(v(50),rm[11]);
			vectoraddsubassign<cmatrix_subv,rvector_slice,50,cvector> c18(m[12],rw(50));

			vectormult<rvector,cmatrix_subv,50,complex> b19(rv,m[11]);
			vectormult<cvector,rmatrix_subv,50,complex> b20(v,rm[11]);
			vectormult<cmatrix_subv,rvector,50,complex> c19(m[11],rv);
			vectormult<rmatrix_subv,cvector,50,complex> c20(rm[11],v);
			vectormult<rmatrix_subv,cmatrix_subv,50,complex> b21(rm[2],m[3]);
			vectormult<cmatrix_subv,rmatrix_subv,50,complex> b22(m[2],rm[3]);
			vectormult<rvector,cvector,50,complex> b23(rv,w);
			vectormult<cvector,rvector,50,complex> b24(v,rw);
			vectormult<rvector,cvector_slice,50,complex> b25(rv,w(50));
			vectormult<cvector,rvector_slice,50,complex> b26(v,rw(50));
			vectormult<rvector_slice,cvector,50,complex> c23(rv(50),w);
			vectormult<cvector_slice,rvector,50,complex> c24(v(50),rw);
			vectormult<rvector_slice,cvector_slice,50,complex> b27(rv(50),w(50));
			vectormult<cvector_slice,rvector_slice,50,complex> b28(v(50),rw(50));
			vectormult<rvector_slice,cmatrix_subv,50,complex> b29(rv(1,50),m[13]);
			vectormult<cvector_slice,rmatrix_subv,50,complex> b30(v(1,50),rm[13]);
			vectormult<cmatrix_subv,rvector_slice,50,complex> c22(m[11],rv(1,50));
			vectormult<rmatrix_subv,cvector_slice,50,complex> c21(rm[11],v(1,50));
			
			vectoraccu<cvector,rvector,50,cdotprecision> a39(v,rw);
			vectoraccu<cvector,rmatrix_subv,50,cdotprecision> a40(v,rm[11]);
			vectoraccu<cmatrix_subv,rmatrix_subv,50,cdotprecision> a41(m[2],rm[3]);
			vectoraccu<cvector,rvector_slice,50,cdotprecision> a42(v,rw(50));
			vectoraccu<cvector_slice,rvector_slice,50,cdotprecision> a43(v(50),rw(50));
			vectoraccu<cvector_slice,rmatrix_subv,50,cdotprecision> a44(v(1,50),rm[13]);
			
			vectoraccu<rvector,cvector,50,cdotprecision> a45(rv,w);
			vectoraccu<rvector,cmatrix_subv,50,cdotprecision> a46(rv,m[11]);
			vectoraccu<rmatrix_subv,cmatrix_subv,50,cdotprecision> a47(rm[2],m[3]);
			vectoraccu<rvector,cvector_slice,50,cdotprecision> a48(rv,w(50));
			vectoraccu<rvector_slice,cvector_slice,50,cdotprecision> a49(rv(50),w(50));
			vectoraccu<rvector_slice,cmatrix_subv,50,cdotprecision> a50(rv(1,50),m[13]);
			
			vectoraccu<cvector,rvector,50,cidotprecision> a51(v,rw);
			vectoraccu<cvector,rmatrix_subv,50,cidotprecision> a52(v,rm[11]);
			vectoraccu<cmatrix_subv,rmatrix_subv,50,cidotprecision> a53(m[2],rm[3]);
			vectoraccu<cvector,rvector_slice,50,cidotprecision> a54(v,rw(50));
			vectoraccu<cvector_slice,rvector_slice,50,cidotprecision> a55(v(50),rw(50));
			vectoraccu<cvector_slice,rmatrix_subv,50,cidotprecision> a56(v(1,50),rm[13]);

			vectoraccu<rvector,cvector,50,cidotprecision> a57(rv,w);
			vectoraccu<rvector,cmatrix_subv,50,cidotprecision> a58(rv,m[11]);
			vectoraccu<rmatrix_subv,cmatrix_subv,50,cidotprecision> a59(rm[2],m[3]);
			vectoraccu<rvector,cvector_slice,50,cidotprecision> a60(rv,w(50));
			vectoraccu<rvector_slice,cvector_slice,50,cidotprecision> a61(rv(50),w(50));
			vectoraccu<rvector_slice,cmatrix_subv,50,cidotprecision> a62(rv(1,50),m[13]);

         setfail(!a || !b || !c || !d || !e || !f || !g || !h || !i||!j||!k || !l || !a1||!a2||!a3||!a4||!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a20||!a21||!a22||!a23||!a24||!a25||!a14||!a15||!a16||!a17||!a18||!a19||!a26||!b1||!b2||!b3||!b4||!b5||!b6||!b7||!b8||!b9||!b10||!b11||!b12||!b13||!b14||!b15||!b16||!b17||!b18||!b19||!b20||!b21||!b22||!b23||!b24||!b25||!b26||!b27||!b28||!b29||!b30||!c1||!c2||!c3||!c4||!c5||!c6||!c7||!c8||!c9||!c10||!c11||!c12||!c13||!c14||!c15||!c16||!c17||!c18||!c19||!c20||!c21||!c22||!c23||!c24||!a27||!a28||!a29||!a30||!a31||!a39||!a40||!a41||!a42||!a43||!a44||!a45||!a46||!a47||!a48||!a49||!a50||!a51||!a52||!a53||!a54||!a55||!a56||!a57||!a58||!a59||!a60||!a61||!a62);
      }
};
#endif
#ifdef TEST_CMATRIX
template <>
class test<cmatrix> : public testclass
{
   // should test everything a complex should be able to do
   public:
      test<cmatrix>(void)
      {
			cvector v(50);
			cmatrix m1(50,50),m2(50,50);
			rvector rv(50);
			rmatrix rm1(50,50),rm2(50,50);

			matrixconstr<cmatrix,cvector,cmatrix,complex,50> a1(m1,v,m2);
			matrixconstr<cmatrix,cmatrix_subv,cmatrix,complex,50> a2(m1,m2[1],m2);
			
         matrixassign<cmatrix,cvector,cmatrix,complex,50> a3(m1,v,m2);
         matrixassign<cmatrix_slice,cvector,cmatrix,complex,50> a4(m1(50,50),v,m2);
         matrixassign<cmatrix,cvector_slice,cmatrix,complex,50> a5(m1,v(50),m2);
         matrixassign<cmatrix_slice,cvector_slice,cmatrix,complex,50> a6(m1(50,50),v(50),m2);
         matrixassign<cmatrix,cmatrix_subv,cmatrix,complex,50> a7(m1,m2[Col(7)],m2);
         matrixassign<cmatrix_slice,cmatrix_subv,cmatrix,complex,50> a8(m1(50,50),m2[Col(7)],m2);
         
			matrixaddsub<cmatrix,cmatrix,cmatrix,50> e(m1,m2);
			
			matrixscalarmult<cmatrix,complex,cmatrix,50> f(m1);
			matrixscalarmult<cmatrix_slice,complex,cmatrix,40> n(m1(1,40,1,40));
			
			matrixmult<cmatrix,cmatrix,cmatrix,cdotprecision,50> o(m1,m2);
			matrixmult<cmatrix,cmatrix_slice,cmatrix,cdotprecision,50> p(m1,m2(50,50));
			matrixmult<cmatrix_slice,cmatrix_slice,cmatrix,cdotprecision,40> q(m1(1,40,1,40),m2(1,40,1,40));
			// Real ----------------------------------
			matrixconstr<cmatrix,rvector,rmatrix,real,50> a9(m1,rv,rm2);
			matrixconstr<cmatrix,rmatrix_subv,rmatrix,real,50> a10(m1,rm2[1],rm2);
         
         matrixassign<cmatrix,rvector,rmatrix,real,50> a11(m1,rv,rm2);
         matrixassign<cmatrix_slice,rvector,rmatrix,real,50> a12(m1(50,50),rv,rm2);
         matrixassign<cmatrix,rvector_slice,rmatrix,real,50> a13(m1,rv(50),rm2);
         matrixassign<cmatrix_slice,rvector_slice,rmatrix,real,50> a14(m1(50,50),rv(50),rm2);
         matrixassign<cmatrix,rmatrix_subv,rmatrix,real,50> a15(m1,rm2[Col(7)],rm2);
         matrixassign<cmatrix_slice,rmatrix_subv,rmatrix,real,50> a16(m1(50,50),rm2[Col(7)],rm2);
         
			setfail( !e || !f || !n || !o|| !p || !q||!a1||!a2||!a3||!a4|!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a14||!a15||!a16);
//         setfail(!e);
      }
};
#endif
#ifdef TEST_CIVECTOR
template <>
class test<civector> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<civector>(void)
      {
			cimatrix m(50,50),n(50,50);
			civector v(50),w(50);
			rmatrix rm(50,50);
			rvector rv(50),rw(50);
			cmatrix cm(50,50);
			cvector cv(50),cw(50);
			imatrix im(50,50);
			ivector iv(50),iw(50);
			
			vectorconstr<civector,civector,cimatrix,cinterval,50> c(v,w,m);
			vectorconstr<civector,cimatrix_subv,cimatrix,cinterval,50> a26(v,m[9],n);
			
         vectorassign<civector,civector,cimatrix,50,cinterval> a(v,w,m);
         vectorassign<civector_slice,civector,cimatrix,50,cinterval> d(v(50),w,m);
         vectorassign<cimatrix_subv,civector,cimatrix,50,cinterval> b(m[4],w,n);

         vectoraddsub<civector,civector,50,civector> e(v,w);
			vectoraddsub<civector,cimatrix_subv,50,civector> g(v,m[7]);
			vectoraddsub<cimatrix_subv,civector,50,civector> h(m[8],v);
			vectoraddsub<cimatrix_subv,cimatrix_subv,50,civector> i(m[9],m[10]);
			vectoraddsub<civector_slice,civector,50,civector> l(v(50),w);
			vectoraddsub<civector,civector_slice,50,civector> a1(v,w(50));
			vectoraddsub<civector_slice,civector_slice,50,civector> a2(v(50),w(50));
			vectoraddsub<civector_slice,cimatrix_subv,50,civector> a3(v(50),m[11]);
			vectoraddsub<cimatrix_subv,civector_slice,50,civector> a4(m[12],w(50));

         vectoraddsubassign<civector,civector,50,civector> c1(v,w);
			vectoraddsubassign<civector,cimatrix_subv,50,civector> c2(v,m[7]);
			vectoraddsubassign<cimatrix_subv,civector,50,civector> c3(m[8],v);
			vectoraddsubassign<cimatrix_subv,cimatrix_subv,50,civector> c4(m[9],m[10]);
			vectoraddsubassign<civector_slice,civector,50,civector> c5(v(50),w);
			vectoraddsubassign<civector,civector_slice,50,civector> c6(v,w(50));
			vectoraddsubassign<civector_slice,civector_slice,50,civector> c7(v(50),w(50));
			vectoraddsubassign<civector_slice,cimatrix_subv,50,civector> c8(v(50),m[11]);
			vectoraddsubassign<cimatrix_subv,civector_slice,50,civector> c9(m[12],w(50));

			vectorscalar<civector,cinterval,50> a20(v);
			vectorscalar<civector_slice,cinterval,50> a21(v(50));
			vectorscalar<cimatrix_subv,cinterval,50> a22(m[13]);

			vectorscalarassign<civector,cinterval,civector,50> a23(v,w);
			vectorscalarassign<civector_slice,cinterval,civector,50> a24(v(50),w);
			vectorscalarassign<cimatrix_subv,cinterval,civector,50> a25(m[14],w);

			vectormult<civector,cimatrix_subv,50,cinterval> j(v,m[11]);
			vectormult<cimatrix_subv,cimatrix_subv,50,cinterval> k(m[2],m[3]);
			vectormult<civector,civector,50,cinterval> f(v,w);
			vectormult<civector,civector_slice,50,cinterval> a5(v,w(50));
			vectormult<civector_slice,civector_slice,50,cinterval> a6(v(50),w(50));
			vectormult<civector_slice,cimatrix_subv,50,cinterval> a7(v(1,50),m[13]);
			
			vectoraccu<civector,civector,50,cidotprecision> a14(v,w);
			vectoraccu<civector,cimatrix_subv,50,cidotprecision> a15(v,m[11]);
			vectoraccu<cimatrix_subv,cimatrix_subv,50,cidotprecision> a16(m[2],m[3]);
			vectoraccu<civector,civector_slice,50,cidotprecision> a17(v,w(50));
			vectoraccu<civector_slice,civector_slice,50,cidotprecision> a18(v(50),w(50));
			vectoraccu<civector_slice,cimatrix_subv,50,cidotprecision> a19(v(1,50),m[13]);

         vectorconv<civector,civector,50,civector> a187(v,w);
			vectorconv<cimatrix_subv,civector,50,civector> a188(m[8],v);
			vectorconv<cimatrix_subv,cimatrix_subv,50,civector> a189(m[9],m[10]);
			vectorconv<civector_slice,civector,50,civector> a190(v(50),w);
			vectorconv<civector_slice,civector_slice,50,civector> a191(v(50),w(50));
			vectorconv<civector_slice,cimatrix_subv,50,civector> a192(v(50),m[11]);

         vectorconvassign<civector,civector,50,civector> a193(v,w);
			vectorconvassign<cimatrix_subv,civector,50,civector> a194(m[8],v);
			vectorconvassign<cimatrix_subv,cimatrix_subv,50,civector> a195(m[9],m[10]);
			vectorconvassign<civector_slice,civector,50,civector> a196(v(50),w);
			vectorconvassign<civector_slice,civector_slice,50,civector> a197(v(50),w(50));
			vectorconvassign<civector_slice,cimatrix_subv,50,civector> a198(v(50),m[11]);

         vectorsect<civector,civector,50,civector,complex> a199(v,w);
			vectorsect<cimatrix_subv,civector,50,civector,complex> a200(m[8],v);
			vectorsect<cimatrix_subv,cimatrix_subv,50,civector,complex> a201(m[9],m[10]);
			vectorsect<civector_slice,civector,50,civector,complex> a202(v(50),w);
			vectorsect<civector_slice,civector_slice,50,civector,complex> a203(v(50),w(50));
			vectorsect<civector_slice,cimatrix_subv,50,civector,complex> a204(v(50),m[11]);

         vectorsectassign<civector,civector,50,civector,complex> a205(v,w);
			vectorsectassign<cimatrix_subv,civector,50,civector,complex> a206(m[8],v);
			vectorsectassign<cimatrix_subv,cimatrix_subv,50,civector,complex> a207(m[9],m[10]);
			vectorsectassign<civector_slice,civector,50,civector,complex> a208(v(50),w);
			vectorsectassign<civector_slice,civector_slice,50,civector,complex> a209(v(50),w(50));
			vectorsectassign<civector_slice,cimatrix_subv,50,civector,complex> a210(v(50),m[11]);

		//-------- real ------------	
			vectorconstr<civector,rvector,rmatrix,real,50> a27(v,rv,rm);
			vectorconstr<civector,rmatrix_subv,rmatrix,real,50> a29(v,rm[3],rm);
			
			vectorassign<civector,rvector,rmatrix,50,real> a28(v,rv,rm);
			vectorassign<civector_slice,rvector,rmatrix,50,real> a30(v(50),rv,rm);
			vectorassign<cimatrix_subv,rvector,rmatrix,50,real> a31(m[5],rv,rm);
			
         vectoraddsub<civector,rvector,50,civector> b1(v,rw);
         vectoraddsub<rvector,civector,50,civector> b2(rv,w);
			vectoraddsub<rvector,cimatrix_subv,50,civector> b3(rv,m[7]);
			vectoraddsub<civector,rmatrix_subv,50,civector> b4(v,rm[7]);
			vectoraddsub<rmatrix_subv,civector,50,civector> b5(rm[8],v);
			vectoraddsub<cimatrix_subv,rvector,50,civector> b6(m[8],rv);
			vectoraddsub<rmatrix_subv,cimatrix_subv,50,civector> b7(rm[9],m[10]);
			vectoraddsub<cimatrix_subv,rmatrix_subv,50,civector> b8(m[9],rm[10]);
			vectoraddsub<rvector_slice,civector,50,civector> b9(rv(50),w);
			vectoraddsub<civector_slice,rvector,50,civector> b10(v(50),rw);
			vectoraddsub<rvector,civector_slice,50,civector> b11(rv,w(50));
			vectoraddsub<civector,rvector_slice,50,civector> b12(v,rw(50));
			vectoraddsub<rvector_slice,civector_slice,50,civector> b13(rv(50),w(50));
			vectoraddsub<civector_slice,rvector_slice,50,civector> b14(v(50),rw(50));
			vectoraddsub<rvector_slice,cimatrix_subv,50,civector> b15(rv(50),m[11]);
			vectoraddsub<civector_slice,rmatrix_subv,50,civector> b16(v(50),rm[11]);
			vectoraddsub<rmatrix_subv,civector_slice,50,civector> b17(rm[12],w(50));
			vectoraddsub<cimatrix_subv,rvector_slice,50,civector> b18(m[12],rw(50));

         vectoraddsubassign<civector,rvector,50,civector> c10(v,rw);
			vectoraddsubassign<civector,rmatrix_subv,50,civector> c11(v,rm[7]);
			vectoraddsubassign<cimatrix_subv,rvector,50,civector> c12(m[8],rv);
			vectoraddsubassign<cimatrix_subv,rmatrix_subv,50,civector> c13(m[9],rm[10]);
			vectoraddsubassign<civector_slice,rvector,50,civector> c14(v(50),rw);
			vectoraddsubassign<civector,rvector_slice,50,civector> c15(v,rw(50));
			vectoraddsubassign<civector_slice,rvector_slice,50,civector> c16(v(50),rw(50));
			vectoraddsubassign<civector_slice,rmatrix_subv,50,civector> c17(v(50),rm[11]);
			vectoraddsubassign<cimatrix_subv,rvector_slice,50,civector> c18(m[12],rw(50));

			vectormult<rvector,cimatrix_subv,50,cinterval> b19(rv,m[11]);
			vectormult<civector,rmatrix_subv,50,cinterval> b20(v,rm[11]);
			vectormult<cimatrix_subv,rvector,50,cinterval> c19(m[11],rv);
			vectormult<rmatrix_subv,civector,50,cinterval> c20(rm[11],v);
			vectormult<rmatrix_subv,cimatrix_subv,50,cinterval> b21(rm[2],m[3]);
			vectormult<cimatrix_subv,rmatrix_subv,50,cinterval> b22(m[2],rm[3]);
			vectormult<rvector,civector,50,cinterval> b23(rv,w);
			vectormult<civector,rvector,50,cinterval> b24(v,rw);
			vectormult<rvector,civector_slice,50,cinterval> b25(rv,w(50));
			vectormult<civector,rvector_slice,50,cinterval> b26(v,rw(50));
			vectormult<rvector_slice,civector,50,cinterval> c23(rv(50),w);
			vectormult<civector_slice,rvector,50,cinterval> c24(v(50),rw);
			vectormult<rvector_slice,civector_slice,50,cinterval> b27(rv(50),w(50));
			vectormult<civector_slice,rvector_slice,50,cinterval> b28(v(50),rw(50));
			vectormult<rvector_slice,cimatrix_subv,50,cinterval> b29(rv(1,50),m[13]);
			vectormult<civector_slice,rmatrix_subv,50,cinterval> b30(v(1,50),rm[13]);
			vectormult<cimatrix_subv,rvector_slice,50,cinterval> c22(m[11],rv(1,50));
			vectormult<rmatrix_subv,civector_slice,50,cinterval> c21(rm[11],v(1,50));
			
			vectoraccu<civector,rvector,50,cidotprecision> a51(v,rw);
			vectoraccu<civector,rmatrix_subv,50,cidotprecision> a52(v,rm[11]);
			vectoraccu<cimatrix_subv,rmatrix_subv,50,cidotprecision> a53(m[2],rm[3]);
			vectoraccu<civector,rvector_slice,50,cidotprecision> a54(v,rw(50));
			vectoraccu<civector_slice,rvector_slice,50,cidotprecision> a55(v(50),rw(50));
			vectoraccu<civector_slice,rmatrix_subv,50,cidotprecision> a56(v(1,50),rm[13]);

			vectoraccu<rvector,civector,50,cidotprecision> a57(rv,w);
			vectoraccu<rvector,cimatrix_subv,50,cidotprecision> a58(rv,m[11]);
			vectoraccu<rmatrix_subv,cimatrix_subv,50,cidotprecision> a59(rm[2],m[3]);
			vectoraccu<rvector,civector_slice,50,cidotprecision> a60(rv,w(50));
			vectoraccu<rvector_slice,civector_slice,50,cidotprecision> a61(rv(50),w(50));
			vectoraccu<rvector_slice,cimatrix_subv,50,cidotprecision> a62(rv(1,50),m[13]);

         vectorconv<civector,rvector,50,civector> a211(v,rw);
			vectorconv<cimatrix_subv,rvector,50,civector> a212(m[8],rv);
			vectorconv<cimatrix_subv,rmatrix_subv,50,civector> a213(m[9],rm[10]);
			vectorconv<civector_slice,rvector,50,civector> a214(v(50),rw);
			vectorconv<civector_slice,rvector_slice,50,civector> a215(v(50),rw(50));
			vectorconv<civector_slice,rmatrix_subv,50,civector> a216(v(50),rm[11]);

         vectorconvassign<civector,rvector,50,civector> a217(v,rw);
			vectorconvassign<cimatrix_subv,rvector,50,civector> a218(m[8],rv);
			vectorconvassign<cimatrix_subv,rmatrix_subv,50,civector> a219(m[9],rm[10]);
			vectorconvassign<civector_slice,rvector,50,civector> a220(v(50),rw);
			vectorconvassign<civector_slice,rvector_slice,50,civector> a221(v(50),rw(50));
			vectorconvassign<civector_slice,rmatrix_subv,50,civector> a222(v(50),rm[11]);

         vectorsect<civector,rvector,50,civector,complex> a223(v,rw);
			vectorsect<cimatrix_subv,rvector,50,civector,complex> a224(m[8],rv);
			vectorsect<cimatrix_subv,rmatrix_subv,50,civector,complex> a225(m[9],rm[10]);
			vectorsect<civector_slice,rvector,50,civector,complex> a226(v(50),rw);
			vectorsect<civector_slice,rvector_slice,50,civector,complex> a227(v(50),rw(50));
			vectorsect<civector_slice,rmatrix_subv,50,civector,complex> a228(v(50),rm[11]);

         vectorsectassign<civector,rvector,50,civector,complex> a229(v,rw);
			vectorsectassign<cimatrix_subv,rvector,50,civector,complex> a230(m[8],rv);
			vectorsectassign<cimatrix_subv,rmatrix_subv,50,civector,complex> a231(m[9],rm[10]);
			vectorsectassign<civector_slice,rvector,50,civector,complex> a232(v(50),rw);
			vectorsectassign<civector_slice,rvector_slice,50,civector,complex> a233(v(50),rw(50));
			vectorsectassign<civector_slice,rmatrix_subv,50,civector,complex> a234(v(50),rm[11]);

		//-------- interval ------------	
			vectorconstr<civector,ivector,imatrix,interval,50> a63(v,iv,im);
			vectorconstr<civector,imatrix_subv,imatrix,interval,50> a64(v,im[3],im);
			
			vectorassign<civector,ivector,imatrix,50,interval> a65(v,iv,im);
			vectorassign<civector_slice,ivector,imatrix,50,interval> a66(v(50),iv,im);
			vectorassign<cimatrix_subv,ivector,imatrix,50,interval> a67(m[5],iv,im);
			
         vectoraddsub<civector,ivector,50,civector> a68(v,iw);
         vectoraddsub<ivector,civector,50,civector> a69(iv,w);
			vectoraddsub<ivector,cimatrix_subv,50,civector> a70(iv,m[7]);
			vectoraddsub<civector,imatrix_subv,50,civector> a71(v,im[7]);
			vectoraddsub<imatrix_subv,civector,50,civector> a72(im[8],v);
			vectoraddsub<cimatrix_subv,ivector,50,civector> a73(m[8],iv);
			vectoraddsub<imatrix_subv,cimatrix_subv,50,civector> a74(im[9],m[10]);
			vectoraddsub<cimatrix_subv,imatrix_subv,50,civector> a75(m[9],im[10]);
			vectoraddsub<ivector_slice,civector,50,civector> a76(iv(50),w);
			vectoraddsub<civector_slice,ivector,50,civector> a77(v(50),iw);
			vectoraddsub<ivector,civector_slice,50,civector> a78(iv,w(50));
			vectoraddsub<civector,ivector_slice,50,civector> a79(v,iw(50));
			vectoraddsub<ivector_slice,civector_slice,50,civector> a80(iv(50),w(50));
			vectoraddsub<civector_slice,ivector_slice,50,civector> a81(v(50),iw(50));
			vectoraddsub<ivector_slice,cimatrix_subv,50,civector> a82(iv(50),m[11]);
			vectoraddsub<civector_slice,imatrix_subv,50,civector> a83(v(50),im[11]);
			vectoraddsub<imatrix_subv,civector_slice,50,civector> a84(im[12],w(50));
			vectoraddsub<cimatrix_subv,ivector_slice,50,civector> a85(m[12],iw(50));

         vectoraddsubassign<civector,ivector,50,civector> a86(v,iw);
			vectoraddsubassign<civector,imatrix_subv,50,civector> a87(v,im[7]);
			vectoraddsubassign<cimatrix_subv,ivector,50,civector> a88(m[8],iv);
			vectoraddsubassign<cimatrix_subv,imatrix_subv,50,civector> a89(m[9],im[10]);
			vectoraddsubassign<civector_slice,ivector,50,civector> a90(v(50),iw);
			vectoraddsubassign<civector,ivector_slice,50,civector> a91(v,iw(50));
			vectoraddsubassign<civector_slice,ivector_slice,50,civector> a92(v(50),iw(50));
			vectoraddsubassign<civector_slice,imatrix_subv,50,civector> a93(v(50),im[11]);
			vectoraddsubassign<cimatrix_subv,ivector_slice,50,civector> a94(m[12],iw(50));

			vectormult<ivector,cimatrix_subv,50,cinterval> a95(iv,m[11]);
			vectormult<civector,imatrix_subv,50,cinterval> a96(v,im[11]);
			vectormult<cimatrix_subv,ivector,50,cinterval> a97(m[11],iv);
			vectormult<imatrix_subv,civector,50,cinterval> a98(im[11],v);
			vectormult<imatrix_subv,cimatrix_subv,50,cinterval> a99(im[2],m[3]);
			vectormult<cimatrix_subv,imatrix_subv,50,cinterval> a100(m[2],im[3]);
			vectormult<ivector,civector,50,cinterval> a101(iv,w);
			vectormult<civector,ivector,50,cinterval> a102(v,iw);
			vectormult<ivector,civector_slice,50,cinterval> a103(iv,w(50));
			vectormult<civector,ivector_slice,50,cinterval> a104(v,iw(50));
			vectormult<ivector_slice,civector,50,cinterval> a105(iv(50),w);
			vectormult<civector_slice,ivector,50,cinterval> a106(v(50),iw);
			vectormult<ivector_slice,civector_slice,50,cinterval> a107(iv(50),w(50));
			vectormult<civector_slice,ivector_slice,50,cinterval> a108(v(50),iw(50));
			vectormult<ivector_slice,cimatrix_subv,50,cinterval> a109(iv(1,50),m[13]);
			vectormult<civector_slice,imatrix_subv,50,cinterval> a110(v(1,50),im[13]);
			vectormult<cimatrix_subv,ivector_slice,50,cinterval> a111(m[11],iv(1,50));
			vectormult<imatrix_subv,civector_slice,50,cinterval> a112(im[11],v(1,50));
			
			vectoraccu<civector,ivector,50,cidotprecision> a113(v,iw);
			vectoraccu<civector,imatrix_subv,50,cidotprecision> a114(v,im[11]);
			vectoraccu<cimatrix_subv,imatrix_subv,50,cidotprecision> a115(m[2],im[3]);
			vectoraccu<civector,ivector_slice,50,cidotprecision> a116(v,iw(50));
			vectoraccu<civector_slice,ivector_slice,50,cidotprecision> a117(v(50),iw(50));
			vectoraccu<civector_slice,imatrix_subv,50,cidotprecision> a118(v(1,50),im[13]);

			vectoraccu<ivector,civector,50,cidotprecision> a119(iv,w);
			vectoraccu<ivector,cimatrix_subv,50,cidotprecision> a120(iv,m[11]);
			vectoraccu<imatrix_subv,cimatrix_subv,50,cidotprecision> a121(im[2],m[3]);
			vectoraccu<ivector,civector_slice,50,cidotprecision> a122(iv,w(50));
			vectoraccu<ivector_slice,civector_slice,50,cidotprecision> a123(iv(50),w(50));
			vectoraccu<ivector_slice,cimatrix_subv,50,cidotprecision> a124(iv(1,50),m[13]);

         vectorconv<civector,ivector,50,civector> a235(v,iw);
			vectorconv<cimatrix_subv,ivector,50,civector> a236(m[8],iv);
			vectorconv<cimatrix_subv,imatrix_subv,50,civector> a237(m[9],im[10]);
			vectorconv<civector_slice,ivector,50,civector> a238(v(50),iw);
			vectorconv<civector_slice,ivector_slice,50,civector> a239(v(50),iw(50));
			vectorconv<civector_slice,imatrix_subv,50,civector> a240(v(50),im[11]);

         vectorconvassign<civector,ivector,50,civector> a241(v,iw);
			vectorconvassign<cimatrix_subv,ivector,50,civector> a242(m[8],iv);
			vectorconvassign<cimatrix_subv,imatrix_subv,50,civector> a243(m[9],im[10]);
			vectorconvassign<civector_slice,ivector,50,civector> a244(v(50),iw);
			vectorconvassign<civector_slice,ivector_slice,50,civector> a245(v(50),iw(50));
			vectorconvassign<civector_slice,imatrix_subv,50,civector> a246(v(50),im[11]);

         vectorsect<civector,ivector,50,civector,complex> a247(v,iw);
			vectorsect<cimatrix_subv,ivector,50,civector,complex> a248(m[8],iv);
			vectorsect<cimatrix_subv,imatrix_subv,50,civector,complex> a249(m[9],im[10]);
			vectorsect<civector_slice,ivector,50,civector,complex> a250(v(50),iw);
			vectorsect<civector_slice,ivector_slice,50,civector,complex> a251(v(50),iw(50));
			vectorsect<civector_slice,imatrix_subv,50,civector,complex> a252(v(50),im[11]);

         vectorsectassign<civector,ivector,50,civector,complex> a253(v,iw);
			vectorsectassign<cimatrix_subv,ivector,50,civector,complex> a254(m[8],iv);
			vectorsectassign<cimatrix_subv,imatrix_subv,50,civector,complex> a255(m[9],im[10]);
			vectorsectassign<civector_slice,ivector,50,civector,complex> a256(v(50),iw);
			vectorsectassign<civector_slice,ivector_slice,50,civector,complex> a257(v(50),iw(50));
			vectorsectassign<civector_slice,imatrix_subv,50,civector,complex> a258(v(50),im[11]);

		//-------- complex ------------	
			vectorconstr<civector,cvector,cmatrix,complex,50> a125(v,cv,cm);
			vectorconstr<civector,cmatrix_subv,cmatrix,complex,50> a126(v,cm[3],cm);
			
			vectorassign<civector,cvector,cmatrix,50,complex> a127(v,cv,cm);
			vectorassign<civector_slice,cvector,cmatrix,50,complex> a128(v(50),cv,cm);
			vectorassign<cimatrix_subv,cvector,cmatrix,50,complex> a129(m[5],cv,cm);
			
         vectoraddsub<civector,cvector,50,civector> a130(v,cw);
         vectoraddsub<cvector,civector,50,civector> a131(cv,w);
			vectoraddsub<cvector,cimatrix_subv,50,civector> a132(cv,m[7]);
			vectoraddsub<civector,cmatrix_subv,50,civector> a133(v,cm[7]);
			vectoraddsub<cmatrix_subv,civector,50,civector> a134(cm[8],v);
			vectoraddsub<cimatrix_subv,cvector,50,civector> a135(m[8],cv);
			vectoraddsub<cmatrix_subv,cimatrix_subv,50,civector> a136(cm[9],m[10]);
			vectoraddsub<cimatrix_subv,cmatrix_subv,50,civector> a137(m[9],cm[10]);
			vectoraddsub<cvector_slice,civector,50,civector> a138(cv(50),w);
			vectoraddsub<civector_slice,cvector,50,civector> a139(v(50),cw);
			vectoraddsub<cvector,civector_slice,50,civector> a140(cv,w(50));
			vectoraddsub<civector,cvector_slice,50,civector> a141(v,cw(50));
			vectoraddsub<cvector_slice,civector_slice,50,civector> a142(cv(50),w(50));
			vectoraddsub<civector_slice,cvector_slice,50,civector> a143(v(50),cw(50));
			vectoraddsub<cvector_slice,cimatrix_subv,50,civector> a144(cv(50),m[11]);
			vectoraddsub<civector_slice,cmatrix_subv,50,civector> a145(v(50),cm[11]);
			vectoraddsub<cmatrix_subv,civector_slice,50,civector> a146(cm[12],w(50));
			vectoraddsub<cimatrix_subv,cvector_slice,50,civector> a147(m[12],cw(50));

         vectoraddsubassign<civector,cvector,50,civector> a148(v,cw);
			vectoraddsubassign<civector,cmatrix_subv,50,civector> a149(v,cm[7]);
			vectoraddsubassign<cimatrix_subv,cvector,50,civector> a150(m[8],cv);
			vectoraddsubassign<cimatrix_subv,cmatrix_subv,50,civector> a151(m[9],cm[10]);
			vectoraddsubassign<civector_slice,cvector,50,civector> a152(v(50),cw);
			vectoraddsubassign<civector,cvector_slice,50,civector> a153(v,cw(50));
			vectoraddsubassign<civector_slice,cvector_slice,50,civector> a154(v(50),cw(50));
			vectoraddsubassign<civector_slice,cmatrix_subv,50,civector> a155(v(50),cm[11]);
			vectoraddsubassign<cimatrix_subv,cvector_slice,50,civector> a156(m[12],cw(50));

			vectormult<cvector,cimatrix_subv,50,cinterval> a157(cv,m[11]);
			vectormult<civector,cmatrix_subv,50,cinterval> a158(v,cm[11]);
			vectormult<cimatrix_subv,cvector,50,cinterval> a159(m[11],cv);
			vectormult<cmatrix_subv,civector,50,cinterval> a160(cm[11],v);
			vectormult<cmatrix_subv,cimatrix_subv,50,cinterval> a161(cm[2],m[3]);
			vectormult<cimatrix_subv,cmatrix_subv,50,cinterval> a162(m[2],cm[3]);
			vectormult<cvector,civector,50,cinterval> a163(cv,w);
			vectormult<civector,cvector,50,cinterval> a164(v,cw);
			vectormult<cvector,civector_slice,50,cinterval> a165(cv,w(50));
			vectormult<civector,cvector_slice,50,cinterval> a166(v,cw(50));
			vectormult<cvector_slice,civector,50,cinterval> a167(cv(50),w);
			vectormult<civector_slice,cvector,50,cinterval> a168(v(50),cw);
			vectormult<cvector_slice,civector_slice,50,cinterval> a169(cv(50),w(50));
			vectormult<civector_slice,cvector_slice,50,cinterval> a170(v(50),cw(50));
			vectormult<cvector_slice,cimatrix_subv,50,cinterval> a171(cv(1,50),m[13]);
			vectormult<civector_slice,cmatrix_subv,50,cinterval> a172(v(1,50),cm[13]);
			vectormult<cimatrix_subv,cvector_slice,50,cinterval> a173(m[11],cv(1,50));
			vectormult<cmatrix_subv,civector_slice,50,cinterval> a174(cm[11],v(1,50));
			
			vectoraccu<civector,cvector,50,cidotprecision> a175(v,cw);
			vectoraccu<civector,cmatrix_subv,50,cidotprecision> a176(v,cm[11]);
			vectoraccu<cimatrix_subv,cmatrix_subv,50,cidotprecision> a177(m[2],cm[3]);
			vectoraccu<civector,cvector_slice,50,cidotprecision> a178(v,cw(50));
			vectoraccu<civector_slice,cvector_slice,50,cidotprecision> a179(v(50),cw(50));
			vectoraccu<civector_slice,cmatrix_subv,50,cidotprecision> a180(v(1,50),cm[13]);

			vectoraccu<cvector,civector,50,cidotprecision> a181(cv,w);
			vectoraccu<cvector,cimatrix_subv,50,cidotprecision> a182(cv,m[11]);
			vectoraccu<cmatrix_subv,cimatrix_subv,50,cidotprecision> a183(cm[2],m[3]);
			vectoraccu<cvector,civector_slice,50,cidotprecision> a184(cv,w(50));
			vectoraccu<cvector_slice,civector_slice,50,cidotprecision> a185(cv(50),w(50));
			vectoraccu<cvector_slice,cimatrix_subv,50,cidotprecision> a186(cv(1,50),m[13]);

         vectorconv<civector,cvector,50,civector> a259(v,cw);
			vectorconv<cimatrix_subv,cvector,50,civector> a260(m[8],cv);
			vectorconv<cimatrix_subv,cmatrix_subv,50,civector> a261(m[9],cm[10]);
			vectorconv<civector_slice,cvector,50,civector> a262(v(50),cw);
			vectorconv<civector_slice,cvector_slice,50,civector> a263(v(50),cw(50));
			vectorconv<civector_slice,cmatrix_subv,50,civector> a264(v(50),cm[11]);

         vectorconvassign<civector,cvector,50,civector> a265(v,cw);
			vectorconvassign<cimatrix_subv,cvector,50,civector> a266(m[8],cv);
			vectorconvassign<cimatrix_subv,cmatrix_subv,50,civector> a267(m[9],cm[10]);
			vectorconvassign<civector_slice,cvector,50,civector> a268(v(50),cw);
			vectorconvassign<civector_slice,cvector_slice,50,civector> a269(v(50),cw(50));
			vectorconvassign<civector_slice,cmatrix_subv,50,civector> a270(v(50),cm[11]);

         vectorsect<civector,cvector,50,civector,complex> a271(v,cw);
			vectorsect<cimatrix_subv,cvector,50,civector,complex> a272(m[8],cv);
			vectorsect<cimatrix_subv,cmatrix_subv,50,civector,complex> a273(m[9],cm[10]);
			vectorsect<civector_slice,cvector,50,civector,complex> a274(v(50),cw);
			vectorsect<civector_slice,cvector_slice,50,civector,complex> a275(v(50),cw(50));
			vectorsect<civector_slice,cmatrix_subv,50,civector,complex> a276(v(50),cm[11]);

         vectorsectassign<civector,cvector,50,civector,complex> a277(v,cw);
			vectorsectassign<cimatrix_subv,cvector,50,civector,complex> a278(m[8],cv);
			vectorsectassign<cimatrix_subv,cmatrix_subv,50,civector,complex> a279(m[9],cm[10]);
			vectorsectassign<civector_slice,cvector,50,civector,complex> a280(v(50),cw);
			vectorsectassign<civector_slice,cvector_slice,50,civector,complex> a281(v(50),cw(50));
			vectorsectassign<civector_slice,cmatrix_subv,50,civector,complex> a282(v(50),cm[11]);

//       complex x complex -------------
         vectorconv<cvector,cvector,50,civector> a283(cv,cw);
			vectorconv<cmatrix_subv,cvector,50,civector> a284(cm[8],cv);
			vectorconv<cmatrix_subv,cmatrix_subv,50,civector> a285(cm[9],cm[10]);
			vectorconv<cvector_slice,cvector,50,civector> a286(cv(50),cw);
			vectorconv<cvector_slice,cvector_slice,50,civector> a287(cv(50),cw(50));
			vectorconv<cvector_slice,cmatrix_subv,50,civector> a288(cv(50),cm[11]);

// 		real x complex ----------------
         vectorconv<cvector,rvector,50,civector> a289(cv,rw);
			vectorconv<cmatrix_subv,rvector,50,civector> a290(cm[8],rv);
			vectorconv<cmatrix_subv,rmatrix_subv,50,civector> a291(cm[9],rm[10]);
			vectorconv<cvector_slice,rvector,50,civector> a292(cv(50),rw);
			vectorconv<cvector_slice,rvector_slice,50,civector> a293(cv(50),rw(50));
			vectorconv<cvector_slice,rmatrix_subv,50,civector> a294(cv(50),rm[11]);

//       interval x complex ------------

         vectoraddsub<ivector,cvector,50,civector> a295(iv,cw);
         vectoraddsub<cvector,ivector,50,civector> a296(cv,iw);
			vectoraddsub<cvector,imatrix_subv,50,civector> a297(cv,im[7]);
			vectoraddsub<ivector,cmatrix_subv,50,civector> a298(iv,cm[7]);
			vectoraddsub<cmatrix_subv,ivector,50,civector> a299(cm[8],iv);
			vectoraddsub<imatrix_subv,cvector,50,civector> a300(im[8],cv);
			vectoraddsub<cmatrix_subv,imatrix_subv,50,civector> a301(cm[9],im[10]);
			vectoraddsub<imatrix_subv,cmatrix_subv,50,civector> a302(im[9],cm[10]);
			vectoraddsub<cvector_slice,ivector,50,civector> a303(cv(50),iw);
			vectoraddsub<ivector_slice,cvector,50,civector> a304(iv(50),cw);
			vectoraddsub<cvector,ivector_slice,50,civector> a305(cv,iw(50));
			vectoraddsub<ivector,cvector_slice,50,civector> a306(iv,cw(50));
			vectoraddsub<cvector_slice,ivector_slice,50,civector> a307(cv(50),iw(50));
			vectoraddsub<ivector_slice,cvector_slice,50,civector> a308(iv(50),cw(50));
			vectoraddsub<cvector_slice,imatrix_subv,50,civector> a309(cv(50),im[11]);
			vectoraddsub<ivector_slice,cmatrix_subv,50,civector> a310(iv(50),cm[11]);
			vectoraddsub<cmatrix_subv,ivector_slice,50,civector> a311(cm[12],iw(50));
			vectoraddsub<imatrix_subv,cvector_slice,50,civector> a312(im[12],cw(50));

			vectormult<cvector,imatrix_subv,50,cinterval> a313(cv,im[11]);
			vectormult<ivector,cmatrix_subv,50,cinterval> a314(iv,cm[11]);
			vectormult<imatrix_subv,cvector,50,cinterval> a315(im[11],cv);
			vectormult<cmatrix_subv,ivector,50,cinterval> a316(cm[11],iv);
			vectormult<cmatrix_subv,imatrix_subv,50,cinterval> a317(cm[2],im[3]);
			vectormult<imatrix_subv,cmatrix_subv,50,cinterval> a318(im[2],cm[3]);
			vectormult<cvector,ivector,50,cinterval> a319(cv,iw);
			vectormult<ivector,cvector,50,cinterval> a320(iv,cw);
			vectormult<cvector,ivector_slice,50,cinterval> a321(cv,iw(50));
			vectormult<ivector,cvector_slice,50,cinterval> a322(iv,cw(50));
			vectormult<cvector_slice,ivector,50,cinterval> a323(cv(50),iw);
			vectormult<ivector_slice,cvector,50,cinterval> a324(iv(50),cw);
			vectormult<cvector_slice,ivector_slice,50,cinterval> a325(cv(50),iw(50));
			vectormult<ivector_slice,cvector_slice,50,cinterval> a326(iv(50),cw(50));
			vectormult<cvector_slice,imatrix_subv,50,cinterval> a327(cv(1,50),im[13]);
			vectormult<ivector_slice,cmatrix_subv,50,cinterval> a328(iv(1,50),cm[13]);
			vectormult<imatrix_subv,cvector_slice,50,cinterval> a329(im[11],cv(1,50));
			vectormult<cmatrix_subv,ivector_slice,50,cinterval> a330(cm[11],iv(1,50));
			
			vectoraccu<ivector,cvector,50,cidotprecision> a331(iv,cw);
			vectoraccu<ivector,cmatrix_subv,50,cidotprecision> a332(iv,cm[11]);
			vectoraccu<imatrix_subv,cmatrix_subv,50,cidotprecision> a333(im[2],cm[3]);
			vectoraccu<ivector,cvector_slice,50,cidotprecision> a334(iv,cw(50));
			vectoraccu<ivector_slice,cvector_slice,50,cidotprecision> a335(iv(50),cw(50));
			vectoraccu<ivector_slice,cmatrix_subv,50,cidotprecision> a336(iv(1,50),cm[13]);

			vectoraccu<cvector,ivector,50,cidotprecision> a337(cv,iw);
			vectoraccu<cvector,imatrix_subv,50,cidotprecision> a338(cv,im[11]);
			vectoraccu<cmatrix_subv,imatrix_subv,50,cidotprecision> a339(cm[2],im[3]);
			vectoraccu<cvector,ivector_slice,50,cidotprecision> a340(cv,iw(50));
			vectoraccu<cvector_slice,ivector_slice,50,cidotprecision> a341(cv(50),iw(50));
			vectoraccu<cvector_slice,imatrix_subv,50,cidotprecision> a342(cv(1,50),im[13]);

         vectorconv<ivector,cvector,50,civector> a343(iv,cw);
			vectorconv<imatrix_subv,cvector,50,civector> a344(im[8],cv);
			vectorconv<imatrix_subv,cmatrix_subv,50,civector> a345(im[9],cm[10]);
			vectorconv<ivector_slice,cvector,50,civector> a346(iv(50),cw);
			vectorconv<ivector_slice,cvector_slice,50,civector> a347(iv(50),cw(50));
			vectorconv<ivector_slice,cmatrix_subv,50,civector> a348(iv(50),cm[11]);


         setfail(!a || !b || !c || !d || !e || !f || !g || !h || !i||!j||!k || !l || !a1||!a2||!a3||!a4||!a5||!a6||!a7||!a20||!a21||!a22||!a23||!a24||!a25||!a14||!a15||!a16||!a17||!a18||!a19||!a26||!b1||!b2||!b3||!b4||!b5||!b6||!b7||!b8||!b9||!b10||!b11||!b12||!b13||!b14||!b15||!b16||!b17||!b18||!b19||!b20||!b21||!b22||!b23||!b24||!b25||!b26||!b27||!b28||!b29||!b30||!c1||!c2||!c3||!c4||!c5||!c6||!c7||!c8||!c9||!c10||!c11||!c12||!c13||!c14||!c15||!c16||!c17||!c18||!c19||!c20||!c21||!c22||!c23||!c24||!a27||!a28||!a29||!a30||!a31||!a51||!a52||!a53||!a54||!a55||!a56||!a57||!a58||!a59||!a60||!a61||!a62||!a63||!a64||!a65||!a66||!a67||!a68||!a69||!a70||!a71||!a72||!a73||!a74||!a75||!a76||!a77||!a78||!a79||!a80||!a81||!a82||!a83||!a84||!a85||!a86||!a87||!a88||!a89||!a90||!a91||!a92||!a93||!a94||!a95||!a96||!a97||!a98||!a99||!a100||!a101||!a102||!a103||!a104||!a105||!a106||!a107||!a108||!a109||!a110||!a111||!a112||!a113||!a114||!a115||!a116||!a117||!a118||!a119||!a120||!a121||!a122||!a123||!a124||!a125||!a126||!a127||!a128||!a129||!a130||!a131||!a132||!a133||!a134||!a135||!a136||!a137||!a138||!a139||!a140||!a141||!a142||!a143||!a144||!a145||!a146||!a147||!a148||!a149||!a150||!a151||!a152||!a153||!a154||!a155||!a156||!a157||!a158||!a159||!a160||!a161||!a162||!a163||!a164||!a165||!a166||!a167||!a168||!a169||!a170||!a171||!a172||!a173||!a174||!a175||!a176||!a177||!a178||!a179||!a180||!a181||!a182||!a183||!a184||!a185||!a186||!a187||!a188||!a189||!a190||!a191||!a192||!a193||!a194||!a195||!a196||!a197||!a198||!a199||!a200||!a201||!a202||!a203||!a204||!a205||!a206||!a207||!a208||!a209||!a210||!a211||!a212||!a213||!a214||!a215||!a216||!a217||!a218||!a219||!a220||!a221||!a222||!a223||!a224||!a225||!a226||!a227||!a228||!a229||!a230||!a231||!a232||!a233||!a234||!a235||!a236||!a237||!a238||!a239||!a240||!a241||!a242||!a243||!a244||!a245||!a246||!a247||!a248||!a249||!a250||!a251||!a252||!a253||!a254||!a255||!a256||!a257||!a258||!a259||!a260||!a261||!a262||!a263||!a264||!a265||!a266||!a267||!a268||!a269||!a270||!a271||!a272||!a273||!a274||!a275||!a276||!a277||!a278||!a279||!a280||!a281||!a282||!a283||!a284||!a285||!a286||!a287||!a288||!a289||!a290||!a291||!a292||!a293||!a294||!a295||!a296||!a297||!a298||!a299||!a300||!a301||!a302||!a303||!a304||!a305||!a306||!a307||!a308||!a309||!a310||!a311||!a312||!a313||!a314||!a315||!a316||!a317||!a318||!a319||!a320||!a321||!a322||!a323||!a324||!a325||!a326||!a327||!a328||!a329||!a330||!a331||!a332||!a333||!a334||!a335||!a336||!a337||!a338||!a339||!a340||!a341||!a342||!a343||!a344||!a345||!a346||!a347||!a348);
      }
};
#endif
#ifdef TEST_CIMATRIX
template <>
class test<cimatrix> : public testclass
{
   // should test everything a cinterval should be able to do
   public:
      test<cimatrix>(void)
      {
			civector v(50);
			cimatrix m1(50,50),m2(50,50);
			rvector rv(50);
			rmatrix rm1(50,50),rm2(50,50);
			ivector iv(50);
			imatrix im1(50,50),im2(50,50);
			cvector cv(50);
			cmatrix cm1(50,50),cm2(50,50);

			matrixconstr<cimatrix,civector,cimatrix,cinterval,50> a1(m1,v,m2);
			matrixconstr<cimatrix,cimatrix_subv,cimatrix,cinterval,50> a2(m1,m2[1],m2);
			
         matrixassign<cimatrix,civector,cimatrix,cinterval,50> a3(m1,v,m2);
         matrixassign<cimatrix_slice,civector,cimatrix,cinterval,50> a4(m1(50,50),v,m2);
         matrixassign<cimatrix,civector_slice,cimatrix,cinterval,50> a5(m1,v(50),m2(50,50));
         matrixassign<cimatrix_slice,civector_slice,cimatrix,cinterval,50> a6(m1(50,50),v(50),m2);
         matrixassign<cimatrix,cimatrix_subv,cimatrix,cinterval,50> a7(m1,m2[Col(7)],m2);
         matrixassign<cimatrix_slice,cimatrix_subv,cimatrix,cinterval,50> a8(m1(50,50),m2[Col(7)],m2);
         
			matrixaddsub<cimatrix,cimatrix,cimatrix,50> e(m1,m2);
			
			matrixscalarmult<cimatrix,cinterval,cimatrix,50> f(m1);
			matrixscalarmult<cimatrix_slice,cinterval,cimatrix,40> n(m1(1,40,1,40));
			
			matrixmult<cimatrix,cimatrix,cimatrix,cidotprecision,50> o(m1,m2);
			matrixmult<cimatrix,cimatrix_slice,cimatrix,cidotprecision,50> p(m1,m2(50,50));
			matrixmult<cimatrix_slice,cimatrix_slice,cimatrix,cidotprecision,40> q(m1(1,40,1,40),m2(1,40,1,40));
			// Real ----------------------------------
			matrixconstr<cimatrix,rvector,rmatrix,real,50> a9(m1,rv,rm2);
			matrixconstr<cimatrix,rmatrix_subv,rmatrix,real,50> a10(m1,rm2[1],rm2);
         
         matrixassign<cimatrix,rvector,rmatrix,real,50> a11(m1,rv,rm2);
         matrixassign<cimatrix_slice,rvector,rmatrix,real,50> a12(m1(50,50),rv,rm2);
         matrixassign<cimatrix,rvector_slice,rmatrix,real,50> a13(m1,rv(50),rm2);
         matrixassign<cimatrix_slice,rvector_slice,rmatrix,real,50> a14(m1(50,50),rv(50),rm2);
         matrixassign<cimatrix,rmatrix_subv,rmatrix,real,50> a15(m1,rm2[Col(7)],rm2);
         matrixassign<cimatrix_slice,rmatrix_subv,rmatrix,real,50> a16(m1(50,50),rm2[Col(7)],rm2);
         
			// complex ----------------------------------
			matrixconstr<cimatrix,cvector,cmatrix,complex,50> a25(m1,cv,cm2);
			matrixconstr<cimatrix,cmatrix_subv,cmatrix,complex,50> a26(m1,cm2[1],cm2);
         
         matrixassign<cimatrix,cvector,cmatrix,complex,50> a27(m1,cv,cm2);
         matrixassign<cimatrix_slice,cvector,cmatrix,complex,50> a28(m1(50,50),cv,cm2);
         matrixassign<cimatrix,cvector_slice,cmatrix,complex,50> a29(m1,cv(50),cm2);
         matrixassign<cimatrix_slice,cvector_slice,cmatrix,complex,50> a30(m1(50,50),cv(50),cm2);
         matrixassign<cimatrix,cmatrix_subv,cmatrix,complex,50> a31(m1,cm2[Col(7)],cm2);
         matrixassign<cimatrix_slice,cmatrix_subv,cmatrix,complex,50> a32(m1(50,50),cm2[Col(7)],cm2);
         
			// interval ----------------------------------
			matrixconstr<cimatrix,ivector,imatrix,interval,50> a17(m1,iv,im2);
			matrixconstr<cimatrix,imatrix_subv,imatrix,interval,50> a18(m1,im2[1],im2);
         
         matrixassign<cimatrix,ivector,imatrix,interval,50> a19(m1,iv,im2);
         matrixassign<cimatrix_slice,ivector,imatrix,interval,50> a20(m1(50,50),iv,im2);
         matrixassign<cimatrix,ivector_slice,imatrix,interval,50> a21(m1,iv(50),im2);
         matrixassign<cimatrix_slice,ivector_slice,imatrix,interval,50> a22(m1(50,50),iv(50),im2);
         matrixassign<cimatrix,imatrix_subv,imatrix,interval,50> a23(m1,im2[Col(7)],im2);
         matrixassign<cimatrix_slice,imatrix_subv,imatrix,interval,50> a24(m1(50,50),im2[Col(7)],im2);
         
			setfail( !e || !f || !n || !o|| !p || !q||!a1||!a2||!a3||!a4|!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a14||!a15||!a16||!a17||!a18||!a19||!a20||!a21||!a22||!a23||!a24||!a25||!a26||!a27||!a28||!a29||!a30||!a31||!a32);
//         setfail(!e);
      }
};
#endif
#ifdef TEST_LRVECTOR
template <>
class test<l_rvector> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<l_rvector>(void)
      {
			l_rmatrix m(50,50),n(50,50);
			rmatrix rm(50,50);
			rvector rv(50),rw(50);
			l_rvector v(50),w(50);
			
			vectorconstr<l_rvector,l_rvector,l_rmatrix,l_real,50> c(v,w,m);
			vectorconstr<l_rvector,l_rmatrix_subv,l_rmatrix,l_real,50> a26(v,m[9],n);
			
         vectorassign<l_rvector,l_rvector,l_rmatrix,50,l_real> a(v,w,m);
         vectorassign<l_rvector_slice,l_rvector,l_rmatrix,50,l_real> d(v(50),w,m);
         vectorassign<l_rmatrix_subv,l_rvector,l_rmatrix,50,l_real> b(m[4],w,n);

         vectoraddsub<l_rvector,l_rvector,50,l_rvector> e(v,w);
			vectoraddsub<l_rvector,l_rmatrix_subv,50,l_rvector> g(v,m[7]);
			vectoraddsub<l_rmatrix_subv,l_rvector,50,l_rvector> h(m[8],v);
			vectoraddsub<l_rmatrix_subv,l_rmatrix_subv,50,l_rvector> i(m[9],m[10]);
			vectoraddsub<l_rvector_slice,l_rvector,50,l_rvector> l(v(50),w);
			vectoraddsub<l_rvector,l_rvector_slice,50,l_rvector> a1(v,w(50));
			vectoraddsub<l_rvector_slice,l_rvector_slice,50,l_rvector> a2(v(50),w(50));
			vectoraddsub<l_rvector_slice,l_rmatrix_subv,50,l_rvector> a3(v(50),m[11]);
			vectoraddsub<l_rmatrix_subv,l_rvector_slice,50,l_rvector> a4(m[12],w(50));

         vectoraddsubassign<l_rvector,l_rvector,50,l_rvector> c1(v,w);
			vectoraddsubassign<l_rvector,l_rmatrix_subv,50,l_rvector> c2(v,m[7]);
			vectoraddsubassign<l_rmatrix_subv,l_rvector,50,l_rvector> c3(m[8],v);
			vectoraddsubassign<l_rmatrix_subv,l_rmatrix_subv,50,l_rvector> c4(m[9],m[10]);
			vectoraddsubassign<l_rvector_slice,l_rvector,50,l_rvector> c5(v(50),w);
			vectoraddsubassign<l_rvector,l_rvector_slice,50,l_rvector> c6(v,w(50));
			vectoraddsubassign<l_rvector_slice,l_rvector_slice,50,l_rvector> c7(v(50),w(50));
			vectoraddsubassign<l_rvector_slice,l_rmatrix_subv,50,l_rvector> c8(v(50),m[11]);
			vectoraddsubassign<l_rmatrix_subv,l_rvector_slice,50,l_rvector> c9(m[12],w(50));

			vectorscalar<l_rvector,l_real,50> a20(v);
			vectorscalar<l_rvector_slice,l_real,50> a21(v(50));
			vectorscalar<l_rmatrix_subv,l_real,50> a22(m[13]);

			vectorscalarassign<l_rvector,l_real,l_rvector,50> a23(v,w);
			vectorscalarassign<l_rvector_slice,l_real,l_rvector,50> a24(v(50),w);
			vectorscalarassign<l_rmatrix_subv,l_real,l_rvector,50> a25(m[14],w);

			vectormult<l_rvector,l_rmatrix_subv,50,l_real> j(v,m[11]);
			vectormult<l_rmatrix_subv,l_rmatrix_subv,50,l_real> k(m[2],m[3]);
			vectormult<l_rvector,l_rvector,50,l_real> f(v,w);
			vectormult<l_rvector,l_rvector_slice,50,l_real> a5(v,w(50));
			vectormult<l_rvector_slice,l_rvector_slice,50,l_real> a6(v(50),w(50));
			vectormult<l_rvector_slice,l_rmatrix_subv,50,l_real> a7(v(1,50),m[13]);
			
			vectoraccu<l_rvector,l_rvector,50,dotprecision> a8(v,w);
			vectoraccu<l_rvector,l_rmatrix_subv,50,dotprecision> a9(v,m[11]);
			vectoraccu<l_rmatrix_subv,l_rmatrix_subv,50,dotprecision> a10(m[2],m[3]);
			vectoraccu<l_rvector,l_rvector_slice,50,dotprecision> a11(v,w(50));
			vectoraccu<l_rvector_slice,l_rvector_slice,50,dotprecision> a12(v(50),w(50));
			vectoraccu<l_rvector_slice,l_rmatrix_subv,50,dotprecision> a13(v(1,50),m[13]);
			
			vectoraccu<l_rvector,l_rvector,50,idotprecision> a14(v,w);
			vectoraccu<l_rvector,l_rmatrix_subv,50,idotprecision> a15(v,m[11]);
			vectoraccu<l_rmatrix_subv,l_rmatrix_subv,50,idotprecision> a16(m[2],m[3]);
			vectoraccu<l_rvector,l_rvector_slice,50,idotprecision> a17(v,w(50));
			vectoraccu<l_rvector_slice,l_rvector_slice,50,idotprecision> a18(v(50),w(50));
			vectoraccu<l_rvector_slice,l_rmatrix_subv,50,idotprecision> a19(v(1,50),m[13]);

		//-------- real ------------	
			vectorconstr<l_rvector,rvector,rmatrix,real,50> a27(v,rv,rm);
			vectorconstr<l_rvector,rmatrix_subv,rmatrix,real,50> a29(v,rm[3],rm);
			
			vectorassign<l_rvector,rvector,rmatrix,50,real> a28(v,rv,rm);
			vectorassign<l_rvector_slice,rvector,rmatrix,50,real> a30(v(50),rv,rm);
			vectorassign<l_rmatrix_subv,rvector,rmatrix,50,real> a31(m[5],rv,rm);
			
         vectoraddsub<l_rvector,rvector,50,l_rvector> b1(v,rw);
         vectoraddsub<rvector,l_rvector,50,l_rvector> b2(rv,w);
			vectoraddsub<rvector,l_rmatrix_subv,50,l_rvector> b3(rv,m[7]);
			vectoraddsub<l_rvector,rmatrix_subv,50,l_rvector> b4(v,rm[7]);
			vectoraddsub<rmatrix_subv,l_rvector,50,l_rvector> b5(rm[8],v);
			vectoraddsub<l_rmatrix_subv,rvector,50,l_rvector> b6(m[8],rv);
			vectoraddsub<rmatrix_subv,l_rmatrix_subv,50,l_rvector> b7(rm[9],m[10]);
			vectoraddsub<l_rmatrix_subv,rmatrix_subv,50,l_rvector> b8(m[9],rm[10]);
			vectoraddsub<rvector_slice,l_rvector,50,l_rvector> b9(rv(50),w);
			vectoraddsub<l_rvector_slice,rvector,50,l_rvector> b10(v(50),rw);
			vectoraddsub<rvector,l_rvector_slice,50,l_rvector> b11(rv,w(50));
			vectoraddsub<l_rvector,rvector_slice,50,l_rvector> b12(v,rw(50));
			vectoraddsub<rvector_slice,l_rvector_slice,50,l_rvector> b13(rv(50),w(50));
			vectoraddsub<l_rvector_slice,rvector_slice,50,l_rvector> b14(v(50),rw(50));
			vectoraddsub<rvector_slice,l_rmatrix_subv,50,l_rvector> b15(rv(50),m[11]);
			vectoraddsub<l_rvector_slice,rmatrix_subv,50,l_rvector> b16(v(50),rm[11]);
			vectoraddsub<rmatrix_subv,l_rvector_slice,50,l_rvector> b17(rm[12],w(50));
			vectoraddsub<l_rmatrix_subv,rvector_slice,50,l_rvector> b18(m[12],rw(50));

         vectoraddsubassign<l_rvector,rvector,50,l_rvector> c10(v,rw);
			vectoraddsubassign<l_rvector,rmatrix_subv,50,l_rvector> c11(v,rm[7]);
			vectoraddsubassign<l_rmatrix_subv,rvector,50,l_rvector> c12(m[8],rv);
			vectoraddsubassign<l_rmatrix_subv,rmatrix_subv,50,l_rvector> c13(m[9],rm[10]);
			vectoraddsubassign<l_rvector_slice,rvector,50,l_rvector> c14(v(50),rw);
			vectoraddsubassign<l_rvector,rvector_slice,50,l_rvector> c15(v,rw(50));
			vectoraddsubassign<l_rvector_slice,rvector_slice,50,l_rvector> c16(v(50),rw(50));
			vectoraddsubassign<l_rvector_slice,rmatrix_subv,50,l_rvector> c17(v(50),rm[11]);
			vectoraddsubassign<l_rmatrix_subv,rvector_slice,50,l_rvector> c18(m[12],rw(50));

			vectormult<rvector,l_rmatrix_subv,50,l_real> b19(rv,m[11]);
			vectormult<l_rvector,rmatrix_subv,50,l_real> b20(v,rm[11]);
			vectormult<l_rmatrix_subv,rvector,50,l_real> c19(m[11],rv);
			vectormult<rmatrix_subv,l_rvector,50,l_real> c20(rm[11],v);
			vectormult<rmatrix_subv,l_rmatrix_subv,50,l_real> b21(rm[2],m[3]);
			vectormult<l_rmatrix_subv,rmatrix_subv,50,l_real> b22(m[2],rm[3]);
			vectormult<rvector,l_rvector,50,l_real> b23(rv,w);
			vectormult<l_rvector,rvector,50,l_real> b24(v,rw);
			vectormult<rvector,l_rvector_slice,50,l_real> b25(rv,w(50));
			vectormult<l_rvector,rvector_slice,50,l_real> b26(v,rw(50));
			vectormult<rvector_slice,l_rvector,50,l_real> c23(rv(50),w);
			vectormult<l_rvector_slice,rvector,50,l_real> c24(v(50),rw);
			vectormult<rvector_slice,l_rvector_slice,50,l_real> b27(rv(50),w(50));
			vectormult<l_rvector_slice,rvector_slice,50,l_real> b28(v(50),rw(50));
			vectormult<rvector_slice,l_rmatrix_subv,50,l_real> b29(rv(1,50),m[13]);
			vectormult<l_rvector_slice,rmatrix_subv,50,l_real> b30(v(1,50),rm[13]);
			vectormult<l_rmatrix_subv,rvector_slice,50,l_real> c22(m[11],rv(1,50));
			vectormult<rmatrix_subv,l_rvector_slice,50,l_real> c21(rm[11],v(1,50));
			
			vectoraccu<l_rvector,rvector,50,dotprecision> a39(v,rw);
			vectoraccu<l_rvector,rmatrix_subv,50,dotprecision> a40(v,rm[11]);
			vectoraccu<l_rmatrix_subv,rmatrix_subv,50,dotprecision> a41(m[2],rm[3]);
			vectoraccu<l_rvector,rvector_slice,50,dotprecision> a42(v,rw(50));
			vectoraccu<l_rvector_slice,rvector_slice,50,dotprecision> a43(v(50),rw(50));
			vectoraccu<l_rvector_slice,rmatrix_subv,50,dotprecision> a44(v(1,50),rm[13]);
			
			vectoraccu<rvector,l_rvector,50,dotprecision> a45(rv,w);
			vectoraccu<rvector,l_rmatrix_subv,50,dotprecision> a46(rv,m[11]);
			vectoraccu<rmatrix_subv,l_rmatrix_subv,50,dotprecision> a47(rm[2],m[3]);
			vectoraccu<rvector,l_rvector_slice,50,dotprecision> a48(rv,w(50));
			vectoraccu<rvector_slice,l_rvector_slice,50,dotprecision> a49(rv(50),w(50));
			vectoraccu<rvector_slice,l_rmatrix_subv,50,dotprecision> a50(rv(1,50),m[13]);
			
			vectoraccu<l_rvector,rvector,50,idotprecision> a51(v,rw);
			vectoraccu<l_rvector,rmatrix_subv,50,idotprecision> a52(v,rm[11]);
			vectoraccu<l_rmatrix_subv,rmatrix_subv,50,idotprecision> a53(m[2],rm[3]);
			vectoraccu<l_rvector,rvector_slice,50,idotprecision> a54(v,rw(50));
			vectoraccu<l_rvector_slice,rvector_slice,50,idotprecision> a55(v(50),rw(50));
			vectoraccu<l_rvector_slice,rmatrix_subv,50,idotprecision> a56(v(1,50),rm[13]);

			vectoraccu<rvector,l_rvector,50,idotprecision> a57(rv,w);
			vectoraccu<rvector,l_rmatrix_subv,50,idotprecision> a58(rv,m[11]);
			vectoraccu<rmatrix_subv,l_rmatrix_subv,50,idotprecision> a59(rm[2],m[3]);
			vectoraccu<rvector,l_rvector_slice,50,idotprecision> a60(rv,w(50));
			vectoraccu<rvector_slice,l_rvector_slice,50,idotprecision> a61(rv(50),w(50));
			vectoraccu<rvector_slice,l_rmatrix_subv,50,idotprecision> a62(rv(1,50),m[13]);

         setfail(!a || !b || !c || !d || !e || !f || !g || !h || !i||!j||!k || !l || !a1||!a2||!a3||!a4||!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a20||!a21||!a22||!a23||!a24||!a25||!a14||!a15||!a16||!a17||!a18||!a19||!a26||!b1||!b2||!b3||!b4||!b5||!b6||!b7||!b8||!b9||!b10||!b11||!b12||!b13||!b14||!b15||!b16||!b17||!b18||!b19||!b20||!b21||!b22||!b23||!b24||!b25||!b26||!b27||!b28||!b29||!b30||!c1||!c2||!c3||!c4||!c5||!c6||!c7||!c8||!c9||!c10||!c11||!c12||!c13||!c14||!c15||!c16||!c17||!c18||!c19||!c20||!c21||!c22||!c23||!c24||!a27||!a28||!a29||!a30||!a31||!a39||!a40||!a41||!a42||!a43||!a44||!a45||!a46||!a47||!a48||!a49||!a50||!a51||!a52||!a53||!a54||!a55||!a56||!a57||!a58||!a59||!a60||!a61||!a62);
      }
};
#endif
#ifdef TEST_LRMATRIX
template <>
class test<l_rmatrix> : public testclass
{
   // should test everything a l_real should be able to do
   public:
      test<l_rmatrix>(void)
      {
			l_rvector v(50);
			l_rmatrix m1(50,50),m2(50,50);
			rvector rv(50);
			rmatrix rm1(50,50),rm2(50,50);

			matrixconstr<l_rmatrix,l_rvector,l_rmatrix,l_real,50> a1(m1,v,m2);
			matrixconstr<l_rmatrix,l_rmatrix_subv,l_rmatrix,l_real,50> a2(m1,m2[1],m2);
			
         matrixassign<l_rmatrix,l_rvector,l_rmatrix,l_real,50> a3(m1,v,m2);
         matrixassign<l_rmatrix_slice,l_rvector,l_rmatrix,l_real,50> a4(m1(50,50),v,m2);
         matrixassign<l_rmatrix,l_rvector_slice,l_rmatrix,l_real,50> a5(m1,v(50),m2);
         matrixassign<l_rmatrix_slice,l_rvector_slice,l_rmatrix,l_real,50> a6(m1(50,50),v(50),m2);
         matrixassign<l_rmatrix,l_rmatrix_subv,l_rmatrix,l_real,50> a7(m1,m2[Col(7)],m2);
         matrixassign<l_rmatrix_slice,l_rmatrix_subv,l_rmatrix,l_real,50> a8(m1(50,50),m2[Col(7)],m2);
         
			matrixaddsub<l_rmatrix,l_rmatrix,l_rmatrix,50> e(m1,m2);
			
			matrixscalarmult<l_rmatrix,l_real,l_rmatrix,50> f(m1);
			matrixscalarmult<l_rmatrix_slice,l_real,l_rmatrix,40> n(m1(1,40,1,40));
			
			matrixmult<l_rmatrix,l_rmatrix,l_rmatrix,dotprecision,50> o(m1,m2);
			matrixmult<l_rmatrix,l_rmatrix_slice,l_rmatrix,dotprecision,50> p(m1,m2(50,50));
			matrixmult<l_rmatrix_slice,l_rmatrix_slice,l_rmatrix,dotprecision,40> q(m1(1,40,1,40),m2(1,40,1,40));
			// Real ----------------------------------
			matrixconstr<l_rmatrix,rvector,rmatrix,real,50> a9(m1,rv,rm2);
			matrixconstr<l_rmatrix,rmatrix_subv,rmatrix,real,50> a10(m1,rm2[1],rm2);
         
         matrixassign<l_rmatrix,rvector,rmatrix,real,50> a11(m1,rv,rm2);
         matrixassign<l_rmatrix_slice,rvector,rmatrix,real,50> a12(m1(50,50),rv,rm2);
         matrixassign<l_rmatrix,rvector_slice,rmatrix,real,50> a13(m1,rv(50),rm2);
         matrixassign<l_rmatrix_slice,rvector_slice,rmatrix,real,50> a14(m1(50,50),rv(50),rm2);
         matrixassign<l_rmatrix,rmatrix_subv,rmatrix,real,50> a15(m1,rm2[Col(7)],rm2);
         matrixassign<l_rmatrix_slice,rmatrix_subv,rmatrix,real,50> a16(m1(50,50),rm2[Col(7)],rm2);
         
			setfail( !e || !f || !n || !o|| !p || !q||!a1||!a2||!a3||!a4|!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a14||!a15||!a16);
//         setfail(!e);
      }
};
#endif
#ifdef TEST_LIVECTOR
template <>
class test<l_ivector> : public testclass
{
   // should test everything a real should be able to do
   public:
      test<l_ivector>(void)
      {
			l_imatrix m(50,50),n(50,50);
			l_ivector v(50),w(50);
			rmatrix rm(50,50);
			rvector rv(50),rw(50);
			l_rmatrix cm(50,50);
			l_rvector cv(50),cw(50);
			imatrix im(50,50);
			ivector iv(50),iw(50);
			
			vectorconstr<l_ivector,l_ivector,l_imatrix,l_interval,50> c(v,w,m);
			vectorconstr<l_ivector,l_imatrix_subv,l_imatrix,l_interval,50> a26(v,m[9],n);
			
         vectorassign<l_ivector,l_ivector,l_imatrix,50,l_interval> a(v,w,m);
         vectorassign<l_ivector_slice,l_ivector,l_imatrix,50,l_interval> d(v(50),w,m);
         vectorassign<l_imatrix_subv,l_ivector,l_imatrix,50,l_interval> b(m[4],w,n);

         vectoraddsub<l_ivector,l_ivector,50,l_ivector> e(v,w);
			vectoraddsub<l_ivector,l_imatrix_subv,50,l_ivector> g(v,m[7]);
			vectoraddsub<l_imatrix_subv,l_ivector,50,l_ivector> h(m[8],v);
			vectoraddsub<l_imatrix_subv,l_imatrix_subv,50,l_ivector> i(m[9],m[10]);
			vectoraddsub<l_ivector_slice,l_ivector,50,l_ivector> l(v(50),w);
			vectoraddsub<l_ivector,l_ivector_slice,50,l_ivector> a1(v,w(50));
			vectoraddsub<l_ivector_slice,l_ivector_slice,50,l_ivector> a2(v(50),w(50));
			vectoraddsub<l_ivector_slice,l_imatrix_subv,50,l_ivector> a3(v(50),m[11]);
			vectoraddsub<l_imatrix_subv,l_ivector_slice,50,l_ivector> a4(m[12],w(50));

         vectoraddsubassign<l_ivector,l_ivector,50,l_ivector> c1(v,w);
			vectoraddsubassign<l_ivector,l_imatrix_subv,50,l_ivector> c2(v,m[7]);
			vectoraddsubassign<l_imatrix_subv,l_ivector,50,l_ivector> c3(m[8],v);
			vectoraddsubassign<l_imatrix_subv,l_imatrix_subv,50,l_ivector> c4(m[9],m[10]);
			vectoraddsubassign<l_ivector_slice,l_ivector,50,l_ivector> c5(v(50),w);
			vectoraddsubassign<l_ivector,l_ivector_slice,50,l_ivector> c6(v,w(50));
			vectoraddsubassign<l_ivector_slice,l_ivector_slice,50,l_ivector> c7(v(50),w(50));
			vectoraddsubassign<l_ivector_slice,l_imatrix_subv,50,l_ivector> c8(v(50),m[11]);
			vectoraddsubassign<l_imatrix_subv,l_ivector_slice,50,l_ivector> c9(m[12],w(50));

			vectorscalar<l_ivector,l_interval,50> a20(v);
			vectorscalar<l_ivector_slice,l_interval,50> a21(v(50));
			vectorscalar<l_imatrix_subv,l_interval,50> a22(m[13]);

			vectorscalarassign<l_ivector,l_interval,l_ivector,50> a23(v,w);
			vectorscalarassign<l_ivector_slice,l_interval,l_ivector,50> a24(v(50),w);
			vectorscalarassign<l_imatrix_subv,l_interval,l_ivector,50> a25(m[14],w);

			vectormult<l_ivector,l_imatrix_subv,50,l_interval> j(v,m[11]);
			vectormult<l_imatrix_subv,l_imatrix_subv,50,l_interval> k(m[2],m[3]);
			vectormult<l_ivector,l_ivector,50,l_interval> f(v,w);
			vectormult<l_ivector,l_ivector_slice,50,l_interval> a5(v,w(50));
			vectormult<l_ivector_slice,l_ivector_slice,50,l_interval> a6(v(50),w(50));
			vectormult<l_ivector_slice,l_imatrix_subv,50,l_interval> a7(v(1,50),m[13]);
			
			vectoraccu<l_ivector,l_ivector,50,idotprecision> a14(v,w);
			vectoraccu<l_ivector,l_imatrix_subv,50,idotprecision> a15(v,m[11]);
			vectoraccu<l_imatrix_subv,l_imatrix_subv,50,idotprecision> a16(m[2],m[3]);
			vectoraccu<l_ivector,l_ivector_slice,50,idotprecision> a17(v,w(50));
			vectoraccu<l_ivector_slice,l_ivector_slice,50,idotprecision> a18(v(50),w(50));
			vectoraccu<l_ivector_slice,l_imatrix_subv,50,idotprecision> a19(v(1,50),m[13]);

         vectorconv<l_ivector,l_ivector,50,l_ivector> a187(v,w);
			vectorconv<l_imatrix_subv,l_ivector,50,l_ivector> a188(m[8],v);
			vectorconv<l_imatrix_subv,l_imatrix_subv,50,l_ivector> a189(m[9],m[10]);
			vectorconv<l_ivector_slice,l_ivector,50,l_ivector> a190(v(50),w);
			vectorconv<l_ivector_slice,l_ivector_slice,50,l_ivector> a191(v(50),w(50));
			vectorconv<l_ivector_slice,l_imatrix_subv,50,l_ivector> a192(v(50),m[11]);

         vectorconvassign<l_ivector,l_ivector,50,l_ivector> a193(v,w);
			vectorconvassign<l_imatrix_subv,l_ivector,50,l_ivector> a194(m[8],v);
			vectorconvassign<l_imatrix_subv,l_imatrix_subv,50,l_ivector> a195(m[9],m[10]);
			vectorconvassign<l_ivector_slice,l_ivector,50,l_ivector> a196(v(50),w);
			vectorconvassign<l_ivector_slice,l_ivector_slice,50,l_ivector> a197(v(50),w(50));
			vectorconvassign<l_ivector_slice,l_imatrix_subv,50,l_ivector> a198(v(50),m[11]);

         vectorsect<l_ivector,l_ivector,50,l_ivector,l_real> a199(v,w);
			vectorsect<l_imatrix_subv,l_ivector,50,l_ivector,l_real> a200(m[8],v);
			vectorsect<l_imatrix_subv,l_imatrix_subv,50,l_ivector,l_real> a201(m[9],m[10]);
			vectorsect<l_ivector_slice,l_ivector,50,l_ivector,l_real> a202(v(50),w);
			vectorsect<l_ivector_slice,l_ivector_slice,50,l_ivector,l_real> a203(v(50),w(50));
			vectorsect<l_ivector_slice,l_imatrix_subv,50,l_ivector,l_real> a204(v(50),m[11]);

         vectorsectassign<l_ivector,l_ivector,50,l_ivector,l_real> a205(v,w);
			vectorsectassign<l_imatrix_subv,l_ivector,50,l_ivector,l_real> a206(m[8],v);
			vectorsectassign<l_imatrix_subv,l_imatrix_subv,50,l_ivector,l_real> a207(m[9],m[10]);
			vectorsectassign<l_ivector_slice,l_ivector,50,l_ivector,l_real> a208(v(50),w);
			vectorsectassign<l_ivector_slice,l_ivector_slice,50,l_ivector,l_real> a209(v(50),w(50));
			vectorsectassign<l_ivector_slice,l_imatrix_subv,50,l_ivector,l_real> a210(v(50),m[11]);

		//-------- real ------------	
			vectorconstr<l_ivector,rvector,rmatrix,real,50> a27(v,rv,rm);
			vectorconstr<l_ivector,rmatrix_subv,rmatrix,real,50> a29(v,rm[3],rm);
			
			vectorassign<l_ivector,rvector,rmatrix,50,real> a28(v,rv,rm);
			vectorassign<l_ivector_slice,rvector,rmatrix,50,real> a30(v(50),rv,rm);
			vectorassign<l_imatrix_subv,rvector,rmatrix,50,real> a31(m[5],rv,rm);
			
         vectoraddsub<l_ivector,rvector,50,l_ivector> b1(v,rw);
         vectoraddsub<rvector,l_ivector,50,l_ivector> b2(rv,w);
			vectoraddsub<rvector,l_imatrix_subv,50,l_ivector> b3(rv,m[7]);
			vectoraddsub<l_ivector,rmatrix_subv,50,l_ivector> b4(v,rm[7]);
			vectoraddsub<rmatrix_subv,l_ivector,50,l_ivector> b5(rm[8],v);
			vectoraddsub<l_imatrix_subv,rvector,50,l_ivector> b6(m[8],rv);
			vectoraddsub<rmatrix_subv,l_imatrix_subv,50,l_ivector> b7(rm[9],m[10]);
			vectoraddsub<l_imatrix_subv,rmatrix_subv,50,l_ivector> b8(m[9],rm[10]);
			vectoraddsub<rvector_slice,l_ivector,50,l_ivector> b9(rv(50),w);
			vectoraddsub<l_ivector_slice,rvector,50,l_ivector> b10(v(50),rw);
			vectoraddsub<rvector,l_ivector_slice,50,l_ivector> b11(rv,w(50));
			vectoraddsub<l_ivector,rvector_slice,50,l_ivector> b12(v,rw(50));
			vectoraddsub<rvector_slice,l_ivector_slice,50,l_ivector> b13(rv(50),w(50));
			vectoraddsub<l_ivector_slice,rvector_slice,50,l_ivector> b14(v(50),rw(50));
			vectoraddsub<rvector_slice,l_imatrix_subv,50,l_ivector> b15(rv(50),m[11]);
			vectoraddsub<l_ivector_slice,rmatrix_subv,50,l_ivector> b16(v(50),rm[11]);
			vectoraddsub<rmatrix_subv,l_ivector_slice,50,l_ivector> b17(rm[12],w(50));
			vectoraddsub<l_imatrix_subv,rvector_slice,50,l_ivector> b18(m[12],rw(50));

         vectoraddsubassign<l_ivector,rvector,50,l_ivector> c10(v,rw);
			vectoraddsubassign<l_ivector,rmatrix_subv,50,l_ivector> c11(v,rm[7]);
			vectoraddsubassign<l_imatrix_subv,rvector,50,l_ivector> c12(m[8],rv);
			vectoraddsubassign<l_imatrix_subv,rmatrix_subv,50,l_ivector> c13(m[9],rm[10]);
			vectoraddsubassign<l_ivector_slice,rvector,50,l_ivector> c14(v(50),rw);
			vectoraddsubassign<l_ivector,rvector_slice,50,l_ivector> c15(v,rw(50));
			vectoraddsubassign<l_ivector_slice,rvector_slice,50,l_ivector> c16(v(50),rw(50));
			vectoraddsubassign<l_ivector_slice,rmatrix_subv,50,l_ivector> c17(v(50),rm[11]);
			vectoraddsubassign<l_imatrix_subv,rvector_slice,50,l_ivector> c18(m[12],rw(50));

			vectormult<rvector,l_imatrix_subv,50,l_interval> b19(rv,m[11]);
			vectormult<l_ivector,rmatrix_subv,50,l_interval> b20(v,rm[11]);
			vectormult<l_imatrix_subv,rvector,50,l_interval> c19(m[11],rv);
			vectormult<rmatrix_subv,l_ivector,50,l_interval> c20(rm[11],v);
			vectormult<rmatrix_subv,l_imatrix_subv,50,l_interval> b21(rm[2],m[3]);
			vectormult<l_imatrix_subv,rmatrix_subv,50,l_interval> b22(m[2],rm[3]);
			vectormult<rvector,l_ivector,50,l_interval> b23(rv,w);
			vectormult<l_ivector,rvector,50,l_interval> b24(v,rw);
			vectormult<rvector,l_ivector_slice,50,l_interval> b25(rv,w(50));
			vectormult<l_ivector,rvector_slice,50,l_interval> b26(v,rw(50));
			vectormult<rvector_slice,l_ivector,50,l_interval> c23(rv(50),w);
			vectormult<l_ivector_slice,rvector,50,l_interval> c24(v(50),rw);
			vectormult<rvector_slice,l_ivector_slice,50,l_interval> b27(rv(50),w(50));
			vectormult<l_ivector_slice,rvector_slice,50,l_interval> b28(v(50),rw(50));
			vectormult<rvector_slice,l_imatrix_subv,50,l_interval> b29(rv(1,50),m[13]);
			vectormult<l_ivector_slice,rmatrix_subv,50,l_interval> b30(v(1,50),rm[13]);
			vectormult<l_imatrix_subv,rvector_slice,50,l_interval> c22(m[11],rv(1,50));
			vectormult<rmatrix_subv,l_ivector_slice,50,l_interval> c21(rm[11],v(1,50));
			
			vectoraccu<l_ivector,rvector,50,idotprecision> a51(v,rw);
			vectoraccu<l_ivector,rmatrix_subv,50,idotprecision> a52(v,rm[11]);
			vectoraccu<l_imatrix_subv,rmatrix_subv,50,idotprecision> a53(m[2],rm[3]);
			vectoraccu<l_ivector,rvector_slice,50,idotprecision> a54(v,rw(50));
			vectoraccu<l_ivector_slice,rvector_slice,50,idotprecision> a55(v(50),rw(50));
			vectoraccu<l_ivector_slice,rmatrix_subv,50,idotprecision> a56(v(1,50),rm[13]);

			vectoraccu<rvector,l_ivector,50,idotprecision> a57(rv,w);
			vectoraccu<rvector,l_imatrix_subv,50,idotprecision> a58(rv,m[11]);
			vectoraccu<rmatrix_subv,l_imatrix_subv,50,idotprecision> a59(rm[2],m[3]);
			vectoraccu<rvector,l_ivector_slice,50,idotprecision> a60(rv,w(50));
			vectoraccu<rvector_slice,l_ivector_slice,50,idotprecision> a61(rv(50),w(50));
			vectoraccu<rvector_slice,l_imatrix_subv,50,idotprecision> a62(rv(1,50),m[13]);

         vectorconv<l_ivector,rvector,50,l_ivector> a211(v,rw);
			vectorconv<l_imatrix_subv,rvector,50,l_ivector> a212(m[8],rv);
			vectorconv<l_imatrix_subv,rmatrix_subv,50,l_ivector> a213(m[9],rm[10]);
			vectorconv<l_ivector_slice,rvector,50,l_ivector> a214(v(50),rw);
			vectorconv<l_ivector_slice,rvector_slice,50,l_ivector> a215(v(50),rw(50));
			vectorconv<l_ivector_slice,rmatrix_subv,50,l_ivector> a216(v(50),rm[11]);

         vectorconvassign<l_ivector,rvector,50,l_ivector> a217(v,rw);
			vectorconvassign<l_imatrix_subv,rvector,50,l_ivector> a218(m[8],rv);
			vectorconvassign<l_imatrix_subv,rmatrix_subv,50,l_ivector> a219(m[9],rm[10]);
			vectorconvassign<l_ivector_slice,rvector,50,l_ivector> a220(v(50),rw);
			vectorconvassign<l_ivector_slice,rvector_slice,50,l_ivector> a221(v(50),rw(50));
			vectorconvassign<l_ivector_slice,rmatrix_subv,50,l_ivector> a222(v(50),rm[11]);

         vectorsect<l_ivector,rvector,50,l_ivector,l_real> a223(v,rw);
			vectorsect<l_imatrix_subv,rvector,50,l_ivector,l_real> a224(m[8],rv);
			vectorsect<l_imatrix_subv,rmatrix_subv,50,l_ivector,l_real> a225(m[9],rm[10]);
			vectorsect<l_ivector_slice,rvector,50,l_ivector,l_real> a226(v(50),rw);
			vectorsect<l_ivector_slice,rvector_slice,50,l_ivector,l_real> a227(v(50),rw(50));
			vectorsect<l_ivector_slice,rmatrix_subv,50,l_ivector,l_real> a228(v(50),rm[11]);

         vectorsectassign<l_ivector,rvector,50,l_ivector,l_real> a229(v,rw);
			vectorsectassign<l_imatrix_subv,rvector,50,l_ivector,l_real> a230(m[8],rv);
			vectorsectassign<l_imatrix_subv,rmatrix_subv,50,l_ivector,l_real> a231(m[9],rm[10]);
			vectorsectassign<l_ivector_slice,rvector,50,l_ivector,l_real> a232(v(50),rw);
			vectorsectassign<l_ivector_slice,rvector_slice,50,l_ivector,l_real> a233(v(50),rw(50));
			vectorsectassign<l_ivector_slice,rmatrix_subv,50,l_ivector,l_real> a234(v(50),rm[11]);

		//-------- interval ------------	
			vectorconstr<l_ivector,ivector,imatrix,interval,50> a63(v,iv,im);
			vectorconstr<l_ivector,imatrix_subv,imatrix,interval,50> a64(v,im[3],im);
			
			vectorassign<l_ivector,ivector,imatrix,50,interval> a65(v,iv,im);
			vectorassign<l_ivector_slice,ivector,imatrix,50,interval> a66(v(50),iv,im);
			vectorassign<l_imatrix_subv,ivector,imatrix,50,interval> a67(m[5],iv,im);
			
         vectoraddsub<l_ivector,ivector,50,l_ivector> a68(v,iw);
         vectoraddsub<ivector,l_ivector,50,l_ivector> a69(iv,w);
			vectoraddsub<ivector,l_imatrix_subv,50,l_ivector> a70(iv,m[7]);
			vectoraddsub<l_ivector,imatrix_subv,50,l_ivector> a71(v,im[7]);
			vectoraddsub<imatrix_subv,l_ivector,50,l_ivector> a72(im[8],v);
			vectoraddsub<l_imatrix_subv,ivector,50,l_ivector> a73(m[8],iv);
			vectoraddsub<imatrix_subv,l_imatrix_subv,50,l_ivector> a74(im[9],m[10]);
			vectoraddsub<l_imatrix_subv,imatrix_subv,50,l_ivector> a75(m[9],im[10]);
			vectoraddsub<ivector_slice,l_ivector,50,l_ivector> a76(iv(50),w);
			vectoraddsub<l_ivector_slice,ivector,50,l_ivector> a77(v(50),iw);
			vectoraddsub<ivector,l_ivector_slice,50,l_ivector> a78(iv,w(50));
			vectoraddsub<l_ivector,ivector_slice,50,l_ivector> a79(v,iw(50));
			vectoraddsub<ivector_slice,l_ivector_slice,50,l_ivector> a80(iv(50),w(50));
			vectoraddsub<l_ivector_slice,ivector_slice,50,l_ivector> a81(v(50),iw(50));
			vectoraddsub<ivector_slice,l_imatrix_subv,50,l_ivector> a82(iv(50),m[11]);
			vectoraddsub<l_ivector_slice,imatrix_subv,50,l_ivector> a83(v(50),im[11]);
			vectoraddsub<imatrix_subv,l_ivector_slice,50,l_ivector> a84(im[12],w(50));
			vectoraddsub<l_imatrix_subv,ivector_slice,50,l_ivector> a85(m[12],iw(50));

         vectoraddsubassign<l_ivector,ivector,50,l_ivector> a86(v,iw);
			vectoraddsubassign<l_ivector,imatrix_subv,50,l_ivector> a87(v,im[7]);
			vectoraddsubassign<l_imatrix_subv,ivector,50,l_ivector> a88(m[8],iv);
			vectoraddsubassign<l_imatrix_subv,imatrix_subv,50,l_ivector> a89(m[9],im[10]);
			vectoraddsubassign<l_ivector_slice,ivector,50,l_ivector> a90(v(50),iw);
			vectoraddsubassign<l_ivector,ivector_slice,50,l_ivector> a91(v,iw(50));
			vectoraddsubassign<l_ivector_slice,ivector_slice,50,l_ivector> a92(v(50),iw(50));
			vectoraddsubassign<l_ivector_slice,imatrix_subv,50,l_ivector> a93(v(50),im[11]);
			vectoraddsubassign<l_imatrix_subv,ivector_slice,50,l_ivector> a94(m[12],iw(50));

			vectormult<ivector,l_imatrix_subv,50,l_interval> a95(iv,m[11]);
			vectormult<l_ivector,imatrix_subv,50,l_interval> a96(v,im[11]);
			vectormult<l_imatrix_subv,ivector,50,l_interval> a97(m[11],iv);
			vectormult<imatrix_subv,l_ivector,50,l_interval> a98(im[11],v);
			vectormult<imatrix_subv,l_imatrix_subv,50,l_interval> a99(im[2],m[3]);
			vectormult<l_imatrix_subv,imatrix_subv,50,l_interval> a100(m[2],im[3]);
			vectormult<ivector,l_ivector,50,l_interval> a101(iv,w);
			vectormult<l_ivector,ivector,50,l_interval> a102(v,iw);
			vectormult<ivector,l_ivector_slice,50,l_interval> a103(iv,w(50));
			vectormult<l_ivector,ivector_slice,50,l_interval> a104(v,iw(50));
			vectormult<ivector_slice,l_ivector,50,l_interval> a105(iv(50),w);
			vectormult<l_ivector_slice,ivector,50,l_interval> a106(v(50),iw);
			vectormult<ivector_slice,l_ivector_slice,50,l_interval> a107(iv(50),w(50));
			vectormult<l_ivector_slice,ivector_slice,50,l_interval> a108(v(50),iw(50));
			vectormult<ivector_slice,l_imatrix_subv,50,l_interval> a109(iv(1,50),m[13]);
			vectormult<l_ivector_slice,imatrix_subv,50,l_interval> a110(v(1,50),im[13]);
			vectormult<l_imatrix_subv,ivector_slice,50,l_interval> a111(m[11],iv(1,50));
			vectormult<imatrix_subv,l_ivector_slice,50,l_interval> a112(im[11],v(1,50));
			
			vectoraccu<l_ivector,ivector,50,idotprecision> a113(v,iw);
			vectoraccu<l_ivector,imatrix_subv,50,idotprecision> a114(v,im[11]);
			vectoraccu<l_imatrix_subv,imatrix_subv,50,idotprecision> a115(m[2],im[3]);
			vectoraccu<l_ivector,ivector_slice,50,idotprecision> a116(v,iw(50));
			vectoraccu<l_ivector_slice,ivector_slice,50,idotprecision> a117(v(50),iw(50));
			vectoraccu<l_ivector_slice,imatrix_subv,50,idotprecision> a118(v(1,50),im[13]);

			vectoraccu<ivector,l_ivector,50,idotprecision> a119(iv,w);
			vectoraccu<ivector,l_imatrix_subv,50,idotprecision> a120(iv,m[11]);
			vectoraccu<imatrix_subv,l_imatrix_subv,50,idotprecision> a121(im[2],m[3]);
			vectoraccu<ivector,l_ivector_slice,50,idotprecision> a122(iv,w(50));
			vectoraccu<ivector_slice,l_ivector_slice,50,idotprecision> a123(iv(50),w(50));
			vectoraccu<ivector_slice,l_imatrix_subv,50,idotprecision> a124(iv(1,50),m[13]);

         vectorconv<l_ivector,ivector,50,l_ivector> a235(v,iw);
			vectorconv<l_imatrix_subv,ivector,50,l_ivector> a236(m[8],iv);
			vectorconv<l_imatrix_subv,imatrix_subv,50,l_ivector> a237(m[9],im[10]);
			vectorconv<l_ivector_slice,ivector,50,l_ivector> a238(v(50),iw);
			vectorconv<l_ivector_slice,ivector_slice,50,l_ivector> a239(v(50),iw(50));
			vectorconv<l_ivector_slice,imatrix_subv,50,l_ivector> a240(v(50),im[11]);

         vectorconvassign<l_ivector,ivector,50,l_ivector> a241(v,iw);
			vectorconvassign<l_imatrix_subv,ivector,50,l_ivector> a242(m[8],iv);
			vectorconvassign<l_imatrix_subv,imatrix_subv,50,l_ivector> a243(m[9],im[10]);
			vectorconvassign<l_ivector_slice,ivector,50,l_ivector> a244(v(50),iw);
			vectorconvassign<l_ivector_slice,ivector_slice,50,l_ivector> a245(v(50),iw(50));
			vectorconvassign<l_ivector_slice,imatrix_subv,50,l_ivector> a246(v(50),im[11]);

         vectorsect<l_ivector,ivector,50,l_ivector,l_real> a247(v,iw);
			vectorsect<l_imatrix_subv,ivector,50,l_ivector,l_real> a248(m[8],iv);
			vectorsect<l_imatrix_subv,imatrix_subv,50,l_ivector,l_real> a249(m[9],im[10]);
			vectorsect<l_ivector_slice,ivector,50,l_ivector,l_real> a250(v(50),iw);
			vectorsect<l_ivector_slice,ivector_slice,50,l_ivector,l_real> a251(v(50),iw(50));
			vectorsect<l_ivector_slice,imatrix_subv,50,l_ivector,l_real> a252(v(50),im[11]);

         vectorsectassign<l_ivector,ivector,50,l_ivector,l_real> a253(v,iw);
			vectorsectassign<l_imatrix_subv,ivector,50,l_ivector,l_real> a254(m[8],iv);
			vectorsectassign<l_imatrix_subv,imatrix_subv,50,l_ivector,l_real> a255(m[9],im[10]);
			vectorsectassign<l_ivector_slice,ivector,50,l_ivector,l_real> a256(v(50),iw);
			vectorsectassign<l_ivector_slice,ivector_slice,50,l_ivector,l_real> a257(v(50),iw(50));
			vectorsectassign<l_ivector_slice,imatrix_subv,50,l_ivector,l_real> a258(v(50),im[11]);

		//-------- l_real ------------	
			vectorconstr<l_ivector,l_rvector,l_rmatrix,l_real,50> a125(v,cv,cm);
			vectorconstr<l_ivector,l_rmatrix_subv,l_rmatrix,l_real,50> a126(v,cm[3],cm);
			
			vectorassign<l_ivector,l_rvector,l_rmatrix,50,l_real> a127(v,cv,cm);
			vectorassign<l_ivector_slice,l_rvector,l_rmatrix,50,l_real> a128(v(50),cv,cm);
			vectorassign<l_imatrix_subv,l_rvector,l_rmatrix,50,l_real> a129(m[5],cv,cm);
			
         vectoraddsub<l_ivector,l_rvector,50,l_ivector> a130(v,cw);
         vectoraddsub<l_rvector,l_ivector,50,l_ivector> a131(cv,w);
			vectoraddsub<l_rvector,l_imatrix_subv,50,l_ivector> a132(cv,m[7]);
			vectoraddsub<l_ivector,l_rmatrix_subv,50,l_ivector> a133(v,cm[7]);
			vectoraddsub<l_rmatrix_subv,l_ivector,50,l_ivector> a134(cm[8],v);
			vectoraddsub<l_imatrix_subv,l_rvector,50,l_ivector> a135(m[8],cv);
			vectoraddsub<l_rmatrix_subv,l_imatrix_subv,50,l_ivector> a136(cm[9],m[10]);
			vectoraddsub<l_imatrix_subv,l_rmatrix_subv,50,l_ivector> a137(m[9],cm[10]);
			vectoraddsub<l_rvector_slice,l_ivector,50,l_ivector> a138(cv(50),w);
			vectoraddsub<l_ivector_slice,l_rvector,50,l_ivector> a139(v(50),cw);
			vectoraddsub<l_rvector,l_ivector_slice,50,l_ivector> a140(cv,w(50));
			vectoraddsub<l_ivector,l_rvector_slice,50,l_ivector> a141(v,cw(50));
			vectoraddsub<l_rvector_slice,l_ivector_slice,50,l_ivector> a142(cv(50),w(50));
			vectoraddsub<l_ivector_slice,l_rvector_slice,50,l_ivector> a143(v(50),cw(50));
			vectoraddsub<l_rvector_slice,l_imatrix_subv,50,l_ivector> a144(cv(50),m[11]);
			vectoraddsub<l_ivector_slice,l_rmatrix_subv,50,l_ivector> a145(v(50),cm[11]);
			vectoraddsub<l_rmatrix_subv,l_ivector_slice,50,l_ivector> a146(cm[12],w(50));
			vectoraddsub<l_imatrix_subv,l_rvector_slice,50,l_ivector> a147(m[12],cw(50));

         vectoraddsubassign<l_ivector,l_rvector,50,l_ivector> a148(v,cw);
			vectoraddsubassign<l_ivector,l_rmatrix_subv,50,l_ivector> a149(v,cm[7]);
			vectoraddsubassign<l_imatrix_subv,l_rvector,50,l_ivector> a150(m[8],cv);
			vectoraddsubassign<l_imatrix_subv,l_rmatrix_subv,50,l_ivector> a151(m[9],cm[10]);
			vectoraddsubassign<l_ivector_slice,l_rvector,50,l_ivector> a152(v(50),cw);
			vectoraddsubassign<l_ivector,l_rvector_slice,50,l_ivector> a153(v,cw(50));
			vectoraddsubassign<l_ivector_slice,l_rvector_slice,50,l_ivector> a154(v(50),cw(50));
			vectoraddsubassign<l_ivector_slice,l_rmatrix_subv,50,l_ivector> a155(v(50),cm[11]);
			vectoraddsubassign<l_imatrix_subv,l_rvector_slice,50,l_ivector> a156(m[12],cw(50));

			vectormult<l_rvector,l_imatrix_subv,50,l_interval> a157(cv,m[11]);
			vectormult<l_ivector,l_rmatrix_subv,50,l_interval> a158(v,cm[11]);
			vectormult<l_imatrix_subv,l_rvector,50,l_interval> a159(m[11],cv);
			vectormult<l_rmatrix_subv,l_ivector,50,l_interval> a160(cm[11],v);
			vectormult<l_rmatrix_subv,l_imatrix_subv,50,l_interval> a161(cm[2],m[3]);
			vectormult<l_imatrix_subv,l_rmatrix_subv,50,l_interval> a162(m[2],cm[3]);
			vectormult<l_rvector,l_ivector,50,l_interval> a163(cv,w);
			vectormult<l_ivector,l_rvector,50,l_interval> a164(v,cw);
			vectormult<l_rvector,l_ivector_slice,50,l_interval> a165(cv,w(50));
			vectormult<l_ivector,l_rvector_slice,50,l_interval> a166(v,cw(50));
			vectormult<l_rvector_slice,l_ivector,50,l_interval> a167(cv(50),w);
			vectormult<l_ivector_slice,l_rvector,50,l_interval> a168(v(50),cw);
			vectormult<l_rvector_slice,l_ivector_slice,50,l_interval> a169(cv(50),w(50));
			vectormult<l_ivector_slice,l_rvector_slice,50,l_interval> a170(v(50),cw(50));
			vectormult<l_rvector_slice,l_imatrix_subv,50,l_interval> a171(cv(1,50),m[13]);
			vectormult<l_ivector_slice,l_rmatrix_subv,50,l_interval> a172(v(1,50),cm[13]);
			vectormult<l_imatrix_subv,l_rvector_slice,50,l_interval> a173(m[11],cv(1,50));
			vectormult<l_rmatrix_subv,l_ivector_slice,50,l_interval> a174(cm[11],v(1,50));
			
			vectoraccu<l_ivector,l_rvector,50,idotprecision> a175(v,cw);
			vectoraccu<l_ivector,l_rmatrix_subv,50,idotprecision> a176(v,cm[11]);
			vectoraccu<l_imatrix_subv,l_rmatrix_subv,50,idotprecision> a177(m[2],cm[3]);
			vectoraccu<l_ivector,l_rvector_slice,50,idotprecision> a178(v,cw(50));
			vectoraccu<l_ivector_slice,l_rvector_slice,50,idotprecision> a179(v(50),cw(50));
			vectoraccu<l_ivector_slice,l_rmatrix_subv,50,idotprecision> a180(v(1,50),cm[13]);

			vectoraccu<l_rvector,l_ivector,50,idotprecision> a181(cv,w);
			vectoraccu<l_rvector,l_imatrix_subv,50,idotprecision> a182(cv,m[11]);
			vectoraccu<l_rmatrix_subv,l_imatrix_subv,50,idotprecision> a183(cm[2],m[3]);
			vectoraccu<l_rvector,l_ivector_slice,50,idotprecision> a184(cv,w(50));
			vectoraccu<l_rvector_slice,l_ivector_slice,50,idotprecision> a185(cv(50),w(50));
			vectoraccu<l_rvector_slice,l_imatrix_subv,50,idotprecision> a186(cv(1,50),m[13]);

         vectorconv<l_ivector,l_rvector,50,l_ivector> a259(v,cw);
			vectorconv<l_imatrix_subv,l_rvector,50,l_ivector> a260(m[8],cv);
			vectorconv<l_imatrix_subv,l_rmatrix_subv,50,l_ivector> a261(m[9],cm[10]);
			vectorconv<l_ivector_slice,l_rvector,50,l_ivector> a262(v(50),cw);
			vectorconv<l_ivector_slice,l_rvector_slice,50,l_ivector> a263(v(50),cw(50));
			vectorconv<l_ivector_slice,l_rmatrix_subv,50,l_ivector> a264(v(50),cm[11]);

         vectorconvassign<l_ivector,l_rvector,50,l_ivector> a265(v,cw);
			vectorconvassign<l_imatrix_subv,l_rvector,50,l_ivector> a266(m[8],cv);
			vectorconvassign<l_imatrix_subv,l_rmatrix_subv,50,l_ivector> a267(m[9],cm[10]);
			vectorconvassign<l_ivector_slice,l_rvector,50,l_ivector> a268(v(50),cw);
			vectorconvassign<l_ivector_slice,l_rvector_slice,50,l_ivector> a269(v(50),cw(50));
			vectorconvassign<l_ivector_slice,l_rmatrix_subv,50,l_ivector> a270(v(50),cm[11]);

         vectorsect<l_ivector,l_rvector,50,l_ivector,l_real> a271(v,cw);
			vectorsect<l_imatrix_subv,l_rvector,50,l_ivector,l_real> a272(m[8],cv);
			vectorsect<l_imatrix_subv,l_rmatrix_subv,50,l_ivector,l_real> a273(m[9],cm[10]);
			vectorsect<l_ivector_slice,l_rvector,50,l_ivector,l_real> a274(v(50),cw);
			vectorsect<l_ivector_slice,l_rvector_slice,50,l_ivector,l_real> a275(v(50),cw(50));
			vectorsect<l_ivector_slice,l_rmatrix_subv,50,l_ivector,l_real> a276(v(50),cm[11]);

         vectorsectassign<l_ivector,l_rvector,50,l_ivector,l_real> a277(v,cw);
			vectorsectassign<l_imatrix_subv,l_rvector,50,l_ivector,l_real> a278(m[8],cv);
			vectorsectassign<l_imatrix_subv,l_rmatrix_subv,50,l_ivector,l_real> a279(m[9],cm[10]);
			vectorsectassign<l_ivector_slice,l_rvector,50,l_ivector,l_real> a280(v(50),cw);
			vectorsectassign<l_ivector_slice,l_rvector_slice,50,l_ivector,l_real> a281(v(50),cw(50));
			vectorsectassign<l_ivector_slice,l_rmatrix_subv,50,l_ivector,l_real> a282(v(50),cm[11]);

//       l_real x l_real -------------
         vectorconv<l_rvector,l_rvector,50,l_ivector> a283(cv,cw);
			vectorconv<l_rmatrix_subv,l_rvector,50,l_ivector> a284(cm[8],cv);
			vectorconv<l_rmatrix_subv,l_rmatrix_subv,50,l_ivector> a285(cm[9],cm[10]);
			vectorconv<l_rvector_slice,l_rvector,50,l_ivector> a286(cv(50),cw);
			vectorconv<l_rvector_slice,l_rvector_slice,50,l_ivector> a287(cv(50),cw(50));
			vectorconv<l_rvector_slice,l_rmatrix_subv,50,l_ivector> a288(cv(50),cm[11]);

// 		real x l_real ----------------
         vectorconv<l_rvector,rvector,50,l_ivector> a289(cv,rw);
			vectorconv<l_rmatrix_subv,rvector,50,l_ivector> a290(cm[8],rv);
			vectorconv<l_rmatrix_subv,rmatrix_subv,50,l_ivector> a291(cm[9],rm[10]);
			vectorconv<l_rvector_slice,rvector,50,l_ivector> a292(cv(50),rw);
			vectorconv<l_rvector_slice,rvector_slice,50,l_ivector> a293(cv(50),rw(50));
			vectorconv<l_rvector_slice,rmatrix_subv,50,l_ivector> a294(cv(50),rm[11]);

//       interval x l_real ------------

         vectoraddsub<ivector,l_rvector,50,l_ivector> a295(iv,cw);
         vectoraddsub<l_rvector,ivector,50,l_ivector> a296(cv,iw);
			vectoraddsub<l_rvector,imatrix_subv,50,l_ivector> a297(cv,im[7]);
			vectoraddsub<ivector,l_rmatrix_subv,50,l_ivector> a298(iv,cm[7]);
			vectoraddsub<l_rmatrix_subv,ivector,50,l_ivector> a299(cm[8],iv);
			vectoraddsub<imatrix_subv,l_rvector,50,l_ivector> a300(im[8],cv);
			vectoraddsub<l_rmatrix_subv,imatrix_subv,50,l_ivector> a301(cm[9],im[10]);
			vectoraddsub<imatrix_subv,l_rmatrix_subv,50,l_ivector> a302(im[9],cm[10]);
			vectoraddsub<l_rvector_slice,ivector,50,l_ivector> a303(cv(50),iw);
			vectoraddsub<ivector_slice,l_rvector,50,l_ivector> a304(iv(50),cw);
			vectoraddsub<l_rvector,ivector_slice,50,l_ivector> a305(cv,iw(50));
			vectoraddsub<ivector,l_rvector_slice,50,l_ivector> a306(iv,cw(50));
			vectoraddsub<l_rvector_slice,ivector_slice,50,l_ivector> a307(cv(50),iw(50));
			vectoraddsub<ivector_slice,l_rvector_slice,50,l_ivector> a308(iv(50),cw(50));
			vectoraddsub<l_rvector_slice,imatrix_subv,50,l_ivector> a309(cv(50),im[11]);
			vectoraddsub<ivector_slice,l_rmatrix_subv,50,l_ivector> a310(iv(50),cm[11]);
			vectoraddsub<l_rmatrix_subv,ivector_slice,50,l_ivector> a311(cm[12],iw(50));
			vectoraddsub<imatrix_subv,l_rvector_slice,50,l_ivector> a312(im[12],cw(50));

			vectormult<l_rvector,imatrix_subv,50,l_interval> a313(cv,im[11]);
			vectormult<ivector,l_rmatrix_subv,50,l_interval> a314(iv,cm[11]);
			vectormult<imatrix_subv,l_rvector,50,l_interval> a315(im[11],cv);
			vectormult<l_rmatrix_subv,ivector,50,l_interval> a316(cm[11],iv);
			vectormult<l_rmatrix_subv,imatrix_subv,50,l_interval> a317(cm[2],im[3]);
			vectormult<imatrix_subv,l_rmatrix_subv,50,l_interval> a318(im[2],cm[3]);
			vectormult<l_rvector,ivector,50,l_interval> a319(cv,iw);
			vectormult<ivector,l_rvector,50,l_interval> a320(iv,cw);
			vectormult<l_rvector,ivector_slice,50,l_interval> a321(cv,iw(50));
			vectormult<ivector,l_rvector_slice,50,l_interval> a322(iv,cw(50));
			vectormult<l_rvector_slice,ivector,50,l_interval> a323(cv(50),iw);
			vectormult<ivector_slice,l_rvector,50,l_interval> a324(iv(50),cw);
			vectormult<l_rvector_slice,ivector_slice,50,l_interval> a325(cv(50),iw(50));
			vectormult<ivector_slice,l_rvector_slice,50,l_interval> a326(iv(50),cw(50));
			vectormult<l_rvector_slice,imatrix_subv,50,l_interval> a327(cv(1,50),im[13]);
			vectormult<ivector_slice,l_rmatrix_subv,50,l_interval> a328(iv(1,50),cm[13]);
			vectormult<imatrix_subv,l_rvector_slice,50,l_interval> a329(im[11],cv(1,50));
			vectormult<l_rmatrix_subv,ivector_slice,50,l_interval> a330(cm[11],iv(1,50));
			
			vectoraccu<ivector,l_rvector,50,idotprecision> a331(iv,cw);
			vectoraccu<ivector,l_rmatrix_subv,50,idotprecision> a332(iv,cm[11]);
			vectoraccu<imatrix_subv,l_rmatrix_subv,50,idotprecision> a333(im[2],cm[3]);
			vectoraccu<ivector,l_rvector_slice,50,idotprecision> a334(iv,cw(50));
			vectoraccu<ivector_slice,l_rvector_slice,50,idotprecision> a335(iv(50),cw(50));
			vectoraccu<ivector_slice,l_rmatrix_subv,50,idotprecision> a336(iv(1,50),cm[13]);

			vectoraccu<l_rvector,ivector,50,idotprecision> a337(cv,iw);
			vectoraccu<l_rvector,imatrix_subv,50,idotprecision> a338(cv,im[11]);
			vectoraccu<l_rmatrix_subv,imatrix_subv,50,idotprecision> a339(cm[2],im[3]);
			vectoraccu<l_rvector,ivector_slice,50,idotprecision> a340(cv,iw(50));
			vectoraccu<l_rvector_slice,ivector_slice,50,idotprecision> a341(cv(50),iw(50));
			vectoraccu<l_rvector_slice,imatrix_subv,50,idotprecision> a342(cv(1,50),im[13]);

         vectorconv<ivector,l_rvector,50,l_ivector> a343(iv,cw);
			vectorconv<imatrix_subv,l_rvector,50,l_ivector> a344(im[8],cv);
			vectorconv<imatrix_subv,l_rmatrix_subv,50,l_ivector> a345(im[9],cm[10]);
			vectorconv<ivector_slice,l_rvector,50,l_ivector> a346(iv(50),cw);
			vectorconv<ivector_slice,l_rvector_slice,50,l_ivector> a347(iv(50),cw(50));
			vectorconv<ivector_slice,l_rmatrix_subv,50,l_ivector> a348(iv(50),cm[11]);


         setfail(!a || !b || !c || !d || !e || !f || !g || !h || !i||!j||!k || !l || !a1||!a2||!a3||!a4||!a5||!a6||!a7||!a20||!a21||!a22||!a23||!a24||!a25||!a14||!a15||!a16||!a17||!a18||!a19||!a26||!b1||!b2||!b3||!b4||!b5||!b6||!b7||!b8||!b9||!b10||!b11||!b12||!b13||!b14||!b15||!b16||!b17||!b18||!b19||!b20||!b21||!b22||!b23||!b24||!b25||!b26||!b27||!b28||!b29||!b30||!c1||!c2||!c3||!c4||!c5||!c6||!c7||!c8||!c9||!c10||!c11||!c12||!c13||!c14||!c15||!c16||!c17||!c18||!c19||!c20||!c21||!c22||!c23||!c24||!a27||!a28||!a29||!a30||!a31||!a51||!a52||!a53||!a54||!a55||!a56||!a57||!a58||!a59||!a60||!a61||!a62||!a63||!a64||!a65||!a66||!a67||!a68||!a69||!a70||!a71||!a72||!a73||!a74||!a75||!a76||!a77||!a78||!a79||!a80||!a81||!a82||!a83||!a84||!a85||!a86||!a87||!a88||!a89||!a90||!a91||!a92||!a93||!a94||!a95||!a96||!a97||!a98||!a99||!a100||!a101||!a102||!a103||!a104||!a105||!a106||!a107||!a108||!a109||!a110||!a111||!a112||!a113||!a114||!a115||!a116||!a117||!a118||!a119||!a120||!a121||!a122||!a123||!a124||!a125||!a126||!a127||!a128||!a129||!a130||!a131||!a132||!a133||!a134||!a135||!a136||!a137||!a138||!a139||!a140||!a141||!a142||!a143||!a144||!a145||!a146||!a147||!a148||!a149||!a150||!a151||!a152||!a153||!a154||!a155||!a156||!a157||!a158||!a159||!a160||!a161||!a162||!a163||!a164||!a165||!a166||!a167||!a168||!a169||!a170||!a171||!a172||!a173||!a174||!a175||!a176||!a177||!a178||!a179||!a180||!a181||!a182||!a183||!a184||!a185||!a186||!a187||!a188||!a189||!a190||!a191||!a192||!a193||!a194||!a195||!a196||!a197||!a198||!a199||!a200||!a201||!a202||!a203||!a204||!a205||!a206||!a207||!a208||!a209||!a210||!a211||!a212||!a213||!a214||!a215||!a216||!a217||!a218||!a219||!a220||!a221||!a222||!a223||!a224||!a225||!a226||!a227||!a228||!a229||!a230||!a231||!a232||!a233||!a234||!a235||!a236||!a237||!a238||!a239||!a240||!a241||!a242||!a243||!a244||!a245||!a246||!a247||!a248||!a249||!a250||!a251||!a252||!a253||!a254||!a255||!a256||!a257||!a258||!a259||!a260||!a261||!a262||!a263||!a264||!a265||!a266||!a267||!a268||!a269||!a270||!a271||!a272||!a273||!a274||!a275||!a276||!a277||!a278||!a279||!a280||!a281||!a282||!a283||!a284||!a285||!a286||!a287||!a288||!a289||!a290||!a291||!a292||!a293||!a294||!a295||!a296||!a297||!a298||!a299||!a300||!a301||!a302||!a303||!a304||!a305||!a306||!a307||!a308||!a309||!a310||!a311||!a312||!a313||!a314||!a315||!a316||!a317||!a318||!a319||!a320||!a321||!a322||!a323||!a324||!a325||!a326||!a327||!a328||!a329||!a330||!a331||!a332||!a333||!a334||!a335||!a336||!a337||!a338||!a339||!a340||!a341||!a342||!a343||!a344||!a345||!a346||!a347||!a348);
      }
};
#endif
#ifdef TEST_LIMATRIX
template <>
class test<l_imatrix> : public testclass
{
   // should test everything a l_interval should be able to do
   public:
      test<l_imatrix>(void)
      {
			l_ivector v(50);
			l_imatrix m1(50,50),m2(50,50);
			rvector rv(50);
			rmatrix rm1(50,50),rm2(50,50);
			ivector iv(50);
			imatrix im1(50,50),im2(50,50);
			l_rvector cv(50);
			l_rmatrix cm1(50,50),cm2(50,50);

			matrixconstr<l_imatrix,l_ivector,l_imatrix,l_interval,50> a1(m1,v,m2);
			matrixconstr<l_imatrix,l_imatrix_subv,l_imatrix,l_interval,50> a2(m1,m2[1],m2);
			
         matrixassign<l_imatrix,l_ivector,l_imatrix,l_interval,50> a3(m1,v,m2);
         matrixassign<l_imatrix_slice,l_ivector,l_imatrix,l_interval,50> a4(m1(50,50),v,m2);
         matrixassign<l_imatrix,l_ivector_slice,l_imatrix,l_interval,50> a5(m1,v(50),m2(50,50));
         matrixassign<l_imatrix_slice,l_ivector_slice,l_imatrix,l_interval,50> a6(m1(50,50),v(50),m2);
         matrixassign<l_imatrix,l_imatrix_subv,l_imatrix,l_interval,50> a7(m1,m2[Col(7)],m2);
         matrixassign<l_imatrix_slice,l_imatrix_subv,l_imatrix,l_interval,50> a8(m1(50,50),m2[Col(7)],m2);
         
			matrixaddsub<l_imatrix,l_imatrix,l_imatrix,50> e(m1,m2);
			
			matrixscalarmult<l_imatrix,l_interval,l_imatrix,50> f(m1);
			matrixscalarmult<l_imatrix_slice,l_interval,l_imatrix,40> n(m1(1,40,1,40));
			
			matrixmult<l_imatrix,l_imatrix,l_imatrix,idotprecision,50> o(m1,m2);
			matrixmult<l_imatrix,l_imatrix_slice,l_imatrix,idotprecision,50> p(m1,m2(50,50));
			matrixmult<l_imatrix_slice,l_imatrix_slice,l_imatrix,idotprecision,40> q(m1(1,40,1,40),m2(1,40,1,40));
			// Real ----------------------------------
			matrixconstr<l_imatrix,rvector,rmatrix,real,50> a9(m1,rv,rm2);
			matrixconstr<l_imatrix,rmatrix_subv,rmatrix,real,50> a10(m1,rm2[1],rm2);
         
         matrixassign<l_imatrix,rvector,rmatrix,real,50> a11(m1,rv,rm2);
         matrixassign<l_imatrix_slice,rvector,rmatrix,real,50> a12(m1(50,50),rv,rm2);
         matrixassign<l_imatrix,rvector_slice,rmatrix,real,50> a13(m1,rv(50),rm2);
         matrixassign<l_imatrix_slice,rvector_slice,rmatrix,real,50> a14(m1(50,50),rv(50),rm2);
         matrixassign<l_imatrix,rmatrix_subv,rmatrix,real,50> a15(m1,rm2[Col(7)],rm2);
         matrixassign<l_imatrix_slice,rmatrix_subv,rmatrix,real,50> a16(m1(50,50),rm2[Col(7)],rm2);
         
			// l_real ----------------------------------
			matrixconstr<l_imatrix,l_rvector,l_rmatrix,l_real,50> a25(m1,cv,cm2);
			matrixconstr<l_imatrix,l_rmatrix_subv,l_rmatrix,l_real,50> a26(m1,cm2[1],cm2);
         
         matrixassign<l_imatrix,l_rvector,l_rmatrix,l_real,50> a27(m1,cv,cm2);
         matrixassign<l_imatrix_slice,l_rvector,l_rmatrix,l_real,50> a28(m1(50,50),cv,cm2);
         matrixassign<l_imatrix,l_rvector_slice,l_rmatrix,l_real,50> a29(m1,cv(50),cm2);
         matrixassign<l_imatrix_slice,l_rvector_slice,l_rmatrix,l_real,50> a30(m1(50,50),cv(50),cm2);
         matrixassign<l_imatrix,l_rmatrix_subv,l_rmatrix,l_real,50> a31(m1,cm2[Col(7)],cm2);
         matrixassign<l_imatrix_slice,l_rmatrix_subv,l_rmatrix,l_real,50> a32(m1(50,50),cm2[Col(7)],cm2);
         
			// interval ----------------------------------
			matrixconstr<l_imatrix,ivector,imatrix,interval,50> a17(m1,iv,im2);
			matrixconstr<l_imatrix,imatrix_subv,imatrix,interval,50> a18(m1,im2[1],im2);
         
         matrixassign<l_imatrix,ivector,imatrix,interval,50> a19(m1,iv,im2);
         matrixassign<l_imatrix_slice,ivector,imatrix,interval,50> a20(m1(50,50),iv,im2);
         matrixassign<l_imatrix,ivector_slice,imatrix,interval,50> a21(m1,iv(50),im2);
         matrixassign<l_imatrix_slice,ivector_slice,imatrix,interval,50> a22(m1(50,50),iv(50),im2);
         matrixassign<l_imatrix,imatrix_subv,imatrix,interval,50> a23(m1,im2[Col(7)],im2);
         matrixassign<l_imatrix_slice,imatrix_subv,imatrix,interval,50> a24(m1(50,50),im2[Col(7)],im2);
         
			setfail( !e || !f || !n || !o|| !p || !q||!a1||!a2||!a3||!a4|!a5||!a6||!a7||!a8||!a9||!a10||!a11||!a12||!a13||!a14||!a15||!a16||!a17||!a18||!a19||!a20||!a21||!a22||!a23||!a24||!a25||!a26||!a27||!a28||!a29||!a30||!a31||!a32);
//         setfail(!e);
      }
};

#endif

} // namespace cxsc 

