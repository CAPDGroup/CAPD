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

/* CVS $Id: test.cpp,v 1.25 2014/01/30 17:23:49 cxsc Exp $ */

#undef CXSC_INDEX_CHECK
#define CXSC_INDEX_CHECK 1

#define TEST_REAL
#define TEST_DOTPRECISION            
#define TEST_L_REAL            
#define TEST_INTERVAL
#undef TEST_IDOTPRECISION            
#define TEST_L_INTERVAL               
#define TEST_COMPLEX
#define TEST_CDOTPRECISION
#define TEST_RVECTOR
#undef TEST_RMATRIX              
#undef TEST_IVECTOR
#undef TEST_IMATRIX              
#undef TEST_CVECTOR
#undef TEST_CMATRIX              
#undef TEST_CIVECTOR
#undef TEST_CIMATRIX              
#undef TEST_LRVECTOR
#undef TEST_LRMATRIX              
#undef TEST_LIVECTOR
#undef TEST_LIMATRIX              

#include <iostream>
#include "test.hpp"

using namespace std;
using namespace cxsc;

map<string,int> testclass::protocol; 
int             testclass::verbose=0;      
int             testclass::testno=0;
int             testclass::testclassno=0;

#include <cstdlib>
namespace cxsc {
int                          usual;
map<long,int>                testoptions;
}

int main(int argc,char *argv[])
{
   cout << "C-XSC-Test" << endl;

   if(argc==1)
   {
      cout << "  Usage: ./test [<loglevel> [-times <times>] [-only <section> <test>] [-dont <section> <test>] ]" << endl;
   }
   
   int verbose=0;
   int times=1;
   if(argc>1)
      verbose=atoi(argv[1]);
  
   int i;
   for(i=2;i<argc;i++)
   {
      if(argv[i]==string("-times"))
      {  
         times=atoi(argv[++i]);
      } else
      if(argv[i]==string("-only"))
      {
         int se=atoi(argv[++i]),no=atoi(argv[++i]);
         usual=no3;
         testoptions[se*1000L+no]=yes3;
         cout << "I will only test " << se << "-" << no << endl; 
      } else
      if(argv[i]==string("-dont"))
      {
         int se=atoi(argv[++i]),no=atoi(argv[++i]);
         testoptions[se*1000L+no]=no3;
         cout << "I won't test " << se << "-" << no << endl;
      } 
      
   }
   
   for(i=0;i<times;i++)
   {
      if(times!=1)
         cout << i+1 << ":";
         
      {

         testclass report(verbose);
         bool      failed=false;
                  
         report.reset();
      
         try
         {
            
   
              
#ifdef TEST_REAL             
            test<real>          real_test;
             
            if(!real_test)
               failed=true;
#endif
#ifdef TEST_DOTPRECISION            
            test<dotprecision>  dotprecision_test;

            if(!dotprecision_test)
               failed=true;
#endif
#ifdef TEST_L_REAL            
            test<l_real>        l_real_test;

            if(!l_real_test)
               failed=true;
#endif
#ifdef TEST_INTERVAL
            test<interval>      interval_test;

            if(!interval_test)
               failed=true;
#endif
#ifdef TEST_IDOTPRECISION            
            test<idotprecision> idotprecision_test;

            if(!idotprecision_test)
               failed=true;
#endif
#ifdef TEST_L_INTERVAL               
            test<l_interval>    l_interval_test;

            if(!l_interval_test)
               failed=true;
#endif
#ifdef TEST_COMPLEX
            test<complex>       complex_test;
            if(!complex_test)
               failed=true;
#endif
#ifdef TEST_CDOTPRECISION
            test<cdotprecision> cdotprecision_test;
            if(!cdotprecision_test)
               failed=true;
#endif
#ifdef TEST_CINTERVAL
            test<cinterval>       cinterval_test;
            if(!cinterval_test)
               failed=true;
#endif
#ifdef TEST_CIDOTPRECISION
            test<cidotprecision> cidotprecision_test;
            if(!cidotprecision_test)
               failed=true;
#endif
#ifdef TEST_RVECTOR
            test<rvector>      rvector_test;

            if(!rvector_test)
               failed=true;
#endif
#ifdef TEST_RMATRIX              
            test<rmatrix>       rmatrix_test;

            if(!rmatrix_test)
               failed=true;
#endif               
#ifdef TEST_IVECTOR
            test<ivector>      ivector_test;

            if(!ivector_test)
               failed=true;
#endif
#ifdef TEST_IMATRIX              
            test<imatrix>       imatrix_test;

            if(!imatrix_test)
               failed=true;
#endif               
#ifdef TEST_CVECTOR
            test<cvector>      cvector_test;

            if(!cvector_test)
               failed=true;
#endif
#ifdef TEST_CMATRIX              
            test<cmatrix>       cmatrix_test;

            if(!cmatrix_test)
               failed=true;
#endif               
#ifdef TEST_CIVECTOR
            test<civector>      civector_test;

            if(!civector_test)
               failed=true;
#endif
#ifdef TEST_CIMATRIX              
            test<cimatrix>       cimatrix_test;

            if(!cimatrix_test)
               failed=true;
#endif               
#ifdef TEST_LRVECTOR
            test<l_rvector>      l_rvector_test;

            if(!l_rvector_test)
               failed=true;
#endif
#ifdef TEST_LRMATRIX              
            test<l_rmatrix>       l_rmatrix_test;

            if(!l_rmatrix_test)
               failed=true;
#endif               
#ifdef TEST_LIVECTOR
            test<l_ivector>      l_ivector_test;

            if(!l_ivector_test)
               failed=true;
#endif
#ifdef TEST_LIMATRIX              
            test<l_imatrix>       l_imatrix_test;

            if(!l_imatrix_test)
               failed=true;
#endif               

            if(failed)
               cout << "test failed!" << endl;
            else
               cout << "test was successful!" << endl;
         }
         catch(ERROR_ALL &e)
         {
            cerr<<"Fehler: "<<e.errtext()<<endl;
         }
   
         if(times==1)
            report.report();

         report.reset();
      }
   }
   return 0;
}

