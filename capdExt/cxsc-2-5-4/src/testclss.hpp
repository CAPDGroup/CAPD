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

/* CVS $Id: testclss.hpp,v 1.25 2014/01/30 17:23:49 cxsc Exp $ */

#ifndef _CXSC_TESTCLASS_HPP_INCLUDED
#define _CXSC_TESTCLASS_HPP_INCLUDED

#ifndef CXSC_INDEX_CHECK
#define CXSC_INDEX_CHECK 1
#endif

#include <string>
#include <iostream>
#include <map>


#define CSTRING_PI   "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811,3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284812"
#define CSTRING_PI_2 "1.5707963267948966192313216916397514420985846996875529104874722961539082031431044993140174126710585339910740432566411533235469223047752911158626797040642405,1.5707963267948966192313216916397514420985846996875529104874722961539082031431044993140174126710585339910740432566411533235469223047752911158626797040642406"
#define CSTRING_PI_4 "0.7853981633974483096156608458198757210492923498437764552437361480769541015715522496570087063355292669955370216283205766617734611523876455579313398520321202,0.7853981633974483096156608458198757210492923498437764552437361480769541015715522496570087063355292669955370216283205766617734611523876455579313398520321203"

#define CSTRING_E    "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841675092,2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841675093"
#define CSTRING_EM1  "1.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841675092,1.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190115738341879307021540891499348841675093"

#include "real.hpp"
#include "l_real.hpp"
#include "interval.hpp"
#include "l_interv.hpp"

using namespace std;
namespace cxsc {


#define asusual3  0
#define yes3      1
#define no3       2

extern int                          usual;
extern map<long,int>                testoptions;

class testclass
{
      string whowhat;
      int  fail;
      enum { not_tested=0,good,bad} teststatus;
      static map<string,int> protocol;
      static             int verbose; 
      static             int testno;
      static             int testclassno;      
      int                    internaltestno;
   public:
      testclass(int v=-1) 
      { 
         testclassno++;
         internaltestno=0;
         fail=0;
         if(v!=-1)
            verbose=v; 
      }
      void reset(void)
      {
         testclassno=0;
         testno=0;
         internaltestno=0;
         protocol.clear();
      }
   protected:
      void  testing(const string &a,const string &b) 
      { 
         testno++;
         internaltestno++;
         whowhat=a+": "+b;
         if(verbose>2)
            cerr << "Checking (" << testclassno << "-" << internaltestno << "): " 
                 << whowhat << endl; 
      }
      int   tested(void) 
      {
         // cout << "ID:" << testno*1000L+internaltestno << endl;
         //cout << "testop:" << testoptions[testno*1000L+internaltestno] << endl;
         if(protocol[whowhat]==not_tested
          && !(testoptions[testclassno*1000L+internaltestno]==asusual3 && usual==no3)
          && !(testoptions[testclassno*1000L+internaltestno]==no3)
           )
         {
            if(verbose>1)
               cerr << "Testing (" << testclassno << "-" << internaltestno << "): " 
                    << whowhat << endl; 
            return(0);
         }
         return(1);
      }
      
      bool compare(const real & arg,const real &ist, const real &soll);
      bool compare(const l_real & arg,const l_real &ist, const l_real &soll);      
      
      template <class T>
      bool compare(const T & arg,const T & arg2,const T &ist, const T &soll)
      {
         // arg: just for information if test failed
         // ist: result of computation
         // soll: result that should have been computed
         
         double maxerr=1e-10;
         
         if(abs(ist-soll)<=maxerr)
         {
            /*if(!fail)
               protocol[whowhat]=good;
            else
               protocol[whowhat]=bad;*/
               
            return true;
         } else
         {
            if(verbose>0)
               cout << whowhat << ": ("<< arg<<"," << arg2 << ")->" 
                    << ist << "-" << soll << "=" << ist-soll 
                    << ">" << maxerr << endl;
            fail=1;
            //protocol[whowhat]=bad;
            return false;
         }
         
      }
      
      bool comparei(const interval &,const interval &, const interval &);
      bool comparei(const l_interval &,const l_interval &, const l_interval &);

      template <class T>
      bool comparei(const T & arg,const T & arg2,const T &ist, const T &soll)
      {
         // arg: just for information if test failed
         // ist: result of computation
         // soll: result that should have been computed
         double maxerr=1e-5; // This functions test only functionability
                              // in spite of numerical exactness!
      
         if(soll+T(-maxerr,maxerr)>=ist)
         {
            /*if(!fail)
               protocol[whowhat]=good;
            else
               protocol[whowhat]=bad;*/
               
            return true;
         } else
         {
            if(verbose>0)
               cout << whowhat << ": ("<< arg<<"," << arg2 << ")->" 
                    << ist << "  " << soll << endl;
            fail=1;
            //protocol[whowhat]=bad;
            return false;
         }
         
      }
      
      void ok()        
      {
         if(verbose>2)
            cerr << whowhat << " ok                   " << endl; 
         protocol[whowhat]=good;
      }
      void error()     
      {
         if(verbose>1)
            cerr << whowhat << " FAILED!              " << endl; 
         fail=1;
         protocol[whowhat]=bad; 
      }
      
      void setfail(int a=1) { fail|=a; }

   public:
      void report(void) 
      {
         typedef map<string,int>::const_iterator CI;
         string  result[bad+1]={"not tested","ok","FAILED"};
         int tested=0,errors=0;
         
         
         for(CI p=protocol.begin();p!=protocol.end();++p)
         {
            if(verbose>0 || (verbose==0 && p->second==bad))
               cout << p->first << "\t" << result[p->second] << endl;
            if(p->second==bad)
               errors++;
            tested++;
         }
         cout << "----------------------------" << endl;
         cout << tested << " functions were tested" << endl;
         cout << errors << " errors found" << endl;
      }
      
      int operator!() const { return fail; }
};

bool testclass::compare(const real & arg,const real &ist, const real &soll)
{
   // arg: just for information if test failed
   // ist: result of computation
   // soll: result that should have been computed

   double maxerr=1e-10; // This functions test only functionability
                        // in spite of numerical exactness!

   if(abs(ist-soll)<=maxerr)
   {
      /*if(!fail)
         protocol[whowhat]=good;
      else
         protocol[whowhat]=bad;*/
         
      return true;
   } else
   {
      if(verbose>0)
         cout << whowhat << ": ("<< arg<<")->" 
              << ist << "-" << soll << "=" << ist-soll 
              << ">" << maxerr << endl;
      fail=1;
      // protocol[whowhat]=bad;
      
      return false;
   }
   
}

bool testclass::compare(const l_real & arg,const l_real &ist, const l_real &soll)
{
   // arg: just for information if test failed
   // ist: result of computation
   // soll: result that should have been computed

   double maxerr=1e-10; // This functions test only functionability
                        // in spite of numerical exactness!

   if(abs(ist-soll)<=maxerr)
   {
      /*if(!fail)
         protocol[whowhat]=good;
      else
         protocol[whowhat]=bad;*/
         
      return true;
   } else
   {
      if(verbose>0)
         cout << whowhat << ": ("<< arg<<")->" 
              << ist << "-" << soll << "=" << ist-soll 
              << ">" << maxerr << endl;
      fail=1;
      // protocol[whowhat]=bad;
      
      return false;
   }
   
}

bool testclass::comparei(const interval & arg,const interval & ist, const interval & soll)
{
   // arg: just for information if test failed
   // ist: result of computation
   // soll: result that should have been computed
   double maxerr=1e-5; // This functions test only functionability
                        // in spite of numerical exactness!

   if(soll+interval(-maxerr,maxerr)>=ist)
   {
      return true;
   } else
   {
      if(verbose>0)
         cout << whowhat << ": ("<< arg<<")->" 
              << ist << "  " << soll << endl;
      fail=1;
      
      return false;
   }
   
}

bool testclass::comparei(const l_interval & arg,const l_interval & ist, const l_interval & soll)
{
   // arg: just for information if test failed
   // ist: result of computation
   // soll: result that should have been computed
   double maxerr=1e-5; // This functions test only functionability
                        // in spite of numerical exactness!

   if(soll+l_interval(-maxerr,maxerr)>=ist)
   {
      return true;
   } else
   {
      if(verbose>0)
         cout << whowhat << ": ("<< arg<<")->" 
              << ist << "  " << soll << endl;
      fail=1;
      
      return false;
   }
   
}

} // namespace cxsc 

#endif // _CXSC_TESTCLASS_HPP_INCLUDED
