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

/* CVS $Id: q_tan.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_TAN_CPP
#define Q_TAN_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_tan(real x){
  real res;
  long int m,n,k;
  real ysq,y,q,s,c;


  /* Special cases */
  if(NANTEST(x))                                             /* Test: x=NaN */
      res=q_abortnan(INV_ARG,&x,12);
  else {
  if ((x<-q_sint[2])||(x>q_sint[2]))
      res=q_abortr1(INV_ARG,&x,12);               /* abs. argument too big */

  /* Argument reduction */
  if(x==0) res=0; else
    {
      y=x*q_pi2i; 
      if (y>0) k=CUTINT(y+0.5); else k=CUTINT(y-0.5);
      y=q_rtrg(x,k); 
      n=k%4; if (n<0) n+=4; m=n%2; 

     /* Approximation */
     if ((-q_sint[4]<y)&&(y<q_sint[4]))
       if (m==0) res=y;
       else      res=-1/y;
     else {
       ysq=y*y;

       /* Computation sine */
       q=ysq*(((((((q_sins[5]*ysq)+q_sins[4])
         *ysq+q_sins[3])*ysq+q_sins[2])*ysq+q_sins[1])*ysq)+q_sins[0]);
      /* if (n==0) */
         s=y+y*q;
      /* else
         s=-(y+y*q);*/  
 
       /* Computation cosine */
       q=ysq*ysq*(((((((q_sinc[5]*ysq)+q_sinc[4])
          *ysq+q_sinc[3])*ysq+q_sinc[2])*ysq+q_sinc[1])*ysq)+q_sinc[0]);
 
       if (ysq >= q_sint[0])
         c=0.625+(0.375-(0.5*ysq)+q);
       else if (ysq >= q_sint[1])
         c=0.8125+((0.1875-(0.5*ysq))+q);
       else
         c=1.0-(0.5*ysq - q);
     /*  if (n==3) c=-c;*/ 
 
       /* Computation tangens */
       if (m==0) res=s/c;
       else      res=-c/s;
     }
   }
  }

  return(res);
 }

} // Namespace

#endif





