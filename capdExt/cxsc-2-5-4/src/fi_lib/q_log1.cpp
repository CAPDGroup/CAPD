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

/* CVS $Id: q_log1.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_LOG1_CPP
#define Q_LOG1_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

/* --------------------------------------------------------------------- */
/* - Computation of log(x), table lookup method                        - */
/* - We use the idea of Mr. P.T.P. Tang                                - */
/* - Version without argument check                                    - */
/* --------------------------------------------------------------------- */


/* --- Computation of log(x): Range I  --------------------------------- */
 real q_p1l1(int m, real fg, real fk, real y){
     int j;
     real l_lead,l_trail,u,q;
     real v;       

     /* Step 1 */
     j=CUTINT((fg-1.0)*128);              /* floor */
     l_lead =m*q_lgld[128]+q_lgld[j];
     l_trail=m*q_lgtl[128]+q_lgtl[j];
   
     /* Step 2: Approximation  */
     u=(fk+fk)/(y+fg);
     v=u*u;
     q=u*v*(q_lgb[0]+v*q_lgb[1]);

     /* Step 3 */
     return(l_lead+(u+(q+l_trail)));
   }


/* ---  Computation of log(x): Range II --------------------------------- */

 real q_p2l1(real fk){
     real g,q,u,v,u1,f1,f2,u2;

     /* Step 1 */
     g=1/(2+fk);
     u=2*fk*g;
     v=u*u;

     /* Step 2 */
     q=u*v*(q_lgc[0]+v*(q_lgc[1]+v*(q_lgc[2]+v*q_lgc[3])));

     /* Step 3 */
     u1=CUT24(_double(u));                        /* u1 = 24 leading bits of u   */
     f1=CUT24(_double(fk));                       /* f1 = 24 leading bits von fk */
     f2=fk-f1;
     u2=((2*(fk-u1)-u1*f1)-u1*f2)*g;

     /* Step 4 */
     return(u1+(u2+q));
    }

/* --- Version without argument check !!! ----------------------------- */

 real q_log1(real x){
  int m;
  real fg,fk,y,res;

  /* main program with different cases */
    if (x<q_minr)                  /* only positive normalised arguments */
       res=q_abortr1(INV_ARG,&x,6);
    else if (x==1) res=0.0;
    else if ((q_lgt1<x) && (x<q_lgt2)) 
      {
        fk=x-1;
        res=q_p2l1(fk);
      }
    else 
      {
        FREXPO(x,m); 
        m-=1023; 
                                  
        y=x;
        POWER2(y,-m);
        fg=CUTINT(128*y+0.5);          /* exp2(+7)=128       */
        fg=0.0078125*fg;               /* exp2(-7)=0.0078125 */  
        fk=y-fg;

        res=q_p1l1(m,fg,fk,y);
      }

  return(res);

 }  /* function q_log1 */

/* --------------------------------------------------------------------- */
/* - Computation of log1p(x)=log(1+x), table lookup method             - */
/* --------------------------------------------------------------------- */

 real q_l1p1(real x){
   int m;
   real fg,fk,y,t,h,res;

   /* main program with different cases */
     if (x<=-1) 
        res=q_abortr1(INV_ARG,&x,7);
    else if (x==0) res=x;
    else if ((-q_lgt5<x) && (x<q_lgt5)) res=x;
                              /* res=(8*x-ldexp(1.0,-1074))*0.125; */
    else if ((q_lgt3<x) && (x<q_lgt4)) 
      {
        fk=x;
        res=q_p2l1(fk);
      }
    else 
      {
        t=q_lgt6;
        if (x<t) y=1+x;
        else y=x;
        FREXPO(y,m);
        m-=1023;
        
        POWER2(y,-m);
        fg=CUTINT(128*y+0.5);       /* exp2(+7)=128       */
        fg=0.0078125*fg;            /* exp2(-7)=0.0078125 */
        
        if (m<=-2) 
          fk=y-fg;
        else if ((-1<=m) && (m<=52))
          {
            fk=1.0;
            POWER2(fk,-m);
            h=x;
            POWER2(h,-m); 
            fk=(fk-fg)+h;
          }
        else
          {
            fk=1.0;
            POWER2(fk,-m);
            h=x;
            POWER2(h,-m);
            fk=(h-fg)+fk;
          }

        res=q_p1l1(m,fg,fk,y);
      }

  return(res);

  }    /* function q_l1p1 */

} // Namespace

#endif





