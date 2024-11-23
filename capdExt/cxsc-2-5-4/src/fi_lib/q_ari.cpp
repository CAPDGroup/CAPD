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

/* CVS $Id: q_ari.cpp,v 1.16 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef Q_ARI_CPP
#define Q_ARI_CPP

#include "fi_lib.hpp"

namespace fi_lib{

 using cxsc::real; 

/* ------------------------------------------------------------------- */
/* --- utilities                                                   --- */
/* ------------------------------------------------------------------- */

 real q_abs(real x)
 { a_diee y;
   y.f=_double(x);
   y.ieee.sign=0;
   x=y.f;
   return x;
 }

 real q_min(real x,real y)
 { return x < y ? x : y; }

 real q_max(real x,real y)
 { return x > y ? x : y; }

 real q_mid(interval x) 
 { return (0.5*(Inf(x)+Sup(x))); }

/* ------------------------------------------------------------------- */
/* --- assignment                                                  --- */
/* ------------------------------------------------------------------- */

 interval eq_ii(interval y)
 { interval res;

  Inf(res) = Inf(y);
  Sup(res) = Sup(y);

  return res;
 }

/* ------------------------------------------------------------------- */

 interval eq_id(real y)
 { interval res;

  Inf(res)=y;
  Sup(res)=y;

  return res;
 }

/* ------------------------------------------------------------------- */
/* --- functions for addition                                      --- */
/* ------------------------------------------------------------------- */

 interval add_ii(interval x, interval y)
 { interval res;

  if (Inf(x)==-Inf(y)) Inf(res)=0; else Inf(res) = q_pred(Inf(x)+Inf(y));
  if (Sup(x)==-Sup(y)) Sup(res)=0; else Sup(res) = q_succ(Sup(x)+Sup(y));

  return res;
 }

/* ------------------------------------------------------------------- */

 interval add_id(interval x, real y)
 { interval res;

  if (Inf(x)==-y) Inf(res)=0; else Inf(res) = q_pred(Inf(x)+y);
  if (Sup(x)==-y) Sup(res)=0; else Sup(res) = q_succ(Sup(x)+y);

  return res;
 }

/* ------------------------------------------------------------------- */

 interval add_di(real x, interval y)
 { interval res;

  if (x==-Inf(y)) Inf(res)=0; else Inf(res) = q_pred(x+Inf(y));
  if (x==-Sup(y)) Sup(res)=0; else Sup(res) = q_succ(x+Sup(y));

  return res;
 }

/* ------------------------------------------------------------------- */
/* --- functions for subtraction                                   --- */
/* ------------------------------------------------------------------- */

 interval sub_ii(interval x, interval y)
 { interval res;

  if (Inf(x)==Sup(y)) Inf(res)=0; else Inf(res) = q_pred(Inf(x)-Sup(y));
  if (Sup(x)==Inf(y)) Sup(res)=0; else Sup(res) = q_succ(Sup(x)-Inf(y));

  return res;
 }

/* ------------------------------------------------------------------- */

 interval sub_id(interval x, real y)
 { interval res;

  if (Inf(x)==y) Inf(res)=0; else Inf(res) = q_pred(Inf(x)-y);
  if (Sup(x)==y) Sup(res)=0; else Sup(res) = q_succ(Sup(x)-y);

  return res;
 }

/* ------------------------------------------------------------------- */

 interval sub_di(real x, interval y)
 { interval res;

  if (x==Sup(y)) Inf(res)=0; else Inf(res) = q_pred(x-Sup(y));
  if (x==Inf(y)) Sup(res)=0; else Sup(res) = q_succ(x-Inf(y));

  return res;
 }

/* ------------------------------------------------------------------- */
/* --- functions for multiplication                                --- */
/* ------------------------------------------------------------------- */

 interval mul_ii(interval x, interval y)         /*  [x] * [y]                */
 { interval res;
  real h;

  if (Inf(x) >=0) {                        /*  0 <= [x]                 */

    if (Inf(y) >=0) {                      /*  0 <= [y]                 */
      h=Inf(x)*Inf(y);
      Inf(res)=(h==0 ? 0 : q_pred(h));  
    } else {                              /*  [y] <= 0  or  0 \in [y]  */
      h=Sup(x)*Inf(y);
      Inf(res)=(Sup(x)==0 ? 0 : q_pred(h));
    } 

    if (Sup(y) <=0) {                      /*  [y] <= 0                 */
      h=Inf(x)*Sup(y);  
      Sup(res)=(h==0 ? 0 : q_succ(h));  
    } else {                              /*  0 <= [y]  or  0 \in [y]  */
      h=Sup(x)*Sup(y);
      Sup(res)=(Sup(x)==0 ? 0 : q_succ(h));  
    }

  } else if (Sup(x)<=0) {                  /*  [x] <= 0                 */

    if (Sup(y)<=0) {                       /*  [y] <= 0                 */
      h=Sup(x)*Sup(y);
      Inf(res)=(h==0 ? 0 : q_pred(h));
    } else                                /*  0 <= [y]  or  0 \in [y]  */
      Inf(res)=q_pred(Inf(x)*Sup(y)); 

    if (Inf(y)>=0) {                       /*  0 <= [y]                 */
      h=Sup(x)*Inf(y);
      Sup(res)=(h==0 ? 0 : q_succ(h));   
    } else                                /*  [y] <= 0  or  0 \in [y]  */
      Sup(res)=q_succ(Inf(x)*Inf(y));

  } else {                                /*  0 \in [x]                */

    if (Inf(y)>=0) {                       /*  0 <= [y]                 */
      Inf(res)=q_pred(Inf(x)*Sup(y));
      Sup(res)=q_succ(Sup(x)*Sup(y));
    } else if (Sup(y)<=0) {                /*  [y] <= 0                 */
      Inf(res)=q_pred(Sup(x)*Inf(y));
      Sup(res)=q_succ(Inf(x)*Inf(y));
    } else {                              /*  0 \in [x], 0 \in [y]     */
      Inf(res)=q_pred( q_min(Inf(x)*Sup(y), Sup(x)*Inf(y)) );
      Sup(res)=q_succ( q_max(Inf(x)*Inf(y), Sup(x)*Sup(y)) );
    }

  }

  return res;
 }

/* ------------------------------------------------------------------- */

 interval mul_id(interval x, real y)
 { interval res;
  real h;

  if (y>0) { 

    h=Inf(x)*y;
    if ((h==0) && (Inf(x)>=0)) 
      Inf(res)=0;
    else 
      Inf(res)=q_pred(h);

    h=Sup(x)*y;
    if ((h==0) && (Sup(x)<=0))
      Sup(res)=0;
    else 
      Sup(res)=q_succ(h); 

  } else if (y<0) {

    h=Sup(x)*y;
    if ((h==0) && (Sup(x)<=0)) 
      Inf(res)=0;
    else 
      Inf(res)=q_pred(h);

    h=Inf(x)*y;
    if ((h==0) && (Inf(x)>=0)) 
      Sup(res)=0;
    else 
      Sup(res)=q_succ(h); 

  } else {  /* y==0 */  
    Inf(res)=0;
    Sup(res)=0;
  }

  return res;
}

/* ------------------------------------------------------------------- */

 interval mul_di(real y, interval x)
 { interval res;
  real h;

  if (y>0) { 

    h=Inf(x)*y;
    if ((h==0) && (Inf(x)>=0)) 
      Inf(res)=0;
    else 
      Inf(res)=q_pred(h);

    h=Sup(x)*y;
    if ((h==0) && (Sup(x)<=0))
      Sup(res)=0;
    else 
      Sup(res)=q_succ(h); 

  } else if (y<0) {

    h=Sup(x)*y;
    if ((h==0) && (Sup(x)<=0)) 
      Inf(res)=0;
    else 
      Inf(res)=q_pred(h);

    h=Inf(x)*y;
    if ((h==0) && (Inf(x)>=0)) 
      Sup(res)=0;
    else 
      Sup(res)=q_succ(h); 

  } else {  /* y==0 */  
    Inf(res)=0;
    Sup(res)=0;
  }

  return res;
}

/* ------------------------------------------------------------------- */
/* --- functions for division                                      --- */
/* ------------------------------------------------------------------- */

 interval div_ii(interval x,interval y)
 { interval res;
  real h;

  if (Inf(y)>0) { 

    if (Inf(x)>=0) {
      h=Inf(x)/Sup(y);     
      Inf(res)=((h==0) ? 0 : q_pred(h));
    } else 
      Inf(res)=q_pred(Inf(x)/Inf(y));

    if (Sup(x)<=0) { 
      h=Sup(x)/Sup(y);     
      Sup(res)=((h==0) ? 0 : q_succ(h));
    } else 
      Sup(res)=q_succ(Sup(x)/Inf(y));

  } else if (Sup(y)<0) {

    if (Sup(x)<=0) {
      h=Sup(x)/Inf(y);     
      Inf(res)=((h==0) ? 0 : q_pred(h));
    } else 
      Inf(res)=q_pred(Sup(x)/Sup(y));

    if (Inf(x)>=0) { 
      h=Inf(x)/Inf(y);     
      Sup(res)=((h==0) ? 0 : q_succ(h));
    } else 
      Sup(res)=q_succ(Inf(x)/Sup(y));
  } else {
    res=x;
    q_abortdivi(DIV_ZERO,&Inf(y),&Sup(y));  /* Error: Division by zero! */
  }

 return res;
 }

/* ------------------------------------------------------------------- */

 interval div_di(real x,interval y)
 { interval res; 
  interval xi;
  
  Inf(xi)=x;
  Sup(xi)=x;
  res=div_ii(xi,y);

  return res; 
 }

/* ------------------------------------------------------------------- */

 interval div_id(interval x, real y)
 { interval res;
  real h;

  if (y>0) {

    h=Inf(x)/y;
    if ((h==0) && (Inf(x)>=0)) 
      Inf(res)=0;
    else 
      Inf(res)=q_pred(h);

    h=Sup(x)/y;
    if ((h==0) && (Sup(x)<=0)) 
      Sup(res)=0;
    else 
      Sup(res)=q_succ(h);

   } else if (y<0) {

    h=Sup(x)/y;
    if ((h==0) && (Sup(x)<=0)) 
      Inf(res)=0;
    else 
      Inf(res)=q_pred(h);

    h=Inf(x)/y;
    if ((h==0) && (Inf(x)>=0)) 
      Sup(res)=0;
    else 
      Sup(res)=q_succ(h);

   } else {   /* y==0 */
     res=x;
     q_abortdivd(DIV_ZERO,&y);             /* Error: Division by zero! */
   }

  return res;
 }

/* ------------------------------------------------------------------- */
/* --- logical operations                                          --- */
/* ------------------------------------------------------------------- */

 int in_di(real x,interval y)
 { 
  if (Inf(y)<=x && x<=Sup(y)) 
    return 1; 
  else 
    return 0;
 }

 int in_ii(interval x,interval y)
 { 
  if (Inf(y)<Inf(x) && Sup(x)<Sup(y)) 
    return 1; 
  else 
    return 0;
 }

 int ieq_ii(interval x,interval y)
 { 
  if (Inf(x)==Inf(y) && Sup(x)==Sup(y)) 
    return 1; 
  else 
    return 0; 
 }

 int is_ii(interval x,interval y)
 { 
  if (Inf(x) < Inf(y) && Sup(x) < Sup(y)) 
    return 1; 
  else 
    return 0; 
 }

 int ig_ii(interval x,interval y)
 { 
  if (Inf(x) > Inf(y) && Sup(x) > Sup(y)) 
    return 1; 
  else 
    return 0; 
 }

 int ise_ii(interval x,interval y)
 { 
  if (Inf(x) <=Inf(y) && Sup(x) <= Sup(y)) 
    return 1; 
  else 
    return 0; 
 }

 int ige_ii(interval x,interval y)
 { 
  if (Inf(x) >= Inf(y) && Sup(x) >= Sup(y)) 
    return 1; 
  else 
    return 0; 
 }

 int dis_ii(interval x,interval y)
 { 
  if (Sup(x)<Inf(y) || Sup(y) < Inf(x)) 
    return 1; 
  else 
    return 0; 
 }


/* ------------------------------------------------------------------- */
/* --- convex hull, intersection, diam                             --- */
/* ------------------------------------------------------------------- */

 interval hull(interval x, interval y) 
 { interval res;

  if (Inf(x) <= Inf(y)) 
    Inf(res) = Inf(x); 
  else 
    Inf(res) = Inf(y);

  if (Sup(x) >= Sup(y)) 
    Sup(res) = Sup(x); 
  else 
    Sup(res) = Sup(y);

  return res;
 }

 interval intsec(interval x, interval y)
 { interval res;
 
  if (Inf(x) >= Inf(y)) 
    Inf(res)=Inf(x); 
  else 
    Inf(res)=Inf(y);

  if (Sup(x) <= Sup(y)) 
    Sup(res)=Sup(x); 
  else 
    Sup(res)=Sup(y);
 
  return res;
 }

 real q_diam(interval x) 
 { real res;

  if (Inf(x)==Sup(x)) 
    res=0; 
  else 
    res=q_succ(Sup(x)-Inf(x));

  return res;
 }

/* ------------------------------------------------------------------- */
/* --- end of file q_ari.c                                         --- */
/* ------------------------------------------------------------------- */

} // Namespace

#endif





