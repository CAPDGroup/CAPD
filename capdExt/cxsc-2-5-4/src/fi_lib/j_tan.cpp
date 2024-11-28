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

/* CVS $Id: j_tan.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_TAN_CPP
#define J_TAN_CPP

#include "fi_lib.hpp" 

namespace fi_lib{
 
 using cxsc::interval;
 using cxsc::real;

 interval j_tan(interval x){
  interval res;
  real h1,h2;
  long int k1,k2,q1;


  if ((Inf(x)<-q_sint[2])||(Sup(x)>q_sint[2]))
    res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),12);  /* abs. argument too big */

  if (Inf(x)==Sup(x))                           /* point interval */
  { 
    if ((Inf(x)>=-q_sint[4])&&(Inf(x)<0))
    {
      Inf(res)=q_pred(Inf(x));
      Sup(res)=Inf(x);
    }
    else if ((Inf(x)>=0)&&(Inf(x)<=q_sint[4]))
    {         
      Inf(res)=Inf(x);
      if (Inf(x)==0)
        Sup(res)=0; 
      else
        Sup(res)=q_succ(Inf(x));
    }
    else
    {
      Inf(res)=q_tan(Inf(x));
      if (Inf(res)<0)
      {
        Sup(res)=Inf(res)*q_tanm;
        Inf(res)*=q_tanp;
      }
      else
      {
        Sup(res)=Inf(res)*q_tanp;
        Inf(res)*=q_tanm;
      }
    }
  }
  else                                      /* x is not a point interval */
  {
    h1=Inf(x)*q_pi2i; 
    k1=CUTINT(h1);
    if (k1<0) q1=(k1-1)%2; else q1=k1%2; if (q1<0) q1+=2;

    h2=Sup(x)*q_pi2i;
    k2=CUTINT(h2); 

    if ((k1==k2) || (q1==1)&(k1==k2-1))
    {
      if ((-q_sint[4]<Inf(x))&&(Inf(x)<0))
        Inf(res)=q_pred(Inf(x));
      else if ((0<=Inf(x))&&(Inf(x)<q_sint[4]))
        Inf(res)=Inf(x);
      else
      { 
        Inf(res)=q_tan(Inf(x));
        if (Inf(res)>=0)
          Inf(res)*=q_tanm;
        else
          Inf(res)*=q_tanp;
      }

      if ((-q_sint[4]<Sup(x))&&(Sup(x)<=0))
        Sup(res)=Sup(x);
      else if ((0<Sup(x))&&(Sup(x)<q_sint[4]))
        Sup(res)=q_succ(Sup(x));
      else
      { 
        Sup(res)=q_tan(Sup(x));
        if (Sup(res)>=0)
          Sup(res)*=q_tanp;
        else
          Sup(res)*=q_tanm;
      }
    }
    else                                           /* invalid argument */
    {
      res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),12);
    }
  }   
  

  return(res);
 }

} // Namespace

#endif





