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

/* CVS $Id: j_cot.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_COT_CPP
#define J_COT_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;
 using cxsc::real;

 interval j_cot(interval x){
  interval res;
  real h1,h2;
  long int k1,k2,q1;


  if ((Inf(x)<-q_sint[2])||(Sup(x)>q_sint[2]))
    res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),13);  /* abs. argument too big */

  if (Inf(x)==Sup(x))                           /* point interval */
    { 
      Inf(res)=q_cot(Inf(x));
      if (Inf(res)<0)
        {
          Sup(res)=Inf(res)*q_cotm;
          Inf(res)*=q_cotp;
        }
      else
        {
          Sup(res)=Inf(res)*q_cotp;
          Inf(res)*=q_cotm;
        }
    }
  else if (((Inf(x)<=0)&(Sup(x)>=0))||((Sup(x)<0)&(Sup(x)>-q_minr))
                          ||((Inf(x)>0)&&(Inf(x)<q_minr)))
    res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),13);  /* Singularitaet */
  else
    {
      h1=Inf(x)*q_pi2i;
      k1=CUTINT(h1);
      if (k1<0) q1=(k1-1)%2; else q1=k1%2; if (q1<0) q1+=2;

      h2=Sup(x)*q_pi2i;
      k2=CUTINT(h2); 

      if ((k1==k2) || ((q1==0)&&(k1==k2-1)))
       {
         Inf(res)=q_cot(Sup(x));
         if (Inf(res)>=0)
           Inf(res)*=q_cotm;
         else
           Inf(res)*=q_cotp;
         Sup(res)=q_cot(Inf(x));
         if (Sup(res)>=0)
           Sup(res)*=q_cotp;
         else
           Sup(res)*=q_cotm;
        }
       else                                          /* invalid argument */
        {
          res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),13);
        }
   }   

  return(res);
 }

} // Namespace

#endif





