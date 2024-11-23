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

/* CVS $Id: j_sin.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_SIN_CPP
#define J_SIN_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;
 using cxsc::real;

 interval j_sin(interval x){
  interval res;
  real erg1,erg2,y1,y2;
  long int k1,k2,q1,q2;

  if (Inf(x)==Sup(x))             /* point interval */
    { 
      if ((Inf(x)<-q_sint[2])||(Sup(x)>q_sint[2]))
        { Inf(res)=-1.0; Sup(res)=1.0; }
      else 
        { 
          if ((Inf(x)>=-q_sint[3])&&(Inf(x)<0))
            {
              Inf(res)=Inf(x);
              Sup(res)=q_succ(Inf(x));
            }
          else if ((Inf(x)>=0)&&(Inf(x)<=q_sint[3]))
            {         
              Sup(res)=Inf(x);
              if (Inf(x)==0)
                Inf(res)=0; 
              else
                Inf(res)=q_pred(Inf(x));
            }
          else
            {
              Inf(res)=q_sin(Inf(x));
              if (Inf(res)<0)
                {
                  Sup(res)=Inf(res)*q_sinm;
                  Inf(res)*=q_sinp;
                }
              else
                {
                  Sup(res)=Inf(res)*q_sinp;
                  Inf(res)*=q_sinm;
                }
            }

         }
    }
  else
    {
      if ((Sup(x)-Inf(x))>=(2.0*q_pi)) 
        { Inf(res)=-1.0; Sup(res)=1.0; }
      else if ((Inf(x)<-q_sint[2])||(Sup(x)>q_sint[2]))
        { Inf(res)=-1.0; Sup(res)=1.0; }
      else 
        {
          /* argument reduction infimum */
          erg1=Inf(x)*q_pi2i; 
          if (erg1>0) {k1=CUTINT(erg1+0.5);
                       q1=(CUTINT(erg1))%4; }
          else        {k1=CUTINT(erg1-0.5);
                       q1=((CUTINT(erg1))-1)%4; }
          y1=q_rtrg(Inf(x),k1); 
            
          /* argument reduction supremum */
          erg2=Sup(x)*q_pi2i; 
          if (erg2>0) {k2=CUTINT(erg2+0.5);
                       q2=(CUTINT(erg2))%4;}
          else        {k2=CUTINT(erg2-0.5);
                       q2=((CUTINT(erg2))-1)%4;}
          y2=q_rtrg(Sup(x),k2);
          if (q1<0) q1+=4;
          if (q2<0) q2+=4;

          if (q1==q2) 
          {
            if ((Sup(x)-Inf(x))>=q_pi) 
              {
                Inf(res)=-1.0; Sup(res)=1.0;
              }
            else if ((q1==1) || (q1==2))
              {
                Inf(res)=q_sin1(y2,k2);
                if (Inf(res)<0) Inf(res)*=q_sinp; else Inf(res)*=q_sinm;
                Sup(res)=q_sin1(y1,k1);
                if (Sup(res)<0) Sup(res)*=q_sinm; else Sup(res)*=q_sinp; 
              }
            else if (q1==0)
              {
                if ((Inf(x)>0) & (Inf(x)<=q_sint[3]))
                  Inf(res)=q_pred(Inf(x));
                else
                  Inf(res)=q_sin1(y1,k1)*q_sinm;
                if ((Sup(x)>0) & (Sup(x)<=q_sint[3]))
                  Sup(res)=Sup(x);
                else
                  Sup(res)=q_sin1(y2,k2)*q_sinp;                    
              }
            else   /* q1 == 3 */
              {
                 if ((Inf(x)>=-q_sint[3]) & (Inf(x)<0))
                   Inf(res)=Inf(x);
                 else
                   Inf(res)=q_sin1(y1,k1)*q_sinp;
                 if ((Sup(x)>=-q_sint[3]) & (Sup(x)<0))
                   Sup(res)=q_succ(Sup(x));
                 else
                   Sup(res)=q_sin1(y2,k2)*q_sinm;                   
              }
           }
           else  /* now we have q1<>q2 */
           { 
             if (q1==0)
             {
               if (q2==1)
               {
                 if ((Inf(x)>0) & (Inf(x)<=q_sint[3]))
                   Inf(res)=q_pred(Inf(x));
                 else
                   {
                     erg1=q_sin1(y1,k1);
                     erg2=q_sin1(y2,k2);
                     if (erg1<erg2) Inf(res)=erg1*q_sinm;
                     else           Inf(res)=erg2*q_sinm;
                   }
                 Sup(res)=1.0;
               }
               else if (q2==2)
               {
                 Inf(res)=q_sin1(y2,k2)*q_sinp;
                 Sup(res)=1.0;
               }
               else  /* q2==3 */
               {
                 Inf(res)=-1.0;
                 Sup(res)=1.0;
               }
             }
             else if (q1==1)
             {
               if (q2==0)
               {
                 Inf(res)=-1.0;
                 erg1=q_sin1(y1,k1);
                 erg2=q_sin1(y2,k2);
                 if (erg1>erg2) Sup(res)=erg1*q_sinp; else Sup(res)=erg2*q_sinp;
               }
               else if (q2==2)
               {
                 Inf(res)=q_sin1(y2,k2)*q_sinp;
                 Sup(res)=q_sin1(y1,k1)*q_sinp;
               }
               else  /* q2==3 */
               {
                 Inf(res)=-1.0;
                 Sup(res)=q_sin1(y1,k1)*q_sinp;
               }
             }
             else if (q1==2)
             {
               if (q2==0)
               {
                 Inf(res)=-1.0;
                 if ((Sup(x)>0) & (Sup(x)<=q_sint[3]))
                   Sup(res)=Sup(x);
                 else
                   Sup(res)=q_sin1(y2,k2)*q_sinp;
               }
               else if (q2==1)
               {
                 Inf(res)=-1.0;
                 Sup(res)=1.0;
               }
               else  /* q2==3 */
               {
                 Inf(res)=-1.0;
                 if ((Sup(x)>=-q_sint[3]) & (Sup(x)<0))
                   Sup(res)=q_succ(Sup(x));
                 else
                   {
                     erg1=q_sin1(y1,k1);
                     erg2=q_sin1(y2,k2);
                     if (erg1>erg2) Sup(res)=erg1*q_sinm;
                     else Sup(res)=erg2*q_sinm;
                   }
               }
             }
             else  /* now we have q1==3 */            
             {
               if (q2==0)
               {
                 if ((Inf(x)>=-q_sint[3]) & (Inf(x)<0))
                   Inf(res)=Inf(x);
                 else
                   Inf(res)=q_sin1(y1,k1)*q_sinp;
                 if ((Sup(x)>0) & (Sup(x)<=q_sint[3]))
                   Sup(res)=Sup(x);
                 else
                   Sup(res)=q_sin1(y2,k2)*q_sinp;
               }
               else if (q2==1)
               {
                 if ((Inf(x)>=-q_sint[3]) & (Inf(x)<0))
                   Inf(res)=Inf(x);
                 else
                   Inf(res)=q_sin1(y1,k1)*q_sinp;
                 Sup(res)=1.0;
               }
               else  /* q2==2 */
               {
                 erg1=q_sin1(y1,k1);
                 erg2=q_sin1(y2,k2);
                 if (erg1<erg2) Inf(res)=erg1*q_sinp; else Inf(res)=erg2*q_sinp;
                 Sup(res)=1.0;
               }
             }
           }
        }
    }

  if (Inf(res)<-1.0) Inf(res)=-1.0;
  if (Sup(res)>1.0) Sup(res)=1.0;
  return(res);
 }

} // Namespace

#endif





