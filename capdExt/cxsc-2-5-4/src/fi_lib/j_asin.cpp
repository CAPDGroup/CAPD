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

/* CVS $Id: j_asin.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_ASIN_CPP
#define J_ASIN_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_asin(interval x){
  interval res;
  if (Inf(x)==Sup(x))             /* point interval */
    { 
       if (Inf(x)<0)
        {
          if (Inf(x)>-q_atnt)
            {
               Inf(res)=q_pred(Inf(x));
               Sup(res)=Inf(x);
            }
          else
            {
              Inf(res)=q_asin(Inf(x));
              Sup(res)=Inf(res)*q_csnm;
              Inf(res)*=q_csnp;
              if (Sup(res)>Inf(x)) Sup(res)=Inf(x);
            } 
        }
       else 
         {
           if (Inf(x)<q_atnt)
             {         
               Inf(res)=Inf(x);
               if (Inf(x)==0)
                  Sup(res)=0; 
               else
                  Sup(res)=q_succ(Inf(x));
             }
           else
             {
                Inf(res)=q_asin(Inf(x));
                Sup(res)=Inf(res)*q_csnp;
                Inf(res)*=q_csnm;
                if (Inf(res)<Inf(x)) Inf(res)=Inf(x);
              }
        }
    }
  else
    {
      if (Inf(x)<0)
        {
          if (Inf(x)>-q_atnt)
            Inf(res)=q_pred(Inf(x)); 
          else
            Inf(res)=q_asin(Inf(x))*q_csnp;          
        }
      else  /* Inf(x)>=0 */
        {
          if (Inf(x)<q_atnt)
            Inf(res)=Inf(x);      
          else 
            {
              Inf(res)=q_asin(Inf(x))*q_csnm;
              if (Inf(res)<Inf(x)) Inf(res)=Inf(x);
            }
        }
      if (Sup(x)<=0)
        {
          if (Sup(x)>-q_atnt)
            Sup(res)=Sup(x);
          else
            {
              Sup(res)=q_asin(Sup(x))*q_csnm;
              if (Sup(res)>Sup(x)) Sup(res)=Sup(x);
            }          
        }
      else  /* Sup(x)>0 */
        {
          if (Sup(x)<q_atnt)
            Sup(res)=q_succ(Sup(x));         /* includes the case Sup(x)=0 */       
          else 
            Sup(res)=q_asin(Sup(x))*q_csnp;
        }
    }   
  return(res);
 }
 
} // Namespace

#endif





