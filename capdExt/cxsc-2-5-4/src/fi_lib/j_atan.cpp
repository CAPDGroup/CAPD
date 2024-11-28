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

/* CVS $Id: j_atan.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_ATAN_CPP
#define J_ATAN_CPP

#include "fi_lib.hpp" 

namespace fi_lib{
using cxsc::interval;

 interval j_atan(interval x){
  interval res;
  if (Inf(x)==Sup(x))             /* point interval */
    { 
       if (Inf(x)<0)
        {
          if (Inf(x)>-q_atnt)
            {
               Inf(res)=Inf(x);
               Sup(res)=q_succ(Inf(x));
            }
          else
            {
              Inf(res)=q_atan(Inf(x));
              Sup(res)=Inf(res)*q_ctnm;
              Inf(res)*=q_ctnp;
              if (Inf(res)<Inf(x)) Inf(res)=Inf(x);
            } 
        }
       else 
         {
           if (Inf(x)<q_atnt)
             {         
               Sup(res)=Inf(x);
               if (Inf(x)==0)
                  Inf(res)=0; 
               else
                  Inf(res)=q_pred(Inf(x));
             }
           else
             {
                Inf(res)=q_atan(Inf(x));
                Sup(res)=Inf(res)*q_ctnp;
                Inf(res)*=q_ctnm;
                if (Sup(res)>Inf(x)) Sup(res)=Inf(x);
             }
        }
    }
  else
    {
      if (Inf(x)<=0)
        {
          if (Inf(x)>-q_atnt)
            Inf(res)=Inf(x);         /* includes the case Inf(x)=0 */
          else
            {
              Inf(res)=q_atan(Inf(x))*q_ctnp;
              if (Inf(res)<Inf(x)) Inf(res)=Inf(x);
            }          
        }
      else  /* now Inf(x)>0 */
        {
          if (Inf(x)<q_atnt)
            Inf(res)=q_pred(Inf(x));      
          else 
            Inf(res)=q_atan(Inf(x))*q_ctnm;
        }
      if (Sup(x)<0)
        {
          if (Sup(x)>-q_atnt)
            Sup(res)=q_succ(Sup(x));
          else
            Sup(res)=q_atan(Sup(x))*q_ctnm;          
        }
      else  /* Sup(x)>=0 */
        {
          if (Sup(x)<q_atnt)
            Sup(res)=Sup(x);          /* includes the case Sup(x)=0  */       
          else
            { 
              Sup(res)=q_atan(Sup(x))*q_ctnp;
              if (Sup(res)>Sup(x)) Sup(res)=Sup(x);
            }
        }
    }   
  return(res);
 }

} // Namespace

#endif





