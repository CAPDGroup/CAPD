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

/* CVS $Id: j_expm.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_EXPM_CPP
#define J_EXPM_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_expm(interval x){
  interval res;
  if (Inf(x)==Sup(x))                     /* point interval */
    { 
       if (Inf(x)<0)
        {
          if (Inf(x)>-q_minr)
            {
               Inf(res)=Inf(x);
               Sup(res)=q_succ(Inf(x));
            }
          else
            {
              Inf(res)=q_expm(Inf(x));
              Sup(res)=Inf(res)*q_exmm;
              Inf(res)*=q_exmp;
            } 
        }
       else 
         {
           if (Inf(x)<q_minr)
             {         
               Inf(res)=Inf(x);
               if (Inf(x)==0)
                  Sup(res)=0; 
               else
                  Sup(res)=q_succ(Inf(x));
             }
           else
             {
                Inf(res)=q_expm(Inf(x));
                Sup(res)=Inf(res)*q_exmp;
                Inf(res)*=q_exmm;
              }
        }
    }
  else
    {
      if (Inf(x)<=0)
        {
          if (Inf(x)>-q_minr)
            Inf(res)=Inf(x);                 /* includes case Inf(x)=0 */
          else
            Inf(res)=q_expm(Inf(x))*q_exmp;          
        }
      else  /* Inf(x)>0 */
        {
          if (Inf(x)<q_minr)
            Inf(res)=Inf(x);      
          else 
            Inf(res)=q_expm(Inf(x))*q_exmm;
        }
      if (Sup(x)<0)
        {
          if (Sup(x)>-q_minr)
            Sup(res)=q_succ(Sup(x));
          else
            Sup(res)=q_expm(Sup(x))*q_exmm;          
        }
      else  /* Sup(x)>=0 */
        {
          if (Sup(x)<q_minr)
            Sup(res)=q_succ(Sup(x));        
          else 
            Sup(res)=q_expm(Sup(x))*q_exmp;
        }
    }   

  if (Inf(res)<-1.0) Inf(res)=-1.0;
  return(res);
 }

} // Namespace

#endif





