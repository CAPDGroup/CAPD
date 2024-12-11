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

/* CVS $Id: j_sinh.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_SINH_CPP
#define J_SINH_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_sinh(interval x){
  interval res;
  if (Inf(x)==Sup(x))             /* point interval */
    { 
       if (Inf(x)<0)
        {
          if (Inf(x)>-q_minr)
            {
               Inf(res)=q_pred(Inf(x));
               Sup(res)=Inf(x);
            }
          else
            {
              Inf(res)=q_sinh(Inf(x));
              Sup(res)=Inf(res)*q_snhm;
              Inf(res)*=q_snhp;
              if (Sup(res)>Inf(x)) Sup(res)=Inf(x);
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
                Inf(res)=q_sinh(Inf(x));
                Sup(res)=Inf(res)*q_snhp;
                Inf(res)*=q_snhm;
                if (Inf(res)<Inf(x)) Inf(res)=Inf(x);
              }
        }
    }
  else
    {
      if (Inf(x)<0)
        {
          if (Inf(x)>-q_minr)
            Inf(res)=q_pred(Inf(x));
          else
            Inf(res)=q_sinh(Inf(x))*q_snhp;          
        }
      else  /* Inf(x)>=0 */
        {
          if (Inf(x)<q_minr)
            Inf(res)=Inf(x);        /* includes the case Inf(x)=0 */
          else 
            {  
              Inf(res)=q_sinh(Inf(x))*q_snhm;
              if (Inf(res)<Inf(x)) Inf(res)=Inf(x);
            }
        }
      if (Sup(x)<=0)
        {
          if (Sup(x)>-q_minr)
            Sup(res)=Sup(x);        /* includes the case Sup(x)=0 */
          else
            {
               Sup(res)=q_sinh(Sup(x))*q_snhm;
               if (Sup(res)>Sup(x)) Sup(res)=Sup(x);
            }          
        }
      else  /* now Sup(x)>=0 */
        {
          if (Sup(x)<q_minr)
            Sup(res)=q_succ(Sup(x));        
          else 
            Sup(res)=q_sinh(Sup(x))*q_snhp;
        }
    }   
  return(res);
 }

}

#endif





