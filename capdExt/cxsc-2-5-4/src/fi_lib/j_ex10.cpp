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

/* CVS $Id: j_ex10.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_EX10_CPP
#define J_EX10_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_ex10(interval x){
  interval res;
  if (Inf(x)==Sup(x))                  /* point interval */
    { 
      if ((Inf(x)>=0) && (Inf(x)<=22) && (CUTINT(Inf(x))==Inf(x)))
        {  
          Inf(res)=q_ex10(Inf(x));
          Sup(res)=Inf(res);
        }
      else if (Inf(x)<=q_extn)
        { 
          Inf(res)=0.0;
          Sup(res)=q_minr;
        }
      else 
        { 
          Inf(res)=q_ex10(Inf(x));
          Sup(res)=Inf(res)*q_e10p;
          Inf(res)*=q_e10m;
        }
    }
    else
    {
      if (Inf(x)<=q_extn) 
        Inf(res)=0.0;
      else if ((CUTINT(Inf(x))==Inf(x)) && (Inf(x)>=0) && (Inf(x)<=22))
        Inf(res)=q_ex10(Inf(x));
      else
        Inf(res)=q_ex10(Inf(x))*q_e10m;
      if (Sup(x)<=q_extn)
        Sup(res)=q_minr;
      else if ((CUTINT(Sup(x))==Sup(x)) && (Sup(x)>=0) && (Sup(x)<=22))
        Sup(res)=q_ex10(Sup(x));
      else  
        Sup(res)=q_ex10(Sup(x))*q_e10p;
    }   
  if (Inf(res)<0.0) Inf(res)=0.0;
  if ((Sup(x)<=0.0) && (Sup(res)>1.0)) Sup(res)=1.0;
  if ((Inf(x)>=0.0) && (Inf(res)<1.0)) Inf(res)=1.0;

  return(res);
 }
 
} // Namespace

#endif





