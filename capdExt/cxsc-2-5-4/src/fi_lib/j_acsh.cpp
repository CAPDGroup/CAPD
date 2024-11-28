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

/* CVS $Id: j_acsh.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_ACSH_CPP
#define J_ACSH_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_acsh(interval x){
  interval res;
  if (Inf(x)<1)                    /* Invalid argument */
  {

      res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),23);
  } 
  else
  {
    if (Inf(x)==Sup(x))             /* point interval */
    { 
        if (Inf(x)==1) 
          Inf(res)=Sup(res)=0;
        else
        { 
          Inf(res)=q_acsh(Inf(x));
          Sup(res)=Inf(res)*q_acsp;
          Inf(res)*=q_acsm;
        }
    }
    else
    {
      Inf(res)=q_acsh(Inf(x))*q_acsm;
      Sup(res)=q_acsh(Sup(x))*q_acsp;
    }   
  }
  return(res);
 }

} // Namespace

#endif





