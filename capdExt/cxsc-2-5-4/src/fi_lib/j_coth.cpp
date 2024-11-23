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

/* CVS $Id: j_coth.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_COTH_CPP
#define J_COTH_CPP

#include "fi_lib.hpp" 

namespace fi_lib{
 using cxsc::interval;

 interval j_coth(interval x){
  interval res;

  if (Sup(x)<0)
  {
   if (Inf(x)==Sup(x))             /* point interval */
    { 
      Inf(res)=q_coth(Inf(x));
      Sup(res)=Inf(res)*q_cthm;
      Inf(res)*=q_cthp;
    }
   else
    {
      Inf(res)=q_coth(Sup(x))*q_cthp;
      Sup(res)=q_coth(Inf(x))*q_cthm;
    } 
   if (Sup(res)>-1) Sup(res)=-1.0;
  }    /* end  if (Sup(x)<0) */
  else if(Inf(x)>0)
  {
   if (Inf(x)==Sup(x))             /* point interval */
    { 
      Inf(res)=q_coth(Inf(x));
      Sup(res)=Inf(res)*q_cthp;
      Inf(res)*=q_cthm;
    }
   else
    {
      Inf(res)=q_coth(Sup(x))*q_cthm;
      Sup(res)=q_coth(Inf(x))*q_cthp;
    }  
   if (Inf(res)<1) Inf(res)=1.0;
  }
  else   /* invalid argument */
  {

   res=q_abortr2(INV_ARG,&Inf(x),&Sup(x),21);

  }
  return(res);
 }

} // Namespace

#endif





