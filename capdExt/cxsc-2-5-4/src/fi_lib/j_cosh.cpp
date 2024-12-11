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

/* CVS $Id: j_cosh.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_COSH_CPP
#define J_COSH_CPP

#include "fi_lib.hpp" 

namespace fi_lib {

 using cxsc::interval;

 interval j_cosh(interval x){
  interval res;

  if (Sup(x)<0)
  {
   if (Inf(x)==Sup(x))             /* point interval */
    { 
      Inf(res)=q_cosh(Inf(x));
      Sup(res)=Inf(res)*q_cshp;
      Inf(res)*=q_cshm;
    }
   else
    {
      Inf(res)=q_cosh(Sup(x))*q_cshm;
      Sup(res)=q_cosh(Inf(x))*q_cshp;
    }
    if (Inf(res)<1.0) Inf(res)=1.0;
  }    /* end  if (Sup(x)<0) */
  else if (Inf(x)>0)
  {
   if (Inf(x)==Sup(x))             /* point interval */
    { 
      Inf(res)=q_cosh(Inf(x));
      Sup(res)=Inf(res)*q_cshp;
      Inf(res)*=q_cshm;
    }
   else
    {
      Inf(res)=q_cosh(Inf(x))*q_cshm;
      Sup(res)=q_cosh(Sup(x))*q_cshp;
    }  
    if (Inf(res)<1.0) Inf(res)=1.0;
  }
  else if (-Inf(x)>Sup(x))
  {
      Inf(res)=1.0;
      Sup(res)=q_cosh(Inf(x))*q_cshp;
  }
  else
  {
      Inf(res)=1.0;
      Sup(res)=q_cosh(Sup(x))*q_cshp;
  }
  return(res);
 }

} // Namespace

#endif





