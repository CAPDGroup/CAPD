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

/* CVS $Id: j_lg1p.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_LG1P_CPP
#define J_LG1P_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_lg1p(interval x){
  interval res;
  if (Inf(x)==Sup(x))                   /* point interval */
    { 
      Inf(res)=q_lg1p(Inf(x));
      if (Inf(res)>=0) 
        {
          Sup(res)=Inf(res)*q_lgpp;
          Inf(res)*=q_lgpm;
        }
      else
        {
          Sup(res)=Inf(res)*q_lgpm;
          Inf(res)*=q_lgpp;
        }
    }
  else
    {
      Inf(res)=q_lg1p(Inf(x));
      if (Inf(res)>=0)
        Inf(res)*=q_lgpm;
      else
        Inf(res)*=q_lgpp;
      Sup(res)=q_lg1p(Sup(x));
      if (Sup(res)>=0)
        Sup(res)*=q_lgpp;
      else
        Sup(res)*=q_lgpm;
    }   
  return(res);
 }

} // Namespace

#endif





