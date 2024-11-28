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

/* CVS $Id: j_sqrt.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef J_SQRT_CPP
#define J_SQRT_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::interval;

 interval j_sqrt(interval x){
  interval res;
  if (Inf(x)==Sup(x))                        /* point interval */
    { 
      if (Inf(x)==0)
        Inf(res)=Sup(res)=0;
      else { Inf(res)=q_sqrt(Inf(x)); 
             Sup(res)=q_succ(Inf(res));
             Inf(res)=q_pred(Inf(res));
           }
    }
  else
    {
      if (Inf(x)==0) Inf(res)=0;
      else Inf(res)=q_pred(q_sqrt(Inf(x)));
      Sup(res)=q_succ(q_sqrt(Sup(x)));
    }   
  return(res);
 }

} // Namespace

#endif





