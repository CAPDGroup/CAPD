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

/* CVS $Id: q_asnh.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef Q_ASNH_CPP
#define Q_ASNH_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_asnh(real x){
  real h, res;
  int neg;

  if(NANTEST(x))                                       /* Test: x=NaN */
      res=q_abortnan(INV_ARG,&x,22);
  else {
  if ((x>-2.5e-8)&&(x<2.5e-8)) res=x;
  else
  {
    if (x<0) { x=-x; neg=1; } else neg=0;
    if (x>1e150) 
      {
        if (neg==1) res=-(q_l2+q_log1(x));
        else        res=  q_l2+q_log1(x);
      }
    else {
      if (x>=1.25)  /* old: x>=0.03125 */
      { 
         if (neg==1) res=-q_log1(x+sqrt(x*x+1));
         else        res= q_log1(x+sqrt(x*x+1)); 
      }
      else
      {
         h=1/x;
         if (neg==1) res=-q_l1p1(x+x/(sqrt(1+h*h)+h));
         else        res= q_l1p1(x+x/(sqrt(1+h*h)+h));
      }
    }  
  }
  }

  return(res);
 }

} // Namespace

#endif





