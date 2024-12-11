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

/* CVS $Id: q_sinh.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_SINH_CPP
#define Q_SINH_CPP

#include "fi_lib.hpp" 

namespace fi_lib {

 using cxsc::real;
 
 real q_sinh(real x){
  real absx, h;
  int sgn;
  real res;

  if(NANTEST(x))                                       /* Test: x=NaN */
      res=q_abortnan(INV_ARG,&x,18);
  else {
  if (x<0) {sgn=-1; absx=-x;}
  else     {sgn=1;  absx=x; }
 
  if (absx>q_ex2a)
      res=q_abortr1(OVER_FLOW,&x,18);                 /* Overflow */

  if (absx<2.5783798e-8) res=x;
  else if (absx>=0.662) 
    {
       h=q_ep1(absx);
       res=sgn*0.5*(h-1.0/h);
    }
  else
    {
       h=q_epm1(absx);
       res=sgn*0.5*(h+h/(h+1.0));
    }

  }
  return(res);
 }
}

#endif





