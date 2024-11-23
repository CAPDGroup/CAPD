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

/* CVS $Id: q_atnh.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef Q_ATNH_CPP
#define Q_ATNH_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_atnh(real x){
  real absx, res;

  if(NANTEST(x))                                       /* Test: x=NaN */
      res=q_abortnan(INV_ARG,&x,24);
  else {
  if ((x<=-1.0)||(1.0<=x))       
      res=q_abortr1(INV_ARG,&x,24);
  if (x<0) absx=-x; else absx=x;
  if (absx>=q_at3i) res=0.5*q_log1((1+absx)/(1-absx));
  else              res=0.5*q_l1p1((2*absx)/(1-absx));
  if (x<0) res=-res;
  }

  return(res);
 }
 
} // Namespace

#endif





