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

/* CVS $Id: q_acth.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef Q_ACTH_CPP
#define Q_ACTH_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_acth(real x){
  real absx, res;
 
  if(NANTEST(x))                                       /* Test: x=NaN */
    res=q_abortnan(INV_ARG,&x,25);
  else {
    if (x<0) absx=-x; else absx=x;  
    if (absx<=1.0) 
      res=q_abortr1(INV_ARG,&x,25);
    res=0.5*q_l1p1(2.0/(absx-1.0));
    if (x!=absx) res=-res;
  }

  return(res);
 }

} // Namespace

#endif





