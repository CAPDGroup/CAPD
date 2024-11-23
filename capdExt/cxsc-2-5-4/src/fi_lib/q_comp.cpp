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

/* CVS $Id: q_comp.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef Q_COMP_CPP
#define Q_COMP_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

/* ------------------------------------------------------------------- */
/* --- utilities                                                   --- */
/* ------------------------------------------------------------------- */

 real q_comp(int s, real m, int e)
 {
  a_diee su;

  if (!((s==1)||(s==-1))) { m=s; q_abortr1(INV_ARG,&m,26); }
  if ((e<-1023)||(e>1024)) { m=e; q_abortr1(INV_ARG,&m,26); }
  if ((m<0)||(m>=2)) q_abortr1(INV_ARG,&m,26);
  if ((e!=-1023)&&(m<1)) q_abortr1(INV_ARG,&m,26);
  if (e==-1023) m+=1;   /* hidden-bit */
  su.f=_double(m);
  if (s==1) su.ieee.sign=0;
  else      su.ieee.sign=1;
  su.ieee.expo=e+1023;  

  return su.f;
 }         /* end function q_comp */

/* ------------------------------------------------------------------- */

 real q_cmps(real m, int e)
 {
  a_diee su;

  if ((e<-1023)||(e>1024)) { m=e; q_abortr1(INV_ARG,&m,26); }
  if ((m<=-2)||(m>=2)) q_abortr1(INV_ARG,&m,26);
  if ((e!=-1023)&&(m<1)&&(m>-1)) q_abortr1(INV_ARG,&m,26);
  if (e==-1023) {if (m>=0) m+=1; else m-=1;}   /* hidden-bit */
  su.f=_double(m);
  su.ieee.expo=e+1023;  

  return su.f;
 }         /* end function q_cmps */

/* ------------------------------------------------------------------- */

 int q_sign(real x)
 {
  a_diee su;

  su.f=_double(x);
  if (su.ieee.sign==0) return +1;
  else                 return -1;
 }         /* end function q_sign */

/* ------------------------------------------------------------------- */

 real q_mant(real x)
 {
  a_diee su;

  su.f=_double(x);
  su.ieee.sign=0;
  su.ieee.expo=1023;
  if ((-q_minr<x) && (x<q_minr)) su.f-=1;  /* no hidden-bit */

  return su.f;
 }         /* end function q_mant */

/* ------------------------------------------------------------------- */

 real q_mnts(real x)
 {
  a_diee su;

  su.f=_double(x);
  su.ieee.expo=1023;
  if ((0<=x) && (x<q_minr)) su.f-=1;  /* no hidden-bit */
  if ((-q_minr<x) && (x<0)) su.f+=1;  /* no hidden-bit */

  return su.f;
 }         /* end function q_mnts */

/* ------------------------------------------------------------------- */

 int q_expo(real x)
 {
  a_diee su;

  su.f=_double(x);

  return (su.ieee.expo-1023);
 }         /* end function q_expo */

/* ------------------------------------------------------------------- */

} // Namespace

#endif





