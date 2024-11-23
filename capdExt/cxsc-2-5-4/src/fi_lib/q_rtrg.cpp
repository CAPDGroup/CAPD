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

/* CVS $Id: q_rtrg.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_RTRG_CPP
#define Q_RTRG_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_rtrg(real x, long int k){
  real red,h;
  a_diee r1,r2;

  if ((-512<k)&&(k<512)) {
     r2.f = _double(x-k*(q_pih[0]+q_pih[1]));
     red = q_r2tr(r2.f,k);
  } else {
     r1.f =  _double(x-k*q_pih[0]);
     h    = k*q_pih[1];
     r2.f = _double(r1.f-h);
     if (r1.ieee.expo == r2.ieee.expo ) 
        red = r1.f - ( ((((k*q_pih[6] + k*q_pih[5]) + k*q_pih[4]) 
                         + k*q_pih[3]) + k*q_pih[2]) + h );
     else
        red = q_r2tr(r2.f,k);  
  } 
  return(red);
 }

 real q_r2tr(real r, long int k){
  real red,h;
  a_diee r1,r2;
  
  r2.f = _double(r);
  h    = k*q_pih[2];
  r1.f = _double(r2.f-h);
  if (r1.ieee.expo == r2.ieee.expo ) 
     red = r2.f - ( (((k*q_pih[6] + k*q_pih[5]) + k*q_pih[4]) 
                         + k*q_pih[3]) + h);
  else {
     h    = k*q_pih[3];
     r2.f = _double(r1.f-h);
     if (r1.ieee.expo == r2.ieee.expo ) 
        red = r1.f - ( ((k*q_pih[6] + k*q_pih[5]) + k*q_pih[4]) + h);
     else {
        h    = k*q_pih[4];
        r1.f = _double(r2.f-h);
        if (r1.ieee.expo == r2.ieee.expo ) 
           red = r2.f - ( (k*q_pih[6] + k*q_pih[5]) + h);
        else {
           h    = k*q_pih[5];
           r2.f = _double(r1.f-h);
           if (r1.ieee.expo == r2.ieee.expo ) 
              red = r1.f - (k*q_pih[6] + h);
           else
              red = r2.f - k*q_pih[6];
        }       
     }
  }

  return(red);
 }

} // Namespace

#endif





