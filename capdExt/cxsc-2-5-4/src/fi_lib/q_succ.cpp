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

/* CVS $Id: q_succ.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_SUCC_CPP
#define Q_SUCC_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

/* ------------------------------------------------------------------- */
/* --- successor of a real value                                 --- */
/* ------------------------------------------------------------------- */

 real q_succ(real y){
  a_diee su;
 
  su.f=_double(y);
  if (su.ieee.sign==0) {   /*  y >= 0 */
    if (su.ieee.expo==2047 &&  su.ieee.mant0==0 && su.ieee.mant1==0 ) 
      su.f=su.f; 
    else { 
      if (su.ieee.mant1==0xffffffff) { 
        su.ieee.mant1=0; 
        if (su.ieee.mant0==1048575) { 
          su.ieee.mant0=0; 
	  su.ieee.expo++;
        } else { 
          su.ieee.mant0++;
        }
      } else { 
        su.ieee.mant1++;
      }
    }
  } else {                  /* y < 0 */
    if (su.ieee.expo==2047 && su.ieee.mant0==0 && su.ieee.mant1==0 ) 
      su.f=su.f; 
    else 
      if (su.ieee.sign==1 && su.ieee.expo==0 && su.ieee.mant0==0 && su.ieee.mant1==0) {
        su.ieee.sign=0;
        su.ieee.mant1=1;
      } else {
        if (su.ieee.mant1==0) { 
          su.ieee.mant1=0xffffffff; 
          if (su.ieee.mant0==0) { 
            su.ieee.mant0=1048575; 
	    su.ieee.expo--;
          } else { 
            su.ieee.mant0--;
          }
        } else { 
          su.ieee.mant1--;
        }
      }
    }

  return su.f;
 }         /* end function q_succ */

} // Namespace

#endif





