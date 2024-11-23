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

/* CVS $Id: q_atan.cpp,v 1.15 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef Q_ATAN_CPP
#define Q_ATAN_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_atan(real x){
  real res;
  real absx,ym,y,ysq;
  int    ind,sgn;


  if(NANTEST(x))
    res=q_abortnan(INV_ARG,&x,16);
  else {
     
    if (x<0) absx=-x; else absx=x;    
    if (absx<=q_atnt) res=x;
    else {
            if (absx<8) {sgn=1; ym=0;}
            else        {sgn=-1; ym=q_piha; absx=1/absx;}
          
            ind=0;
            while (absx>=q_atnb[ind+1]) ind+=1;
            y=(absx-q_atnc[ind])/(1+absx*q_atnc[ind]);    
            ysq=y*y; 
            res = (y+y*(ysq*(((((q_atnd[5]*ysq+q_atnd[4])
                  *ysq+q_atnd[3])*ysq+q_atnd[2])
                  *ysq+q_atnd[1])*ysq+q_atnd[0])))+q_atna[ind]; 
            if (x<0) res=-(res*sgn+ym);
            else     res= (res*sgn+ym);  
          }
  }

  return(res);
 }

} // Namespace

#endif





