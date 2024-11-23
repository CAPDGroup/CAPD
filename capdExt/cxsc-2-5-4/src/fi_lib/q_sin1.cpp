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

/* CVS $Id: q_sin1.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_SIN1_CPP
#define Q_SIN1_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

 real q_sin1(real x, long int k){
  real res;
  long int m,n;
  real ysq,q;

  /* Special cases  */
  if(NANTEST(x))                                       /* Test: x=NaN */
      res=q_abortnan(INV_ARG,&x,10);
  else {
  if ((x<-q_sint[2])||(x>q_sint[2]))
    res=q_abortr1(INV_ARG,&x,10);           /* abs. argument too big */

    /* Argument reduction */
    n=k%4; if(n<0) n+=4; m=n%2;

    /* Approximation */
    ysq=x*x;
    if (m==0) {        /* Approximation sine-function, scheme of Horner */

      if ((-q_sint[3]<x)&&(x<q_sint[3]))
      {
        if (n==0) res=x; 
        else res=-x;
      }
      else {
        q=ysq*(((((((q_sins[5]*ysq)+q_sins[4])
          *ysq+q_sins[3])*ysq+q_sins[2])*ysq+q_sins[1])*ysq)+q_sins[0]);
        if (n==0)
          res=x+x*q;
        else
          res=-(x+x*q);
      }
    } else {           /* Approximation cosine-function, scheme of Horner */
      q=ysq*ysq*(((((((q_sinc[5]*ysq)+q_sinc[4])
        *ysq+q_sinc[3])*ysq+q_sinc[2])*ysq+q_sinc[1])*ysq)+q_sinc[0]);

      if (ysq >= q_sint[0])
        res=0.625+(0.375-(0.5*ysq)+q);
      else if (ysq >= q_sint[1])
        res=0.8125+((0.1875-(0.5*ysq))+q);
      else
        res=1.0-(0.5*ysq - q);
      if (n==3) res=-res; 
    }  /* end if */

  }

  return(res);
 }

} // Namespace

#endif





