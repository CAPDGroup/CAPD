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

/* CVS $Id: q_epm1.cpp,v 1.15 2014/01/30 17:23:55 cxsc Exp $ */

#ifndef Q_EPM1_CPP
#define Q_EPM1_CPP

#include "fi_lib.hpp" 

namespace fi_lib{

 using cxsc::real;

/* --------------------------------------------------------------------- */
/* - Computation of exp(x)-1, table lookup method                      - */
/* - We use the idea of Mr. P.T.P. Tang                                - */
/* --------------------------------------------------------------------- */

/* ------   Prozedure for Range 1   ------------------------------------ */

 real q_p1e1(real x){
    int j;
    long int n,m;
    real r,r1,r2,q,s;
    real res;

    /* Step 1 */
    if (x>0) n=CUTINT((x*q_exil)+0.5);
    else     n=CUTINT((x*q_exil)-0.5);      /* round (x)       */
    j=n % 32;                               /* n2=n mod 32     */
    if (j<0) j+=32;                         /* We force n2>=0  */
    m=(n-j)/32;
    r1=x-n*q_exl1;
    r2=-(n*q_exl2);

    /* Step 2 */
    r=r1+r2;
    q=(((q_exa[4]*r+q_exa[3])*r+q_exa[2])*r+q_exa[1])*r+q_exa[0];
    q=r*r*q;
    q=r1+(r2+q);

    /* Step 3 */
    s=q_exld[j]+q_extl[j];
      if (m>=53)
	{
          if (m<1023) { res=1.0; POWER2(res,-m); } else res=0.0;
	  res=(q_exld[j]+(s*q+(q_extl[j]-res)));
	  POWER2(res,m);
	}
      else
	{ if (m<=-8)
	    {
	     res=(q_exld[j]+(s*q+q_extl[j]));
	     POWER2(res,m);
	     res-=1;
	    }
	  else
	    {
	      res=1.0;  
	      POWER2(res,-m); 
	      res=((q_exld[j]-res)+(q_exld[j]*q+q_extl[j]*(1+q)));
	      POWER2(res,m);
	    }
      }
    return(res);
  }

/* ------   Prozedure for Range 2   ------------------------------------ */
 real q_p2e1(real x){
    real  u,v,y,z,q;

    /* Step 1 */
    u=(real) (CUT24(_double(x)));
    v=x-u;
    y=u*u*0.5;
    z=v*(x+u)*0.5;

    /* Step 2 */   
    q=(((((((q_exb[8]*x+q_exb[7])*x+q_exb[6])*x+q_exb[5])
	     *x+q_exb[4])*x+q_exb[3])*x+q_exb[2])*x+q_exb[1])*x+q_exb[0];
    q=x*x*x*q;

    /* Schritt 3 */
    if (y>=7.8125e-3)              /* = 2^-7 */
      return ((u+y)+(q+(v+z)) );
    else
      return (x+(y+(q+z)) );
  }

/* ------   Main program with different cases   ----------------------- */

 real q_epm1(real x){
  real fabsx,res;
  
  if (x<0) fabsx=-x; else fabsx=x;
  if (fabsx<q_ext1)
      {
	res = (q_p2h * x + fabsx) * q_p2mh; 
      }
    else
      { if (q_ex2a<x)
	  res=q_abortr1(OVER_FLOW,&x,3);         /* Overflow */
	else
	  { if (x<q_ext3)
	      {
		 res=-1.0+q_p2mh;
	      }
	    else
	      { if (x==0)
		  res=x;
		else
		  { if ((q_ext4<x) && (x<q_ext5))
		      res=q_p2e1(x);
		    else
		      res=q_p1e1(x);
		  }
	      }
	  }
      }
  return(res);
 }

} // Namespace

#endif





