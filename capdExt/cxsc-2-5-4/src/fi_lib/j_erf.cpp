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

/* CVS $Id: */

/***********************************************************************/  
/* Stand: 18.04.2000                                                   */
/* Autor: cand.math.oec Stefan Traub, IAM, Universitaet Karlsruhe (TH) */    
/***********************************************************************/

#ifndef J_ERF_CPP
#define J_ERF_CPP

#include "fi_lib.hpp" 

namespace fi_lib {

 using cxsc::interval;
 using cxsc::real;

#ifdef erf_kettenbruch

interval j_erf(interval x)
// Intervallfunktion von erf(x);
{
    real a,b,erf_b(0),y1(0),y2(0);
    bool gl;

    a = Inf(x);  b = Sup(x);
    gl = a==b;

// Berechnung einer Oberschranke:
    if (b>0) 
	if (b>=q_erfa1) 
	{
	    erf_b = erf_intv(b);
	    y2 = erf_b * q_erf_p;
	}
	else 
	{ 
	    y2 = q_erfx0 * q_erf_p;
	    y2 = comp(y2,-1021); 
	}
    else 
	if (b<=-q_erfa1) 
	{
	    erf_b = erf_intv(b);
	    y2 = erf_b * q_erf_m;	    
	}
    if (y2>1.0) y2 = 1.0;

// Berechnung einer Unterschranke:
    if (a<0) {
	if (a<=-q_erfa1) {
	    if (gl) y1 = erf_b * q_erf_p;
	    else y1 = erf_intv(a) * q_erf_p;
	} else {
	    y1 = -q_erfx0 * q_erf_p;
	    y1 = comp(y1,-1021);
	}
    } else { 
	if (a>=q_erfa1) {
	    if (gl) y1 = erf_b * q_erf_m;
	    else y1 = erf_intv(a) * q_erf_m;
	}
    }
    if (y1<-1) y1 = -1;	

    return interval(y1,y2);
}

interval j_erfc(interval x)
// Intervallfunktion von erfc(x);
{
    real y1(0.0),y2,a,b,erfc_b(0);
    a = Inf(x);  b = Sup(x);
    bool gl(a==b);
// Berechnung einer Unterschranke y1:
    if (b<=-6) 
       y1 = pred(2.0);
    else {
	if (b<=26.5) {
	    if (b==0) 
	        y1 = 1.0;
	    else 
	    {
		erfc_b = erfc_intv(b);
		y1 = erfc_b * q_erfc_m; // Abrundung
	    }
	}    
    }	    
// Berechnung einer Oberschranke y2:	 
    if (a<=-6) 
        y2 = 2.0;
    else { 
	if (a<=26.5) { 
	    if (a==0) y2 = 1.0;
	    else 
	    {
		if (gl) y2 = erfc_b * q_erfc_p;
		else y2 = erfc_intv(a) * q_erfc_p;
	    }
	} else 
	    y2 = comp(q_erfc_1,-1018);	
    }
    
    return interval(y1,y2);
}

#else
/* ------------------------------------------------------------------- */
/* ----                    the function j_erf                   ------ */
/* ------------------------------------------------------------------- */

 interval j_erf(interval x){
  interval res;
  real y;

  if (Inf(x)==Sup(x))             /* point interval */
    { 
      if (Inf(x)==q_erft[0])
        { Inf(res)=0.0; Sup(res)=0.0; }
      else 
        { 
          if ((Inf(x)>-q_erft[1])&&(Inf(x)<q_erft[0]))
            {
              Inf(res)=-q_minr;
              Sup(res)=0.0;
            }
          else 
	    {
	      if ((Inf(x)>q_erft[0])&&(Inf(x)<q_erft[1]))
		{         
		  Inf(res)=0.0;
		  Sup(res)=q_minr;
		}
	      else 
		{
		  if (Inf(x)<=-q_erft[5])
		    {
		      Inf(res)=-1.0;
		      Sup(res)=-1.0+1e-15;
		    }
		  else
		    {
		      if (Inf(x)>=q_erft[5])
			{
			  Inf(res)=1.0-1e-15;
			  Sup(res)=1.0;
			}
		      else 
			{
			  if (Inf(x)<=-q_erft[1])
			    {
			      y=q_erf(Inf(x));
			      Inf(res)=y*q_erfp;
			      Sup(res)=y*q_erfm;
			    }
			  else
			    {
			      y=q_erf(Inf(x));
			      Inf(res)=y*q_erfm;
			      Sup(res)=y*q_erfp;
			    }  
			}
		    }
		}
	    }
	}
    }
  else
    {
      if (Inf(x)<=-q_erft[5])
	Inf(res)=-1.0;
      else
	{
	  if (Inf(x)<=-q_erft[1])
	    Inf(res)=q_erf(Inf(x))*q_erfp;
	  else
	    {
	      if (Inf(x)<q_erft[0]) 
		Inf(res)=-q_minr;
	      else
		{
		  if (Inf(x)<q_erft[1])
		    Inf(res)=0.0;
		  else
		    {
		      if (Inf(x)<q_erft[5])
			Inf(res)=q_erf(Inf(x))*q_erfm;
		      else
			Inf(res)=1.0-1e-15;
		    }
		}
	    }
	}
      if (Sup(x)<=-q_erft[5])
	Sup(res)=-1.0+1e-15;
      else
	{
	  if (Sup(x)<=-q_erft[1])
	    Sup(res)=q_erf(Sup(x))*q_erfm;
	  else
	    {
	      if (Sup(x)<q_erft[0]) 
		Sup(res)=0.0;
	      else
		{
		  if (Sup(x)<q_erft[1])
		    Sup(res)=q_minr;
		  else 
		    {
		      if (Sup(x)<q_erft[5])
			Sup(res)=q_erf(Sup(x))*q_erfp;
		      else
			Sup(res)=1.0;
		    }
		}
	    }
	}
    }     

  if (Inf(res)<-1.0) Inf(res)=-1.0;
  if (Sup(res)<=-1.0) Sup(res)=-1.0+1e-15;
  if (Sup(res)>1.0) Sup(res)=1.0;
  if (Inf(res)>=1.0) Inf(res)=1.0-1e-15;

  return(res);
 }


/* ------------------------------------------------------------------- */
/* ----                    the function j_erfc                  ------ */
/* ------------------------------------------------------------------- */

 interval j_erfc(interval x){
  interval res;
  real y;

  if (Inf(x)==Sup(x))             /* point interval */
    { 
      if (Inf(x)==q_erft[0])
        { Inf(res)=1.0; Sup(res)=1.0; }
      else 
        { 
          if ((Inf(x)>-q_erft[1])&&(Inf(x)<q_erft[0]))
            {
              Inf(res)=1.0;
              Sup(res)=1.0+1e-15;   
            }
          else 
	    {
	      if ((Inf(x)>q_erft[0])&&(Inf(x)<q_erft[1]))
		{         
		  Inf(res)=1.0-1e-15;
		  Sup(res)=1.0;
		}
	      else 
		{
		  if (Inf(x)>=q_erft[6])
		    {
		      Inf(res)=0.0;
		      Sup(res)=q_minr;
		    }
		  else 
		    {
		      if (Inf(x)<=-q_erft[6])
			{
			  Inf(res)=2.0-1e-15;
			  Sup(res)=2.0;
			}
		      else
			{
			  y=q_erfc(Inf(x));
			  Inf(res)=y*q_efcm;
			  Sup(res)=y*q_efcp;
			}
		    }
		}  
	    }
	}
    }
  else
    {
      if (Inf(x)<=-q_erft[6])
	Sup(res)=2.0;
      else
	{
	  if (Inf(x)<=-q_erft[1])
	    Sup(res)=q_erfc(Inf(x))*q_efcp;
	  else
	    {
	      if (Inf(x)<q_erft[0]) 
		Sup(res)=1.0+1e-15;
	      else
		{
		  if (Inf(x)<q_erft[1])
		    Sup(res)=1.0;
		  else
		    {
		      if (Inf(x)<q_erft[6])
			Sup(res)=q_erfc(Inf(x))*q_efcp;
		      else
			Sup(res)=q_minr;
		    }
		}
	    }
	}
      if (Sup(x)<=-q_erft[6])
	Inf(res)=2.0-1e-15;
      else
	{
	  if (Sup(x)<=-q_erft[1])
	    Inf(res)=q_erfc(Sup(x))*q_efcm;
	  else
	    {
	      if (Sup(x)<q_erft[0]) 
		Inf(res)=1.0;
	      else
		{
		  if (Sup(x)<q_erft[1])
		    Inf(res)=1.0-1e-15;
		  else
		    {
		      if (Sup(x)<q_erft[6])
			Inf(res)=q_erfc(Sup(x))*q_efcm;
		      else
			Inf(res)=0.0;
		    }
		}
	    }
	}
    }     

  if (Inf(res)<0.0) Inf(res)=0.0;
  if (Sup(res)<=0.0) Sup(res)=q_minr;
  if (Sup(res)>2.0) Sup(res)=2.0;  
  if (Inf(res)>=2.0) Inf(res)=2.0-1e-15;

  return(res);
 }
#endif  
  
  
} // Namespace

#endif





