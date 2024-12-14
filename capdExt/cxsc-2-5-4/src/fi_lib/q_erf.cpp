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

/* CVS $Id: q_erf.cpp,v 1.24 2014/01/30 17:23:55 cxsc Exp $ */

/***********************************************************************/  
/* Stand: 18.04.2000                                                   */
/* Autor: cand.math.oec Stefan Traub, IAM, Universitaet Karlsruhe (TH) */    
/***********************************************************************/

#ifndef Q_ERF_CPP
#define Q_ERF_CPP

#include "fi_lib.hpp" 

namespace fi_lib {

 using cxsc::real;

/* ------------------------------------------------------------------- */
/* ----           the function q_expx2 = e^(- x^2)              ------ */
/* ------------------------------------------------------------------- */

 real q_expx2(real x){
  real m, res;
  long int z;

  if (x<0.0) x = -x;
  z = CUTINT(x);    /* ganzzahliger Anteil */
  m = x - z;
  if (m>0.5)
    {
       z = z + 1;
       m = m - 1;   /* m <= 0.5  und  x = z + m */
    }
  res = q_expz[z] * q_exp(-(z+z)*m) * q_exp(-m*m);
  if (z==27) res = res * q_exp2(-64);
  
  return(res);
 }
 
int int_no(real *a, 
           const int n,    // n: Anzahl der Elemente des Feldes a
           const real& x)  // x: Eine real-Zahl
// Ein Intervall [A,B] wird durch ein Feld a mit n Elementen in
// n-1 Teilintervalle unterteilt. Für ein x aus [A,B] wird eine
// Intervall-Nummer zurückgegeben, wobei das erste Teilintervall
// die Nr. 0 und das letzte Teilintervall die Nr. n-2 erhält.
// Für  x<A  ist Nr. = -1, und für  x>=B  gilt  Nr. = n-1; 
{
    int i,j,k;
    i = 0;
    j = n-1;
    do 
    {
	k = (i+j)/2;
	if (x<a[k]) j = k-1;
	else i = k+1;
    } 
    while (i <= j);

    return j;
} 

#ifdef erf_kettenbruch
real erf_a(const real& x)
// Berechnet Naeherungen fuer erf(x) in:  0 <= x <= 0.5;
// Rel. Fehlerschranke in 1.97193e-308 <= x <= 0.5 ist:  4.23454e-16;
{
    real y(0),v;
    int ex = expo(x);
    if ( (x<=q_erfa1) && (ex > -2147483647) ) 
    { 
		 std::cerr << "erf(x) in denormalized range." << std::endl;
		 exit(1);
    }
    else 
	 if (ex<=-30) y = nu_ko_1 * x;
	 else 
         { // Kettenbruch:  x*K_4(v), v = 1/x^2;
	    v = 1/(x*x);
	    y = q_erfa3_a[4] / (  v + q_erfa3_b[4]);
	    y = q_erfa3_a[3] / ( (v + q_erfa3_b[3]) + y );
	    y = q_erfa3_a[2] / ( (v + q_erfa3_b[2]) + y );
	    y = q_erfa3_a[1] / ( (v + q_erfa3_b[1]) + y ) + q_erfa3_b[0];
	    y *= x;
		 y += x;
         }
    return y;
}

real erfa_intv(const real& x)
//  Berechnet fuer die Intervallfunktion Naeherungen fuer 
//  erf(x) im Bereich:     0 <= x <= 0.5;
//  Im Bereich   0 <= x <= 1.97193E-308  wird der Funktionswert auf
//  Null gesetzt, wodurch in diesem Bereich eine Unterschranke 
//  vorgegeben ist.
{
    real y(0.0),v;
    int ex = expo(x);
    if (x>=q_erfa1) 
    {
	 	if (ex<=-30) y = nu_ko_1 * x;
	 	else 
                { // Kettenbruch:  x*K_4(v), v = 1/x^2;
	    	v = 1/(x*x);
	    	y = q_erfa3_a[4] / (  v + q_erfa3_b[4]);
	    	y = q_erfa3_a[3] / ( (v + q_erfa3_b[3]) + y );
	    	y = q_erfa3_a[2] / ( (v + q_erfa3_b[2]) + y );
	    	y = q_erfa3_a[1] / ( (v + q_erfa3_b[1]) + y ) + q_erfa3_b[0];
	    	y *= x;
			y += x;
                }
   }
   return y;
}

real erfc_a(const real& x)
// Berechnung von erfc(x) für x in [0,0.5]
// Fehlerschranke in diesem Bereich:  7.477839E-16;
{
    real y;
    if (x < 1e-18) y = 1.0;  // Dies ist OK !
    else y = 1.0 - erf_a(x);

    return y;
}

real erf_b(const real& x)
// Berechnet Naeherungen fuer erf(x) in:  0.5 <= x <= 1.125;
// Rel. Fehlerschranke in diesem Bereich: 4.4262e-16;
{
    real y,v;
    // Kettenbruch:  K_7(v), v = 1/(x-x0);  x0 = 45/64 = 0.703125
    if (x==q_erfB_x0) y = q_erfb_b[0]; // Auswertefehler = 0.0
    else
    {
		v = 1/(x-q_erfB_x0);
		y = q_erfb_a[7] / (  v + q_erfb_b[7]);
		y = q_erfb_a[6] / ( (v + q_erfb_b[6]) + y );
		y = q_erfb_a[5] / ( (v + q_erfb_b[5]) + y );
		y = q_erfb_a[4] / ( (v + q_erfb_b[4]) + y );
		y = q_erfb_a[3] / ( (v + q_erfb_b[3]) + y );
		y = q_erfb_a[2] / ( (v + q_erfb_b[2]) + y ); 
		y = q_erfb_a[1] / ( (v + q_erfb_b[1]) + y ) + q_erfb_b[0];
    }
    y += 1.0;
    times2pown(y,-1); // Division durch 2
	 
    return y;
}

real erf_c(const real& x)
// Berechnet Näherungen für erf(x) in:  1.125 <= x <= 1.75;
// Rel. Fehlerschranke in diesem Bereich: 3.0055e-16;
{
    real y,v;
    // Kettenbruch:  K_6(v), v = 1/(x-x0);  x0 = 11/8 = 1.375
    if (x==q_erfC_x0) y = q_erfc_b[0];  // Auswertefehler = 0.0
    else
    {
	v = 1/(x-q_erfC_x0);
	y = q_erfc_a[6] / (  v + q_erfc_b[6]);
	y = q_erfc_a[5] / ( (v + q_erfc_b[5]) + y );
	y = q_erfc_a[4] / ( (v + q_erfc_b[4]) + y );
	y = q_erfc_a[3] / ( (v + q_erfc_b[3]) + y );
	y = q_erfc_a[2] / ( (v + q_erfc_b[2]) + y ); 
	y = q_erfc_a[1] / ( (v + q_erfc_b[1]) + y ) + q_erfc_b[0];
    }
	
    return y;
}


real erf_d(const real& x)
// Berechnet Näherungen für erf(x) in:  1.75 <= x <= 2.375;
// Rel. Fehlerschranke in diesem Bereich: 3.2284e-16;
{
    real y,v;
    // Kettenbruch:  K_6(v), v = 1/(x-x0);  x0 = 2.0
    if (x==q_erfD_x0) y = q_erfd_b[0];  // Auswertefehler = 0.0
    else
    {
	v = 1/(x-q_erfD_x0);
	y = q_erfd_a[6] / (  v + q_erfd_b[6]);
	y = q_erfd_a[5] / ( (v + q_erfd_b[5]) + y );
	y = q_erfd_a[4] / ( (v + q_erfd_b[4]) + y );
	y = q_erfd_a[3] / ( (v + q_erfd_b[3]) + y );
	y = q_erfd_a[2] / ( (v + q_erfd_b[2]) + y ); 
	y = q_erfd_a[1] / ( (v + q_erfd_b[1]) + y ) + q_erfd_b[0];
    }
	
    return y;
}

real erf_e(const real& x)
// Berechnet Näherungen für erf(x) in:  2.375 <= x <= 3.0;
// Rel. Fehlerschranke in diesem Bereich: 2.8806e-16;
{
    real y,v;
    // Kettenbruch:  K_6(v), v = 1/(x-x0);  x0 = 2.5;
    if (x==q_erfE_x0) y = q_erfe_b[0];  // Auswertefehler = 0.0
    else
    {
	v = 1/(x-q_erfE_x0);
	y = q_erfe_a[6] / (  v + q_erfe_b[6]);
	y = q_erfe_a[5] / ( (v + q_erfe_b[5]) + y );
	y = q_erfe_a[4] / ( (v + q_erfe_b[4]) + y );
	y = q_erfe_a[3] / ( (v + q_erfe_b[3]) + y );
	y = q_erfe_a[2] / ( (v + q_erfe_b[2]) + y ); 
	y = q_erfe_a[1] / ( (v + q_erfe_b[1]) + y ) + q_erfe_b[0];
    }
	
    return y;
}

real erf_f(const real& x)
// Berechnet Näherungen für erf(x) in:  3 <= x <= 3.75;
// Rel. Fehlerschranke in diesem Bereich: 3.0008e-16;
{
    real y,v;
    // Kettenbruch:  K_6(v), v = 1/(x-x0);  x0 = 3.5;
    if (x==q_erfF_x0) y = q_erff_b[0];  // Auswertefehler = 0.0
    else
    {
	v = 1/(x-q_erfF_x0);
	y = q_erff_a[6] / (  v + q_erff_b[6]);
	y = q_erff_a[5] / ( (v + q_erff_b[5]) + y );
	y = q_erff_a[4] / ( (v + q_erff_b[4]) + y );
	y = q_erff_a[3] / ( (v + q_erff_b[3]) + y );
	y = q_erff_a[2] / ( (v + q_erff_b[2]) + y ); 
	y = q_erff_a[1] / ( (v + q_erff_b[1]) + y ) + q_erff_b[0];
    }
	
    return y;
}

real erf_g(const real& x)
// Berechnet Näherungen für erf(x) in:  3.75 <= x <= 4.75;
// Rel. Fehlerschranke in diesem Bereich: 2.4453e-16;
{
    real y,v;
    // Kettenbruch:  K_6(v), v = 1/(x-x0);  x0 = 4.25;
    if (x==q_erfG_x0) y = q_erfg_b[0];  // Auswertefehler = 0.0
    else
    {
	v = 1/(x-q_erfG_x0);
	y = q_erfg_a[6] / (  v + q_erfg_b[6]);
	y = q_erfg_a[5] / ( (v + q_erfg_b[5]) + y );
	y = q_erfg_a[4] / ( (v + q_erfg_b[4]) + y );
	y = q_erfg_a[3] / ( (v + q_erfg_b[3]) + y );
	y = q_erfg_a[2] / ( (v + q_erfg_b[2]) + y ); 
	y = q_erfg_a[1] / ( (v + q_erfg_b[1]) + y ) + q_erfg_b[0];
    }
	
    return y;
}

real erf_h(const real& x)
// Berechnet Näherungen für erf(x) in:  4.75 <= x <= 6;
// Rel. Fehlerschranke in diesem Bereich: 2.4303e-16;
{
    real y,v;
    // Kettenbruch:  K_5(v), v = 1/(x-x0);  x0 = 5.375;
    if (x==q_erfH_x0) y = q_erfh_b[0];  // Auswertefehler = 0.0
    else
    {
	v = 1/(x-q_erfH_x0);
	y = q_erfh_a[5] / (  v + q_erfh_b[5]);
	y = q_erfh_a[4] / ( (v + q_erfh_b[4]) + y );
	y = q_erfh_a[3] / ( (v + q_erfh_b[3]) + y );
	y = q_erfh_a[2] / ( (v + q_erfh_b[2]) + y ); 
	y = q_erfh_a[1] / ( (v + q_erfh_b[1]) + y ) + q_erfh_b[0];
    }
	
    return y;
}

real erfc_b(const real& x)
// Berechnet Näherungen für erf(x) in:  0.5 <= x <= 1.125;
// Rel. Fehlerschranke in diesem Bereich: 8.4078e-16;
{
    real y,v;
    // Kettenbruch:  K_7(v), v = 1/(x-x0) = 1/(x-0.90625);

    if (x==q_erfcB_x0) y = q_erfcb_b[0]; // Auswertefehler = 0.0
    else 
    {
	v = 1/(x-q_erfcB_x0);
	y = q_erfcb_a[7] / (  v + q_erfcb_b[7]);
	y = q_erfcb_a[6] / ( (v + q_erfcb_b[6]) + y );
	y = q_erfcb_a[5] / ( (v + q_erfcb_b[5]) + y );
	y = q_erfcb_a[4] / ( (v + q_erfcb_b[4]) + y ); 
	y = q_erfcb_a[3] / ( (v + q_erfcb_b[3]) + y ); 
	y = q_erfcb_a[2] / ( (v + q_erfcb_b[2]) + y ); 
	y = q_erfcb_a[1] / ( (v + q_erfcb_b[1]) + y ) + q_erfcb_b[0];
    }

    return y;
}

real erfc_c(const real& x)
// Berechnet Näherungen für erfc(x) in:  1.125 <= x <= 1.5;
// Rel. Fehlerschranke in diesem Bereich: 7.9198e-16;
{
    real y,v;
    // Kettenbruch:  K_6(v), v = 1/(x-x0) = 1/(x-1.375);

    if (x==q_erfcC_x0) y = q_erfcc_b[0]; 
    else 
    {
	v = 1/(x-q_erfcC_x0);
	y = q_erfcc_a[6] /   (v + q_erfcc_b[6]); 
	y = q_erfcc_a[5] / ( (v + q_erfcc_b[5]) + y );
	y = q_erfcc_a[4] / ( (v + q_erfcc_b[4]) + y );
	y = q_erfcc_a[3] / ( (v + q_erfcc_b[3]) + y ); 
	y = q_erfcc_a[2] / ( (v + q_erfcc_b[2]) + y ); 
	y = q_erfcc_a[1] / ( (v + q_erfcc_b[1]) + y ) + q_erfcc_b[0];
    }

    return y;
}

real erfc_d(const real& x)
// Berechnet Näherungen für erfc(x) in:  1.5 <= x <= 2;
// Rel. Fehlerschranke in diesem Bereich: 1.620e-15;
{
    real y,v;
    // Kettenbruch:  K_7(v), v = 1/(x-x0) = 1/(x-1.8125);

    if (x==q_erfcD_x0) y = q_erfcd_b[0]; 
    else 
    {
	v = 1/(x-q_erfcD_x0);
	y = q_erfcd_a[7] /   (v + q_erfcd_b[7]);
	y = q_erfcd_a[6] / ( (v + q_erfcd_b[6]) + y ); 
	y = q_erfcd_a[5] / ( (v + q_erfcd_b[5]) + y );
	y = q_erfcd_a[4] / ( (v + q_erfcd_b[4]) + y );
	y = q_erfcd_a[3] / ( (v + q_erfcd_b[3]) + y ); 
	y = q_erfcd_a[2] / ( (v + q_erfcd_b[2]) + y );
	y = q_erfcd_a[1] / ( (v + q_erfcd_b[1]) + y ) + q_erfcd_b[0];
    }

    return y;
}

real erfc_e(const real& x)
// Berechnet Näherungen für erfc(x) in:  2 <= x <= 3.5;
// Rel. Fehlerschranke in diesem Bereich E: 1.2109e-15;
{
    real y,v;
    // Kettenbruch:  K_5(v), v = 1/(x-x0) = 1/(x-2.75);

    if (x == q_erfcE_x0) y = q_erfce_b[0]; 
    else 
    {
	v = 1/(x-q_erfcE_x0);
	y = q_erfce_a[5] /   (v + q_erfce_b[5]);
	y = q_erfce_a[4] / ( (v + q_erfce_b[4]) + y );
	y = q_erfce_a[3] / ( (v + q_erfce_b[3]) + y ); 
	y = q_erfce_a[2] / ( (v + q_erfce_b[2]) + y );
	y = q_erfce_a[1] / ( (v + q_erfce_b[1]) + y ) + q_erfce_b[0];
    }
    y = y * expmx2(x);
    return y;
}

real erfc_f(const real& x)
// Berechnet Näherungen für erfc(x) in:  2 <= x <= 3.5;
// Rel. Fehlerschranke in diesem Bereich E: 1.2109e-15;
{
    real y,v;
    // Kettenbruch:  K_5(v), v = 1/(x-x0) = 1/(x-2.75);

    if (x == q_erfcF_x0) y = q_erfcf_b[0]; 
    else 
    {
	v = 1/(x-q_erfcF_x0);
	y = q_erfcf_a[5] /   (v + q_erfcf_b[5]);
	y = q_erfcf_a[4] / ( (v + q_erfcf_b[4]) + y );
	y = q_erfcf_a[3] / ( (v + q_erfcf_b[3]) + y ); 
	y = q_erfcf_a[2] / ( (v + q_erfcf_b[2]) + y );
	y = q_erfcf_a[1] / ( (v + q_erfcf_b[1]) + y ) + q_erfcf_b[0];
    }
    y = y * expmx2(x);
    return y;
}

real erfc_g(const real& x)
// Berechnet Näherungen für erfc(x) in:  2 <= x <= 3.5;
// Rel. Fehlerschranke in diesem Bereich E: 1.2109e-15;
{
    real y,v;
    // Kettenbruch:  K_5(v), v = 1/(x-x0) = 1/(x-2.75);

    if (x == q_erfcG_x0) y = q_erfcg_b[0]; 
    else 
    {
	v = 1/(x-q_erfcG_x0);
	y = q_erfcg_a[5] /   (v + q_erfcg_b[5]);
	y = q_erfcg_a[4] / ( (v + q_erfcg_b[4]) + y );
	y = q_erfcg_a[3] / ( (v + q_erfcg_b[3]) + y ); 
	y = q_erfcg_a[2] / ( (v + q_erfcg_b[2]) + y );
	y = q_erfcg_a[1] / ( (v + q_erfcg_b[1]) + y ) + q_erfcg_b[0];
    }
    y = y * expmx2(x);
    return y;
}


real erfc_h(const real& x)
// Berechnet Näherungen für erfc(x) in:  14 <= x <= 26.5;
// Rel. Fehlerschranke in diesem Bereich H: 1.3109e-15;
{
    real y,v;
    // Kettenbruch:  K_4(v), v = 1/(x-x0) = 1/(x-20.5);

    if (x == q_erfcH_x0) y = q_erfch_b[0]; 
    else 
    {
	v = 1/(x-q_erfcH_x0);
	y = q_erfch_a[4] /   (v + q_erfch_b[4]);
	y = q_erfch_a[3] / ( (v + q_erfch_b[3]) + y );
	y = q_erfch_a[2] / ( (v + q_erfch_b[2]) + y );
	y = q_erfch_a[1] / ( (v + q_erfch_b[1]) + y ) + q_erfch_b[0];
    }
    y *= expmx2(x);
    return y;
}

real a_erf[9] = { 0,
		  0.5,
		  1.125,
		  1.75,
		  2.375,
		  3.0,
		  3.75,
		  4.75,
		  6 };  // Die a_erf[j] sind die Intervallrandpunkte.

real erf_pos(const real& x)
// Berechnet erf(x) für x >= 0;
{
    int Nr=int_no(a_erf,9,x);
    real y;

    switch(Nr)
    {
	case 0: y = erf_a(x); break;
	case 1: y = erf_b(x); break;
	case 2: y = erf_c(x); break;
	case 3: y = erf_d(x); break;
	case 4: y = erf_e(x); break;
	case 5: y = erf_f(x); break;
	case 6: y = erf_g(x); break;
	case 7: y = erf_h(x); break;
	// case 8
	default: y = 1.0;
    }

    return y;
}

real q_erf(real x)
// Berechnung von erf(x) für alle Punktargumente x in IR;
{
    real y;
    int singx( sign(x) );

    if (singx>=0) y = erf_pos(x);
    else y = -erf_pos(-x);

    return y;
}

real erf_pos_intv(const real& x)
// Berechnet erf(x) für 0 <= x < +MaxReal zur Implementierung
// einer Intervallfunktion bez. erf(x).
{
    int Nr = int_no(a_erf,9,x);
    real y;

    switch(Nr)
    {
	case 0 : y = erfa_intv(x); break;
	case 1 : y = erf_b(x);      break;
	case 2 : y = erf_c(x);      break;
	case 3 : y = erf_d(x);      break;
	case 4 : y = erf_e(x);      break;
	case 5 : y = erf_f(x);      break;
	case 6 : y = erf_g(x);      break;
	case 7 : y = erf_h(x);      break;
	// case 8 
	default: y = 1.0;
    }

    return y;
}

real erf_intv(const real& x)
// Berechnet erf(x) für Punktargumente zur Implementierung
// einer Intervallfunktion bez. erf(x).
{
    real y;
    int singx( sign(x) );

    if (singx>=0) y = erf_pos_intv(x);
    else y = -erf_pos_intv(-x);
    return y;
}

real erfc_j(const real& x)
// Berechnet Näherungen für erfc(x) in:  -6 <= x <= 0;
// Rel. Fehlerschranke in diesem Bereich J:  8.5881E-16;
{
    real y(1.0);

    if (x <= -q_erfa1) 
	y = 1 - erf(x); 
    return y;
}

real a_erfc[10] = {-6,0,0.5,1.125,1.5,2.0,3.5,7.0,14.0,succ(26.5)};  
// Die a_erfc[j] sind die Intervallrandpunkte.

real q_erfc(real x)
// Berechnet erfc(x) für -MaxReal <= x <= +26.5;
{
    int Nr = int_no(a_erfc,10,x);
    real y;

    switch(Nr)
    {
	case -1 : y = 2.0;        break;
	case 0  : y = erfc_j(x);  break;
	case 1  : y = erfc_a(x);  break;
	case 2  : y = erfc_b(x);  break;
	case 3  : y = erfc_c(x);  break;
	case 4  : y = erfc_d(x);  break;
	case 5  : y = erfc_e(x);  break;
	case 6  : y = erfc_f(x);  break;
	case 7  : y = erfc_g(x);  break;
	case 8  : y = erfc_h(x);  break;
	// case 9  
	default : y=x; 
	          std::cerr << "erfc(x) probably in denormalized range." << std::endl;
	          exit(1);
    }

    return y;
}

real erfc_intv(const real& x)
// Berechnet erfc(x) für  -MaxReal <= x <= +MaxReal
// für die Intervallfunktion bez.  erfc(x);
{
    int Nr = int_no(a_erfc,10,x);
    real y;

    switch(Nr)
    {
	case -1 : y = 2.0;        break;
	case 0  : y = erfc_j(x);  break;
	case 1  : y = erfc_a(x);  break;
	case 2  : y = erfc_b(x);  break;
	case 3  : y = erfc_c(x);  break;
	case 4  : y = erfc_d(x);  break;
	case 5  : y = erfc_e(x);  break;
	case 6  : y = erfc_f(x);  break;
	case 7  : y = erfc_g(x);  break;
	case 8  : y = erfc_h(x);  break;
	// case 9  
	default : y = 0.0;
    }

    return y;
}

#else

/* ------------------------------------------------------------------- */
/* ----                    the function q_erf                   ------ */
/* ------------------------------------------------------------------- */

 real q_erf(real x){
  real p, q, res, x2;
  
  if(NANTEST(x)) res = q_abortnan(INV_ARG,&x,27);        /* Test: x=NaN */ 
  else 
  {if (x==q_erft[0]) res = 0.0;
   else 
   {if (x<q_erft[0]) res = -q_erf(-x);
    else 
    {if (x<q_erft[1]) res = q_abortnan(INV_ARG,&x,27); /* UNDERFLOW */  
     else 
     {if (x<q_erft[2]) res = x * q_epA2[0];
      else 
      {if (x<q_erft[3]) 
         { x2 = x*x;
           p = (((q_epA2[4]*x2+q_epA2[3])*x2+q_epA2[2])*x2+q_epA2[1])*x2+q_epA2[0];  
           q = (((q_eqA2[4]*x2+q_eqA2[3])*x2+q_eqA2[2])*x2+q_eqA2[1])*x2+q_eqA2[0];  
           res = x * (p / q);
	 }
       else 
       {if (x<q_erft[4])
          { p = (((((q_epB1[6]*x+q_epB1[5])*x+q_epB1[4])*x+q_epB1[3])*x+q_epB1[2])*x+q_epB1[1])*x+q_epB1[0];
            q = (((((q_eqB1[6]*x+q_eqB1[5])*x+q_eqB1[4])*x+q_eqB1[3])*x+q_eqB1[2])*x+q_eqB1[1])*x+q_eqB1[0];
            res = 1.0 - (q_expx2(x) * (p / q));
          } 
        else 
        {if (x<q_erft[5])
           { p = ((((q_epB2[5]*x+q_epB2[4])*x+q_epB2[3])*x+q_epB2[2])*x+q_epB2[1])*x+q_epB2[0];
             q = (((((q_eqB2[6]*x+q_eqB2[5])*x+q_eqB2[4])*x+q_eqB2[3])*x+q_eqB2[2])*x+q_eqB2[1])*x+q_eqB2[0];
             res = 1.0 - (q_expx2(x) * (p / q));
           }
         else res = 1.0;
	}                     
       }
      }
     }
    }
   }
  }
  return(res);
 }


/* ------------------------------------------------------------------- */
/* ----                   the function q_erfc                   ------ */
/* ------------------------------------------------------------------- */

 real q_erfc(real x){
  real p, q, res, x2;
  
  if(NANTEST(x)) res = q_abortnan(INV_ARG,&x,28);        /* Test: x=NaN */ 
  else 
  {if (x<-q_erft[1]) res = 1.0 + q_erf(-x);
   else 
   {if (x<q_erft[1]) res = 1.0; 
    else 
    {if (x<q_erft[2]) res = 1.0 - (x * q_epA2[0]);
     else 
     {if (x<q_erft[3]) 
        { x2 = x*x;
          p = (((q_epA2[4]*x2+q_epA2[3])*x2+q_epA2[2])*x2+q_epA2[1])*x2+q_epA2[0];  
          q = (((q_eqA2[4]*x2+q_eqA2[3])*x2+q_eqA2[2])*x2+q_eqA2[1])*x2+q_eqA2[0];  
          res = 1.0 - (x * (p / q));
	 }
      else 
      {if (x<q_erft[4])
         { p = (((((q_epB1[6]*x+q_epB1[5])*x+q_epB1[4])*x+q_epB1[3])*x+q_epB1[2])*x+q_epB1[1])*x+q_epB1[0];
           q = (((((q_eqB1[6]*x+q_eqB1[5])*x+q_eqB1[4])*x+q_eqB1[3])*x+q_eqB1[2])*x+q_eqB1[1])*x+q_eqB1[0];
           res = q_expx2(x) * (p / q);
          } 
       else 
       {if (x<q_erft[5])
          { p = ((((q_epB2[5]*x+q_epB2[4])*x+q_epB2[3])*x+q_epB2[2])*x+q_epB2[1])*x+q_epB2[0];
            q = (((((q_eqB2[6]*x+q_eqB2[5])*x+q_eqB2[4])*x+q_eqB2[3])*x+q_eqB2[2])*x+q_eqB2[1])*x+q_eqB2[0];
            res = q_expx2(x) * (p / q);
           } 
        else 
        {if (x<q_erft[6])
           { x2 = x*x;
             p = (((q_epB3[4]/x2+q_epB3[3])/x2+q_epB3[2])/x2+q_epB3[1])/x2+q_epB3[0];
             q = (((q_eqB3[4]/x2+q_eqB3[3])/x2+q_eqB3[2])/x2+q_eqB3[1])/x2+q_eqB3[0];
             res = ((q_expx2(x) * p) / (x*q));
            } 
	 else res = q_abortnan(INV_ARG,&x,28);        /* UNDERFLOW */ 
	}             
       }        
      }
     }
    }
   }
  }
  return(res);
 }
#endif
}
#endif





