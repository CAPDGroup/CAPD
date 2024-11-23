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

/* CVS $Id: fi_lib_consts.hpp,v 1.17 2014/01/30 17:23:54 cxsc Exp $ */

#ifndef FI_LIB_CONSTS_HPP
#define FI_LIB_CONSTS_HPP

#include "fi_lib.hpp" 

namespace fi_lib{

using cxsc::real;
 
/*********************************************************************/
/* replacement for function ldexp from math.h (libm.a)               */
/*********************************************************************/

#define POWER2(x,k) {if (x!=0) {a_diee *test=(a_diee *) &x; test->ieee.expo+=k;}}               /* only for normalised numbers ! */

/*********************************************************************/
/* replacement for function frexp from math.h (libm.a)               */
/*********************************************************************/

#define FREXPO(x,k) {if (x!=0) {a_diee *test=(a_diee *) &x; k=test->ieee.expo;} else k=0;}                          /* only for normalised numbers ! */

/*********************************************************************/
/* first 24 bits of a real value ( cut / round )                   */
/*********************************************************************/

#define CUT24(x) (float)(x)

/*********************************************************************/
/* data type union to access sign, mantisse, exponent                */
/*********************************************************************/

#ifndef CXSC_INCLUDE

typedef union __a_diee
        {
                double f;
#ifdef IS_64_BIT
                long l;
#endif
                struct
                {
#if INTEL
                        unsigned int mant1 :32;
                        unsigned int mant0 :20;
                        unsigned int expo  :11;
                        unsigned int sign  : 1;
#else
                        unsigned int sign  : 1;
                        unsigned int expo  :11;
                        unsigned int mant0 :20;
                        unsigned int mant1 :32;
#endif
                } ieee;
        } a_diee;

#endif


/*********************************************************************/
/* assignment real -> int  (32 bits)                               */
/*********************************************************************/

#define CUTINT(x) (int)(_double(x))


/*********************************************************************/
/* error handling                                                    */
/*********************************************************************/

const int INV_ARG=1;  /* for error handling: invalid argument */
const int OVER_FLOW=2; /* for error handling: overflow         */
const int DIV_ZERO=3;  /* for error handling: division by zero */

/*********************************************************************/
/* constants for the elementary functions                            */
/*********************************************************************/

extern real q_ln2h;
extern real q_l10i;
extern real q_l2i;
extern real q_l10;
extern real q_l2;
extern real q_p2h;
extern real q_p2mh;
extern real q_mine;
extern real q_minr;
extern real q_pi;
extern real q_piha;
extern real q_pih[7];
extern real q_pi2i;
extern real q_pi2p[3];
extern real q_sqra;
extern real q_ctht;

extern real q_ext1;
extern real q_ext2;
extern real q_ex2a;
extern real q_ex2b;
extern real q_ex2c;
extern real q_ext3;
extern real q_ext4;
extern real q_ext5;
extern real q_extm;
extern real q_extn;
extern real q_exil;
extern real q_exl1;
extern real q_exl2;
extern real q_e10i;
extern real q_e1l1;
extern real q_e1l2;
extern real q_exa[5];
extern real q_exb[9];
extern real q_exc[7];
extern real q_exd[7];
extern real q_exld[32];
extern real q_extl[32];

extern real q_exem;
extern real q_exep;
extern real q_exmm;
extern real q_exmp;

extern real q_snhm;
extern real q_snhp;
extern real q_cshm;
extern real q_cshp;
extern real q_cthm;
extern real q_cthp;
extern real q_tnhm;
extern real q_tnhp;

extern real q_lgt1;
extern real q_lgt2;
extern real q_lgt3;
extern real q_lgt4;
extern real q_lgt5;
extern real q_lgt6;
extern real q_lgb[2];
extern real q_lgc[4];
extern real q_lgld[129];
extern real q_lgtl[129];

extern real q_logm;
extern real q_logp;
extern real q_lgpm;
extern real q_lgpp;
extern real q_sqtm;
extern real q_sqtp;

extern real q_atna[7];
extern real q_atnb[8];
extern real q_atnc[7];
extern real q_atnd[6];
extern real q_atnt;

extern real q_asnm;
extern real q_asnp;
extern real q_acsm;
extern real q_acsp;
extern real q_actm;
extern real q_actp;
extern real q_atnm;
extern real q_atnp;

extern real q_csnm;
extern real q_csnp;
extern real q_ccsm;
extern real q_ccsp;
extern real q_cctm;
extern real q_cctp;
extern real q_ctnm;
extern real q_ctnp;

extern real q_sinc[6];
extern real q_sins[6];
extern real q_sint[5];

extern real q_sinm;
extern real q_sinp;
extern real q_cosm;
extern real q_cosp;
extern real q_cotm;
extern real q_cotp;
extern real q_tanm;
extern real q_tanp;

extern real q_lg2m;
extern real q_lg2p;
extern real q_l10m;
extern real q_l10p;
extern real q_e2em;
extern real q_e2ep;
extern real q_e10m;
extern real q_e10p;

extern real q_at3i;

// nu_ko_1 = 2/sqrt(pi):
extern real nu_ko_1;

#ifdef erf_kettenbruch
extern real q_erf_m;
extern real q_erf_p;
extern real q_erfc_m;
extern real q_erfc_p;
extern real q_erfa1;
extern real q_erfx0;
extern real q_erfc_1;
extern real q_erfa3_a[5];
extern real q_erfa3_b[5];
extern real q_erfb_a[8];
extern real q_erfb_b[8];
extern real q_erfB_x0;
extern real q_erfc_a[7];
extern real q_erfc_b[7];
extern real q_erfC_x0;
extern real q_erfd_a[7];
extern real q_erfd_b[7];
extern real q_erfD_x0;
extern real q_erfe_a[7];
extern real q_erfe_b[7];
extern real q_erfE_x0;
extern real q_erff_a[7];
extern real q_erff_b[7];
extern real q_erfF_x0;
extern real q_erfg_a[7];
extern real q_erfg_b[7];
extern real q_erfG_x0;
extern real q_erfh_a[6];
extern real q_erfh_b[6];
extern real q_erfH_x0;
extern real q_erfcb_a[8];
extern real q_erfcb_b[8];
extern real q_erfcB_x0;
extern real q_erfcc_a[7];
extern real q_erfcc_b[7];
extern real q_erfcC_x0;
extern real q_erfcd_a[8];
extern real q_erfcd_b[8];
extern real q_erfcD_x0;
extern real q_erfce_a[6];
extern real q_erfce_b[6];
extern real q_erfcE_x0;
extern real q_erfcf_a[6];
extern real q_erfcf_b[6];
extern real q_erfcF_x0;
extern real q_erfcg_a[6];
extern real q_erfcg_b[6];
extern real q_erfcG_x0;
extern real q_erfch_a[5];
extern real q_erfch_b[5];
extern real q_erfcH_x0;
#else
extern real q_erfm;
extern real q_erfp;
extern real q_efcm;
extern real q_efcp;
extern real q_epA2[5];
extern real q_eqA2[5];
extern real q_epB1[7];
extern real q_eqB1[7];
extern real q_eqB2[7];
extern real q_epB2[6];
extern real q_epB3[5];
extern real q_eqB3[5];
extern real q_erft[7];
#endif
extern real q_expz[28];
} // Namespace

#endif





