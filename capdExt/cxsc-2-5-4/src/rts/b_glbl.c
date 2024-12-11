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

/* CVS $Id: b_glbl.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_glbl.c                              */
/*                                                              */
/*      Description     : Definition of global variables.       */
/*                                                              */
/*      Note            : All variables must be initialized.    */
/*                                                              */
/*                   (n,u,d),fsrq to e_srq_ =b=                 */
/*                   hardware control word variable b_ctrl      */
/*                   f_name reserves working space for filenames*/
/*                   in r_pred(), r_succ().                     */
/*                   DIRNAME_LENGTH used.                       */
/*                   final record corrected                     */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

local unsigned long int b_ctrl = 0;

local struct a_pptr *a_ptop = NULL;

local char *tp_000 = NULL;
local char *tp_001 = NULL;

local bentry *e_bstk = (bentry *)NULL;
local bentry *e_btop = (bentry *)NULL;

local int e_line = 0;

local a_bool e_efio = FALSE;
local a_bool e_efdz = FALSE;
local a_bool e_efof = FALSE;
local a_bool e_efuf = FALSE;
local a_bool e_efie = FALSE;
local a_bool e_ofio = FALSE;
local a_bool e_ofdz = FALSE;
local a_bool e_ofof = FALSE;
local a_bool e_ofuf = FALSE;
local a_bool e_ofie = FALSE;

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C+SUN4_OS5_GNU_C
local unsigned int e_srdn = 0x2d000000;		/* SPARC FSR double, near */
local unsigned int e_srdu = 0xad000000;		/* SPARC FSR double, up   */
local unsigned int e_srdd = 0xed000000;		/* SPARC FSR double, down */
/* default traps: division by zero, invalid operation, overflow */

local unsigned int e_srqn = 0x30000000;		/* SPARC FSR quad(80), near */
local unsigned int e_srqu = 0xb0000000;		/* SPARC FSR quad(80), up   */
local unsigned int e_srqd = 0xf0000000;		/* SPARC FSR quad(80), down */
local unsigned int e_srqc = 0x70000000;		/* SPARC FSR quad(80), chop */
local unsigned int e_srq_;                   /* SPARC FSR quad, aktuell  */
#endif
#if HP_9000_C
local unsigned int e_srdn = 0x00000000;		/* HP9000 FSR, near */
local unsigned int e_srdu = 0x00000400;		/* HP9000 FSR, up   */
local unsigned int e_srdd = 0x00000600;		/* HP9000 FSR, down */
/* default traps: division by zero, invalid operation, overflow */
local unsigned int e_temd  = 0xe0000000;
#endif
#if IBM_LINUX_C+IBM_EMX_C
local unsigned short e_cwdn = 0x0272;		/* i486 cw double, near */
local unsigned short e_cwdu = 0x0a72;		/* i486 cw double, up   */
local unsigned short e_cwdd = 0x0672;		/* i486 cw double, down */
/* default traps: division by zero, invalid operation, overflow */

local unsigned short e_cwen = 0x037f;		/* i486 cw extended, near */
local unsigned short e_cweu = 0x0b7f;		/* i486 cw extended, up   */
local unsigned short e_cwed = 0x077f;		/* i486 cw extended, down */
local unsigned short e_cwec = 0x0f7f;		/* i486 cw extended, chop */
/* default traps: none */
local unsigned short e_cwe_;			/* i486 cw extended, aktuell */
#endif
#if IBM_RS6000_C
local double e_srdn = 1.0;			/* RS/6000 fpscr near */
local double e_srdu = 1.0;			/* RS/6000 fpscr up   */
local double e_srdd = 1.0;			/* RS/6000 fpscr down */
/* default traps: division by zero, invalid operation, overflow */
/* Initialisierung in r_fini */
#endif
#endif

local a_btyp b_maxl = 3;

local char *e_head = NULL;

local char f_name[f_fnsz];

local int f_argc = 0;
local int f_orgc = 0;
local char **f_argv = NULL;
local char **f_orgv = NULL;
local int f_apos = 0;
local a_bool f_ppcc = FALSE;
#ifdef NORMALIZE_ENABLE
local a_bool f_ppdn = FALSE;
#endif
#ifdef COPYIB_ENABLE
local a_bool f_ppib = FALSE;
#endif
local a_bool f_ppmt = FALSE;

/****** f_pppd: indicates that program parameter list is displayed ******/
local a_bool f_pppd = FALSE;

/****** f_pppl: indicates that program parameter list is processed ******/
local a_bool f_pppl = TRUE;

/****** f_ppsz: indicates that signed IEEE zero is displayed as -0 ******/
#ifdef SIGNED_ENABLE
local a_bool f_ppsz = FALSE;
#endif
#ifdef STATTB_ENABLE
local a_bool f_pptb = FALSE;
#endif

/****** f_ftop: top of linked list of file descriptors ******************/
local a_VOID f_ftop = NULL;

/****** f_pptf: indicates that no temporary files are generated *********/
local a_bool f_pptf = FALSE;

#ifdef TRACEC_ENABLE
local a_bool f_pptr = FALSE;
local int e_tlvl = 0;
#endif

local f_text f_inpu;
local f_text f_outp;
local f_text f_errr = { NULL,   FALSE, FALSE, TRUE, FALSE, TRUE,
                                FALSE, TRUE, TRUE, FALSE, FALSE, 1 };
local f_text f_pmti = { NULL,   FALSE, TRUE, TRUE, TRUE, FALSE,
                                TRUE, FALSE, TRUE, FALSE, FALSE, 1 };
local f_text f_pmto = { NULL,   FALSE, FALSE, TRUE, FALSE, TRUE,
                                FALSE, TRUE, TRUE, FALSE, FALSE, 1 };

#if SHORTABTYP
static a_btyp o_1o2_[] = B_TYPINI( 0x3fe00000, 0x00000000 );
static a_btyp o_eps_[] = B_TYPINI( 0x00000000, 0x00000001 );
static a_btyp o_flnb[] = B_TYPINI( 0x40362e42, 0xfefa39ef );
static a_btyp o_fln2[] = B_TYPINI( 0x3fe62e42, 0xfefa39ef );
static a_btyp o_max_[] = B_TYPINI( 0x7fefffff, 0xffffffff );
static a_btyp o_meps[] = B_TYPINI( 0x80000000, 0x00000001 );
static a_btyp o_minf[] = B_TYPINI( 0xfff00000, 0x00000000 );
static a_btyp o_mmax[] = B_TYPINI( 0xffefffff, 0xffffffff );
static a_btyp o_mone[] = B_TYPINI( 0xbff00000, 0x00000000 );
static a_btyp o_one_[] = B_TYPINI( 0x3ff00000, 0x00000000 );
static a_btyp o_pinf[] = B_TYPINI( 0x7ff00000, 0x00000000 );
static a_btyp o_pio2[] = B_TYPINI( 0x3ff921fb, 0x54442d18 );  /* chopped */
static a_btyp o_pio4[] = B_TYPINI( 0x3fe921fb, 0x54442d18 );  /* chopped */
static a_btyp o_sero[] = B_TYPINI( 0x80000000, 0x00000000 );
static a_btyp o_ten_[] = B_TYPINI( 0x40240000, 0x00000000 );
static a_btyp o_two_[] = B_TYPINI( 0x40000000, 0x00000000 );
static a_btyp o_zero[] = B_TYPINI( 0x00000000, 0x00000000 );
#else
static a_btyp o_1o2_[] = B_TYPINI( 0x3fe00000L, 0x00000000L );
static a_btyp o_eps_[] = B_TYPINI( 0x00000000L, 0x00000001L );
static a_btyp o_flnb[] = B_TYPINI( 0x40362e42L, 0xfefa39efL );
static a_btyp o_fln2[] = B_TYPINI( 0x3fe62e42L, 0xfefa39efL );
static a_btyp o_max_[] = B_TYPINI( 0x7fefffffL, 0xffffffffL );
static a_btyp o_meps[] = B_TYPINI( 0x80000000L, 0x00000001L );
static a_btyp o_minf[] = B_TYPINI( 0xfff00000L, 0x00000000L );
static a_btyp o_mmax[] = B_TYPINI( 0xffefffffL, 0xffffffffL );
static a_btyp o_mone[] = B_TYPINI( 0xbff00000L, 0x00000000L );
static a_btyp o_one_[] = B_TYPINI( 0x3ff00000L, 0x00000000L );
static a_btyp o_pinf[] = B_TYPINI( 0x7ff00000L, 0x00000000L );
static a_btyp o_pio2[] = B_TYPINI( 0x3ff921fbL, 0x54442d18L );  /* chopped */
static a_btyp o_pio4[] = B_TYPINI( 0x3fe921fbL, 0x54442d18L );  /* chopped */
static a_btyp o_sero[] = B_TYPINI( 0x80000000L, 0x00000000L );
static a_btyp o_ten_[] = B_TYPINI( 0x40240000L, 0x00000000L );
static a_btyp o_two_[] = B_TYPINI( 0x40000000L, 0x00000000L );
static a_btyp o_zero[] = B_TYPINI( 0x00000000L, 0x00000000L );
#endif

/* possible allignment problems if a_real and a_btyp are not    */
/* alligned the same way                                        */
local a_real *r_1o2_ = (a_real *)&o_1o2_[0];
local a_real *r_eps_ = (a_real *)&o_eps_[0];
local a_real *r_fln2 = (a_real *)&o_fln2[0];
local a_real *r_flnb = (a_real *)&o_flnb[0];
local a_real *r_max_ = (a_real *)&o_max_[0];
local a_real *r_meps = (a_real *)&o_meps[0];
local a_real *r_minf = (a_real *)&o_minf[0];
local a_real *r_mmax = (a_real *)&o_mmax[0];
local a_real *r_mone = (a_real *)&o_mone[0];
local a_real *r_one_ = (a_real *)&o_one_[0];
local a_real *r_pinf = (a_real *)&o_pinf[0];
local a_real *r_pio2 = (a_real *)&o_pio2[0];
local a_real *r_pio4 = (a_real *)&o_pio4[0];
local a_real *r_sero = (a_real *)&o_sero[0];
local a_real *r_ten_ = (a_real *)&o_ten_[0];
local a_real *r_two_ = (a_real *)&o_two_[0];
local a_real *r_zero = (a_real *)&o_zero[0];

local aentry *e_astk = (aentry *)NULL;

local a_VOID e_rptr = NULL;
local int e_rtyp = 0;

local char *f_pplt[FPPLT] = { NULL };

local int b_rflg = NEAREST;

#if VAX_VMS_C+CRAY_UNIX_C
#ifdef LINT_ARGS
local void b_glbl(void)
#else
local void b_glbl()
#endif
        {
        }
#endif





