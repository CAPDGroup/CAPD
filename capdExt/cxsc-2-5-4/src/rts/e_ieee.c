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

/* CVS $Id: e_ieee.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_ieee.c                              */
/*                                                              */
/*      Entries         : a_bool e_io_e( );                     */
/*                        a_bool e_dz_e( );                     */
/*                        a_bool e_of_e( );                     */
/*                        a_bool e_uf_e( );                     */
/*                        a_bool e_ie_e( );                     */
/*                        a_bool e_io_o( );                     */
/*                        a_bool e_dz_o( );                     */
/*                        a_bool e_of_o( );                     */
/*                        a_bool e_uf_o( );                     */
/*                        a_bool e_ie_o( );                     */
/*                        void e_sioo( );                       */
/*                        void e_sdzo( );                       */
/*                        void e_sofo( );                       */
/*                        void e_sufo( );                       */
/*                        void e_sieo( );                       */
/*                        void e_rioo( );                       */
/*                        void e_rdzo( );                       */
/*                        void e_rofo( );                       */
/*                        void e_rufo( );                       */
/*                        void e_rieo( );                       */
/*                        void e_sioe( );                       */
/*                        void e_sdze( );                       */
/*                        void e_sofe( );                       */
/*                        void e_sufe( );                       */
/*                        void e_siee( );                       */
/*                        void e_rioe( );                       */
/*                        void e_rdze( );                       */
/*                        void e_rofe( );                       */
/*                        void e_rufe( );                       */
/*                        void e_riee( );                       */
/*                        void e_sall( );                       */
/*                        void e_rall( );                       */
/*                        void e_rest(d);                       */
/*                        void e_save(p);                       */
/*                        a_intg *p;                            */
/*                        a_intg d;                             */
/*                                                              */
/*      Arguments       : p - place to store flags              */
/*                        d - flags to be restored              */
/*                                                              */
/*                   and e_temd =b=                             */
/*                   changed r_add funct. in LOAD_fsr to r_lfsr */
/*                   init hardware FSR, when fsr has changed    */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_bool e_efio;
extern a_bool e_efdz;
extern a_bool e_efof;
extern a_bool e_efuf;
extern a_bool e_efie;
extern a_bool e_ofio;
extern a_bool e_ofdz;
extern a_bool e_ofof;
extern a_bool e_ofuf;
extern a_bool e_ofie;
#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C+SUN4_OS5_GNU_C
extern long e_srdn, e_srdu, e_srdd;
#endif
#if HP_9000_C
extern long e_srdn, e_srdu, e_srdd, e_temd;
#endif
#if IBM_LINUX_C+IBM_EMX_C
extern short e_cwdn, e_cwdu, e_cwdd;
#endif
#if IBM_RS6000_C
extern double e_srdn, e_srdu, e_srdd;
#endif
#endif
#endif

#ifdef IEEE_HARDWARE

#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C+SUN4_OS5_GNU_C

/* trap enable mask */
#define TEM_NVM		0x08000000	/* invalid operation */
#define TEM_OFM		0x04000000	/* overflow */
#define TEM_UFM		0x02000000	/* underflow */
#define TEM_DZM		0x01000000	/* div. by zero */
#define TEM_NXM		0x00800000	/* inexact */

/* accrued exceptions */
#define EXC_NVA		0x00000200	/* invalid operation */
#define EXC_OFA		0x00000100	/* overflow */
#define EXC_UFA		0x00000080	/* underflow */
#define EXC_DZA		0x00000040	/* underflow */
#define EXC_NXA		0x00000020	/* inexact */

/* set/clear fsr bits */
#define SET_TEM(m)	{ e_srdu |= (m); e_srdn |= (m); e_srdd |= (m); }
#define CLR_TEM(m)	{ e_srdu &= ~(m); e_srdn &= ~(m); e_srdd &= ~(m); }
#define SET_FSR(m)	{ e_srdu |= (m); e_srdn |= (m); e_srdd |= (m); }
#define CLR_FSR(m)	{ e_srdu &= ~(m); e_srdn &= ~(m); e_srdd &= ~(m); }
/*#define LOAD_FSR()  { r_addn(0.0,0.0); }*/
#define LOAD_FSR()  { r_lfsr(); }
#define GF

/* load e_of* bits from fsr */
#define LOAD_FLAGS { unsigned long fsr;\
	    fsr = e_srdu|e_srdn|e_srdd;\
	    if (fsr&EXC_NVA) e_ofio = TRUE;\
	    if (fsr&EXC_OFA) e_ofof = TRUE;\
	    if (fsr&EXC_UFA) e_ofuf = TRUE;\
	    if (fsr&EXC_DZA) e_ofdz = TRUE;\
	    if (fsr&EXC_NXA) e_ofie = TRUE;\
		}

#endif	/* SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C+SUN4_OS5_GNU_C */

#if HP_9000_C

/* trap enable mask */
#define TEM_NVM		0x80000000	/* invalid operation */
#define TEM_DZM		0x40000000	/* div. by zero */
#define TEM_OFM		0x20000000	/* overflow */
#define TEM_UFM		0x10000000	/* underflow */
#define TEM_NXM		0x08000000	/* inexact */

/* accrued exceptions */
#define EXC_NVA		0x80000000	/* invalid operation */
#define EXC_DZA		0x40000000	/* div. by zero */
#define EXC_OFA		0x20000000	/* overflow */
#define EXC_UFA		0x10000000	/* underflow */
#define EXC_NXA		0x08000000	/* inexact */

/* set/clear fsr bits */
#define SET_TEM(m)	{ e_temd |= (m); }
#define CLR_TEM(m)	{ e_temd &= ~(m); }
#define SET_FSR(m)	{ e_srdu |= (m); e_srdn |= (m); e_srdd |= (m); }
#define CLR_FSR(m)	{ e_srdu &= ~(m); e_srdn &= ~(m); e_srdd &= ~(m); }
/*#define LOAD_FSR()  { r_addn(0.0,0.0); }*/
#define LOAD_FSR()  { r_lfsr(); }
#define GF

/* load e_of* bits from fsr */
#define LOAD_FLAGS { unsigned long fsr;\
	    fsr = e_srdu|e_srdn|e_srdd;\
	    if (fsr&EXC_NVA) e_ofio = TRUE;\
	    if (fsr&EXC_OFA) e_ofof = TRUE;\
	    if (fsr&EXC_UFA) e_ofuf = TRUE;\
	    if (fsr&EXC_DZA) e_ofdz = TRUE;\
	    if (fsr&EXC_NXA) e_ofie = TRUE;\
		}

#endif	/* HP_9000_C */

#if IBM_LINUX_C+IBM_EMX_C

/* trap _masks_, i.e. disable */
#define TEM_NVM		0x0001		/* invalid operation */
#define TEM_OFM		0x0008		/* overflow */
#define TEM_UFM		0x0010		/* underflow */
#define TEM_DZM		0x0004		/* div. by zero */
#define TEM_NXM		0x0020		/* inexact */

/* exception flags */
#define EXC_NVA		0x0001		/* invalid operation */
#define EXC_OFA		0x0008		/* overflow */
#define EXC_UFA		0x0010		/* underflow */
#define EXC_DZA		0x0004		/* div. by zero */
#define EXC_NXA		0x0020		/* inexact */

/* set/clear fsr bits */
#define SET_TEM(m)	{ e_cwdu &= ~(m); e_cwdn &= ~(m); e_cwdd &= ~(m); }
#define CLR_TEM(m)	{ e_cwdu |=  (m); e_cwdn |=  (m); e_cwdd |=  (m); }
#define SET_FSR(m)	{ ; }
#define CLR_FSR(m)	{ ; }
#define LOAD_FSR()	{ r_lfsr(); }

/* LOAD_FSR zerstoert Exception-Bits - daher LOAD_FLAGS noetig */
#define GF		LOAD_FLAGS	/* get flags */

/* load e_of* bits from fsr */
#define LOAD_FLAGS { static unsigned short fsr;\
	    __asm__(" fnstsw %0\n fnclex" : "=m" (fsr));\
	    if (fsr&EXC_NVA) e_ofio = TRUE;\
	    if (fsr&EXC_OFA) e_ofof = TRUE;\
	    if (fsr&EXC_UFA) e_ofuf = TRUE;\
	    if (fsr&EXC_DZA) e_ofdz = TRUE;\
	    if (fsr&EXC_NXA) e_ofie = TRUE;\
		}

#endif	/* IBM_LINUX_C */

#if IBM_RS6000_C

/* trap enable mask */
#define TEM_NVM		0x00000080	/* invalid operation */
#define TEM_OFM		0x00000040	/* overflow */
#define TEM_UFM		0x00000020	/* underflow */
#define TEM_DZM		0x00000010	/* div. by zero */
#define TEM_NXM		0x00000008	/* inexact */

/* exception flags */
#define EXC_NVA		0x21f80300	/* invalid operation */
#define EXC_OFA		0x10000000	/* overflow */
#define EXC_UFA		0x08000000	/* underflow */
#define EXC_DZA		0x04000000	/* div. by zero */
#define EXC_NXA		0x02000000	/* inexact */

/* set/clear fp scr bits */
#define S(a)		(*(((unsigned int *)(&(a)))+1))
#define D(a)		((a)&0x3f000000)	/* NV setzen nur VXSNAN */
#define SET_TEM(m)	{ S(e_srdu)|= (m); S(e_srdn)|= (m); S(e_srdd)|= (m); }
#define CLR_TEM(m)	{ S(e_srdu)&=~(m); S(e_srdn)&=~(m); S(e_srdd)&=~(m); }
#define SET_FSR(m)	{ S(e_srdu)|=D(m); S(e_srdn)|=D(m); S(e_srdd)|=D(m); }
#define CLR_FSR(m)	{ S(e_srdu)&=~(m); S(e_srdn)&=~(m); S(e_srdd)&=~(m); }
#define LOAD_FSR()	{ r_lfsr(); }
#define GF

/* load e_of* bits from fsr */
#define LOAD_FLAGS { unsigned int fsr;\
	    fsr = S(e_srdu)|S(e_srdn)|S(e_srdd);\
	    if (fsr&EXC_NVA) e_ofio = TRUE;\
	    if (fsr&EXC_OFA) e_ofof = TRUE;\
	    if (fsr&EXC_UFA) e_ofuf = TRUE;\
	    if (fsr&EXC_DZA) e_ofdz = TRUE;\
	    if (fsr&EXC_NXA) e_ofie = TRUE;\
		}

#endif  /* IBM_RS6000_C */

#else	/* IEEE_HARDWARE */

#define SET_TEM(m)       {;}
#define CLR_TEM(m)        
#define SET_FSR(m)       {;}
#define CLR_FSR(m)        
#define LOAD_FLAGS
#define GF
#define TEM_NVM		0
#define TEM_OFM		0
#define TEM_UFM		0
#define TEM_DZM		0
#define TEM_NXM		0
#define EXC_NVA		0
#define EXC_OFA		0
#define EXC_UFA		0
#define EXC_DZA		0
#define EXC_NXA		0
#define LOAD_FSR()       {;}

#endif	/* IEEE_HARDWARE */

/*--------------------------------------------------------------*/
/* Functions                                                    */
/*--------------------------------------------------------------*/

#ifdef LINT_ARGS
local a_bool e_io_e(void) { return(e_efio); }
local a_bool e_dz_e(void) { return(e_efdz); }
local a_bool e_of_e(void) { return(e_efof); }
local a_bool e_uf_e(void) { return(e_efuf); }
local a_bool e_ie_e(void) { return(e_efie); }

local a_bool e_io_o(void) { LOAD_FLAGS; return(e_ofio); }
local a_bool e_dz_o(void) { LOAD_FLAGS; return(e_ofdz); }
local a_bool e_of_o(void) { LOAD_FLAGS; return(e_ofof); }
local a_bool e_uf_o(void) { LOAD_FLAGS; return(e_ofuf); }
local a_bool e_ie_o(void) { LOAD_FLAGS; return(e_ofie); }

local void e_sioo(void) { GF; e_ofio = TRUE; SET_FSR(EXC_NVA); LOAD_FSR() ;}
local void e_sdzo(void) { GF; e_ofdz = TRUE; SET_FSR(EXC_DZA); LOAD_FSR() ;}
local void e_sofo(void) { GF; e_ofof = TRUE; SET_FSR(EXC_OFA); LOAD_FSR() ;}
local void e_sufo(void) { GF; e_ofuf = TRUE; SET_FSR(EXC_UFA); LOAD_FSR() ;}
local void e_sieo(void) { GF; e_ofie = TRUE; SET_FSR(EXC_NXA); LOAD_FSR() ;}

local void e_rioo(void) { GF; e_ofio = FALSE; CLR_FSR(EXC_NVA); LOAD_FSR(); }
local void e_rdzo(void) { GF; e_ofdz = FALSE; CLR_FSR(EXC_DZA); LOAD_FSR(); }
local void e_rofo(void) { GF; e_ofof = FALSE; CLR_FSR(EXC_OFA); LOAD_FSR(); }
local void e_rufo(void) { GF; e_ofuf = FALSE; CLR_FSR(EXC_UFA); LOAD_FSR(); }
local void e_rieo(void) { GF; e_ofie = FALSE; CLR_FSR(EXC_NXA); LOAD_FSR(); }

local void e_sioe(void) { GF; e_efio = TRUE; SET_TEM(TEM_NVM); LOAD_FSR(); }
local void e_sdze(void) { GF; e_efdz = TRUE; SET_TEM(TEM_DZM); LOAD_FSR(); }
local void e_sofe(void) { GF; e_efof = TRUE; SET_TEM(TEM_OFM); LOAD_FSR(); }
local void e_sufe(void) { GF; e_efuf = TRUE; SET_TEM(TEM_UFM); LOAD_FSR(); }
local void e_siee(void) { GF; e_efie = TRUE; SET_TEM(TEM_NXM); LOAD_FSR(); }

local void e_rioe(void) { GF; e_efio = FALSE; CLR_TEM(TEM_NVM); LOAD_FSR(); }
local void e_rdze(void) { GF; e_efdz = FALSE; CLR_TEM(TEM_DZM); LOAD_FSR(); }
local void e_rofe(void) { GF; e_efof = FALSE; CLR_TEM(TEM_OFM); LOAD_FSR(); }
local void e_rufe(void) { GF; e_efuf = FALSE; CLR_TEM(TEM_UFM); LOAD_FSR(); }
local void e_riee(void) { GF; e_efie = FALSE; CLR_TEM(TEM_NXM); LOAD_FSR(); }
#else
local a_bool e_io_e() { return(e_efio); }
local a_bool e_dz_e() { return(e_efdz); }
local a_bool e_of_e() { return(e_efof); }
local a_bool e_uf_e() { return(e_efuf); }
local a_bool e_ie_e() { return(e_efie); }

local a_bool e_io_o() { LOAD_FLAGS; return(e_ofio); }
local a_bool e_dz_o() { LOAD_FLAGS; return(e_ofdz); }
local a_bool e_of_o() { LOAD_FLAGS; return(e_ofof); }
local a_bool e_uf_o() { LOAD_FLAGS; return(e_ofuf); }
local a_bool e_ie_o() { LOAD_FLAGS; return(e_ofie); }

local void e_sioo() { GF; e_ofio = TRUE; SET_FSR(EXC_NVA); LOAD_FSR(); }
local void e_sdzo() { GF; e_ofdz = TRUE; SET_FSR(EXC_DZA); LOAD_FSR(); }
local void e_sofo() { GF; e_ofof = TRUE; SET_FSR(EXC_OFA); LOAD_FSR(); }
local void e_sufo() { GF; e_ofuf = TRUE; SET_FSR(EXC_UFA); LOAD_FSR(); }
local void e_sieo() { GF; e_ofie = TRUE; SET_FSR(EXC_NXA); LOAD_FSR(); }

local void e_rioo() { GF; e_ofio = FALSE; CLR_FSR(EXC_NVA); LOAD_FSR(); }
local void e_rdzo() { GF; e_ofdz = FALSE; CLR_FSR(EXC_DZA); LOAD_FSR(); }
local void e_rofo() { GF; e_ofof = FALSE; CLR_FSR(EXC_OFA); LOAD_FSR(); }
local void e_rufo() { GF; e_ofuf = FALSE; CLR_FSR(EXC_UFA); LOAD_FSR(); }
local void e_rieo() { GF; e_ofie = FALSE; CLR_FSR(EXC_NXA); LOAD_FSR(); }

local void e_sioe() { GF; e_efio = TRUE; SET_TEM(TEM_NVM); LOAD_FSR(); }
local void e_sdze() { GF; e_efdz = TRUE; SET_TEM(TEM_DZM); LOAD_FSR(); }
local void e_sofe() { GF; e_efof = TRUE; SET_TEM(TEM_OFM); LOAD_FSR(); }
local void e_sufe() { GF; e_efuf = TRUE; SET_TEM(TEM_UFM); LOAD_FSR(); }
local void e_siee() { GF; e_efie = TRUE; SET_TEM(TEM_NXM); LOAD_FSR(); }

local void e_rioe() { GF; e_efio = FALSE; CLR_TEM(TEM_NVM); LOAD_FSR(); }
local void e_rdze() { GF; e_efdz = FALSE; CLR_TEM(TEM_DZM); LOAD_FSR(); }
local void e_rofe() { GF; e_efof = FALSE; CLR_TEM(TEM_OFM); LOAD_FSR(); }
local void e_rufe() { GF; e_efuf = FALSE; CLR_TEM(TEM_UFM); LOAD_FSR(); }
local void e_riee() { GF; e_efie = FALSE; CLR_TEM(TEM_NXM); LOAD_FSR(); }
#endif

#ifdef LINT_ARGS
local void e_sall(void)
#else
local void e_sall()
#endif
{
	GF;
        e_ofio = e_ofdz = e_ofof = e_ofuf = e_ofie = TRUE;
	SET_FSR(EXC_NVA|EXC_DZA|EXC_OFA|EXC_UFA|EXC_NXA);
	LOAD_FSR() ;
}

#ifdef LINT_ARGS
local void e_rall(void)
#else
local void e_rall()
#endif
{
	GF;
        e_ofio = e_ofdz = e_ofof = e_ofuf = e_ofie = FALSE;
	CLR_FSR(EXC_NVA|EXC_DZA|EXC_OFA|EXC_UFA|EXC_NXA);
	LOAD_FSR() ;
}

#ifdef LINT_ARGS
local void e_save(a_intg *p)
#else
local void e_save(p)
a_intg *p;
#endif
        {
	LOAD_FLAGS;
        *p = 0;
        if (e_ofio) *p += 1;
        if (e_ofdz) *p += 2;
        if (e_ofof) *p += 4;
        if (e_ofuf) *p += 8;
        if (e_ofie) *p += 16;
        if (e_efio) *p += 32;
        if (e_efdz) *p += 64;
        if (e_efof) *p += 128;
        if (e_efuf) *p += 256;
        if (e_efie) *p += 512;
        }

#ifdef LINT_ARGS
local void e_rest(a_intg d)
#else
local void e_rest(d)

a_intg d;
#endif
        {
        e_ofio = (d & 1) ? TRUE : FALSE;
        e_ofdz = (d & 2) ? TRUE : FALSE;
        e_ofof = (d & 4) ? TRUE : FALSE;
        e_ofuf = (d & 8) ? TRUE : FALSE;
        e_ofie = (d & 16) ? TRUE : FALSE;
        e_efio = (d & 32) ? TRUE : FALSE;
        e_efdz = (d & 64) ? TRUE : FALSE;
        e_efof = (d & 128) ? TRUE : FALSE;
        e_efuf = (d & 256) ? TRUE : FALSE;
        e_efie = (d & 512) ? TRUE : FALSE;
        if (d &   1) SET_FSR(EXC_NVA) else CLR_FSR(EXC_NVA);
        if (d &   2) SET_FSR(EXC_DZA) else CLR_FSR(EXC_DZA);
        if (d &   4) SET_FSR(EXC_OFA) else CLR_FSR(EXC_OFA);
        if (d &   8) SET_FSR(EXC_UFA) else CLR_FSR(EXC_UFA);
        if (d &  16) SET_FSR(EXC_NXA) else CLR_FSR(EXC_NXA);
        if (d &  32) SET_TEM(TEM_NVM) else CLR_TEM(TEM_NVM);
        if (d &  64) SET_TEM(TEM_DZM) else CLR_TEM(TEM_DZM);
        if (d & 128) SET_TEM(TEM_OFM) else CLR_TEM(TEM_OFM);
        if (d & 256) SET_TEM(TEM_UFM) else CLR_TEM(TEM_UFM);
        if (d & 512) SET_TEM(TEM_NXM) else CLR_TEM(TEM_NXM);
        LOAD_FSR() ;
        }





