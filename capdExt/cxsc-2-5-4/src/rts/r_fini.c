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

/* CVS $Id: r_fini.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_fini.c				              */
/*                                                              */
/*      Entries         : void r_fini()                         */
/*                                                              */
/*      Description     : Initializes the FPUnit of the target. */
/*                                                              */
/*      Note            : Version for SUN4_OS4_C                */
/*                                                              */
/*                        add prototype for ...ieee... at       */
/*                        SUN4_CPP_c                            */
/*                        cast to sigfpe_... at SIGFPE_IGNORE   */
/*                        changed pxsc_trap_handler to r_ftrp   */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C

#ifdef LINT_ARGS
/* 	Sun FUNCTIONS for C Programmers for IEEE floating point. */
/* These functions are located in math.h =b= */
extern int ieee_flags(char *clear, char *except, char *all, char **out);
extern int ieee_handler(char *set, char *trap, sigfpe_handler_type hdl);
#else
extern int ieee_flags();
extern int ieee_handler();
#endif

#ifdef LINT_ARGS
void r_fini(void)
#else
void r_fini()
#endif
{
        int status;
        char *out;
        sigfpe_handler_type hdl;

        hdl = (sigfpe_handler_type) r_ftrp;

        status = ieee_handler("clear", "all", 
				(sigfpe_handler_type) SIGFPE_IGNORE);

        /* set trap enable flags and assign trap handler */
        /* trap enable flags are later set in fsrd[ndu] (e_ieee.c) */
        status = ieee_handler("set", "invalid", hdl);
        status = ieee_handler("set", "division", hdl);
        status = ieee_handler("set", "overflow", hdl);
        status = ieee_handler("set", "underflow", hdl);
        status = ieee_handler("set", "inexact", hdl);
        
        /* clear status in processor */
        status = ieee_flags("clear", "exception", "all", &out);

        /* set rounding mode */
        status = ieee_flags("set", "direction", "nearest", &out) ; 

        /* set precision to double */
        status = ieee_flags("set", "precision", "double", &out);

        return;
}

#endif	/* SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C */
#if SUN4_OS5_GNU_C

#ifdef LINT_ARGS
void r_fini(void)
#else
void r_fini()
#endif
{
	r_lfsr();
	signal(SIGFPE, (void (*)(int)) r_ftrp);
}

#endif /* SUN4_OS5_GNU_C */
#if IBM_LINUX_C+IBM_EMX_C

#ifdef LINT_ARGS
void r_fini(void)
#else
void r_fini()
#endif
{
	__asm__("fninit");
        signal(SIGFPE, (void (*)(int)) r_ftrp);
}

#endif	/* IBM_LINUX_C */
#if IBM_RS6000_C

#ifdef LINT_ARGS
void r_fini(void)
#else
void r_fini()
#endif
{
	extern double e_srdn, e_srdu, e_srdd;
	struct sigaction fpehandler;

	/* clear all existing fp exception flags */
	fp_disable_all();
	fp_clr_flag(FP_ALL_XCP);

	/* setup exception handler */
	(void) sigemptyset(&fpehandler.sa_mask);
	fpehandler.sa_flags = FALSE;
	fpehandler.sa_handler = (void (*)(int)) r_ftrp;
	(void) sigaction(SIGFPE, &fpehandler, NULL);

#ifdef RS6000_SYNC_TRAP		/* define to support fp traps */
	/* set trap mode: we need precise exceptions */
	/* execution time increases by a factor of 4-5 */
	(void) fp_trap(FP_TRAP_SYNC);
#endif

	/* set rounding mode, default exceptions are set in p_init */
	*((((unsigned int *)(&e_srdn))+1)) = 0x00000000;
	*((((unsigned int *)(&e_srdu))+1)) = 0x00000002;
	*((((unsigned int *)(&e_srdd))+1)) = 0x00000003;

	/* enable exception bits in FP scr */
	r_lfsr();
}

#endif  /* IBM_RS6000_C */

#if HP_9000_C
#ifdef LINT_ARGS
void r_fini(void)
#else
void r_fini()
#endif
{
}
#endif

#else	/* IEEE_HARDWARE */

#ifdef LINT_ARGS
void r_fini(void)
#else
void r_fini()
#endif
{
}

#endif	/* IEEE_HARDWARE */





