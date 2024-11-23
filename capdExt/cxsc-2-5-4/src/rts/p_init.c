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

/* CVS $Id: p_init.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : p_init.c                              */
/*                                                              */
/*      Entry           : void f_init(argc,argv)                */
/*                        int argc;                             */
/*                        char **argv;                          */
/*                                                              */
/*      Arguments       : argc   - number of arguments          */
/*                        argv   - list of argument strings     */
/*                                                              */
/*      Description     : Initialization of command line scanner*/
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#include "/u/p88c/runtime/o_revs.h"
#else
#include "o_defs.h"
#include "o_revs.h"
#endif
#define local
extern int f_argc;
extern int f_orgc;
extern char **f_orgv;
extern char **f_argv;
extern int f_apos;
extern a_bool f_pppd;
extern a_bool f_pppl;
extern f_text f_pmti;
extern f_text f_pmto;
extern f_text f_errr;
extern FILE *o_pmti;
extern FILE *o_pmto;
extern FILE *o_errr;
extern char *o_text[];
extern char *e_head;
extern int b_rflg;
#endif

#ifdef LINT_ARGS
local void p_init(int argc,char **argv)
#else
local void p_init(argc,argv)

int argc;
char **argv;
#endif

        {
        int i,k;

#if IBM_OS2_ICC_C
        /* special settings for os2 in order t use double for a_real 
	    * and direct output behaviour as used to
	    * position of code is IMPORTANT here !!
	    */
	   setbuf(stdout,(char*) 0);                /* unbuffered output */
	   _control87(EM_UNDERFLOW, EM_UNDERFLOW);  /* direct double trap*/
#endif

        /* user initialization of global variables */
        o_user();

        /* initialization of variables with values from o_user() */
        e_head = o_text[0];
        f_pmti.fp = o_pmti;
        f_pmto.fp = o_pmto;
        f_errr.fp = o_errr;

        f_pmti.win.ch[0] = ' ';

        /* force linkage of module containing global variables only */
#if VAX_VMS_C
        b_cmcp();
        b_glbl();
#ifndef DEC_ARITH
        t_glbl();
        t_cnst();
        t_ecst();
#endif
#endif


        /* initialisation of FP HW trap handler (Sun SPARC) */
        r_fini();

        /* initialization of PASCAL environment         */
        e_sofe();      /* enable overflow trap         */
        e_sioe();      /* enable invalid operation trap*/
        e_sdze();      /* enable divide-by-zero trap   */
        e_riee();      /* disable inexact reult trap   */
        e_rufe();      /* disable underflow trap       */
#ifndef DEC_ARITH
        t_srnd(b_rflg);      /* set rounding mode for ten-byte arithmetic */
#endif

#if IBM_RS6000_C
#ifndef RS6000_SYNC_TRAP	/* synchronous trapping mode not available */
	e_rofe();	/* disable overflow trap         */
	e_rioe();	/* disable invalid operation trap*/
	e_rdze();	/* disable divide-by-zero trap   */
#endif /* RS6000_SYNC_TRAP */
#endif /* IBM_RS6000_C */

        /* initializations for command line scanning    */
        f_orgc = f_argc = argc;
        f_orgv = argv;
        f_argv = (char **)malloc((argc+1)*sizeof(*argv));
        memcpy(f_argv,argv,(argc+1)*sizeof(*argv));
        f_apos = 1;

        /* indicates that program parameter list is handled */
        f_pppl = TRUE;

#ifdef DEMO_VERSION
        /* display version number */
        fprintf(stderr,"%s",PRODUCT_TEXT);
        fprintf(stderr,"%s %s %s %s %s\n",VERSION_TEXT, REVISION_TEXT,
			 VERSION_SYS_TEXT, VERSION_EXT_TEXT, VERSION_ARI_TEXT );
        fprintf(stderr,"%s\n",COPYRIGHT_TEXT);
#endif
        /* scan program parameter list for special options  */
        for (i=1;i<f_argc;)
           {
           switch (b_popt(f_pmto.fp,f_argv[i]))
              {
              case 2: f_pppd = TRUE;
              case 1: for (k=i+1;k<f_argc;k++)
                         f_argv[k-1] = f_argv[k];
                      f_argc--;
                      break;
              case 0: i++;
                      break;
              default: ;
              }
           }

        return;
        }





