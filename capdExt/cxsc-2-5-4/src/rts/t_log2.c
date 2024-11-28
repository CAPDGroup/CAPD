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

/* CVS $Id: t_log2.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_log2.c                              */
/*                                                              */
/*      Entries         : a_real t_log2(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of log10              */
/*                                                              */
/*      Description     : Logarithm to base 10                  */
/*                                                              */
/****************************************************************/
 
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
 
#define log2e      /* log2(e) */           \
           EXTREAL(0x3f, 0xff,              \
                   0xB8, 0xAA, 0x3B, 0x29,  \
                   0x5C, 0x17, 0xF0, 0xBB)  /* +1.44269504088896340735*/\
 
/* ->29 3B AA B8   BB F0 17 5C   D0 FE 87 BE INTEL Darst>!!*/
 
/*--------------------------------------------------------------*
 | Punkt log2                                                   |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int tt_log2(const ExtReal *arg, ExtReal *res)
#else
static int tt_log2(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int tt_log2(ExtReal *arg, ExtReal *res)
#else
static int tt_log2(arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         rnd;        /* RundungsMode                       */
   int         rr;         /* Rueckgabe RundungsArt, hier dummy  */
   int         ret;        /* Rueckgabe                          */
   int         check;      /* Rueckgabe von Makro ArgCheck       */
#ifdef ANSI_C
#if SUN4_CPP_C
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = log2e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#else
static const ExtReal logcnst = log2e;
#endif
#else
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = log2e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#endif
   /* --- pruefe Argument --- */
   ArgCheck1(Ln, arg, res);
 
   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);
 
   /* --- Log2 --- */
   ret = _s_ln(arg, res, &rr);
   mulee(res, &logcnst, res);
   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);
 
   return ret;
} /* log2() */


/*-----------------------------------------------------------------*
 | log2(x) = ln(x)*log2(e)                                         |
 *-----------------------------------------------------------------*/
#ifdef LINT_ARGS
a_real t_log2(a_real arg)
#else
a_real t_log2(arg)
 
a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;
 
        E_SPUSH("t_log2")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);
 
        if ((rc = tt_log2(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);
 
        setrndmode(rnd);

        E_SPOPP("t_log2")
 
        return res;
        }
/* ------------------------------------------------------------ */
 





