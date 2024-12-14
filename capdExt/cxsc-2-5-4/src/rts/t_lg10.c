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

/* CVS $Id: t_lg10.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_lg10.c                              */
/*                                                              */
/*      Entries         : a_real t_lg10(arg)                    */
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
 
#define log10e      /* log10(e) */          \
           EXTREAL(0x3F, 0xFD,              \
                   0xDE, 0x5B, 0xD8, 0xA9,  \
                   0x37, 0x28, 0x71, 0x95)  /* +0.43429448190325182765 */ \
/* ->A9 D8 5B DE  95 71 28 37  AF AA 5B 35   INTEL Darst.!!      */
 
#define log10    llog10 

/*--------------------------------------------------------------*
 | Punkt log10                                                  |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int log10(const ExtReal *arg, ExtReal *res)
#else
static int log10(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int log10(ExtReal *arg, ExtReal *res)
#else
static int log10(arg, res)
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
static char    logvalue[sizeof(ExtReal)] = log10e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#else
static const ExtReal logcnst = log10e;
#endif
#else
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = log10e;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#endif
   /* --- pruefe Argument --- */
   ArgCheck1(Ln, arg, res);
 
   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);
 
   /* --- Log10 --- */
   ret = _s_ln(arg, res, &rr);
   mulee(res, &logcnst, res);
   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);
 
   return ret;
} /* log10() */


/*-----------------------------------------------------------------*
 | log10(x) = ln(x)*log10(e)                                       |
 *-----------------------------------------------------------------*/
#ifdef LINT_ARGS
a_real t_lg10(a_real arg)
#else
a_real t_lg10(arg)
 
a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_lg10")
 
        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);
 
        if ((rc = log10(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);
 
        setrndmode(rnd);

        E_SPOPP("t_lg10")
 
        return res;
        }
/* ------------------------------------------------------------ */
 





