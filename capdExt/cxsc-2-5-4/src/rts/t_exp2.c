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

/* CVS $Id: t_exp2.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_exp2.c                              */
/*                                                              */
/*      Entries         : a_real t_exp2(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of exponential        */
/*                                                              */
/*      Description     : 2**x                                  */
/*                                                              */
/****************************************************************/
 
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
 
#define ln2                  /* ln(2)   */ \
           EXTREAL(0x3F, 0xFE,              \
                   0xB1, 0x72, 0x17, 0xF7,  \
                   0xD1, 0xCF, 0x79, 0xAB)  
/* --- Ende ln(2) Konstante --- */
/* ln(2) = 0.6931471805599453094172321214581765680755001...   */
/*  ->F7 17 72 B1   AB 79 CF D1   98 B3 E3 C9 INTEL Darst.!!! */
 
/*--------------------------------------------------------------*
 | exp2(x) = 2**x  = exp(x*ln(2))      Punktroutine             |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int t_2exp(const ExtReal *arg, ExtReal *res)
#else
int t_2exp(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int t_2exp(ExtReal *arg, ExtReal *res)
#else
int t_2exp(arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         rnd;        /* RundungsMode                       */
   int         ret;        /* Rueckgabe                          */
   int         check;      /* Rueckgabe von Makro ArgCheck       */
#ifdef ANSI_C
#if SUN4_CPP_C
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = ln2;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#else
static const ExtReal logcnst = ln2;
#endif
#else
static ExtReal logcnst;
static char    logvalue[sizeof(ExtReal)] = ln2;
 
memcpy(&logcnst,logvalue,sizeof(logcnst));
#endif
   /* --- pruefe Argument --- */
   ArgCheck1(Exp, arg, res);
 
   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);
 
   /* --- exp2() --- */
   mulee(arg, &logcnst, arg);   /* x = x*ln(2) */
   ret = (expee(arg, res) != 0);
 
   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);
 
   return ret;
} /* exp2() */


#ifdef LINT_ARGS
a_real t_exp2(a_real arg)
#else
a_real t_exp2(arg)
a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;
 
        E_SPUSH("t_exp2")
 
        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);
 
        if ((rc = t_2exp(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);
 
        setrndmode(rnd);
 
        E_SPOPP("t_exp2")
        return res;
        }
/* ------------------------------------------------------------ */
 





