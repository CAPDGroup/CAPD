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

/* CVS $Id: t_sin.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_sin.c                               */
/*                                                              */
/*      Entries         : a_real t_sin(arg)                     */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of sine               */
/*                                                              */
/*      Description     : Sine function.                        */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctReal(t_sin,sinee) */
#ifdef LINT_ARGS
a_real t_sin(a_real arg)
#else
a_real t_sin(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_sin")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = sinee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_sin")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | Punkt Sinus                                                  |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int sinee(const ExtReal *arg, ExtReal *res)
#else
int sinee(arg, res)
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int sinee(ExtReal *arg, ExtReal *res)
#else
int sinee(arg, res)
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DReal       v;          /* Produkt des Arguments mit 2/pi    */
   ExtReal     j;          /* Ganzzahliger Anteil von v         */
   int         jmod4;      /* j modulo 4                        */
   ExtReal     t;          /* reduziertes Argument              */
   int         ret;        /* Rueckgabe                         */
   int         retr;       /* Rueckgabe Reduktion               */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument --- */
   ArgCheck1(Sin, arg, res);

   /* --- ganzzahliger Anteil von arg --- */
   gza_trg(arg, J_Init_Sin, Period_PiHalf, &v, &j, &jmod4);

   /* --- Reduktion --- */
   retr = red_trg(&v, &j, jmod4, &t);

   /* --- Sinus --- */
   ret  = sincos(&t, res);

   /* --- Rueckgabe ohne RundungsInfo --- */
   return (ret!=NoErr?ret:retr);
} /* sinee() */





