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

/* CVS $Id: t_acsh.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_acsh.c                              */
/*                                                              */
/*      Entries         : a_real t_acsh(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of area cosine hyp.   */
/*                                                              */
/*      Description     : Area cosine hyp. function             */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "acosh.h"      */
/* Teilungspunkt           */
/* +1.025000000000000E+000 */
#define Tp  EXTREAL(0x3F, 0xFF, 0x83, 0x33, 0x33, 0x33,  \
                    0x33, 0x33, 0x33, 0x34)

/* StdFctReal(t_acsh,acoshee) */
#ifdef LINT_ARGS
a_real t_acsh(a_real arg)
#else
a_real t_acsh(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_acsh")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = acoshee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_acsh")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | acoshee                                                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int acoshee (const ExtReal *arg, ExtReal *res)
#else
int acoshee (arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int acoshee (ExtReal *arg, ExtReal *res)
#else
int acoshee (arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#else
   static const ExtReal tp=Tp; /* Teilungspunkt 1.025           */
#endif
#else
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#endif
   ExtReal     am1;        /* Argument - 1                      */
   ExtReal     ap1;        /* Argument + 1                      */
   ExtReal     a2m1;       /* Argumetn**2 - 1                   */
   ExtReal     rsqrt;      /* Ergebnis Sqrt                     */
   ExtReal     radd;       /* Ergebnis Add                      */
   int         rnd;        /* RundungsMode                      */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

#ifndef ANSI_C
   memcpy(&tp, tpc, sizeof(tp));
#else
#if SUN4_CPP_C
   memcpy(&tp, tpc, sizeof(tp));
#endif
#endif

   /* --- pruefe Argument --- */
   ArgCheck1(Acosh, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- sqrt(arg**2-1) --- */
   subee(arg, &One, &am1);
   addee(arg, &One, &ap1);
   mulee(&am1, &ap1, &a2m1);
   sqrtee(&a2m1, &rsqrt);

   /* --- arg >= Teilungspunkt --- */
   if(-1!=cmpee(arg, &tp)) {
      /* --- acosh = ln(arg+rsqrt) --- */
      addee(arg, &rsqrt, &radd);
      ret = lnee(&radd, res);
   }
   /* --- |arg| < Teilungspunkt --- */
   else {
      /* --- acosh = lnp1(am1+rsqrt) --- */
      addee(&am1, &rsqrt, &radd);
      ret = lnp1ee(&radd, res);
   }

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret ;
} /* acoshee() */





