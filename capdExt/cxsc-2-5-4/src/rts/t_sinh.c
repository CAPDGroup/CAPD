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

/* CVS $Id: t_sinh.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_sinh.c                              */
/*                                                              */
/*      Entries         : a_real t_sinh(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of hyp. sine          */
/*                                                              */
/*      Description     : Hyperbolic sine function.             */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "sinh.h"       */
/* Teilpunkt               */
/* +1.447715759277344E-001 */
#define Tp                                  \
           EXTREAL(0x3F, 0xFC,              \
                   0x94, 0x3F, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00)

/* StdFctReal(t_sinh,sinhee) */
#ifdef LINT_ARGS
a_real t_sinh(a_real arg)
#else
a_real t_sinh(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_sinh")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = sinhee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_sinh")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | sinh                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int sinhee(const ExtReal *arg, ExtReal *res)
#else
int sinhee(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int sinhee(ExtReal *arg, ExtReal *res)
#else
int sinhee(arg, res)
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
   static const ExtReal tp = Tp; /* Teilpunkt   */
#endif
#else
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#endif
   ExtReal     a;          /* Absolutes Argument                */
   ExtReal     e;          /* Ergebnis Exp                      */
   ExtReal     r;          /* Ergebnis Summe                    */
   ExtReal     ep;         /* Ergebnis positiver Exponent       */
   ExtReal     en;         /* Ergebnis negativer Exponent       */
   ExtReal     n;          /* ZwischenErgebnis Nenner           */
   ExtReal     q;          /* ZwischenErgebnis Quotient         */
   int         sgn;        /* Vorzeichen Argument               */
   int         rnd;        /* RundungsMode                      */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

#ifndef ANSI_C
   memcpy(&tp, tpc, sizeof(tp));
#else
#if SUN4_CPP_C
   memcpy(&tp, tpc, sizeof(tp));
#endif
#endif

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheck1(Sinh, arg, res);
   arg_check = Off;

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- a := |arg| --- */
   sgn = SGNE(arg);
   absee(arg, &a);

   /* --- arg > Teilungspunkt --- */
   if(1==cmpee(&a, &tp)) {
      expee(&a, &ep);
      divee(&One, &ep, &en);
      subee(&ep, &en, &r);
   }
   else {
      expm1ee(&a, &e);
      addee(&e, &One, &n);
      divee(&e, &n, &q);
      addee(&e, &q, &r);
   }
   scaliee(&r, -1, res);         /* /2                          */

   /* --- if arg < 0 ==> res := -res --- */
   if(NEG==sgn)
      chsee(res, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);
   arg_check = On;

   return NoErr;
} /* sinhee() */





