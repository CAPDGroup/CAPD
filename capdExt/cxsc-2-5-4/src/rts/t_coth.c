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

/* CVS $Id: t_coth.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_coth.c                              */
/*                                                              */
/*      Entries         : a_real t_coth(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of hyp. cotangent     */
/*                                                              */
/*      Description     : Hyperbolic cotangent function.        */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "coth.h"       */
/* Teilpunkt               */
/* +2.287500000000000E+001 */
#define Tp                                  \
           EXTREAL(0x40, 0x03,              \
                   0xB7, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00)

/* StdFctReal(t_coth,cothee) */
#ifdef LINT_ARGS
a_real t_coth(a_real arg)
#else
a_real t_coth(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_coth")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = cothee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_coth")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | coth                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int cothee(const ExtReal *arg, ExtReal *res)
#else
int cothee(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int cothee(ExtReal *arg, ExtReal *res)
#else
int cothee(arg, res)
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
   static const ExtReal tp = Tp; /* Teilpunkt                   */
#endif
#else
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#endif
   ExtReal     a;          /* Absolutes Argument                */
   ExtReal     e;          /* ZwischenErgebnis 2*a              */
   ExtReal     exp;        /* exp(e)                            */
   ExtReal     n;          /* ZwischenErgebnis Nenner           */
   ExtReal     q;          /* ZwischenErgebnis Quotient         */
   int         sgn;        /* Vorzeichen Argument               */
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
   ArgCheck1(Coth, arg, res);

   /* --- CheckStatus sichern --- */
   check = arg_check;
   arg_check = Off;

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- a := |arg| --- */
   sgn = SGNE(arg);
   absee(arg, &a);

   /* --- arg > Teilungspunkt --- */
   if(1==cmpee(&a, &tp))
      ret = copyee(&One, res);                  /* |coth| = 1   */

   /* --- arg <= Teilungspunkt --- */
   else {
      /* --- e = 2*a --- */
      scaliee(&a, 1, &e);

      if(-1==cmpee(&e, &Ln2))
         ret = expm1ee(&e, &n);                 /* e**(2*a)-1   */
      else {
         ret = expee(&e, &exp);                 /* e**(2*a)     */
         subee(&exp, &One, &n);
      }
      divee(&Two, &n, &q);
      addee(&One, &q, res);
   }

   /* --- if arg < 0 ==> res := -res --- */
   if(NEG==sgn)
      chsee(res, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   /* --- ArgumentCheck zuruecksetzen --- */
   arg_check = check;

   return ret;
} /* cothee() */





