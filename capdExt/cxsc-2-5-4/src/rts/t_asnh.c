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

/* CVS $Id: t_asnh.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_asnh.c                              */
/*                                                              */
/*      Entries         : a_real t_asnh(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of area sine          */
/*                                                              */
/*      Description     : Area sine function.                   */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "asinh.h"      */
/* Teilpunkt               */
/* +3.125000000000000E-002 */
#define Tp                                  \
           EXTREAL(0x3F, 0xFA,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00)

/* StdFctReal(t_asnh,asinhee) */
#ifdef LINT_ARGS
a_real t_asnh(a_real arg)
#else
a_real t_asnh(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_asnh")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = asinhee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_asnh")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | asinh                                                        |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int asinhee (const ExtReal *arg, ExtReal *res)
#else
int asinhee (arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int asinhee (ExtReal *arg, ExtReal *res)
#else
int asinhee (arg, res)
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
   static const ExtReal tp=Tp; /* Teilungspunkt 0.03125         */
#endif
#else
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#endif
   ExtReal     absa;       /* absolutes Argument                */
   ExtReal     sqra;       /* arg**2                            */
   ExtReal     sqrap1;     /* sqra + 1                          */
   ExtReal     rsqrt;      /* Ergebnis sqrt                     */
   ExtReal     radd;       /* Ergebnis add nach sqrt            */
   int         sgn;        /* Vorzeichen des Argumentes arg     */
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
   ArgCheck1(Asinh, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- absa := |a| --- */
   sgn = SGNE(arg);
   absee(arg, &absa);

   /* --- sqrt(arg**2+1) --- */
   mulee(arg, arg, &sqra);

   /* --- arg**2 >= Teilungspunkt --- */
   if(-1!=cmpee(&sqra, &tp)) {
      addee(&sqra, &One, &sqrap1);
      sqrtee(&sqrap1, &rsqrt);

      /* --- asinh = ln(arg+rsqrt) --- */
      addee(&absa, &rsqrt, &radd);
      ret = lnee(&radd, res);
   }
   /* --- arg < Teilungspunkt --- */
   else {
      sqrtm1ee(&sqra, &rsqrt);
      addee(&absa, &rsqrt, &radd);
      ret = lnp1ee(&radd, res);
   }

   /* --- if arg < 0 ==> res := -res --- */
   if(NEG==sgn)
      chsee(res, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* asinhee() */





