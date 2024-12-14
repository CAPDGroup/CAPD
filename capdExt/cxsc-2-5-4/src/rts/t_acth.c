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

/* CVS $Id: t_acth.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_acth.c                              */
/*                                                              */
/*      Entries         : a_real t_acth(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of area cotangent     */
/*                                                              */
/*      Description     : Area cotangent function               */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "acoth.h"      */
/* Teilpunkt               */
/* +8.000000000000000E+000 */
#define Tp                                 \
           EXTREAL(0x40, 0x02,             \
                   0x80, 0x00, 0x00, 0x00, \
                   0x00, 0x00, 0x00, 0x00)

/* StdFctReal(t_acth,acothee) */
#ifdef LINT_ARGS
a_real t_acth(a_real arg)
#else
a_real t_acth(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_acth")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = acothee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_acth")
        return res;
        }
/* ------------------------------------------------------------ */

/*----------------------------------------------------------------*
 | acothee                                                        |
 *----------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int acothee (const ExtReal *arg, ExtReal *res)
#else
int acothee (arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int acothee (ExtReal *arg, ExtReal *res)
#else
int acothee (arg, res)
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
   static const ExtReal tp=Tp; /* Teilungspunkt = 8     */
#endif
#else
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#endif
   ExtReal     radd;       /* Ergebnis arg + 1          */
   ExtReal     rsub;       /* Ergebnis arg - 1                  */
   ExtReal     rdiv;       /* Ergebnis radd / rsub              */
   ExtReal     rln;        /* Ergebnis ln                       */
   int         rnd;        /* RundungsMode                      */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

#ifndef ANSI_C
   memcpy(&tp, tpc, sizeof(tp));
#else
#if SUN4_CPP_C
   memcpy(&tp, tpc, sizeof(tp));
#endif
#endif

   /* --- pruefe Argument --- */
   ArgCheck1(Acoth, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- rsub = arg - 1 --- */
   subee(arg, &One, &rsub);

   /* --- |arg| <= Teilungspunkt --- */
   if(1!=cmpabsee(&rsub, &tp)) {
      addee(arg, &One, &radd);
      divee(&radd, &rsub, &rdiv);
      lnee(&rdiv, &rln);
   }
   /* --- |arg| > Teilungspunkt --- */
   else {
      divee(&Two, &rsub, &rdiv);
      lnp1ee(&rdiv, &rln);
   }

   /* --- res = 0.5*rln --- */
   scaliee(&rln, -1, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return NoErr;
} /* acothee() */





