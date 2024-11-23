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

/* CVS $Id: t_atnh.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_atnh.c                              */
/*                                                              */
/*      Entries         : a_real t_atnh(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of area tangent       */
/*                                                              */
/*      Description     : Area tangent function.                */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "atanh.h"      */
/* Teilpunkt               */
/* +3.125000000000000E-002 */
#define Tp                                  \
           EXTREAL(0x3F, 0xFA,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00)

/* StdFctReal(t_atnh,atanhee) */
#ifdef LINT_ARGS
a_real t_atnh(a_real arg)
#else
a_real t_atnh(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_atnh")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = atanhee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_atnh")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | atanhee                                                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int atanhee (const ExtReal *arg, ExtReal *res)
#else
int atanhee (arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int atanhee (ExtReal *arg, ExtReal *res)
#else
int atanhee (arg, res)
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
   ExtReal     an;         /* negatives Argument                */
   ExtReal     radd;       /* Ergebnis arg + 1                  */
   ExtReal     rsub;       /* Ergebnis arg - 1                  */
   ExtReal     rdiv;       /* Ergebnis radd / rsub              */
   ExtReal     rln;        /* Ergebnis ln                       */
   ExtReal     rp;         /* Ergebnis lnp1 pos                 */
   ExtReal     rn;         /* Ergebnis lnp1 neg                 */
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
   ArgCheck1(Atanh, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- |arg| >= Teilungspunkt --- */
   if(-1!=cmpabsee(arg, &tp)) {
      addee(&One, arg, &radd);
      subee(&One, arg, &rsub);
      divee(&radd, &rsub, &rdiv);
      ret = lnee(&rdiv, &rln);
   }
   /* --- |arg| < Teilungspunkt --- */
   else {
      chsee(arg, &an);
      ret = lnp1ee(arg, &rp);
      ret = lnp1ee(&an, &rn);
      subee(&rp, &rn, &rln);
   }
   /* --- res = 0.5*rln --- */
   scaliee(&rln, -1, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* atanhee() */





