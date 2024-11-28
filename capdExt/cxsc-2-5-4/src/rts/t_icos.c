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

/* CVS $Id: t_icos.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_icos.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval(t_icos,icosee) */
#ifdef LINT_ARGS
a_intv t_icos(a_intv ai)
#else
a_intv t_icos(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_icos")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = icosee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_icos")
        return res;
        }

/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | icos                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int icosee(const IExtReal *arg, IExtReal *res)
#else
int icosee(arg, res)
const IExtReal *arg;
      IExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int icosee(IExtReal *arg, IExtReal *res)
#else
int icosee(arg, res)
IExtReal *arg;
IExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   IDReal      v;          /* Produkt des Arguments mit 2/pi    */
   IExtReal    j;          /* Ganzzahliger Anteil von v         */
   int         jmod4u;     /* j.u modulo 4                      */
   int         jmod4l;     /* j.l modulo 4                      */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI1(ICos, arg, res);
   arg_check = Off;

   /* --- ganzzahliger Anteil von arg --- */
   gza_trg(&(arg->u), J_Init_Cos, Period_PiHalf, &v.u, &j.u, &jmod4u);
   gza_trg(&(arg->l), J_Init_Cos, Period_PiHalf, &v.l, &j.l, &jmod4l);

   /* --- Reduktion und Cosinus --- */
   if(( ret = isincos(&v, &j, jmod4u, jmod4l, res)) & Except)
      if((ret&Except)!=PLOSS){
         arg_check = On;
         return ret & Except;
      }

   /* nicht runden, falls Argument Null ist */
   if (!cmpee(&(arg->l), &Zero)) ret |= UB_Exact;
   if (!cmpee(&(arg->u), &Zero)) ret |= UB_Exact;

   /* --- Rundungs-Fehler --- */
   if(!(ret&UB_Exact)) round_rel(UP,   &(res->u), &EpsCos, &(res->u));
   if(!(ret&LB_Exact)) round_rel(DOWN, &(res->l), &EpsCos, &(res->l));

   /* Ergebnisverbesserung in Spezialfaellen */
#if INT_HPREC
   {
     if (1==cmpee(&(res->u), &One))             /* cos x > 1  */
        copyee(&One, &(res->u));
     if (1!=cmpee(&(res->l), &MinusOne))        /* cos x <= -1 */
        { copyee(&MinusOne, &(res->l)); }
   }
#endif /* INT_HPREC */

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   /* --- kein Fehler mehr moeglich --- */
   return NoErr;
} /* icos() */





