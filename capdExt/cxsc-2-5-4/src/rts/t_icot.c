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

/* CVS $Id: t_icot.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval(t_icot,icotee) */
#ifdef LINT_ARGS
a_intv t_icot(a_intv ai)
#else
a_intv t_icot(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_icot")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = icotee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_icot")
        return res;
        }

/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | icotee                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int icotee (const IExtReal *arg, IExtReal *res)
#else
int icotee (arg, res)
const IExtReal   *arg;
      IExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int icotee (IExtReal *arg, IExtReal *res)
#else
int icotee (arg, res)
IExtReal   *arg;
IExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DReal       vu;         /* Produkt des Arguments mit 2/pi  upper Bound */
   DReal       vl;         /* Produkt des Arguments mit 2/pi  lower Bound */
   ExtReal     ju;         /* Ganzzahliger Anteil von v       upper Bound */
   ExtReal     jl;         /* Ganzzahliger Anteil von v       lower Bound */
   ExtReal     tu;         /* reduziertes Argument upper Bound            */
   ExtReal     tl;         /* reduziertes Argument lower Bound            */
   ExtReal     jd;         /* Differenz j.u - j.l, max(jdiff):=4          */
   int         jdiff;      /* =jd                                         */
   int         jmod4u;     /* j.u modulo 4                                */
   int         jmod4l;     /* j.l modulo 4                                */
   int         retru;      /* Rueckgabe Reduktion                         */
   int         retrl;      /* Rueckgabe Reduktion                         */
   int         retu;       /* Rueckgabe tancot                            */
   int         retl;       /* Rueckgabe tancot                            */
   int         ret;        /* Rueckgabe                                   */
   int         check;      /* Rueckgabe von Makro ArgCheck                */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI1(ICot, arg, res);
   arg_check = Off;

   /* --- ganzzahliger Anteil von arg --- */
   gza_trg(&(arg->u), J_Init_Cot, Period_PiQuart, &vu, &ju, &jmod4u);
   gza_trg(&(arg->l), J_Init_Cot, Period_PiQuart, &vl, &jl, &jmod4l);

   /* --- jdiff --- */
   subee(&ju, &jl, &jd);
   if (1 == cmpee(&jd, &Four)) copyee(&Four, &jd);
   extreal_to_int(&jd, &jdiff);

   /* --- ArgumentPruefung auf Pol, Abbruch --- */
   if ((jdiff == 4) || (jmod4l <= 1 && jmod4u >= 2)) {
      ret = exc_handle_i1(ICot, ExcISing, arg, res);
      arg_check = On;
      return ret;
   }

   /* --- Reduktion  --- */
   retru = red_trg(&vu, &ju, jmod4u, &tu);
   retrl = red_trg(&vl, &jl, jmod4l, &tl);

   /* --- tan(t*pi/4) --- */
   retu = tancot(&tu, jmod4u, &(res->u));
   retl = tancot(&tl, jmod4l, &(res->l));

   /* --- Abbruch bei Fehler --- */
   if(NoErr!=retu||NoErr!=retl) {
      arg_check = On;
      return max(retu, retl);
   }

   /* --- cot=-tan --- */
   ichsee(res, res);

   /* --- Rundungs-Fehler --- */
   iround_rel(res, &EpsCot, res);

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   /* --- kein Fehler mehr moeglich, nur PLOSS von red_trg --- */
   return max(retru, retrl);
} /* icotee() */





