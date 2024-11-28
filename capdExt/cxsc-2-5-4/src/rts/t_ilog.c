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

/* CVS $Id: t_ilog.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval(t_ilog,ilnee) */
#ifdef LINT_ARGS
a_intv t_ilog(a_intv ai)
#else
a_intv t_ilog(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_ilog")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = ilnee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_ilog")
        return res;
        }

/* ------------------------------------------------------------ */

#ifdef LINT_ARGS
#ifdef ANSI_C
int round_ln(int rnd, int rr, const ExtReal *arg, ExtReal *res);
#else  /* NOT ANSI_C */
int round_ln(int rnd, int rr, ExtReal *arg, ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
int round_ln();
#endif /* LINT_ARGS */

/*--------------------------------------------------------------*
 | Intervall Ln                                                 |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int ilnee (const IExtReal *arg, IExtReal *res)
#else
int ilnee (arg, res)
const IExtReal   *arg;
      IExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int ilnee (IExtReal *arg, IExtReal *res)
#else
int ilnee (arg, res)
IExtReal   *arg;
IExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   IExtReal    r;          /* Ergebnis vor Rundung              */
   int         rnd;        /* RundungsMode                      */
   int         rru;        /* Rueckgabe RundungsArt             */
   int         rrl;        /* Rueckgabe RundungsArt             */
   int         retu;       /* Rueckgabe                         */
   int         retl;       /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI1(ILn, arg, res);
   arg_check = Off;

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- ln --- */
   retu = _s_ln(&(arg->u), &r.u, &rru);
   retl = _s_ln(&(arg->l), &r.l, &rrl);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   /* --- Abbruch bei Fehler --- */
   if(retu!=NoErr || retl!=NoErr) {
      icopyee(&r, res);
      arg_check = On;
      return max(retu, retl);
   }

   /* --- Rundungs-Fehler --- */
   round_ln(UP,   rru, &(r.u), &(res->u));
   round_ln(DOWN, rrl, &(r.l), &(res->l));

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   /* --- kein Fehler mehr moeglich --- */
   return NoErr;
} /* ilnee() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int round_ln(int rnd, int rr, const ExtReal *arg, ExtReal *res)
#else
int round_ln(rnd, rr, arg, res)
      int      rnd;
      int      rr;
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int round_ln(int rnd, int rr, ExtReal *arg, ExtReal *res)
#else
int round_ln(rnd, rr, arg, res)
int      rnd;
int      rr;
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   switch(rr){
   case 1:                              /* im ApproxIntervall   */
      round_rel(rnd, arg, &EpsLn, res);
      break;
   case 2:                              /* ausserhalb ApproxIntervall   */
      sround(rnd, arg, &EpsLnRel, &EpsLnAbs, res);
      break;
   }
   return NoErr;
} /* round_ln() */





