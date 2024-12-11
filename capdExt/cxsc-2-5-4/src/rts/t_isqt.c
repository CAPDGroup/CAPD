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

/* CVS $Id: t_isqt.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval(t_isqt,isqrtee) */
#ifdef LINT_ARGS
a_intv t_isqt(a_intv ai)
#else
a_intv t_isqt(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_isqt")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = isqrtee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_isqt")
        return res;
        }

/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | isqrtee                                                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int isqrtee (const IExtReal *arg, IExtReal *res)
#else
int isqrtee (arg, res)
const IExtReal   *arg;
      IExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int isqrtee (IExtReal *arg, IExtReal *res)
#else
int isqrtee (arg, res)
IExtReal   *arg;
IExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         retu;       /* Rueckgabe                             */
   int         retl;       /* Rueckgabe                             */
   int         rnd;        /* RundungsMode                          */
   int         check;      /* Rueckgabe von Makro ArgCheck          */
   ExtReal     tmp1, tmp2; /* Test auf exaktes Ergebnis             */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI1(ISqrt, arg, res);
   arg_check = Off;

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* Punktargument ? */
   if (!cmpee(&(arg->u), &(arg->l)))
   {
      if (NoErr!=(retu = _s_sqrt(&(arg->u), &(res->u))))
      {
         arg_check = On;
         return(retu);
      }
      copyee(&(res->u), &(res->l));
      setrndmode(UP);
      mulee(&(res->u), &(res->u), &tmp1);
      setrndmode(DOWN);
      mulee(&(res->u), &(res->u), &tmp2);
      setrndmode(NEAR);
      /* ist das Ergebnis exakt ? */
      if (0!=cmpee(&tmp1, &tmp2) || 0!=cmpee(&tmp1, &(arg->u)))
      {
         /* nein, runden */
         if (cmpee(&tmp1, &(arg->u))<=0)
            round_rel(UP, &(res->u), &EpsSqrt, &(res->u));
         else if (cmpee(&tmp2, &(arg->u))>=0)
            round_rel(DOWN, &(res->l), &EpsSqrt, &(res->l));
         else
            iround_rel(res, &EpsSqrt, res);
      }
   }
   else   /* echtes Intervall als Argument */
   {

      /* --- Wurzel --- */
      retu  = _s_sqrt(&(arg->u), &(res->u));
      retl  = _s_sqrt(&(arg->l), &(res->l));

      /* --- Abbruch bei Fehler --- */
      if(NoErr!=retu||NoErr!=retl) {
         arg_check = On;
         return max(retu, retl);
      }

      /* sind linke oder rechte Grenze exakt ? */
      setrndmode(UP);
      mulee(&(res->l), &(res->l), &tmp1);
      setrndmode(DOWN);
      mulee(&(res->l), &(res->l), &tmp2);
      setrndmode(NEAR);
      /* ist das Ergebnis exakt ? */
      if (0!=cmpee(&tmp1, &tmp2) || 0!=cmpee(&tmp1, &(arg->l)))
      {
         /* nein, runden */
         round_rel(DOWN, &(res->l), &EpsSqrt, &(res->l));
      }
      setrndmode(UP);
      mulee(&(res->u), &(res->u), &tmp1);
      setrndmode(DOWN);
      mulee(&(res->u), &(res->u), &tmp2);
      setrndmode(NEAR);
      /* ist das Ergebnis exakt ? */
      if (0!=cmpee(&tmp1, &tmp2) || 0!=cmpee(&tmp1, &(arg->u)))
      {
         /* nein, runden */
         round_rel(UP, &(res->u), &EpsSqrt, &(res->u));
      }
   }

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   /* --- kein Fehler mehr moeglich --- */
   return NoErr;
} /* isqrtee() */





