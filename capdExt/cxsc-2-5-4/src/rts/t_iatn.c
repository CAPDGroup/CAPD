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

/* CVS $Id: t_iatn.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval(t_iatn,iatanee) */
#ifdef LINT_ARGS
a_intv t_iatn(a_intv ai)
#else
a_intv t_iatn(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_iatn")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = iatanee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_iatn")
        return res;
        }

/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | iatan                                                        |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int iatanee (const IExtReal *arg, IExtReal *res)
#else
int iatanee (arg, res)
const IExtReal   *arg;
      IExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int iatanee (IExtReal *arg, IExtReal *res)
#else
int iatanee (arg, res)
IExtReal   *arg;
IExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   IExtReal    r;          /* Ergebnis                          */
   int         retu;       /* Rueckgabe                         */
   int         retl;       /* Rueckgabe                         */
   int         rnd;        /* RundungsMode                      */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI1(IAtan, arg, res);
   arg_check = Off;

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- Arcus Tangens --- */
   retu = _s_atan(&(arg->u), &r.u);
   retl = _s_atan(&(arg->l), &r.l);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   /* --- Abbruch bei Fehler --- */
   if(retu!=NoErr || retl!=NoErr) {
      icopyee(&r, res);
      arg_check = On;
      return max(retu, retl);
   }

   /* --- Rundungs-Fehler --- */
   iround_rel(&r, &EpsATan, res);

   /* Ergebnisverbesserung in Spezialfaellen */
#if INT_HPREC
   { int vz;
     if (1==cmpee(&(res->u), &(IPiHalf.u)))       /* arctan x < pi/2 */
        copyee(&(IPiHalf.u), &(res->u));
     if (1==cmpabsee(&(res->l), &(IPiHalf.u)))    /* arctan x > -pi/2 */
        { copyee(&(IPiHalf.u), &(res->u));
          chsee(&(res->u), &(res->u)); }
     vz = SGNE(&(arg->u));
     if (vz==POS&&-1==cmpee(&(arg->u), &(res->u)))/* arctan x <= x */
        copyee(&(arg->u), &(res->u));
     vz = SGNE(&(arg->l));
     if (vz==NEG&&1==cmpee(&(arg->l), &(res->l))) /* arctan x >= x */
        copyee(&(arg->l), &(res->l));
   }
#endif /* INT_HPREC */

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   /* --- kein Fehler mehr moeglich --- */
   return NoErr;
} /* iatanee() */






