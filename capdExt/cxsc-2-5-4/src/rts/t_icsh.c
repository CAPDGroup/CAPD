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

/* CVS $Id: t_icsh.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctInterval(t_icsh,icoshee) */
#ifdef LINT_ARGS
a_intv t_icsh(a_intv ai)
#else
a_intv t_icsh(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_icsh")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = icoshee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_icsh")
        return res;
        }

/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | coshee                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int icoshee (const IExtReal *arg, IExtReal *res)
#else
int icoshee (arg, res)
const IExtReal   *arg;
      IExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int icoshee (IExtReal *arg, IExtReal *res)
#else
int icoshee (arg, res)
IExtReal   *arg;
IExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   IExtReal    r;          /* Ergebnis vor Rundung              */
   int         retu;       /* Rueckgabe                         */
   int         retl;       /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument, dann Pruefung aus --- */
   ArgCheckI1(ICosh, arg, res);
   arg_check = Off;

   /* --- cosh --- */
   if(-1==cmpee(&(arg->u), &Zero)) { /* --- u<0 --- */
      retu = coshee(&(arg->l), &r.u);
      retl = coshee(&(arg->u), &r.l);
   } else
   if(1==cmpee(&(arg->l), &Zero)) { /* --- l>0 --- */
      retu = coshee(&(arg->u), &r.u);
      retl = coshee(&(arg->l), &r.l);
   } else {
      copyee(&One, &(res->l));             /* --- 0 in argument --- */
      retl = LB_Exact;
      if(!cmpee(&arg->u,&arg->l))          /* --- u = l = 0 --- */
      {
         copyee(&One, &(res->u));
         retu = UB_Exact;
      }
      else
      {
         retu = (-1!=cmpabsee(&(arg->u), &(arg->l)) ?
            coshee(&(arg->u), &r.u) :      /* --- |u|>=|l| --- */
            coshee(&(arg->l), &r.u));      /* --- |u|< |l| --- */
      }
   }

   /* --- Abbruch bei Fehler --- */
   if((retu&Except)!=NoErr || (retl&Except)!=NoErr) {
      icopyee(&r, res);
      arg_check = On;
      return max((retu&Except), (retl&Except));
   }

   /* --- Rundungs-Fehler --- */
   if(!(retu&UB_Exact))
      round_rel(UP, &r.u, &EpsCosh, &(res->u));

   if(!(retl&LB_Exact))
   {
      if(cmpee(&r.l,&One))
         round_rel(DOWN, &r.l, &EpsCosh, &(res->l));
      else
         copyee(&One, &(res->l));
   }

   /* --- ArgPruefung wieder an --- */
   arg_check = On;

   return NoErr;
} /* icoshee() */





