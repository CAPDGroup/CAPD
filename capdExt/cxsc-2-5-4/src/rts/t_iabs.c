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

/* CVS $Id: t_iabs.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_iabs.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctReal(t_abs, absee) */
#ifdef LINT_ARGS
a_real t_abs(a_real arg)
#else
a_real t_abs(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_abs")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = absee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_abs")
        return res;
        }

/* StdFctInterval(t_iabs,iabsee) */
#ifdef LINT_ARGS
a_intv t_iabs(a_intv ai)
#else
a_intv t_iabs(ai)

a_intv ai;
#endif
        {

        a_intv   res;
        int      rnd, rc;
        IExtReal  a, r;

        E_SPUSH("t_iabs")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&ai.INF, &a.l);
        longreal_to_extreal((LongReal *)&ai.SUP, &a.u);

        if ((rc = iabsee(&a, &r))!=0)
            ieee_aborti1(rc, &ai);

        setrndmode(DOWN);
        if ((rc = extreal_to_longreal(&r.l, (LongReal *)&res.INF))!=0)
             ieee_aborti1(rc, &ai);
        setrndmode(UP);
        if ((rc = extreal_to_longreal(&r.u, (LongReal *)&res.SUP))!=0)
            ieee_aborti1(rc, &ai);
        setrndmode(rnd);

        E_SPOPP("t_iabs")
        return res;
        }

/*--------------------------------------------------------------*
 | iabsee                                                       |
 *--------------------------------------------------------------*/
/*
**  geaendert 030990 Baumhof:
**  iabsee berechnete das Intervall von |arg.l| bis |arg.u|,
**  geaendert, so dass {|x|: x in arg} zurueckgegeben wird.
*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int iabsee (const IExtReal *arg, IExtReal *res)
#else
int iabsee (arg, res)
const IExtReal  *arg;
      IExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int iabsee (IExtReal *arg, IExtReal *res)
#else
int iabsee (arg, res)
IExtReal  *arg;
IExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int   ret;

   if (cmpee(&Zero, &(arg->l)) <= 0)            /* arg.l >= 0 */
   {
      ret  = absee(&(arg->u), &(res->u));
      ret += absee(&(arg->l), &(res->l));
   }
   else if (cmpee(&Zero, &(arg->u)) >= 0)       /* arg.u <= 0 */
   {
      ret  = absee(&(arg->u), &(res->l));
      ret += absee(&(arg->l), &(res->u));
   }
   else                                  /* arg.l < 0 < arg.u */
   {
      if (cmpabsee(&(arg->u), &(arg->l)) >= 0)  /* |u| >= |l| */
         ret = absee(&(arg->u), &(res->u));
      else
         ret = absee(&(arg->l), &(res->u));
      copyee(&Zero, &(res->l));
   }

   return ret;
} /* iabsee() */





