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

/* CVS $Id: t_tanh.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_tanh.c                              */
/*                                                              */
/*      Entries         : a_real t_tanh(arg)                    */
/*                        a_real arg;                           */
/*                                                              */
/*      Arguments       : arg  = argument of hyp. tangent       */
/*                                                              */
/*      Description     : Hyperbolic tangent function.          */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctReal(t_tanh,tanhee) */
#ifdef LINT_ARGS
a_real t_tanh(a_real arg)
#else
a_real t_tanh(arg)

a_real arg;
#endif
        {
        int      rnd, rc;
        a_real   res;
        ExtReal  a, r;

        E_SPUSH("t_tanh")

        rnd = getrndmode();
        longreal_to_extreal((LongReal *)&arg, &a);

        if ((rc = tanhee(&a, &r))!=0
            || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
           ieee_abortr1(rc, &arg);

        setrndmode(rnd);

        E_SPOPP("t_tanh")
        return res;
        }
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | tanh                                                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int tanhee(const ExtReal *arg, ExtReal *res)
#else
int tanhee(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int tanhee(ExtReal *arg, ExtReal *res)
#else
int tanhee(arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     r;          /* Ergebnis Coth                     */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument --- */
   ArgCheck1(Tanh, arg, res);

   /* --- Null extra --- */
   if(0==cmpee(arg, &Zero))
      ret = copyee(&Zero, res);

   /* --- Tanh --- */
   else {
      ret = cothee(arg, &r);
      divee(&One, &r, res);
   }

   return ret;
} /* tanhee() */





