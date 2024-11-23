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

/* CVS $Id: t_sqrt.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_sqrt.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* StdFctReal(t_sqrt,sqrtee) */     /* siehe t_sqrt unten, */
                                    /* wegen Rundungsproblemen */
/* ------------------------------------------------------------ */

/*--------------------------------------------------------------*
 | Punkt Wurzel                                                 |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int sqrtee (const ExtReal *arg, ExtReal *res)
#else
int sqrtee (arg, res)
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int sqrtee (ExtReal *arg, ExtReal *res)
#else
int sqrtee (arg, res)
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         rnd;        /* RundungsMode                      */
   int         ret;        /* Rueckgabe                         */
   int         check;      /* Rueckgabe von Makro ArgCheck      */


   /* --- pruefe Argument, evtl. Abbruch --- */
   ArgCheck1(Sqrt, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- Wurzel --- */
   ret = _s_sqrt(arg, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);
   return ret;
} /* sqrtee() */

/*--------------------------------------------------------------*/

#ifdef LINT_ARGS
a_real t_sqrt(a_real arg)
#else
a_real t_sqrt(arg)
a_real arg;
#endif
{
    int rnd, rc;
    a_real res;
    ExtReal a, r, rsucc, rpred, tmp;
    a_real lmin, rs, rp;

    E_SPUSH("t_sqrt")

    rnd = getrndmode();
    setrndmode(NEAR);
    longreal_to_extreal((LongReal *)&arg, &a);

    if ((rc = sqrtee(&a, &r))!=0
        || (rc = extreal_to_longreal(&r, (LongReal *)&res))!=0)
         ieee_abortr1(rc, &arg);

    /* Rundung nachpruefen */
    if (cmpee(&r, &Zero) != 0)
    {
       extreal_to_longreal(&LongRealMin, (LongReal *)&lmin);
      
       R_ASSIGN(rs,r_addu(res, lmin));
       R_ASSIGN(rp,r_subd(res, lmin));
       
       longreal_to_extreal((LongReal *)&res, &r);
       longreal_to_extreal((LongReal *)&rs,  &rsucc);
       longreal_to_extreal((LongReal *)&rp,  &rpred);
       setrndmode(DOWN);
       mulee(&r, &rpred, &tmp);
       if (cmpee(&tmp, &a) < 0)
       {
          mulee(&r, &rsucc, &tmp);
          if (cmpee(&a, &tmp) > 0)
             memcpy(&res, &rs, sizeof(res));
       }
       else
          memcpy(&res, &rp, sizeof(res));
    } 

    setrndmode(rnd);

    E_SPOPP("t_sqrt")
    return(res);
}





