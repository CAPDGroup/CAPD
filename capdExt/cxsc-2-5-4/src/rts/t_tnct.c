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

/* CVS $Id: t_tnct.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*/
/* tancot                                                       */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int tancot(const ExtReal *arg, int jmod4, ExtReal *res)
#else
int tancot(arg, jmod4, res)
const ExtReal  *arg;    /* reduziertes Argument                 */
      int      jmod4;   /* j modulo 4                           */
      ExtReal  *res;    /* Ergebnis                             */
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int tancot(ExtReal *arg, int jmod4, ExtReal *res)
#else
int tancot(arg, jmod4, res)
ExtReal  *arg;          /* reduziertes Argument                 */
int      jmod4;         /* j modulo 4                           */
ExtReal  *res;          /* Ergebnis                             */
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     t;       /* -pi/4 <= Argument t <= pi/4          */
   ExtReal     rn;      /* Ergebnis Numerator                   */
   ExtReal     rd;      /* Ergebnis Denominator                 */
   int         rnd;     /* RundungsMode                         */
   int         ret;     /* Rueckgabe                            */

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- arg*Pi/4 --- */
   mulee(arg, &PiQuart, &t);
   ret = _s_tan(&t, &rn, &rd);

   /* --- cot = rd/rn --- */
   if(jmod4==1||jmod4==2)
      divee(&rd, &rn, res);
   /* --- tan = rn/rd --- */
   else
      divee(&rn, &rd, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return ret;
} /* tancot() */





