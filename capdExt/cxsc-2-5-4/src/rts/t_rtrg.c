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

/* CVS $Id: t_rtrg.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | redtrg reduziert Argument trigonometrischer Funktionen       |
 *--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int red_trg(const DReal *v, const ExtReal *j, int jmod4, ExtReal *t)
#else
int red_trg(v, j, jmod4, t)
const DReal    *v;         /* Produkt des Arguments mit 1/Periode */
const ExtReal  *j;         /* Ganzzahliger Anteil von v           */
      int      jmod4;      /* j modulo 4                          */
      ExtReal  *t;         /* reduziertes Argument */
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int red_trg(DReal *v, ExtReal *j, int jmod4, ExtReal *t)
#else
int red_trg(v, j, jmod4, t)
DReal    *v;               /* Produkt des Arguments mit 1/Periode */
ExtReal  *j;               /* Ganzzahliger Anteil von v           */
int      jmod4;            /* j modulo 4                          */
ExtReal  *t;               /* reduziertes Argument */
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal        tmp;     /* temporaer                           */
   DReal          jd;      /* = j                                 */
   DReal          td;      /* = t                                 */
   int            rnd;     /* RundungsMode                        */

   /* --- RundungsMode sichern, CHOP setzen --- */
   rnd  = getrndmode();
   setrndmode(CHOP);

   /* --- Differenz in Abhaengigkeit von jmod4 --- */
   switch (jmod4) {
   case 0:
      extreal_to_dreal(j, &jd);
      subdd(v, &jd, &td);
      break;
   case 1:
      addee(j, &One, &tmp);
      extreal_to_dreal(&tmp, &jd);
      subdd(&jd, v, &td);
      break;
   case 2:
      extreal_to_dreal(j, &jd);
      subdd(&jd, v, &td);
      break;
   case 3:
      addee(j, &One, &tmp);
      extreal_to_dreal(&tmp, &jd);
      subdd(v, &jd, &td);
      break;
   }

   setrndmode(NEAR);
   dreal_to_extreal(&td, t);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   /* --- partial loss of significance --- */
   if(v->e-td.e>LenOfRedConst-ExtMantLenTenbyte*BitsPerDigit){
   /* printf("!!! t.e:%d, v.e:%d\n",td.e,v->e);*/
      return PLOSS;
   }

   /* --- kein Fehler moeglich --- */
   return NoErr;
} /* red_trg() */





