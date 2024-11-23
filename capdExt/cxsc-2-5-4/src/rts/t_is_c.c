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

/* CVS $Id: t_is_c.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | sincos                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int isincos(const IDReal *v, const IExtReal *j, int jmod4u, int jmod4l,
            IExtReal *res)
#else
int isincos(v, j, jmod4u, jmod4l, res)
const IDReal   *v;         /* Produkt des Arguments mit 2/pi    */
const IExtReal *j;         /* Ganzzahliger Anteil von v         */
      int      jmod4u;     /* j.u modulo 4                      */
      int      jmod4l;     /* j.l modulo 4                      */
      IExtReal *res;       /* Ergebnis                          */
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int isincos(IDReal *v, IExtReal *j, int jmod4u, int jmod4l, IExtReal *res)
#else
int isincos(v, j, jmod4u, jmod4l, res)
IDReal   *v;         /* Produkt des Arguments mit 2/pi                 */
IExtReal *j;         /* Ganzzahliger Anteil von v                      */
int      jmod4u;     /* j.u modulo 4                                   */
int      jmod4l;     /* j.l modulo 4                                   */
IExtReal *res;       /* Ergebnis                                       */
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     tu;         /* reduziertes Argument upper Bound         */
   ExtReal     tl;         /* reduziertes Argument lower Bound         */
   ExtReal     jd;         /* Differenz j->u - j->l, max(jdiff):=4     */
   int         jdiff;      /* =jd                                      */
   int         ret = 0;    /* Rueckgabe Exakt                          */
   int         retru = 0;  /* Rueckgabe Reduktion                      */
   int         retrl = 0;  /* Rueckgabe Reduktion                      */
   int         retu = 0;   /* Rueckgabe sincos                         */
   int         retl = 0;   /* Rueckgabe sincos                         */

   /* --- jdiff --- */
   subee(&(j->u), &(j->l), &jd);
   if (1 == cmpee(&jd, &Four)) copyee(&Four, &jd);
   extreal_to_int(&jd, &jdiff);

   switch (jdiff) {
   case 0: /* ------------------------------------------ jdiff = 0 --- */
      /* --- ganzzahlige Anteile reduzierte Argumente stimmen ueberein --- */
      /* jmod4u = jmod4l; */
      /* --- reduziere Arg --- */
      retru = red_trg(&(v->u), &(j->u), jmod4u, &tu);
      retrl = red_trg(&(v->l), &(j->l), jmod4l, &tl);

      /* --- sincos(t*pi/2) --- */
      switch (jmod4l) {
      case 0:
      case 3:
         retu = sincos(&tu, &(res->u));
         retl = sincos(&tl, &(res->l));
         break;
      case 1:
      case 2:
         retu = sincos(&tl, &(res->u));
         retl = sincos(&tu, &(res->l));
         break;
      }
      break; /* jdiff == 0 */

   case 1: /* ------------------------------------------ jdiff = 1 --- */
      /* --- reduziere Arg --- */
      /*if (jmod4l == 3) jmod4u = 0; */
      /*else jmod4u = jmod4l+1;      */
      retru = red_trg(&(v->u), &(j->u), jmod4u, &tu);
      retrl = red_trg(&(v->l), &(j->l), jmod4l, &tl);

      switch (jmod4l) {
      case 0:
         copyee(&One, &(res->u));
         ret = UB_Exact;
         retl = (-1==cmpee(&tu, &tl) ?
            sincos(&tu, &(res->l)) :
            sincos(&tl, &(res->l)));
         break;
      case 1: /* mon fallend */
         retu = sincos(&tl, &(res->u));
         retl = sincos(&tu, &(res->l));
         break;
      case 2:
         retu = (1==cmpee(&tu, &tl)    ?
            sincos(&tu, &(res->u)) :
            sincos(&tl, &(res->u)));
         copyee(&MinusOne, &(res->l));
         ret = LB_Exact;
         break;
      case 3: /* mon steigend */
         retu = sincos(&tu, &(res->u));
         retl = sincos(&tl, &(res->l));
         break;
      }
      break; /* jdiff == 1 */

   case 2: /* ------------------------------------------ jdiff = 2 --- */
      switch (jmod4l) {
      case 0:
         copyee(&One, &(res->u));
         ret = UB_Exact;
         /* jmod4u = 2; (jmod4l + jdiff) */
         retru = red_trg(&(v->u), &(j->u), jmod4u, &tu);
         retl = sincos(&tu, &(res->l));
         break;
      case 1:
         retrl = red_trg(&(v->l), &(j->l), jmod4l, &tl);
         retu = sincos(&tl, &(res->u));
         copyee(&MinusOne, &(res->l));
         ret = LB_Exact;
         break;
      case 2:
         /* jmod4u = 0; (jmod4l + jdiff) modulo 4*/
         retru = red_trg(&(v->u), &(j->u), jmod4u, &tu);
         retu = sincos(&tu, &(res->u));
         copyee(&MinusOne, &(res->l));
         ret = LB_Exact;
         break;
      case 3:
         copyee(&One, &(res->u));
         ret = UB_Exact;
         retrl = red_trg(&(v->l), &(j->l), jmod4l, &tl);
         retl = sincos(&tl, &(res->l));
         break;
      }
      break; /* jdiff == 2 */

   case 3: /* ------------------------------------------ jdiff = 3 --- */
      switch (jmod4l) {
      case 0:
      case 2:
         copyee(&One, &(res->u));
         copyee(&MinusOne, &(res->l));
         ret = UB_Exact | LB_Exact;
         break;
      case 1:
         /* jmod4u = 0; (jmod4l + jdiff) modulo 4*/
         retru = red_trg(&(v->u), &(j->u), jmod4u, &tu);
         retrl = red_trg(&(v->l), &(j->l), jmod4l, &tl);
         retu = (1==cmpee(&tl, &tu)    ?
            sincos(&tl, &(res->u)) :
            sincos(&tu, &(res->u)));
         copyee(&MinusOne, &(res->l));
         ret = LB_Exact;
         break;
      case 3:
         copyee(&One, &(res->u));
         ret = UB_Exact;
         /* jmod4u = 2; (jmod4l + jdiff) modulo 4*/
         retru = red_trg(&(v->u), &(j->u), jmod4u, &tu);
         retrl = red_trg(&(v->l), &(j->l), jmod4l, &tl);
         retl = (-1==cmpee(&tl, &tu) ?
            sincos(&tl, &(res->l)) :
            sincos(&tu, &(res->l)));
         break;
      }
      break; /* jdiff == 3 */

   case 4: /* ------------------------------------------ jdiff = 4 --- */
      copyee(&One, &(res->u));
      copyee(&MinusOne, &(res->l));
      ret = UB_Exact | LB_Exact;
      break;
   } /* switch(jdiff) */

   if (!cmpee(&(res->l), &MinusOne)) ret |= LB_Exact;
   if (!cmpee(&(res->u), &One     )) ret |= UB_Exact;

   /* --- Rueckgabe: groesste Except & alle RndCtrl --- */
   if(((retu&Except)!=NoErr)||((retl&Except)!=NoErr))
      return (retu&Except)>(retl&Except) ? retu|ret: retl|ret;
   else
      return retru>retrl? retru|ret: retrl|ret;
} /* isincos() */





