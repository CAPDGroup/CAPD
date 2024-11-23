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

/* CVS $Id: t_stan.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/*****************************************************************/
/* Filename : t_stan.c                                           */
/* Date : 1992-07-23                                             */
/*                                                               */
/* 1992-07-23 : add SUN4_CPP_c =b=                               */
/*                                                               */
/*****************************************************************/
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
/* #include "etan.h"       */
/* ------------------------------------------------------------ */
/* Header Tan Emulation Konstanten                              */
/* ------------------------------------------------------------ */
#ifdef LINT_ARGS
#ifdef ANSI_C
static int red(const ExtReal *arg, ExtReal *res);
static int approx(const ExtReal *arg, ExtReal *res);
static int anpassung(const ExtReal *arg, int sgn, int intervall,
                           ExtReal *nres, ExtReal *dres);
#else  /* NOT ANSI_C */
static int red(ExtReal *arg, ExtReal *res);
static int approx(ExtReal *arg, ExtReal *res);
static int anpassung(ExtReal *arg, int sgn, int intervall,
                     ExtReal *nres, ExtReal *dres);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int red(), approx(), anpassung();
#endif /* LINT_ARGS */

#define N 11                      /* (0..N) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
           EXTREAL(0x3F, 0xFF,               \
                   0x80, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +1.000000000000000E+000 */ \
           EXTREAL(0x3F, 0xFD,               \
                   0xAA, 0xAA, 0xAA, 0xAA,   \
                   0xAA, 0xAA, 0xAA, 0xAB),  /* +3.333333333333333E-001 */ \
           EXTREAL(0x3F, 0xFC,               \
                   0x88, 0x88, 0x88, 0x88,   \
                   0x88, 0x88, 0x88, 0x89),  /* +1.333333333333333E-001 */ \
           EXTREAL(0x3F, 0xFA,               \
                   0xDD, 0x0D, 0xD0, 0xDD,   \
                   0x0D, 0xD0, 0xDD, 0x0E),  /* +5.396825396825397E-002 */ \
           EXTREAL(0x3F, 0xF9,               \
                   0xB3, 0x27, 0xA4, 0x41,   \
                   0x60, 0x87, 0xCF, 0x99),  /* +2.186948853615520E-002 */ \
           EXTREAL(0x3F, 0xF8,               \
                   0x91, 0x37, 0x1A, 0xAF,   \
                   0x36, 0x11, 0xE4, 0x7B),  /* +8.863235529902197E-003 */ \
           EXTREAL(0x3F, 0xF6,               \
                   0xEB, 0x69, 0xE8, 0x70,   \
                   0xAB, 0xEE, 0xFD, 0xB0),  /* +3.592128036572481E-003 */ \
           EXTREAL(0x3F, 0xF5,               \
                   0xBE, 0xD1, 0xB2, 0x29,   \
                   0x5B, 0xAF, 0x15, 0xB5),  /* +1.455834387051318E-003 */ \
           EXTREAL(0x3F, 0xF4,               \
                   0x9A, 0xAC, 0x12, 0x40,   \
                   0x1B, 0x3A, 0x22, 0x91),  /* +5.900274409455859E-004 */ \
           EXTREAL(0x3F, 0xF2,               \
                   0xFA, 0xBE, 0xBB, 0x9A,   \
                   0x68, 0xB3, 0x21, 0x0D),  /* +2.391291142435525E-004 */ \
           EXTREAL(0x3F, 0xF1,               \
                   0xCB, 0x3F, 0x0C, 0x57,   \
                   0xE5, 0x7D, 0x64, 0x51),  /* +9.691537956929451E-005 */ \
           EXTREAL(0x3F, 0xF0,               \
                   0xA4, 0xBE, 0xC7, 0x75,   \
                   0x12, 0x92, 0xC9, 0x9F)   /* +3.927832388331683E-005 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

#define AnpassKonst {           /* Anpassungs-Konstanten = tan(c[]) */ \
           EXTREAL(0x3F, 0xFE,               \
                   0xAB, 0x0D, 0xC1, 0x55,   \
                   0xBF, 0xCC, 0x82, 0xF4),  /* +6.681786379192989E-001 */ \
           EXTREAL(0x3F, 0xFD,               \
                   0xD4, 0x13, 0xCC, 0xCF,   \
                   0xE7, 0x79, 0x92, 0x10),  /* +4.142135623730950E-001 */ \
           EXTREAL(0x3F, 0xFC,               \
                   0xCB, 0xAF, 0xAF, 0x02,   \
                   0xA9, 0x8A, 0xC0, 0x3D)   /* +1.989123673796580E-001 */ \
} /* --- Ende Anpassungs-Konstanten = tan(c[]) --- */

/*--------------------------------------------------------------*
 | Emulation Tangens                                            |
 *--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_tan(const ExtReal *arg, ExtReal *nres, ExtReal *dres)
#else
int _s_tan(arg, nres, dres)
const ExtReal  *arg;    /* reduziertes Argument (t)  -pi/4<t<pi/4 */
      ExtReal  *nres;   /* Ergebnis Numerator                     */
      ExtReal  *dres;   /* Ergebnis Denominator                   */
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_tan(ExtReal *arg, ExtReal *nres, ExtReal *dres)
#else
int _s_tan(arg, nres, dres)
ExtReal  *arg;          /* reduziertes Argument (t)  -pi/4<t<pi/4  */
ExtReal  *nres;         /* Ergebnis Numerator                      */
ExtReal  *dres;         /* Ergebnis Denominator                    */
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     absa;    /* absolutes Argument                      */
   ExtReal     red_arg; /* auf [0,pi/16) reduziertes Argument      */
   ExtReal     app;     /* Ergebnis Approx                         */
   int         sgn;     /* Vorzeichen des Argumentes               */
   int         i;       /* Teilintervall 0..3 der zusaetzlichen Reduktion */

   /* --- absa = |arg| --- */
   sgn = SGNE(arg);
   absee(arg, &absa);

   /* --- Teilintervall und Reduktion --- */
   i = red(&absa, &red_arg);

   /* --- Approximation --- */
   approx(&red_arg, &app);

   /* --- Ergebnisanpassung --- */
   anpassung(&app, sgn, i, nres, dres);

   /* --- kein Fehler moeglich --- */
   return NoErr;
} /* _s_tanee() */

/*--------------------------------------------------------------*
 | Teilintervall suchen                                         |
 | Argument-Reduktion                                           |
 *--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int red(const ExtReal *arg, ExtReal *res)
#else
static int red(arg, res)
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int red(ExtReal *arg, ExtReal *res)
#else
static int red(arg, res)
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   if(1==cmpee(arg, &PiDiv8))
      if(1==cmpee(arg, &Pi3Div16)){
         subee(arg, &Pi3Div16, res);
         return 0;
      }
      else {
         subee(arg, &PiDiv8, res);
         return 1;
      }
   else
      if(1==cmpee(arg, &PiDiv16)){
         subee(arg, &PiDiv16, res);
         return 2;
      }
      else {
         copyee(arg,res);
         return 3;
      }

} /* red() */

/*--------------------------------------------------------------*
 | approx                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int approx(const ExtReal *arg, ExtReal *res)
#else
static int approx(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int approx(ExtReal *arg, ExtReal *res)
#else
static int approx(arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
  /* Koeffizienten Approximationspolynom */
#if SUN4_CPP_C
   static ExtReal a[N+1];
   static char ac[N+1][sizeof(ExtReal)] = ApproxKoeff;
#else
   static const ExtReal a[N+1] = ApproxKoeff;
#endif
#else
   static ExtReal a[N+1];
   static char ac[N+1][sizeof(ExtReal)] = ApproxKoeff;
#endif
   ExtReal  u;
   ExtReal  rn, e;
   register  int i;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<N+1; i++)
         memcpy(&(a[i]), ac[i], sizeof(a[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<N+1; i++)
         memcpy(&(a[i]), ac[i], sizeof(a[0]));
   }
#endif
#endif

   /* --- u = arg*arg --- */
   mulee(arg, arg, &u);

   /* --- HornerSchema --- */
   copyee(&a[N], &rn);
   for(i=N-1; i>=0; i--) {
      mulee(&u, &rn, &e);
      addee(&e, &a[i], &rn);
   }
   mulee(arg, &rn, res);

   return NoErr;
} /* approx() */

/*--------------------------------------------------------------*
 | Anpassung                                                    |
 *--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int anpassung(const ExtReal *arg, int sgn, int intervall,
                     ExtReal *nres, ExtReal *dres)
#else
static int anpassung(arg, sgn, intervall, nres, dres)
const ExtReal   *arg;
      int       sgn;
      int       intervall;
      ExtReal   *nres;  /* Ergebnis Numerator                   */
      ExtReal   *dres;  /* Ergebnis Denominator                 */
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int anpassung(ExtReal *arg, int sgn, int intervall,
                     ExtReal *nres, ExtReal *dres)
#else
static int anpassung(arg, sgn, intervall, nres, dres)
ExtReal   *arg;
int       sgn;
int       intervall;
ExtReal   *nres;        /* Ergebnis Numerator                   */
ExtReal   *dres;        /* Ergebnis Denominator                 */
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
   /* Anpassungs-Konstanten = tan(c[]) */
#if SUN4_CPP_C
   static ExtReal alpha[3];
   static char alphac[3][sizeof(ExtReal)] = AnpassKonst;

   {
      int i;
      for (i=0; i<3; i++)
         memcpy(&(alpha[i]), alphac[i], sizeof(alpha[0]));
   }
#else
   static const ExtReal alpha[3] = AnpassKonst;
#endif
#else
   static ExtReal alpha[3];
   static char alphac[3][sizeof(ExtReal)] = AnpassKonst;

   {
      int i;
      for (i=0; i<3; i++)
         memcpy(&(alpha[i]), alphac[i], sizeof(alpha[0]));
   }
#endif

   if(intervall==3){
      copyee(arg, nres);
      copyee(&One, dres);
   }
   else {
      addee(arg, &alpha[intervall], nres);
      mulee(arg, &alpha[intervall], dres);
      subee(&One, dres, dres);
   }

   if(sgn == NEG)
      chsee(nres, nres);

   return NoErr;
} /* anpassung() */





