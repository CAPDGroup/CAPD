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

/* CVS $Id: t_satn.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*      Filename: t_satn.c                                      */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "eatan.h" */
/* ------------------------------------------------------------- */
/* Header Atan Emulation Konstanten                              */
/* ------------------------------------------------------------- */
#ifdef LINT_ARGS
#ifdef ANSI_C
static int teilintervall(const ExtReal *arg);
static int red(const ExtReal *arg, int intervall, ExtReal *res);
static int approx(const ExtReal *arg, ExtReal *res);
static int anpassung(const ExtReal *arg,int sgn,int intervall,
                     ExtReal *res);
#else  /* NOT ANSI_C */
static int teilintervall(ExtReal *arg);
static int red(ExtReal *arg, int intervall, ExtReal *res);
static int approx(ExtReal *arg, ExtReal *res);
static int anpassung(ExtReal *arg, int sgn, int intervall,
                     ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int teilintervall(), red(), approx(), anpassung();
#endif /* LINT_ARGS */

#define N 9                       /* (0..N) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
EXTREAL(0x3F, 0xFF,              \
        0x80, 0x00, 0x00, 0x00,  \
        0x00, 0x00, 0x00, 0x00), /* +1.000000000000000E+000 */ \
EXTREAL(0xBF, 0xFD,              \
        0xAA, 0xAA, 0xAA, 0xAA,  \
        0xAA, 0xAA, 0xAA, 0xAB), /* -3.333333333333333E-001 */ \
EXTREAL(0x3F, 0xFC,              \
        0xCC, 0xCC, 0xCC, 0xCC,  \
        0xCC, 0xCC, 0xCC, 0xCD), /* +2.000000000000000E-001 */ \
EXTREAL(0xBF, 0xFC,              \
        0x92, 0x49, 0x24, 0x92,  \
        0x49, 0x24, 0x92, 0x49), /* -1.428571428571428E-001 */ \
EXTREAL(0x3F, 0xFB,              \
        0xE3, 0x8E, 0x38, 0xE3,  \
        0x8E, 0x38, 0xE3, 0x8E), /* +1.111111111111111E-001 */ \
EXTREAL(0xBF, 0xFB,              \
        0xBA, 0x2E, 0x8B, 0xA2,  \
        0xE8, 0xBA, 0x2E, 0x8C), /* -9.090909090909091E-002 */ \
EXTREAL(0x3F, 0xFB,              \
        0x9D, 0x89, 0xD8, 0x9D,  \
        0x89, 0xD8, 0x9D, 0x8A), /* +7.692307692307693E-002 */ \
EXTREAL(0xBF, 0xFB,              \
        0x88, 0x88, 0x88, 0x88,  \
        0x88, 0x88, 0x88, 0x89), /* -6.666666666666667E-002 */ \
EXTREAL(0x3F, 0xFA,              \
        0xF0, 0xF0, 0xF0, 0xF0,  \
        0xF0, 0xF0, 0xF0, 0xF1), /* +5.882352941176471E-002 */ \
EXTREAL(0xBF, 0xFA,              \
        0xD7, 0x94, 0x35, 0xE5,  \
        0x0D, 0x79, 0x43, 0x5E)  /* -5.263157894736842E-002 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

#define Teilpunkte {             /* Teilpunkte fuer Argumentreduktion */ \
EXTREAL(0x3f, 0xfb,         \
        0xff,  0,  0,  0,   \
         0,  0,  0,  0),    /*  +1.245117187500000E-001 */ \
EXTREAL(0x3f, 0xfd,         \
        0xc7,  0,  0,  0,   \
         0,  0,  0,  0),    /*  +3.886718750000000E-001 */  \
EXTREAL(0x3f, 0xfe,         \
        0xb6,  0,  0,  0,   \
         0,  0,  0,  0),    /*  +7.109375000000000E-001 */  \
EXTREAL(0x3f, 0xff,         \
        0x96,  0,  0,  0,   \
         0,  0,  0,  0),    /*  +1.171875000000000E+000 */  \
EXTREAL(0x40, 00,           \
        0x81,  0,  0,  0,   \
         0,  0,  0,  0),    /*  +2.015625000000000E+000 */  \
EXTREAL(0x40, 01,           \
        0x94,  0,  0,  0,   \
         0,  0,  0,  0),    /*  +4.625000000000000E+000 */  \
EXTREAL(0x40, 02,           \
        0x81,  0,  0,  0,   \
         0,  0,  0,  0)     /*  +8.062500000000000E+000 */  \
} /* --- Ende Teilpunkte fuer Argumentreduktion --- */

#define RedConst {              /* Konstanten fuer Argumentreduktion */ \
EXTREAL(0x3F, 0xFD,         \
        0x81, 00, 00, 00,   \
        00, 00, 00, 00),    /* +2.519531250000000E-001 */ \
EXTREAL(0x3F, 0xFE,         \
        0x8A, 00, 00, 00,   \
        00, 00, 00, 00),    /* +5.390625000000000E-001 */ \
EXTREAL(0x3F, 0xFE,         \
        0xEA, 00, 00, 00,   \
        00, 00, 00, 00),    /* +9.140625000000000E-001 */ \
EXTREAL(0x3F, 0xFF,         \
        0xC2, 00, 00, 00,   \
        00, 00, 00, 00),    /* +1.515625000000000E+000 */ \
EXTREAL(0x40, 00,           \
        0xB7, 00, 00, 00,   \
        00, 00, 00, 00),    /* +2.859375000000000E+000 */ \
EXTREAL(0x40, 02,           \
        0xB4, 00, 00, 00,   \
        00, 00, 00, 00)     /* +1.125000000000000E+001 */ \
} /* --- Ende Konstanten fuer Argumentreduktion --- */

#define AnpassKonst {           /* Anpassungs-Konstanten = tan(c[]) */ \
EXTREAL(0x3F, 0xFC,             \
        0xFC, 0xBD, 0x58, 0xDD, \
        0x41, 0x2F, 0xA8, 0x3E),/* +2.468160519641029E-001 */ \
EXTREAL(0x3F, 0xFD,             \
        0xFD, 0x22, 0xEE, 0x98, \
        0x14, 0x92, 0xC4, 0x6B),/* +4.944071350712753E-001 */ \
EXTREAL(0x3F, 0xFE,             \
        0xBD, 0x93, 0x65, 0x69, \
        0x72, 0x87, 0xEC, 0x63),/* +7.405303366126926E-001 */ \
EXTREAL(0x3F, 0xFE,             \
        0xFC, 0xD1, 0x30, 0x25, \
        0x93, 0x9E, 0xBA, 0xB4),/* +9.875669566860051E-001 */ \
EXTREAL(0x3F, 0xFF,             \
        0x9D, 0xFF, 0xAB, 0x91, \
        0x47, 0xDA, 0xE8, 0x15),/* +1.234364934861979E+000 */ \
EXTREAL(0x3F, 0xFF,             \
        0xBD, 0xB6, 0xC7, 0x31, \
        0x85, 0x6A, 0xF1, 0x8A),/* +1.482140444927459E+000 */ \
EXTREAL(0x3F, 0xFF,             /* NEAR(PiHalf)            */ \
        0xC9, 0x0F, 0xDA, 0xA2, \
        0x21, 0x68, 0xC2, 0x35) /* +1.570796326794897E+000 */ \
} /* --- Ende Anpassungs-Konstanten = tan(c[]) --- */

/*----------------------------------------------------------------*
 | Emulation Arcus Tangens                                        |
 *----------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_atan(const ExtReal *arg, ExtReal *res)
#else
int _s_atan(arg, res)
const ExtReal  *arg;    /* reduziertes Argument (t)  -pi/4<t<pi/4 */
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_atan(ExtReal *arg, ExtReal *res)
#else
int _s_atan(arg, res)
ExtReal  *arg;          /* reduziertes Argument (t)  -pi/4<t<pi/4 */
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal   absa;              /* Absolutes Argument */
   ExtReal   red_arg;           /* Ergebnis der Reduktion */
   ExtReal   app;               /* Ergebnis der Approximation */
   int       sgn;               /* Vorzeichen Argument */
   int       i;                 /* Nummer Teilintervall */

   /* --- absa = |arg| --- */
   sgn = SGNE(arg);
   absee(arg, &absa);

   i = teilintervall(&absa);
   red(&absa, i, &red_arg);
   approx(&red_arg, &app);
   anpassung(&app, sgn, i, res);

   return NoErr;
} /* _s_atan() */


/*------------------------------------------------------------*
 | Teilintervall suchen                                       |
 *------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int teilintervall(const ExtReal *arg)
#else
static int teilintervall(arg)
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int teilintervall(ExtReal *arg)
#else
static int teilintervall(arg)
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
/* Teilpunkte Argumentreduktion */
#if SUN4_CPP_C
   static ExtReal b[7];
   static char bc[7][sizeof(ExtReal)] = Teilpunkte;
#else
   static const ExtReal b[7] = Teilpunkte;
#endif
#else
   static ExtReal b[7];
   static char bc[7][sizeof(ExtReal)] = Teilpunkte;
#endif
   register int i, j;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<7; i++)
         memcpy(&(b[i]), bc[i], sizeof(b[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<7; i++)
         memcpy(&(b[i]), bc[i], sizeof(b[0]));
   }
#endif
#endif

   j=3;
   for(i=2; i>0; i--)
      j+=i*(-1==cmpee(arg, &b[j])?-1:1);
   j+=(-1==cmpee(arg, &b[j])?0:1);

   return j;
} /* teilintervall() */

/*---------------------------------------------------------*
 | Argument-Reduktion                                      |
 *---------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int red(const ExtReal *arg, int intervall, ExtReal *res)
#else
static int red(arg, intervall, res)
const ExtReal   *arg;
      int       intervall;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int red(ExtReal *arg, int intervall, ExtReal *res)
#else
static int red(arg, intervall, res)
ExtReal   *arg;
int       intervall;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
   /* Konstanten fuer Argumentreduktion */
#if SUN4_CPP_C
   static ExtReal c[6];
   static char cc[6][sizeof(ExtReal)] = RedConst;
#else
   static const ExtReal c[6] = RedConst;
#endif
#else
   static ExtReal c[6];
   static char cc[6][sizeof(ExtReal)] = RedConst;
#endif
   ExtReal e1, e2, e3;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<6; i++)
         memcpy(&(c[i]), cc[i], sizeof(c[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<6; i++)
         memcpy(&(c[i]), cc[i], sizeof(c[0]));
   }
#endif
#endif

   if(0==intervall) {
      copyee(arg, res);
      return NoErr;
   }
   if(7==intervall) {
      divee(&MinusOne, arg, res);
      return NoErr;
   }

   subee(arg, &c[intervall-1], &e1);
   mulee(arg, &c[intervall-1], &e2);
   addee(&One, &e2, &e3);
   divee(&e1, &e3, res);

   return NoErr;
} /* red() */

/*-----------------------------------------------------------*
 | Approximation                                             |
 *-----------------------------------------------------------*/

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
   ExtReal   u;         /* x*x */
   ExtReal   rn;        /* n-th result */
   ExtReal   e;         /* Zwischenerg., ExtReal */
   register int i;

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

   mulee(arg, arg, &u);

   copyee(&a[N], &rn);
   for(i=N-1; i>=0; i--) {
      mulee(&u, &rn, &e);
      addee(&e, &a[i], &rn);
   }
   mulee(arg, &rn, res);

   return NoErr;
} /* approx() */

/*---------------------------------------------------------*
 | Anpassung                                               |
 *---------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int anpassung(const ExtReal *arg,int sgn,int intervall,ExtReal *res)
#else
static int anpassung(arg, sgn, intervall, res)
const ExtReal   *arg;
      int       sgn;
      int       intervall;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int anpassung(ExtReal *arg, int sgn, int intervall, ExtReal *res)
#else
static int anpassung(arg, sgn, intervall, res)
ExtReal   *arg;
int       sgn;
int       intervall;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
                      /* Anpassungs-Konstanten = tan(c[]) */
#if SUN4_CPP_C
   static ExtReal alpha[7];
   static char alphac[7][sizeof(ExtReal)] = AnpassKonst;

   {
      int i;
      for (i=0; i<7; i++)
         memcpy(&(alpha[i]), alphac[i], sizeof(alpha[0]));
   }
#else
   static const ExtReal alpha[7] = AnpassKonst;
#endif
#else
   static ExtReal alpha[7];
   static char alphac[7][sizeof(ExtReal)] = AnpassKonst;

   {
      int i;
      for (i=0; i<7; i++)
         memcpy(&(alpha[i]), alphac[i], sizeof(alpha[0]));
   }
#endif

   if(intervall==0)
      copyee(arg, res);
   else
      addee(arg, &alpha[intervall-1], res);

   if(sgn == NEG)
      chsee(res, res);

   return NoErr;
} /* anpassung() */





