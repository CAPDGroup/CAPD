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

/* CVS $Id: t_s_lp.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/******************************************************************/
/*  Filename : t_s_lp.c                                           */
/*                                                                */
/*  Date     : 1992-07-23                                         */
/*  1992-07-23 : add SUN4_CPP_c =b=                               */
/******************************************************************/
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "elnp1.h" */
/* ----------------------------------------------------------- */
/* Header Lnp1 Emulation                                       */
/* ----------------------------------------------------------- */
#define Tp                               /* Teilpunkt               */ \
           EXTREAL(0x3F, 0xFA,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00)  /* +3.125000000000000E-002 */

#define N 12        /* (0..N) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
           EXTREAL(0x3F, 0xFF,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* +1.000000000000000E+000 */ \
           EXTREAL(0xBF, 0xFE,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* -5.000000000000000E-001 */ \
           EXTREAL(0x3F, 0xFD,              \
                   0xAA, 0xAA, 0xAA, 0xAA,  \
                   0xAA, 0xAA, 0xAA, 0xAB), /* +3.333333333333333E-001 */ \
           EXTREAL(0xBF, 0xFD,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* -2.500000000000000E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0xCC, 0xCC, 0xCC, 0xCC,  \
                   0xCC, 0xCC, 0xCC, 0xCD), /* +2.000000000000000E-001 */ \
           EXTREAL(0xBF, 0xFC,              \
                   0xAA, 0xAA, 0xAA, 0xAA,  \
                   0xAA, 0xAA, 0xAA, 0xAB), /* -1.666666666666667E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0x92, 0x49, 0x24, 0x92,  \
                   0x49, 0x24, 0x92, 0x49), /* +1.428571428571428E-001 */ \
           EXTREAL(0xBF, 0xFC,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* -1.250000000000000E-001 */ \
           EXTREAL(0x3F, 0xFB,              \
                   0xE3, 0x8E, 0x38, 0xE3,  \
                   0x8E, 0x38, 0xE3, 0x8E), /* +1.111111111111111E-001 */ \
           EXTREAL(0xBF, 0xFB,              \
                   0xCC, 0xCC, 0xCC, 0xCC,  \
                   0xCC, 0xCC, 0xCC, 0xCD), /* -1.000000000000000E-001 */ \
           EXTREAL(0x3F, 0xFB,              \
                   0xBA, 0x2E, 0x8B, 0xA2,  \
                   0xE8, 0xBA, 0x2E, 0x8C), /* +9.090909090909091E-002 */ \
           EXTREAL(0xBF, 0xFB,              \
                   0xAA, 0xAA, 0xAA, 0xAA,  \
                   0xAA, 0xAA, 0xAA, 0xAB), /* -8.333333333333333E-002 */ \
           EXTREAL(0x3F, 0xFB,              \
                   0x9D, 0x89, 0xD8, 0x9D,  \
                   0x89, 0xD8, 0x9D, 0x8A)  /* +7.692307692307693E-002 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

#ifdef LINT_ARGS
#ifdef ANSI_C
static int approx(const ExtReal *arg, ExtReal *res);
#else  /* NOT ANSI_C */
static int approx(ExtReal *arg, ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int approx();
#endif /* LINT_ARGS */

/*-----------------------------------------------------------------*
 | lnp1 Emulation                                                  |
 | rndmode muss NEAR sein!                                         |
 *-----------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_lnp1(const ExtReal *arg, ExtReal *res)
#else
int _s_lnp1(arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_lnp1(ExtReal *arg, ExtReal *res)
#else
int _s_lnp1(arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#else
   static const ExtReal tp=Tp; /* Teilungspunkt 0.03125  */
#endif
#else
   static ExtReal tp;
   static char tpc[sizeof(ExtReal)] = Tp;
#endif
   ExtReal     a;          /* arg + 1 = Argument fuer ln */
   int         rr;         /* Rueckgabe RundungsArt, hier dummy */

#ifndef ANSI_C
   memcpy(&tp, tpc, sizeof(tp));
#else
#if SUN4_CPP_C
   memcpy(&tp, tpc, sizeof(tp));
#endif
#endif

   /* --- ln(arg+1) falls |arg|>=tp=0.03125 --- */
   if (-1!=cmpabsee(arg, &tp)) {
      addee(arg, &One, &a);
      return _s_ln(&a, res, &rr);
   }

   /* --- lnp1 --- */
   approx(arg, res);

   return NoErr;
} /* $lnp1() */

/*-----------------------------------------------------------------------*
 | Approximation                                                         |
 *-----------------------------------------------------------------------*/

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
   ExtReal     rn;      /* n-th result                              */
   ExtReal     e;       /* Zwischenergebnis                         */
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

   copyee(&a[N], &rn);
   for(i=N-1; i>=0; i--) {
      mulee(arg, &rn, &e);
      addee(&e, &a[i], &rn);
   }
   mulee(arg, &rn, res);

   return NoErr;
} /* approx() */





