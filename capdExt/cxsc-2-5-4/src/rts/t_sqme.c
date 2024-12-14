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

/* CVS $Id: t_sqme.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_sqme.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "sqrtm1.h"     */
/* ----------------------------------------------------------- */
/* Header sqrt(1+x)-1 Emulation                                */
/* ----------------------------------------------------------- */
#ifdef LINT_ARGS
#ifdef ANSI_C
static int approx(const ExtReal *arg, ExtReal *res);
#else  /* NOT ANSI_C */
static int approx(ExtReal *arg, ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int approx();
#endif /* LINT_ARGS */

#define N 12        /* (0..N) N = Grad Approx Polynom */
#define Sqrtm1ApproxKoeff {                  \
           EXTREAL(0x3F, 0xFE,               \
                   0x80, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +5.000000000000000E-001 */ \
           EXTREAL(0xBF, 0xFC,               \
                   0x80, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* -1.250000000000000E-001 */ \
           EXTREAL(0x3F, 0xFB,               \
                   0x80, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +6.250000000000000E-002 */ \
           EXTREAL(0xBF, 0xFA,               \
                   0xA0, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* -3.906250000000000E-002 */ \
           EXTREAL(0x3F, 0xF9,               \
                   0xE0, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +2.734375000000000E-002 */ \
           EXTREAL(0xBF, 0xF9,               \
                   0xA8, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* -2.050781250000000E-002 */ \
           EXTREAL(0x3F, 0xF9,               \
                   0x84, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +1.611328125000000E-002 */ \
           EXTREAL(0xBF, 0xF8,               \
                   0xD6, 0x80, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* -1.309204101562500E-002 */ \
           EXTREAL(0x3F, 0xF8,               \
                   0xB2, 0xC0, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +1.091003417968750E-002 */ \
           EXTREAL(0xBF, 0xF8,               \
                   0x97, 0xF0, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* -9.273529052734375E-003 */ \
           EXTREAL(0x3F, 0xF8,               \
                   0x83, 0x38, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +8.008956909179687E-003 */ \
           EXTREAL(0xBF, 0xF7,               \
                   0xE5, 0xA2, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* -7.007837295532227E-003 */ \
           EXTREAL(0x3F, 0xF7,               \
                   0xCB, 0x23, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00)   /* +6.199240684509277E-003 */ \
} /* --- Ende Sqrtm1Konstanten --- */

/*------------------------------------------------------------*
 | Punkt Sqrtm1                                               |
 *------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int sqrtm1ee (const ExtReal *arg, ExtReal *res)
#else
int sqrtm1ee (arg, res)
const ExtReal   *arg;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int sqrtm1ee (ExtReal *arg, ExtReal *res)
#else
int sqrtm1ee (arg, res)
ExtReal   *arg;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         rnd;        /* RundungsMode                      */
   int         check;      /* Rueckgabe von Makro ArgCheck      */

   /* --- pruefe Argument --- */
   ArgCheck1(Sqrtm1, arg, res);

   /* --- RundungsMode sichern, NEAR setzen --- */
   rnd  = getrndmode();
   setrndmode(NEAR);

   /* --- Sqrtm1 --- */
   approx(arg, res);

   /* --- RundungsMode zuruecksetzen --- */
   setrndmode(rnd);

   return NoErr;
} /* sqrtm1ee() */

/*--------------------------------------------------------------*
 | Approximation                                                |
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
   static char ac[N+1][sizeof(ExtReal)] = Sqrtm1ApproxKoeff;
#else
   static const ExtReal a[N+1] = Sqrtm1ApproxKoeff;
#endif
#else
   static ExtReal a[N+1];
   static char ac[N+1][sizeof(ExtReal)] = Sqrtm1ApproxKoeff;
#endif
   ExtReal     rn;      /* n-th result                          */
   ExtReal     e;       /* Zwischenergebnis                     */
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





