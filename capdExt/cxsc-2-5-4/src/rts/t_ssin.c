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

/* CVS $Id: t_ssin.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*    Filename : t_ssin.c                                       */
/*                                                              */
/*    Date     : 1992-07-23                                     */
/*    1992-07-23 : add SUN4_CPP_c =b=                           */
/*                                                              */
/****************************************************************/
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
/* #include "esin.h"       */
/*---------------------------------------------------------------*/
/* Header Emulation Sinus                                        */
/*---------------------------------------------------------------*/

/* --- Auswahl zwischen einfacher und doppelter Mantisse ----- */
/* ACHTUNG ! bei Aenderung EpsSin in File econst.c verfolgen ! */
#define HPrec   /* falls       definiert: eps(sin) = 2.725e-19 */
                /* falls nicht definiert: eps(sin) = 3.993e-19 */

#define N 11                      /* (0..N) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
           EXTREAL(0x3F, 0xFF,               \
                   0x80, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +1.000000000000000E+000 */ \
           EXTREAL(0xBF, 0xFC,               \
                   0xAA, 0xAA, 0xAA, 0xAA,   \
                   0xAA, 0xAA, 0xAA, 0xAB),  /* -1.666666666666667E-001 */ \
           EXTREAL(0x3F, 0xF8,               \
                   0x88, 0x88, 0x88, 0x88,   \
                   0x88, 0x88, 0x88, 0x89),  /* +8.333333333333333E-003 */ \
           EXTREAL(0xBF, 0xF2,               \
                   0xD0, 0x0D, 0x00, 0xD0,   \
                   0x0D, 0x00, 0xD0, 0x0D),  /* -1.984126984126984E-004 */ \
           EXTREAL(0x3F, 0xEC,               \
                   0xB8, 0xEF, 0x1D, 0x2A,   \
                   0xB6, 0x39, 0x9C, 0x7D),  /* +2.755731922398589E-006 */ \
           EXTREAL(0xBF, 0xE5,               \
                   0xD7, 0x32, 0x2B, 0x3F,   \
                   0xAA, 0x27, 0x1C, 0x7F),  /* -2.505210838544172E-008 */ \
           EXTREAL(0x3F, 0xDE,               \
                   0xB0, 0x92, 0x30, 0x9D,   \
                   0x43, 0x68, 0x4B, 0xE5),  /* +1.605904383682161E-010 */ \
           EXTREAL(0xBF, 0xD6,               \
                   0xD7, 0x3F, 0x9F, 0x39,   \
                   0x9D, 0xC0, 0xF8, 0x8F),  /* -7.647163731819816E-013 */ \
           EXTREAL(0x3F, 0xCE,               \
                   0xCA, 0x96, 0x3B, 0x81,   \
                   0x85, 0x6A, 0x53, 0x59),  /* +2.811457254345521E-015 */ \
           EXTREAL(0xBF, 0xC6,               \
                   0x97, 0xA4, 0xDA, 0x34,   \
                   0x0A, 0x0A, 0xB9, 0x26),  /* -8.220635246624329E-018 */ \
           EXTREAL(0x3F, 0xBD,               \
                   0xB8, 0xDC, 0x77, 0xB6,   \
                   0xE7, 0xAB, 0x8C, 0x5F),  /* +1.957294106339126E-020 */ \
           EXTREAL(0xBF, 0xB4,               \
                   0xBB, 0x0D, 0xA0, 0x98,   \
                   0xB1, 0xC0, 0xCE, 0xCC)   /* -3.868170170630684E-023 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

/* --- Auswahl zwischen einfacher und doppelter Mantisse ----- */
/* --- HPrec jetzt im Header esin.h                      ----- */
/* ACHTUNG ! bei Aenderung EpsSin in File econst.c verfolgen ! */
/* #define HPrec */       /* falls       definiert: eps(sin) = 2.725e-19 */
                /* falls nicht definiert: eps(sin) = 3.993e-19 */

#ifdef LINT_ARGS
#ifdef ANSI_C
#ifdef HPrec
static int approx(const ExtReal *arg, DReal *res);
#else
static int approx(const ExtReal *arg, ExtReal *res);
#endif
#else  /* NOT ANSI_C */
#ifdef HPrec
static int approx(ExtReal *arg, DReal *res);
#else
static int approx(ExtReal *arg, ExtReal *res);
#endif
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int approx();
#endif /* LINT_ARGS */

/*--------------------------------------------------------------*
 | Emulation Sinus                                              |
 *--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_sin(const ExtReal *arg, ExtReal *res)
#else
int _s_sin(arg, res)
const ExtReal  *arg;    /* reduziertes Argument (t)             */
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_sin(ExtReal *arg, ExtReal *res)
#else
int _s_sin(arg, res)
ExtReal  *arg;          /* reduziertes Argument (t)             */
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     u;       /* = arg*arg                            */
#ifdef HPrec
   DReal       app;     /* Ergebnis Approx                      */
   DReal       d, r;    /* Konvertierungsgroessen               */
#else
   ExtReal     app;     /* Ergebnis Approx                      */
#endif

   mulee(arg, arg, &u);
   approx(&u, &app);

#ifdef HPrec
   extreal_to_dreal(arg, &d);
   muldd(&d, &app, &r);
   dreal_to_extreal(&r, res);
#else
   mulee(arg, &app, res);
#endif

   /* kein Fehler moeglich */
   return NoErr;
} /* _s_sinee() */

/*--------------------------------------------------------------*
 | approx                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
#ifdef HPrec
static int approx(const ExtReal *arg, DReal *res)
#else
static int approx(const ExtReal *arg, ExtReal *res)
#endif
#else
#ifdef HPrec
static int approx(arg, res)
const ExtReal *arg;
      DReal   *res;
#else
static int approx(arg, res)
const ExtReal *arg;
      ExtReal *res;
#endif
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
#ifdef HPrec
static int approx(ExtReal *arg, DReal *res)
#else
static int approx(ExtReal *arg, ExtReal *res)
#endif
#else
#ifdef HPrec
static int approx(arg, res)
ExtReal *arg;
DReal   *res;
#else
static int approx(arg, res)
ExtReal *arg;
ExtReal *res;
#endif
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
   ExtReal  rn, e;
#ifdef HPrec
   DReal    r;
#else
   ExtReal  r;
#endif
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

   /* --- HornerSchema --- */
   copyee(&a[N], &rn);
   for(i=N-1; i>0; i--) {
      mulee(arg, &rn, &e);
      addee(&e, &a[i], &rn);
   }
#ifdef HPrec
   muled(arg, &rn, &r);
   adddd(&r, &DOne, res);
#else
   mulee(arg, &rn, &r);
   addee(&r, &One, res);
#endif

   return NoErr;
} /* approx() */





