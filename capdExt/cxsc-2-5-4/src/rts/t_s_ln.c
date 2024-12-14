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

/* CVS $Id: t_s_ln.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_s_ln.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "eln.h" */
/* ------------------------------------------------------------- */
/* Header Ln Emulation                                           */
/* ------------------------------------------------------------- */
#ifdef LINT_ARGS
#ifdef ANSI_C
static int teilintervall(const ExtReal *arg);
static int red(const ExtReal *arg, int intervall, ExtReal *res);
static int approx(const ExtReal *arg, ExtReal *res);
static int anpassung(const ExtReal *arg, ExtReal *exp, int intervall,
                     ExtReal *res);
#else  /* NOT ANSI_C */
static int teilintervall(ExtReal *arg);
static int red(ExtReal *arg, int intervall, ExtReal *res);
static int approx(ExtReal *arg, ExtReal *res);
static int anpassung(ExtReal *arg, ExtReal *exp, int intervall,
                     ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int teilintervall(), red(), approx(), anpassung();
#endif /* LINT_ARGS */

#define N 9         /* (0..N) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
           EXTREAL(0x40, 0x00,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* +2.000000000000000E+000 */ \
           EXTREAL(0x3F, 0xFE,              \
                   0xAA, 0xAA, 0xAA, 0xAA,  \
                   0xAA, 0xAA, 0xAA, 0xAB), /* +6.666666666666666E-001 */ \
           EXTREAL(0x3F, 0xFD,              \
                   0xCC, 0xCC, 0xCC, 0xCC,  \
                   0xCC, 0xCC, 0xCC, 0xCD), /* +4.000000000000000E-001 */ \
           EXTREAL(0x3F, 0xFD,              \
                   0x92, 0x49, 0x24, 0x92,  \
                   0x49, 0x24, 0x92, 0x49), /* +2.857142857142857E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0xE3, 0x8E, 0x38, 0xE3,  \
                   0x8E, 0x38, 0xE3, 0x8E), /* +2.222222222222222E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0xBA, 0x2E, 0x8B, 0xA2,  \
                   0xE8, 0xBA, 0x2E, 0x8C), /* +1.818181818181818E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0x9D, 0x89, 0xD8, 0x9D,  \
                   0x89, 0xD8, 0x9D, 0x8A), /* +1.538461538461539E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0x88, 0x88, 0x88, 0x88,  \
                   0x88, 0x88, 0x88, 0x89), /* +1.333333333333333E-001 */ \
           EXTREAL(0x3F, 0xFB,              \
                   0xF0, 0xF0, 0xF0, 0xF0,  \
                   0xF0, 0xF0, 0xF0, 0xF1), /* +1.176470588235294E-001 */ \
           EXTREAL(0x3F, 0xFB,              \
                   0xD7, 0x94, 0x35, 0xE5,  \
                   0x0D, 0x79, 0x43, 0x5E)  /* +1.052631578947368E-001 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

#define Exactu                  /* obere  Grenze exaktes Intervall  */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xA1, 0xFB, 0x78, 0x12,  \
                   0x1F, 0xB7, 0x81, 0x22)  /* +1.265486725663717E+000 */ \

#define Exactl                  /* untere Grenze exaktes Intervall  */ \
           EXTREAL(0x3F, 0xFE,              \
                   0xCA, 0x4B, 0x30, 0x55,  \
                   0xEE, 0x19, 0x10, 0x1D)  /* +7.902097902097902E-001 */ \
/* --- Ende Intervallgrenzen --- */

#define T  2        /* Anzahl Teilintervalle */
#define RedConst {              /* Konstanten fuer Argumentreduktion */ \
           EXTREAL(0x3F, 0xFF,              \
                   0x80, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* +1.000000000000000E+000 */ \
           EXTREAL(0x3F, 0xFF,              \
                   0xCC, 0xFC, 0x88, 0x16,  \
                   0xEF, 0x80, 0x02, 0x91)  /* +1.601456652831075E+000 */ \
} /* --- Ende Konstanten fuer Argumentreduktion --- */

#define AnpassKonst {            /* Anpassungs-Konstanten = ln(c[]) */ \
           EXTREAL(0x00, 0x00,              \
                   0x00, 0x00, 0x00, 0x00,  \
                   0x00, 0x00, 0x00, 0x00), /* +0.000000000000000E+000 */ \
           EXTREAL(0x3F, 0xFD,              \
                   0xF1, 0x1B, 0x97, 0x24,  \
                   0xDE, 0x72, 0xA2, 0xAB)  /* +4.709136230951334E-001 */ \
} /* --- Ende Anpassungs-Konstanten = ln(c[]) --- */

/*-----------------------------------------------------------------*
 | ln Emulation                                                    |
 *-----------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_ln(const ExtReal *arg, ExtReal *res, int *rr)
#else
int _s_ln(arg, res, rr)
const ExtReal  *arg;
      ExtReal  *res;
      int      *rr;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_ln(ExtReal *arg, ExtReal *res, int *rr)
#else
int _s_ln(arg, res, rr)
ExtReal  *arg;
ExtReal  *res;
int      *rr;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal     earg;          /* Exponent des Argumentes          */
   ExtReal     marg;          /* Mantisse des Argumentes          */
   ExtReal     red_arg;       /* Ergebnis der Reduktion           */
   ExtReal     app;           /* Ergebnis der Approximation       */
   int         i;             /* Nummer Teilintervall             */
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal exactu, exactl;
   static char exactuc[sizeof(ExtReal)] = Exactu;
   static char exactlc[sizeof(ExtReal)] = Exactl;

   memcpy(&exactu, exactuc, sizeof(exactu));
   memcpy(&exactl, exactlc, sizeof(exactl));
#else
   static const ExtReal exactu=Exactu; /* obere  Grenze exaktes Intervall */
   static const ExtReal exactl=Exactl; /* untere Grenze exaktes Intervall */
#endif
#else
   static ExtReal exactu, exactl;
   static char exactuc[sizeof(ExtReal)] = Exactu;
   static char exactlc[sizeof(ExtReal)] = Exactl;

   memcpy(&exactu, exactuc, sizeof(exactu));
   memcpy(&exactl, exactlc, sizeof(exactl));
#endif

   /* --- arg im ApproxIntervall ? --- */
   if(1!=cmpee(arg, &exactu) && -1!=cmpee(arg, &exactl)){

      /* --- Reduktion --- */
      red(arg, 0, &red_arg);

      /* --- Approximation --- */
      approx(&red_arg, res);

      /* --- setze RundungsArt --- */
      *rr=1;
   }

   /* --- Nicht im ApproxIntervall --- */
   else{
      /* --- Aufspaltung in Mantisse und Exponent --- */
      xtracte(arg, &marg, &earg);
      scaliee(&marg, -1, &marg);
      addee(&earg, &One, &earg);

      /* --- Teilintervall ___ */
      i = teilintervall(&marg);

     /* --- Reduktion --- */
     red(&marg, i, &red_arg);

     /* --- Approximation --- */
     approx(&red_arg, &app);

     /* --- Ergebnisanpassung --- */
     anpassung(&app, &earg, i, res);

      /* --- setze RundungsArt --- */
      *rr=2;
   }

   return NoErr;
} /* $ln() */

/*---------------------------------------------------------------*
 | Teilintervall suchen                                          |
 *---------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int teilintervall(const ExtReal *arg)
#else
static int teilintervall(arg)
const ExtReal   *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int teilintervall(ExtReal *arg)
#else
static int teilintervall(arg)
ExtReal   *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal b;
   static char bc[sizeof(ExtReal)] = Exactl;

   memcpy(&b, bc, sizeof(b));
#else
   static const ExtReal b = Exactl;     /* Teilpunkt Argumentreduktion */
#endif
#else
   static ExtReal b;
   static char bc[sizeof(ExtReal)] = Exactl;

   memcpy(&b, bc, sizeof(b));
#endif

   if(-1==cmpee(arg, &b)) return 1;
      else return 0;

} /* teilintervall() */

/*---------------------------------------------------------------------*
 | Argument-Reduktion                                                  |
 *---------------------------------------------------------------------*/

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
   static ExtReal c[T];
   static char cc[T][sizeof(ExtReal)] = RedConst;
#else
   static const ExtReal c[T] = RedConst;
#endif
#else
   static ExtReal c[T];
   static char cc[T][sizeof(ExtReal)] = RedConst;
#endif
   ExtReal t;
   ExtReal z, n;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<T; i++)
         memcpy(&(c[i]), cc[i], sizeof(c[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<T; i++)
         memcpy(&(c[i]), cc[i], sizeof(c[0]));
   }
#endif
#endif

   mulee(arg, &c[intervall], &t);
   subee(&t, &One, &z);
   addee(&t, &One, &n);
   divee(&z, &n, res);

   return NoErr;
} /* red() */

/*-------------------------------------------------------------*
 | Approximation                                               |
 *-------------------------------------------------------------*/

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

/*-------------------------------------------------------------*
 | Anpassung                                                   |
 *-------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int anpassung(const ExtReal *arg, ExtReal *exp, int intervall,
                     ExtReal *res)
#else
static int anpassung(arg, exp, intervall, res)
const ExtReal   *arg;
      ExtReal   *exp;
      int       intervall;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int anpassung(ExtReal *arg, ExtReal *exp, int intervall, ExtReal *res)
#else
static int anpassung(arg, exp, intervall, res)
ExtReal   *arg;
ExtReal   *exp;
int       intervall;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
   /* Anpassungs-Konstanten = tan(c[]) */
#if SUN4_CPP_C
   static ExtReal alpha[T];
   static char alphac[T][sizeof(ExtReal)] = AnpassKonst;
#else
   static const ExtReal alpha[T] = AnpassKonst;
#endif
#else
   static ExtReal alpha[T];
   static char alphac[T][sizeof(ExtReal)] = AnpassKonst;
#endif
   ExtReal e,r;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<T; i++)
         memcpy(&(alpha[i]), alphac[i], sizeof(alpha[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<T; i++)
         memcpy(&(alpha[i]), alphac[i], sizeof(alpha[0]));
   }
#endif
#endif

   subee(arg, &alpha[intervall], &r);
   mulee(&Ln2, exp, &e);
   addee(&r, &e, res);

   return NoErr;
} /* anpassung() */





