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

/* CVS $Id: t_sqt_.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*      Filename: t_sqt_.c                                      */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/* #include "esqrt.h" */
/*--------------------------------------------------------------*/
/* Header Emulation QuadratWurzel                               */
/*--------------------------------------------------------------*/

#define N 4                     /* (0..N) N = Anzahl HeronIterationen */
#define AApprox                             \
EXTREAL(0x3F, 0xFD,              \
        0xD5, 0xA9, 0x57, 0x79,  \
        0xA4, 0x72, 0x91, 0x48)  /* +4.173075996388650E-001 */
#define BApprox                             \
EXTREAL(0x3F, 0xFE,              \
        0x97, 0x14, 0xDC, 0x79,  \
        0x7E, 0x7B, 0xC9, 0xAA)  /* +5.901620670906446E-001 */

#ifdef LINT_ARGS
#ifdef ANSI_C
static int approx(const ExtReal *arg, int ex, ExtReal *res);
static int heron(const ExtReal *arg, const ExtReal *lin_approx,
                 ExtReal *res);
#else  /* NOT ANSI_C */
static int approx(ExtReal *arg, int ex, ExtReal *res);
static int heron(ExtReal *arg, ExtReal *lin_approx, ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int approx(), heron();
#endif /* LINT_ARGS */

/*------------------------------------------------------------------------*
 | Emulation Wurzel                                                       |
 *------------------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_sqrt(const ExtReal *arg, ExtReal *res)
#else
int _s_sqrt(arg, res)
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_sqrt(ExtReal *arg, ExtReal *res)
#else
int _s_sqrt(arg, res)
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal  lin_approx;       /* Ergebnis der linearen Naeherung   */
   ExtReal  earg;             /* Exponent des Arguments            */
   ExtReal  marg;             /* Mantisse des Arguments            */
   ExtReal  linaprx;          /* lineare approximation             */
   int      ex;

   if(0==cmpee(arg, &Zero)) {              /* falls arg==0 gleich raus */
      copyee(&Zero, res);
      return NoErr;
   }

   /* --- spalte arg auf in Mantisse und Exponenten ---*/
   xtracte(arg, &marg, &earg);

   /* --- ex = Exponent modulo 2 ---*/
   ex = mod2e(&earg);

   /* --- lineare Ausgangsnaeherung --- */
   approx(&marg, ex, &lin_approx);

   /* --- Exponenten angleichen ---*/
   extreal_to_int(&earg, &ex);
   if(ex<0) ex--;
   ex /= 2;
   /* lin_approx.s.exp += ex;*/
   scaliee(&lin_approx, ex, &linaprx);

   /* --- Newton Iterationen --- */
   heron(arg, &linaprx, res);

   /* --- kein Fehler moeglich --- */
   return NoErr;
} /* _s_sqrt() */

/*----------------------------------------------------------*
 | Heron                                                    |
 *----------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int heron(const ExtReal *arg, const ExtReal *lin_approx, ExtReal *res)
#else
static int heron(arg, lin_approx, res)
const ExtReal   *arg;
const ExtReal   *lin_approx;
      ExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int heron(ExtReal *arg, ExtReal *lin_approx, ExtReal *res)
#else
static int heron(arg, lin_approx, res)
ExtReal   *arg;
ExtReal   *lin_approx;
ExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal  e1, e2;     /* Zwischenergebnisse                  */
   register int i;

   copyee(lin_approx, res);
   for (i=0; i<N; i++) {
      divee(arg, res, &e1);
      addee(&e1, res, &e2);
      mulee(&e2, &Half, res);
   }

   return NoErr;
} /* heron() */

/*-----------------------------------------------------------*
 | Lineare Approximation                                     |
 *-----------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
static int approx(const ExtReal *arg, int ex, ExtReal *res)
#else
static int approx(arg, ex, res)
const ExtReal  *arg;
      int      ex;  /* Exponent mod 2                       */
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int approx(ExtReal *arg, int ex, ExtReal *res)
#else
static int approx(arg, ex, res)
ExtReal  *arg;
int      ex;        /* Exponent mod 2                   */
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal a, b;
   static char ac[sizeof(ExtReal)] = AApprox;
   static char bc[sizeof(ExtReal)] = BApprox;
#else
   static const ExtReal a = AApprox;
   static const ExtReal b = BApprox;
#endif
#else
   static ExtReal a, b;
   static char ac[sizeof(ExtReal)] = AApprox;
   static char bc[sizeof(ExtReal)] = BApprox;
#endif
   ExtReal  e;          /* Zwischenergebnisse           */

#ifndef ANSI_C
   memcpy(&a, ac, sizeof(a));
   memcpy(&b, bc, sizeof(b));
#else
#if SUN4_CPP_C
   memcpy(&a, ac, sizeof(a));
   memcpy(&b, bc, sizeof(b));
#endif
#endif

   mulee(&a, arg, &e);
   addee(&b, &e, res);

   if(ex) mulee(&SqrtTwo, res, res);

   return NoErr;
} /* approx() */





