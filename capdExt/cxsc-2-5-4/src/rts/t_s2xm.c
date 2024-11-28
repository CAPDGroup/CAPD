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

/* CVS $Id: t_s2xm.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/*****************************************************************/
/*    Filename: t_s2xm.c                                         */
/*                                                               */
/*    Date    : 1992-07-23                                       */
/*    1992-07-23 : add SUN4_CPP_c =b=                            */
/*****************************************************************/
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
/* #include "eexp.h"       */
/*-------------------------------------------------------------*/
/* Header Emulation Exp                                        */
/*-------------------------------------------------------------*/

#define N 11         /* (0..N-1) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
           EXTREAL(0x3F, 0xFE,              \
                   0xB1, 0x72, 0x17, 0xF7,  \
                   0xD1, 0xCF, 0x79, 0xAC), /* +6.931471805599453E-001 */ \
           EXTREAL(0x3F, 0xFC,              \
                   0xF5, 0xFD, 0xEF, 0xFC,  \
                   0x16, 0x2C, 0x75, 0x43), /* +2.402265069591007E-001 */ \
           EXTREAL(0x3F, 0xFA,              \
                   0xE3, 0x58, 0x46, 0xB8,  \
                   0x25, 0x05, 0xFC, 0x5A), /* +5.550410866482158E-002 */ \
           EXTREAL(0x3F, 0xF8,              \
                   0x9D, 0x95, 0x5B, 0x7D,  \
                   0xD2, 0x73, 0xB9, 0x4E), /* +9.618129107628477E-003 */ \
           EXTREAL(0x3F, 0xF5,              \
                   0xAE, 0xC3, 0xFF, 0x3C,  \
                   0x53, 0x39, 0x88, 0x84), /* +1.333355814642844E-003 */ \
           EXTREAL(0x3F, 0xF2,              \
                   0xA1, 0x84, 0x89, 0x7C,  \
                   0x36, 0x3C, 0x3B, 0x7A), /* +1.540353039338161E-004 */ \
           EXTREAL(0x3F, 0xEE,              \
                   0xFF, 0xE5, 0xFE, 0x2C,  \
                   0x45, 0x86, 0x34, 0x36), /* +1.525273380405984E-005 */ \
           EXTREAL(0x3F, 0xEB,              \
                   0xB1, 0x60, 0x11, 0x1D,  \
                   0x2E, 0x41, 0x1F, 0xEC), /* +1.321548679014431E-006 */ \
           EXTREAL(0x3F, 0xE7,              \
                   0xDA, 0x92, 0x9E, 0x9C,  \
                   0xAF, 0x3E, 0x1E, 0xD2), /* +1.017808600923970E-007 */ \
           EXTREAL(0x3F, 0xE3,              \
                   0xF2, 0x67, 0xA8, 0xAC,  \
                   0x5C, 0x76, 0x4F, 0xB8), /* +7.054911620801123E-009 */ \
           EXTREAL(0x3F, 0xDF,              \
                   0xF4, 0x65, 0x63, 0x9A,  \
                   0x8D, 0xD9, 0x26, 0x08)  /* +4.445538271870812E-010 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

/*--------------------------------------------------------------*
 | Emulation Exp 2**x-1                                         |
 | Fehler berechnet fuer 0<=arg<0.1                             |
 *--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int _s_2xm1(const ExtReal *arg, ExtReal *res)
#else
int _s_2xm1(arg, res)
const ExtReal  *arg;    /* 0<=arg<0.1                           */
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int _s_2xm1(ExtReal *arg, ExtReal *res)
#else
int _s_2xm1(arg, res)
ExtReal  *arg;          /* 0<=arg<0.1                           */
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
/* Koeffizienten Approximationspoly */
#if SUN4_CPP_C
   static ExtReal a[N];
   static char ac[N][sizeof(ExtReal)] = ApproxKoeff;
#else
   static const ExtReal a[N]=ApproxKoeff;
#endif
#else
   static ExtReal a[N];
   static char ac[N][sizeof(ExtReal)] = ApproxKoeff;
#endif
   ExtReal   rn;        /* n-th Result Horner                   */
   ExtReal   e;         /* Zwischenergebnis                     */
   register int i;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<N; i++)
         memcpy(&(a[i]), ac[i], sizeof(a[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<N; i++)
         memcpy(&(a[i]), ac[i], sizeof(a[0]));
   }
#endif
#endif

   /* --- Horner --- */
   copyee(&a[N-1], &rn);
   for(i=N-2; i>=0; i--) {
      mulee(arg, &rn, &e);
      addee(&e, &a[i], &rn);
   }
   mulee(arg, &rn, res);

   /* --- kein Fehler moeglich --- */
   return NoErr;
} /* _s_2xm1() */





