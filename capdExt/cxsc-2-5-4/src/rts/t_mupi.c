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

/* CVS $Id: t_mupi.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_mupi.c                              */
/*                                                              */
/****************************************************************/

#include <stdio.h>
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

/*--------------------------------------------------------------*
 | mul ExtReal mit doppelt langem 2/pi, res doppelt lang        |
 | mul au al                                                    |
 |  b0 r0 r1                                                    |
 |  b1 r2 r3                                                    |
 |  b2 r4 r5        bi = 2/Pi                                   |
 |  b3 r6 r7        res = sum(ri)                               |
 |  b4 r8                                                       |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int mulendpi(const ExtReal *arg, int two_or_four, DReal *res)
#else
int mulendpi(arg, two_or_four, res)
const ExtReal   *arg;
      int       two_or_four;
      DReal     *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int mulendpi(ExtReal *arg, int two_or_four, DReal *res)
#else
int mulendpi(arg, two_or_four, res)
ExtReal   *arg;
int       two_or_four;
DReal     *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
/* 2/pi     :> a2f9836e.4E441529FC2757D1F534DDC0DB62959944B2C1E2$-1<*/
/* 2/pi     :> a2 f9 83 6e.4E 44 15 29
               FC 27 57 D1 F5 34 DD C0
               DB 62 95 99 44 B2 C1 E2$-1<*/
#define LenConstPi 5
#define LenRes     LenConstPi*2-1

#define CONSTPI  \
     {           \
   EXTREAL(0x3f, 0xfe, 0xa2, 0xf9, 0x83, 0x6e, 0x00, 0x00, 0x00, 0x00),  \
   EXTREAL(0x3f, 0xdd, 0x9c, 0x88, 0x2a, 0x52, 0x00, 0x00, 0x00, 0x00),  \
   EXTREAL(0x3f, 0xbe, 0xfc, 0x27, 0x57, 0xd1, 0x00, 0x00, 0x00, 0x00),  \
   EXTREAL(0x3f, 0x9e, 0xf5, 0x34, 0xdd, 0xc0, 0x00, 0x00, 0x00, 0x00),  \
   EXTREAL(0x3f, 0x7e, 0xdb, 0x62, 0x95, 0x99, 0x00, 0x00, 0x00, 0x00)   \
     }

#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal TwoDivPi[LenConstPi];
   static char TwoDivPic[LenConstPi][sizeof(ExtReal)] = CONSTPI;
#else
   static ExtReal TwoDivPi[LenConstPi] = CONSTPI;
#endif
#else
   static ExtReal TwoDivPi[LenConstPi];
   static char TwoDivPic[LenConstPi][sizeof(ExtReal)] = CONSTPI;
#endif

/*static ExtReal TwoDivPi1[3] = {{ 0x29, 0x15, 0x44, 0x4e, 0x6e,
                                   0x83, 0xf9, 0xa2, 0xfe, 0x3f },
                                 { 0xc0, 0xdd, 0x34, 0xf5, 0xd1,
                                   0x57, 0x27, 0xfc, 0xbe, 0x3f },
                                 { 0xe2, 0xc1, 0xb2, 0x44, 0x99,
                                   0x95, 0x62, 0xdb, 0x7e, 0x3f }}; */
   ExtReal   r[LenRes];
   ExtReal   au;
   ExtReal   al;
   DReal     d;
   register int   i;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<LenConstPi; i++)
         memcpy(&(TwoDivPi[i]), TwoDivPic[i],
                sizeof(TwoDivPi[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<LenConstPi; i++)
         memcpy(&(TwoDivPi[i]), TwoDivPic[i],
                sizeof(TwoDivPi[0]));
   }
#endif
#endif

/* Pruefung entfaellt */
/* if (two_or_four != 2 && two_or_four != 4) return 1;*/

   msplitee(arg, &au, &al);

   for (i=0; i<LenConstPi-1; i++) {
     mulee (&au, &TwoDivPi[i], &r[2*i]);
     mulee (&al, &TwoDivPi[i], &r[2*i+1]);
   }
   mulee (&au, &TwoDivPi[LenConstPi-1], &r[LenRes-1]);

   initd(res);

   for (i=0; i<LenRes; i++) {
      if (0!=cmpee(&r[i], &Zero)) {
         extreal_to_dreal(&r[i], &d);
         adddd(&d, res, res);
      }
   }

   if (two_or_four == 4) res->e++;

   return NoErr;
} /* mulendpi() */





