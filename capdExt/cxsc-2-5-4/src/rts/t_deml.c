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

/* CVS $Id: t_deml.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/************************************************************************
 * Fuer Reduktion exp()                                                 * 
 * mul ExtReal mit doppelt langem 1/ln(2)=ld(e)=LdE, res DReal          * 
 * res DReal, nur ca. 90 Bit genau:                                     * 
 *  Laenge mant(ExtReal) + Laenge mant(int(maximales Argument exp2())) =* 
 *           64         +           13                = 77              * 
 * mul au al                                                            * 
 *  b0 r0 r1                                                            * 
 *  b1 r2 r3        bi = ld(e)                                          * 
 *  b2 r4           res = sum(ri)                                       * 
 *                                                                      *
 *    Date : 1992-07-23                                                 *
 *    1992-07-23 : add sun4_cpp_c =b=                                   *
 ************************************************************************
 */
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

#ifdef ANSI_C
#ifdef LINT_ARGS
int mulLdE(const ExtReal *arg, DReal *res)
#else
int mulLdE(arg, res)
const ExtReal   *arg;
      DReal     *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int mulLdE(ExtReal *arg, DReal *res)
#else
int mulLdE(arg, res)
ExtReal   *arg;
DReal     *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#define Log2e {                 /* 1/ln(2)=ld(e):                   */ \
           EXTREAL(0x3F, 0xFF,                                         \
                   0xB8, 0xAA, 0x3B, 0x29,                             \
                   0x00, 0x00, 0x00, 0x00), /* +1.442695040721446E+000 */ \
           EXTREAL(0x3F, 0xDE,                                         \
                   0xB8, 0x2F, 0xE1, 0x76,                             \
                   0x00, 0x00, 0x00, 0x00), /* +1.675171315681720E-010 */ \
           EXTREAL(0x3F, 0xBF,                                         \
                   0xBE, 0x87, 0xFE, 0xD0,                             \
                   0x00, 0x00, 0x00, 0x00)  /* +8.069311544173407E-020 */ \
} /* Log2e */

/* "1.71547 652B8 2FE17 77D0F
      FDA0D 23A7D 11D6A EF551 BAD2B 4B116 4A2CD 9A342 648FB C3887 EE$+0";*/
#ifdef ANSI_C
#if SUN4_CPP_C
   static ExtReal LdE[3];
   static char LdEc[3][sizeof(ExtReal)] = Log2e;
#else
   static ExtReal LdE[3]   = Log2e;
#endif
#else
   static ExtReal LdE[3];
   static char LdEc[3][sizeof(ExtReal)] = Log2e;
#endif
   ExtReal   r[5];
   ExtReal   au;
   ExtReal   al;
   DReal     d;
   register int   i;

#ifndef ANSI_C
   {
      int i;
      for (i=0; i<3; i++)
         memcpy(&(LdE[i]), LdEc[i], sizeof(LdE[0]));
   }
#else
#if SUN4_CPP_C
   {
      int i;
      for (i=0; i<3; i++)
         memcpy(&(LdE[i]), LdEc[i], sizeof(LdE[0]));
   }
#endif
#endif

   msplitee(arg, &au, &al);

   for (i=0; i<2; i++) {
     mulee (&au, &LdE[i], &r[2*i]);
     mulee (&al, &LdE[i], &r[2*i+1]);
   }
   mulee (&au, &LdE[2], &r[4]);

   initd(res);

   for (i=0; i<5; i++) {
      if (0!=cmpee(&r[i], &Zero)) {
         extreal_to_dreal(&r[i], &d);
         adddd(&d, res, res);
      }
   }

   return NoErr;
} /* mulLdE() */





