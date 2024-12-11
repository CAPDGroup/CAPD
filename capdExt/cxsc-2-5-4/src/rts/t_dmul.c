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

/* CVS $Id: t_dmul.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

/*--------------------------------------------------------------*
 | mul zweier ExtReal, res zwei ExtReal                         |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dmulee(const ExtReal *arg1, const ExtReal *arg2,
           ExtReal *resh, ExtReal *resl)
#else
int dmulee(arg1, arg2, resh, resl)
const ExtReal   *arg1;
const ExtReal   *arg2;
      ExtReal   *resh;
      ExtReal   *resl;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dmulee(ExtReal *arg1, ExtReal *arg2, ExtReal *resh, ExtReal *resl)
#else
int dmulee(arg1, arg2, resh, resl)
ExtReal   *arg1;
ExtReal   *arg2;
ExtReal   *resh;
ExtReal   *resl;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DReal     dres;

   muled(arg1, arg2, &dres);
   dreal_to_2extreal(&dres, resh, resl);

   return 0;
} /* dmulee() */

/*--------------------------------------------------------------*
 | mul zweier ExtReal, res doppelt lang                         |
 | mul a1u a1l                                                  |
 | a2u  r0  r2                                                  |
 | a2l  r1  r3      res = sum(ri)                               |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int muled(const ExtReal *arg1, const ExtReal *arg2, DReal *res)
#else
int muled(arg1, arg2, res)
const ExtReal *arg1;
const ExtReal *arg2;
      DReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int muled(ExtReal *arg1, ExtReal *arg2, DReal *res)
#else
int muled(arg1, arg2, res)
ExtReal *arg1;
ExtReal *arg2;
DReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal   r[4];
   ExtReal   a1u;
   ExtReal   a1l;
   ExtReal   a2u;
   ExtReal   a2l;
   DReal     d;
   register  int   i;


   msplitee(arg1, &a1u, &a1l);
   msplitee(arg2, &a2u, &a2l);

   mulee (&a1u, &a2u, &r[0]);
   mulee (&a1u, &a2l, &r[1]);
   mulee (&a1l, &a2u, &r[2]);
   mulee (&a1l, &a2l, &r[3]);

   initd(res);

   for (i=0; i<4; i++) {
      if (0!=cmpee(&r[i], &Zero)) {
         extreal_to_dreal(&r[i], &d);
         adddd(&d, res, res);
      }
   }

   return 0;
} /* muled() */

/*--------------------------------------------------------------*
 | mul zweier DReal, res doppelt lang                           |
 | mul a0 a1 a2 a3                                              |
 |  b0 r0 r4 r7 r9                                              |
 |  b1 r1 r5 r8                                                 |
 |  b2 r2 r6                                                    |
 |  b3 r3           res = sum(ri)                               |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int muldd(const DReal *arg1, const DReal *arg2, DReal *res)
#else
int muldd(arg1, arg2, res)
const DReal   *arg1;
const DReal   *arg2;
      DReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int muldd(DReal *arg1, DReal *arg2, DReal *res)
#else
int muldd(arg1, arg2, res)
DReal   *arg1;
DReal   *arg2;
DReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal   r[10];
   ExtReal   a[4];
   ExtReal   b[4];
   ExtReal   ah;
   ExtReal   al;
   ExtReal   bh;
   ExtReal   bl;
   DReal     d;
   int       k;
   register  int   i, j;

   dreal_to_2extreal(arg1, &ah, &al);
   dreal_to_2extreal(arg2, &bh, &bl);

   msplitee(&ah, &a[0], &a[1]);
   msplitee(&al, &a[2], &a[3]);
   msplitee(&bh, &b[0], &b[1]);
   msplitee(&bl, &b[2], &b[3]);

   for (i=0, k=0; i<4; i++)
      for (j=0; j<4-i; j++)
         mulee(&a[i], &b[j], &r[k++]);

   initd(res);

   for (i=0; i<10; i++) {
      if (0!=cmpee(&r[i], &Zero)) {
         extreal_to_dreal(&r[i], &d);
         adddd(&d, res, res);
      }
   }

   return 0;
} /* muled() */





