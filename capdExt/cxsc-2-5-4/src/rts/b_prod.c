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

/* CVS $Id: b_prod.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_prod.c                              */
/*                                                              */
/*      Entry           : void b_prod(a,b,lang)                 */
/*                        a_btyp *a,*b,*lang;                   */
/*                                                              */
/*      Arguments       : a = first operand for product         */
/*                        b = second operand for product        */
/*                        lang = product                        */
/*                                                              */
/*      Description     : Multiply 2-digit a_btyp numbers       */
/*                                                              */
/*      Note            : length of "lang" at least 4           */
/*                        length of "a","b" is 2                */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local void b_prod(a_btyp *a,a_btyp *b,a_btyp *lang)
#else
local void b_prod(a,b,lang)

a_btyp *a;
a_btyp *b;
a_btyp *lang;
#endif
        {
        a_btyp a0,a1,a2,a3,b0,b1,b2,b3;
        a_btyp n2;
        register a_btyp n0,n1,n3;

        /* split mantissas */
        a0 = GETHIGH(a[0]);
        a1 = GETLOW(a[0]);
        a2 = GETHIGH(a[1]);
        a3 = GETLOW(a[1]);
        b0 = GETHIGH(b[0]);
        b1 = GETLOW(b[0]);
        b2 = GETHIGH(b[1]);
        b3 = GETLOW(b[1]);

        n1 = a3*b0+a0*b3;
        if (~(lang[2] = a2*b3)<(n0 = a3*b2)) n1++;
        lang[3] = (lang[2] += n0)<<16;
        n3 = a1*b0+a0*b1;
        if (~(n0 = a1*b2)<n1) n3++;
        if (~(n1 += n0)<(n2 = a2*b1)) n3++;
        lang[2] = (lang[2]>>16) | ((n1 += n2)<<16);
        lang[1] = (n1>>16) | (n3<<16);
        lang[0] = (n3>>16)+a0*b0;
        n2 = a2*b2;
        if (~(n0 = a3*b3)<lang[3]) n2++;
        lang[3] += n0;
        n1 = a0*b2+a2*b0;
        if (~(n0 = a1*b3)<n2) n1++;
        if (~(n2 += n0)<(n3 = a3*b1)) n1++;
        if (~(n2 += n3)<lang[2]) n1++;
        lang[2] += n2;
        if (~(n0 = a1*b1)<n1) lang[0]++;
        if (~(n1 += n0)<lang[1]) lang[0]++;
        lang[1] += n1;

        return;
        }





