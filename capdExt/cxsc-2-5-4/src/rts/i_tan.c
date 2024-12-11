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

/* CVS $Id: i_tan.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_tan.c                               */
/*                                                              */
/*      Entries         : a_intv i_tan(a)                       */
/*                        a_intv a;                             */
/*                                                              */
/*      Arguments       : a = interval argument of function     */
/*                                                              */
/*      Description     : Interval tan function                 */
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

extern a_btyp LhI;    /* entier(x/(pi/4)) mod 16 after call tan(x) */

#ifdef LINT_ARGS
local a_intv i_tan( a_intv a )
#else
local a_intv i_tan(a)
a_intv a;
#endif

{
        /* --- Pi.l --- */
        static a_btyp cpi[] = B_TYPINI( 0x400921fbL, 0x54442d18L);
     /* static unsigned long cpi[2] = {0x182d4454L, 0xfb210940L}; */
        a_real            *pil= (a_real*)&cpi[0];

        a_btyp        rc;
        a_intv          res;
        a_real            dummy;
        unsigned        jub, jlb, d;

        E_SPUSH("i_tan")

        /* --- point argument --- */
        if (i_point(a)) {
            /* printf("point argument!\n"); */
            rc = i_invp(Ltan, a.INF, &res.INF, &res.SUP);
            /* printf("entier part Lhi: %d\n", LhI); */
        }

        /* --- interval argument --- */
        else if (i_iv(a)) {
            /* printf("interval argument!\n");   */

            if (r_lt(r_addd(a.INF,*pil),a.SUP))
                 rc = 1; /* error singularity */
            else {

            rc =  i_invp(Ltan, a.INF, &res.INF, &dummy);
            jlb = (unsigned) LhI;        /* entier(lb/(pi/4)) mod 4 */

            /* printf("entier part Lhi lb mod 16: %u\n", LhI); */

            rc += i_invp(Ltan, a.SUP, &dummy, &res.SUP);
            jub = (unsigned) LhI;        /* entier(lb/(pi/4)) mod 4 */

            /* printf("entier part Lhi ub mod 16: %u\n", LhI); */


            /* --- d --- */
            if(jub<jlb)
                d = (jub+16)-jlb;       /* d in {0,..,15} */
            else
                d = jub - jlb;          /* d in {0,..,15} */
            if (d>=4) rc = 1;           /* error singularity */
            jlb%=4;
            jub%=4;

            if(jub<jlb)
                d = (jub+4)-jlb;        /* d in {0,..,3} */
            else
                d = jub - jlb;          /* d in {0,..,3} */

            if (jlb <= 1 && jub >= 2)   /* singularity */
                rc = 1;

            /* printf("mod 4: d = jub - jlb: %u\n", d);
            printf("jlb: %u jub: %u\n", jlb, jub); */
        }} /* interval argument */
        else
            /* --- no correct interval --- */
            rc = 1;

        /* --- error --- */
        if (rc)
            e_trap(INV_ARG,4,E_TDBL+E_TEXT(5),&a.INF,
                             E_TDBL+E_TEXT(6),&a.SUP);

        E_SPOPP("i_tan")
        return(res);
}





