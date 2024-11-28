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

/* CVS $Id: i_sin.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_sin.c                               */
/*                                                              */
/*      Entries         : a_intv i_sin(a)                       */
/*                        a_intv a;                             */
/*                                                              */
/*      Arguments       : a = interval argument of function     */
/*                                                              */
/*      Description     : Interval sine function                */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_mone,*r_one_;
#endif

extern a_btyp LhI;    /* entier(x/(pi/4)) mod 16 after call sin(x) */

#ifdef LINT_ARGS
local a_intv i_sin( a_intv a )
#else
local a_intv i_sin( a )
a_intv a;
#endif

{
        /* --- 2*Pi rounded downwards --- */
static a_btyp o_2pil[] = B_TYPINI( 0x401921fbL, 0x54442d18L );  /* chopped */

        a_real *twopil = (a_real *)&o_2pil[0];
/*      static unsigned long c2pi[2] = {0x401921fb, 0x54442d18};
        a_real            *twopil= (a_real*)&c2pi[0];    */

        a_btyp        rc;
        a_intv          res;
        a_intv          clb, cub;
        unsigned        jub, jlb, d;

        E_SPUSH("i_sin")

        /* --- point argument --- */
        if (i_point(a)) {
            /* printf("point argument!\n"); */
            if (r_sign(a.INF)==0) {
              res = a;
              rc = 0;
            }
            else {
              rc = i_invp(Lsin,a.INF,&res.INF,&res.SUP);
              /* printf("entier part Lhi: %d\n", LhI); */
            }
        }

        /* --- interval argument --- */
        else if (i_iv(a)) {
/*          printf("interval argument!\n");
*/
            rc = i_invp(Lsin,a.INF,&clb.INF,&clb.SUP);
/* korrigiert 6.8.90
            jlb = (LhI/2)%4;
*/
            jlb = (unsigned) (LhI/2);      /* entier(lb/(pi/2)) */

/*          printf("entier part Lhi lb mod 16: %u\n", LhI);
*/
            rc += i_invp(Lsin,a.SUP,&cub.INF,&cub.SUP);
            jub = (unsigned) (LhI/2);      /* entier(lb/(pi/2)) mod 4 */

/*          printf("entier part Lhi ub mod 16: %u\n", LhI);
*/
            if (rc==0) {
            /* --- d --- */
            if (jlb>jub) jub += 8;
            d = jub - jlb;        /* d in {0,..,7} */

            if (r_lt(r_addd(a.INF,*twopil),a.SUP))  d = 4;

            if (d>4) d=4;

            jlb%=4;
            jub%=4;

/*          printf("d = jub - jlb            : %u\n", d);
            printf("entier part Lhi lb mod 4 : %u\n", jlb);
            printf("entier part Lhi ub mod 4 : %u\n", jub);
*/
            switch (d) {
            case 0: /* -------------------------- d = 0 --- */
                /* jub = jlb */
                switch (jlb) {
                case 0:
                case 3: /* increasing */
                    R_ASSIGN(res.INF,clb.INF);
                    R_ASSIGN(res.SUP,cub.SUP);
                    break;
                case 1:
                case 2: /* decreasing */
                    R_ASSIGN(res.INF,cub.INF);
                    R_ASSIGN(res.SUP,clb.SUP);
                    break;
                    }
                break; /* d == 0 */

            case 1: /* -------------------------- d = 1 --- */
                /* jub = jlb + 1 */
                switch (jlb) {
                case 0:
                    R_ASSIGN(res.INF,
                             (r_lt(cub.INF,clb.INF))?cub.INF:clb.INF);
                    R_ASSIGN(res.SUP,*r_one_);
                    break;
                case 1: /* decreasing */
                    R_ASSIGN(res.INF,cub.INF);
                    R_ASSIGN(res.SUP,clb.SUP);
                    break;
                case 2:
                    R_ASSIGN(res.INF,*r_mone);
                    R_ASSIGN(res.SUP,
                             (r_gt(cub.SUP,clb.SUP))?cub.SUP:clb.SUP);
                    break;
                case 3: /* increasing */
                    R_ASSIGN(res.INF,clb.INF);
                    R_ASSIGN(res.SUP,cub.SUP);
                    break;
                }
                break; /* d == 1 */

            case 2: /* -------------------------- d = 2 --- */
                /* jub = jlb + 2 */
                switch (jlb) {
                case 0:
                    /* jub = 2 */
                    R_ASSIGN(res.INF,cub.INF);
                    R_ASSIGN(res.SUP,*r_one_);
                    break;
                case 1:
                    /* jub = 3 */
                    R_ASSIGN(res.INF,*r_mone);
                    R_ASSIGN(res.SUP,clb.SUP);
                    break;
                case 2:
                    /* jub = 0 */
                    R_ASSIGN(res.INF,*r_mone);
                    R_ASSIGN(res.SUP,cub.SUP);
                    break;
                case 3:
                    /* jub = 1 */
                    R_ASSIGN(res.INF,clb.INF);
                    R_ASSIGN(res.SUP,*r_one_);
                    break;
                }
                break; /* d == 2 */

            case 3: /* -------------------------- d = 3 --- */
                /* jub = jlb + 3 */
                switch (jlb) {
                case 0:
                case 2:
                    R_ASSIGN(res.INF,*r_mone);
                    R_ASSIGN(res.SUP,*r_one_);
                    break;
                case 1:
                    R_ASSIGN(res.INF,*r_mone);
                    R_ASSIGN(res.SUP,
                             (r_gt(cub.SUP,clb.SUP))?cub.SUP:clb.SUP);
                    break;
                case 3:
                    R_ASSIGN(res.INF,
                             (r_lt(cub.INF,clb.INF))?cub.INF:clb.INF);
                    R_ASSIGN(res.SUP,*r_one_);
                    break;
                }
                break; /* d == 3 */

            case 4: /* -------------------------- d = 4 --- */
                /* jub = jlb + 4 */
                R_ASSIGN(res.INF,*r_mone);
                R_ASSIGN(res.SUP,*r_one_);
                break; /* d == 4 */
            } /* switch(jdiff) */

          }
        if (r_lt(res.INF,*r_mone)) R_ASSIGN(res.INF,*r_mone);
        if (r_gt(res.SUP,*r_one_)) R_ASSIGN(res.SUP,*r_one_);
        }/* interval argument */
        else
            /* --- no correct interval --- */
            rc = 1;

        /* --- error --- */
        if (rc)
            e_trap(INV_ARG,4,E_TDBL+E_TEXT(5),&a.INF,
                             E_TDBL+E_TEXT(6),&a.SUP);

        E_SPOPP("i_sin")
        return(res);
}





