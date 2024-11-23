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

/* CVS $Id: r_muld.c,v 1.22 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_muld.c                              */
/*                                                              */
/*      Entries         : a_real r_muld(a,b)                    */
/*                        a_real a,b;                           */
/*                                                              */
/*      Arguments       : a = first IEEE operand of product     */
/*                        b = second IEEE operand of product    */
/*                                                              */
/*      Description     : Multiplication of IEEE numbers.       */
/*                        Result is rounded downwards.          */
/*                                                              */
/*      Include         : mulbody.h                             */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_zero;
#ifndef IEEE_HARDWARE
extern a_bool e_efie;
extern a_bool e_efof;
extern a_bool e_efuf;
extern a_bool e_ofie;
extern a_bool e_ofof;
extern a_bool e_ofuf;
#endif
#endif

#ifdef LINT_ARGS
local a_real r_muld(a_real a,a_real b)
#else
local a_real r_muld(a,b)

a_real a;
a_real b;
#endif
        {
        a_real res;
#ifdef IEEE_HARDWARE
        a_btyp rc;
        int msg;

        E_TPUSH("r_muld")

        a_sets(ROUND_DOWN);

        res = a*b;

        if (a_gets(&rc)!=ZERO)
           {
           if (rc==INV_OP)
              {
              msg = (r_clss(a)==E_CLS0 || r_clss(b)==E_CLS0) ? 5 : 10;
              e_trap(rc,8,E_TMSG,msg,
                          E_TDBL+E_TEXT(1),&a,
                          E_TDBL+E_TEXT(2),&b,
                          E_TDBL+E_TRES,&res);
              }
           else
              e_trap(rc,8,E_TDBL+E_TEXT(1),&a,E_TDBL+E_TEXT(2),&b,
                          E_TDBL+E_TEXT(3),&res,E_TDBL+E_TRES,&res);
           }
#else
        a_btyp manta[D_U_RATIO];
        a_btyp mantb[D_U_RATIO];
        a_btyp lang[BSIZE];
        a_intg expoa,expob,expoc;
        a_bool vza,vzb,vzc,zeroa,zerob;
        a_bool trap;
        a_btyp rc;
#if C_P_3
#else
        a_intg i,j;
#endif

        E_TPUSH("r_muld")

        /* decomposition                                        */
        zeroa = b_deko( a, &expoa, manta, &vza );
        zerob = b_deko( b, &expob, mantb, &vzb );

        /* sign of product                                      */
        vzc = vza^vzb;

        /* NaN or infinity                                      */

/* multiplication of IEEE numbers                               */
#ifdef AIX
#include "/u/p88c/runtime/real/mulbody.h"
#else
#include "mulbody.h"
#endif

        /* rounding                                             */
        if (rc==NO_ERROR) rc = b_rndd( lang, &expoc, vzc );
        else (void)b_rndd( lang, &expoc, vzc );

        /* composition                                          */
        b_comp( &res, expoc, lang, vzc );

        if (rc) {
            switch ((int)rc)
               {
               case INEXACT   : if ((trap = e_efie)==FALSE)
                                   e_ofie = TRUE;
                                break;
               case OVERFLOW  : if ((trap = e_efof)==FALSE)
                                   e_ofof = TRUE;
                                break;
               case UNDERFLOW : if ((trap = e_efuf)==FALSE)
                                   e_ofuf = TRUE;
                                break;
               default : trap = TRUE;
               }
            if (trap) e_trap(rc+E_IEEE,8,E_TDBL+E_TEXT(1),&a,
                                         E_TDBL+E_TEXT(2),&b,
                                         E_TDBL+E_TEXT(3),&res,
                                         E_TDBL+E_TRES,&res);
            }
#endif
        E_TPOPP("r_muld")
        return(res);
        }





