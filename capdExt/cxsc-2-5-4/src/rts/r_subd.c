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

/* CVS $Id: r_subd.c,v 1.22 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_subd.c                              */
/*                                                              */
/*      Entries         : a_real r_subd(a,b)                    */
/*                        a_real a,b;                           */
/*                                                              */
/*      Arguments       : a = first IEEE operand of subtraction */
/*                        b = second IEEE operand of subtraction*/
/*                                                              */
/*      Description     : Subtraction of IEEE numbers.          */
/*                        Result is rounded downwards.          */
/*                                                              */
/*      Include         : subbody.h                             */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
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
local a_real r_subd(a_real a,a_real b)
#else
local a_real r_subd(a,b)

a_real a;
a_real b;
#endif
        {
        a_real res;
#ifdef IEEE_HARDWARE
        a_btyp rc;
        int msg;

        E_TPUSH("r_subd")

        a_sets(ROUND_DOWN);

        res = a-b;

        if (a_gets(&rc)!=ZERO)
           {
           if (rc==INV_OP)
              {
              msg = (r_clss(a)==E_CLS0 || r_clss(b)==E_CLS0) ? 5 : 9;
              e_trap(rc,8,E_TMSG,msg,
                          E_TDBL+E_TEXT(1),&a,
                          E_TDBL+E_TEXT(2),&b,E_TDBL+E_TRES,&res);
              }
           else
              e_trap(rc,8,E_TDBL+E_TEXT(1),&a,E_TDBL+E_TEXT(2),&b,
                          E_TDBL+E_TEXT(3),&res,E_TDBL+E_TRES,&res);
           }
#else
        a_real oldb;
        a_btyp manta[BSIZE];
        a_btyp mantb[BSIZE];
        a_intg expoa,expob;
        a_intg i;
        a_btyp hu,rc;
        a_bool vza,vzb,vzz = TRUE;
        a_bool trap,zeroa,zerob;

        E_TPUSH("r_subd")

        /* decompose IEEE number                                */
        zeroa = b_deko(a,&expoa,manta,&vza);
        zerob = b_deko(b,&expob,mantb,&vzb);

        /* change sign of b                                     */
        vzb = NOT(vzb);

        /* NaN or infinity                                      */

/* subtract IEEE numbers                                        */
#ifdef AIX
#include "/u/p88c/runtime/real/subbody.h"
#else
#include "subbody.h"
#endif

        /* round downwards */
        if (rc==NO_ERROR) rc = b_rndd( manta, &expoa, vza );
        else (void)b_rndd( manta, &expoa, vza );

        /* compose result                                       */
        b_comp( &res, expoa, manta, vza );

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
                                         E_TDBL+E_TEXT(2),&oldb,
                                         E_TDBL+E_TEXT(3),&res,
                                         E_TDBL+E_TRES,&res);
            }
#endif

        E_TPOPP("r_subd")
        return(res);
        }





