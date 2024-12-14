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

/* CVS $Id: b_tdiv.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tdiv.c                              */
/*                                                              */
/*      Entries         : int  b_tdiv(a,b,res)                  */
/*                        tenbyte *a,*b,*res;                   */
/*                                                              */
/*      Arguments       : a = numerator                         */
/*                        b = denominator                       */
/*                        res = quotient                        */
/*                                                              */
/*      Description     : Division of numbers.                  */
/*                        Result is rounded according to        */
/*                        flag "b_rflg".                        */
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
local int b_tdiv(tenbyte *a,tenbyte *b,tenbyte *res)
#else
local int b_tdiv(a,b,res)

tenbyte *a;
tenbyte *b;
tenbyte *res;
#endif
        {
        a_btyp manta[BSIZE];
        a_btyp mantb[BSIZE];
        a_btyp mantc[BSIZE];
        a_intg expoa,expob,expoc,i;
        a_bool  vza,vzb,vzc,zeroa,zerob;

        E_TPUSH("b_tdiv")

        /* initialize unused result positions                   */
        for (i=D_U_RATIO+1;i<BSIZE;i++) mantc[i] = ZERO;

        /* decomposition                                        */
        zeroa = b_tdek(a,&expoa,manta,&vza);
        zerob = b_tdek(b,&expob,mantb,&vzb);

        /* sign of quotient                                     */
        vzc = vza^vzb;

        /* numerator is zero                                    */
        if ( zeroa )
        {
           /* zero/zero                        */
           if (zerob) {
              e_trap(INV_OP+E_IEEE,6,E_TMSG,2,E_TDBL,&a,E_TDBL,&b);

              E_TPOPP("b_tdiv");
              return(1);
              }

           expoc = -tCHARAC;
           b_tcom(res,expoa,manta,vzc);

           E_TPOPP("b_tdiv")
           return(0);
           }

        /* denominator is zero                                  */
        if ( zerob ) {
            e_trap(DIV_BY_ZERO+E_IEEE,4,E_TDBL,&a,E_TDBL,&b);

            E_TPOPP("b_tdiv");
            return(1);
            }

        /* denormalized numbers                                 */
        while (!(tFIRST_BIT & manta[0]))
                {
                b_shl1(manta,BSIZE);
                expoa--;
                }
        while (!(tFIRST_BIT & mantb[0]))
                {
                b_shl1(mantb,BSIZE);
                expob--;
                }

        /* exponent                                             */
        expoc = expoa-expob;

        /* division of mantissas                                */
        b_tmdv(manta,mantb,mantc,&expoc);

        /* adjust denormalized number                           */
        b_tadj(mantc,&expoc);

        /* rounding                                             */
        b_trnd(mantc,&expoc,vzc);

        /* composition                                          */
        b_tcom(res,expoc,mantc,vzc);

        E_TPOPP("b_tdiv")
        return(0);
        }





