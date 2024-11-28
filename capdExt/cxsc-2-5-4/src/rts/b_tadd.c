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

/* CVS $Id: b_tadd.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tadd.c                              */
/*                                                              */
/*      Entries         : int  b_tadd(a,b,res)                  */
/*                        tenbyte *a,*b,*res;                   */
/*                                                              */
/*      Arguments       : a = first operand of addition         */
/*                        b = second operand of addition        */
/*                        res = sum of operands                 */
/*                                                              */
/*      Description     : Addition of numbers.                  */
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
local int b_tadd(tenbyte *a,tenbyte *b,tenbyte *res)
#else
local int b_tadd(a,b,res)

tenbyte *a;
tenbyte *b;
tenbyte *res;
#endif
        {
        a_btyp manta[BSIZE];
        a_btyp mantb[BSIZE];
        a_intg expoa,expob;
        a_bool  vza,vzb;
        a_btyp hu;
        a_intg i;

        E_TPUSH("b_tadd")

        /* decompose IEEE numbers, return if zero               */
        if (b_tdek(a,&expoa,manta,&vza))
           {
           C_COPY(res->c,b->c,t_size)

           E_TPOPP("b_tadd")
           return(0);
           }
        if (b_tdek(b,&expob,mantb,&vzb))
           {
           C_COPY(res->c,a->c,t_size)

           E_TPOPP("b_tadd")
           return(0);
           }

        /* order operands in order to get |a|>=|b|              */
        if ( expob > expoa ) {
                i = expoa; expoa = expob; expob = i;
                i = (a_intg) vza;   vza = vzb;  vzb = (a_bool) i;
                hu = manta[0];
                manta[0] = mantb[0]; mantb[0] = hu;
                hu = manta[1];
                manta[1] = mantb[1]; mantb[1] = hu;
                hu = manta[2];
                manta[2] = mantb[2]; mantb[2] = hu;
                }
        else if ( expoa == expob ) {
                for ( i=0; i<3; i++ )
                    if ( manta[i]<mantb[i] ) {
                        i = (a_intg) vza; vza = vzb; vzb = (a_bool) i;
                        hu = manta[0];
                        manta[0] = mantb[0]; mantb[0] = hu;
                        hu = manta[1];
                        manta[1] = mantb[1]; mantb[1] = hu;
                        hu = manta[2];
                        manta[2] = mantb[2]; mantb[2] = hu;
                        break;
                        }
                    else if ( manta[i]>mantb[i] ) break;
                }

        /* 1. case                                              */
        /* difference of exponents >= length of mantissa+2      */
        if ( expoa-expob>=tMANTL+2 ) {

                /* propagate a borrow to mantissa               */
                if (vza!=vzb) b_subc(manta+2);

                /* set lsb of of result mantissa                */
                else manta[2] |= LSB;
                }

        /* 2. case                                              */
        /* difference of exponents < length of mantissa+2       */
        else {

                /* adjust mantissa of mantb to manta            */
                b_shru( mantb, BSIZE, expoa-expob );

                /* add basetype mantissas                       */
                if ( vza==vzb ) {
                        (void)b_addm(BSIZE,manta,mantb);
                        }

                /* subtract basetype mantissas                  */
                else {
                        (void)b_subm(BSIZE,manta,mantb);

                        /* return zero if no bits are set            */
                        if (b_test(BSIZE,manta)) {
                            for (i=0;i<t_size;i++) res->c[i] = 0;

                            E_TPOPP("b_tadd")
                            return(0);
                            }
                        }
                }

        /* normalization of mantissa                            */
        expoa += BITS_PER_CHAR;
        while ( (MSB & *manta)==ZERO ) {
                b_shl1( manta, BSIZE );
                expoa--;
                }

        /* adjust denormalized number                           */
        b_tadj(manta,&expoa);

        /* round                                                */
        b_trnd(manta,&expoa,vza);

        /* compose result                                       */
        b_tcom(res,expoa,manta,vza);

        E_TPOPP("b_tadd")
        return(0);
        }





