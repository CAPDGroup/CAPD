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

/* CVS $Id: b_tmul.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tmul.c                              */
/*                                                              */
/*      Entries         : int  b_tmul(a,b,res)                  */
/*                        tenbyte *a,*b,*res;                   */
/*                                                              */
/*      Arguments       : a = first operand of product          */
/*                        b = second operand of product         */
/*                        res = product of operands             */
/*                                                              */
/*      Description     : Multiplication of numbers.            */
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
local int b_tmul(tenbyte *a,tenbyte *b,tenbyte *res)
#else
local int b_tmul(a,b,res)

tenbyte *a;
tenbyte *b;
tenbyte *res;
#endif
        {
        a_btyp manta[BSIZE];
        a_btyp mantb[BSIZE];
        a_btyp lang[BSIZE+1];
        a_intg expoa,expob,expoc;
        a_bool  vza,vzb,vzc,zeroa,zerob;

        E_TPUSH("b_tmul")

        /* decomposition                                        */
        zeroa = b_tdek(a,&expoa,manta,&vza);
        zerob = b_tdek(b,&expob,mantb,&vzb);

        /* sign of product                                      */
        vzc = vza^vzb;

        /* initialize exact mantissa array lang                 */
        B_CLEAR(lang,BSIZE+1)

        /* product is zero                                      */
        if (zeroa || zerob) {
            b_tcom(res,-tCHARAC,lang,vzc);

            E_TPOPP("b_tmul")
            return(0);
            }

        /* add exponents                                        */
        expoc = expoa+expob;

        /* determine exact product                              */
        b_muad(manta[2],mantb[5-2-1],lang+5);
        b_muad(manta[1],mantb[4-1-1],lang+4);
        b_muad(manta[0],mantb[3-0-1],lang+3);

        b_muad(manta[2],mantb[4-2-1],lang+4);
        b_muad(manta[1],mantb[3-1-1],lang+3);
        b_muad(manta[0],mantb[2-0-1],lang+2);

        b_muad(manta[2],mantb[3-2-1],lang+3);
        b_muad(manta[1],mantb[2-1-1],lang+2);
        b_muad(manta[0],mantb[1-0-1],lang+1);

        /* adjust mantissa                                      */
        b_shlu(lang,BSIZE+1,2*tSHIFT);

        /* normalization of mantissa                            */
        expoc++;
        while ( (MSB & *lang)==ZERO ) {
                b_shl1(lang,BSIZE+1);
                expoc--;
                }

        /* adjust denormalized number                           */
        b_tadj(lang,&expoc);

        /* round                                                */
        b_trnd(lang,&expoc,vzc);

        /* composition                                          */
        b_tcom(res,expoc,lang,vzc);

        E_TPOPP("b_tmul")
        return(0);
        }





