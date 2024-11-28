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

/* CVS $Id: b_mdiv.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_mdiv.c                              */
/*                                                              */
/*      Entries         : void b_mdiv(ma,mb,mc,expo)            */
/*                        a_btyp *ma,*mb,*mc;                   */
/*                        a_intg *expo;                         */
/*                                                              */
/*      Arguments       : ma = mantissa of length D_U_RATIO     */
/*                        mb = mantissa of length D_U_RATIO     */
/*                        mc = mantissa of length D_U_RATIO+1   */
/*                        expo = exponent of quotient           */
/*                                                              */
/*      Description     : Division of a_btyp mantissas          */
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
local void b_mdiv(a_btyp *ma,a_btyp *mb,a_btyp *mc,a_intg *expo)
#else
local void b_mdiv(ma,mb,mc,expo)

a_btyp *ma;
a_btyp *mb;
a_btyp *mc;
a_intg *expo;
#endif
        {

        a_intg i,j,block;
        a_btyp m[D_U_RATIO],mm[D_U_RATIO],f,h;
        a_bool exact = FALSE;

        /* initialization of mc[]                               */
#if C_P_3
        mc[1] = mc[2] = ZERO;
#else
        B_CLEAR(mc,D_U_RATIO+1)
#endif
        *mc = HIDDEN_BIT;
        block = EXPO_SHIFT;

        /* shift mantissa of numerator left if mb[]>ma[]        */
        for (i=0; i<D_U_RATIO; i++)
                if (ma[i]>mb[i]) break;
                else if (ma[i]<mb[i]) {
#if C_P_3
                    ma[0] = (ma[0]<<1) | (ma[1]>>(B_LENGTH-1));
                    ma[1] <<= 1;
#else
                    b_shl1( ma, D_U_RATIO );
#endif
                    (*expo)--;
                    break;
                    }

        /* subtract denominator once for initialization         */
        (void)b_subm((a_intg)D_U_RATIO,ma,mb);

        /* division loop                                        */
        for (j=1;j<MANTL+1;j += MIN_FAC) {

                /* shift mantissa of numerator                  */
#if C_P_3
                ma[0] = (ma[0]<<MIN_FAC) |
                        (ma[1]>>(B_LENGTH-MIN_FAC));
                ma[1] <<= MIN_FAC;
#else
                b_shlu( ma, D_U_RATIO, MIN_FAC );
#endif

                /* subtract multiple if non-zero factor         */
                if ( (f = *ma/(*mb+1)) >1) {
                    mm[D_U_RATIO-1] = ZERO;
#if C_P_3
                    m[1] = GETLOW(mb[1])*f;
                    h = GETHIGH(mb[1])*f;
                    mm[1] |= MOVEHIGH(h);
                    mm[0] = GETHIGH(h);
#else
                    for (i=D_U_RATIO-1;i>0;i--) {
                        m[i] = GETLOW(mb[i])*f;
                        h = GETHIGH(mb[i])*f;
                        mm[i] |= MOVEHIGH(h);
                        mm[i-1] = GETHIGH(h);
                        }
#endif
                    m[0] = mb[0]*f;
                    (void)b_subm(D_U_RATIO,ma,&m[0]);
                    (void)b_subm(D_U_RATIO,ma,&mm[0]);
                    }
                else if (f==1) {
                    (void)b_subm(D_U_RATIO,ma,mb);
                    }

                /* compare mantissas (at most one subtraction)  */
                for (i=0;i<D_U_RATIO;i++)
                    if ( ma[i]<mb[i] ) break;
                    else if ( ma[i]>mb[i] ) {
                        (void)b_subm(D_U_RATIO,ma,mb);
                        f++;
                        break;
                        }
                if (i==D_U_RATIO) {
                    (void)b_subm(D_U_RATIO,ma,mb);
                    f++;
                    exact = TRUE;
                    }

                /* determine bits of quotient                   */
                if ( (block -= MIN_FAC) <0) {
                        *mc |= (f>>-block);
                        block += B_LENGTH;
                        mc++;
                        }
                *mc |= (f<<block);
                }

        /* set lsb of quotient if remainder is non-zero         */
        if (NOT(exact)) {
#if C_P_3
            if (ma[0] || ma[1]) *mc |= LSB;
#else
            for ( i=0; i<D_U_RATIO; i++ )
            if (ma[i]) {
                *mc |= LSB;
                return;
                }
#endif
            }

        return;
        }





