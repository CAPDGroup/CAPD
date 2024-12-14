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

/* CVS $Id: b_tmdv.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tmdv.c                              */
/*                                                              */
/*      Entries         : void b_tmdv(ma,mb,mc,expo)            */
/*                        a_btyp *ma,*mb,*mc;                   */
/*                        a_intg *expo;                         */
/*                                                              */
/*      Arguments       : ma = mantissa of length D_U_RATIO     */
/*                        mb = mantissa of length D_U_RATIO     */
/*                        mc = mantissa of length D_U_RATIO+1   */
/*                        expo = exponent of quotient           */
/*                                                              */
/*      Description     : Division of basetype mantissas        */
/*                                                              */
/*      External        :                                       */
/*                        b_shlu - shift basetype left          */
/*                        b_shl1 - shift left 1 position        */
/*                        b_subm - subtract basetype arrays     */
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
local void b_tmdv(a_btyp *ma,a_btyp *mb,a_btyp *mc,a_intg *expo)
#else
local void b_tmdv(ma,mb,mc,expo)

a_btyp ma[BSIZE];
a_btyp mb[BSIZE];
a_btyp mc[BSIZE];
a_intg *expo;
#endif
        {

        a_intg i,j,block;
        a_btyp m[BSIZE],mm[BSIZE],f,h;
        a_bool  exact = FALSE;

        E_TPUSH("b_tmdv")

        /* initialization of mc[]                               */
        B_CLEAR(mc,BSIZE)
        B_CLEAR(m, BSIZE)
        B_CLEAR(mm,BSIZE) 

        /* shift mantissa of numerator left if mb[]>ma[]        */
        for (i=0; i<BSIZE; i++)
           {
           if (ma[i]>mb[i]) break;
           else if (ma[i]<mb[i])
              {
              b_shl1(ma,BSIZE);
              (*expo)--;
              break;
              }
           }

        *mc = MSB;
        (void)b_subm(BSIZE,ma,mb);

        block = B_LENGTH-1;

        /* division loop                                        */
        for (j=1;j<tMANTL+1;j += tSHIFT) {

                /* shift mantissa of numerator                  */
                b_shlu(ma,BSIZE,tSHIFT);

                /* subtract multiple if non-zero factor         */
                if ( (f = *ma/(*mb+1)) >1) {
                    mm[2] = m[2] = ZERO;
                    h = GETHIGH(mb[2])*f;
                    mm[2] |= MOVEHIGH(h);
                    mm[1] = GETHIGH(h);
                    m[1] = GETLOW(mb[1])*f;
                    h = GETHIGH(mb[1])*f;
                    mm[1] |= MOVEHIGH(h);
                    mm[0] = GETHIGH(h);
                    m[0] = mb[0]*f;
                    (void)b_subm(BSIZE,ma,&m[0]);
                    (void)b_subm(BSIZE,ma,&mm[0]);
                    }
                else if (f==1) {
                    (void)b_subm(BSIZE,ma,mb);
                    }

                /* compare mantissas (at most one subtraction)  */
                for (i=0;i<BSIZE;i++)
                    if (ma[i]<mb[i]) break;
                    else if (ma[i]>mb[i]) {
                        (void)b_subm(BSIZE,ma,mb);
                        f++;
                        break;
                        }
                if (i==BSIZE) {
                    (void)b_subm(BSIZE,ma,mb);
                    f++;
                    exact = TRUE;
                    }

                /* determine bits of quotient                   */
                if ( (block -= tSHIFT) <0) {
                        *mc |= (f>>-block);
                        block += B_LENGTH;
                        mc++;
                        }
                *mc |= (f<<block);
                }

        /* set lsb of quotient if remainder is non-zero         */
        if (NOT(exact))
            {
            if (NOT(b_test(BSIZE,ma)))
                {
                *mc |= LSB;
                E_TPOPP("b_tmdv")
                return;
                }
            }

        E_TPOPP("b_tmdv")
        return;
        }





