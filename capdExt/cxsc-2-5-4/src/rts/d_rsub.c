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

/* CVS $Id: d_rsub.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_rsub.c                              */
/*                                                              */
/*      Entries         : void d_rsub(c,a)                      */
/*                        d_otpr *c;                            */
/*                        a_real a;                             */
/*                                                              */
/*      Arguments       : c = dotprecision variable             */
/*                        a = real value                        */
/*                                                              */
/*      Description     : subtract real value from dotprecision */
/*                        variable                              */
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
local void d_rsub(d_otpr *c,a_real a)
#else
local void d_rsub(c,a)

d_otpr *c;
a_real a;
#endif
        {
        a_btyp manta[D_U_RATIO+1];
        a_intg expoa,c_ge_a;
        a_intg cout,i,k;
        a_btyp *d;
        a_bool vza;

        E_TPUSH("d_rsub")

        /* decomposition                                        */
        if (b_deko( a, &expoa, manta, &vza ))
           {
           (*c)[A_STATUS] |= (vza) ? A_MZERO : A_PZERO;
           E_TPOPP("d_rsub")
           return;
           }

        (*c)[A_STATUS] |= A_PZERO+A_MZERO;

        manta[D_U_RATIO] = ZERO;

        /* NaN or infinity                                      */
        if ( expoa>EXPO_MAX ) {

                /* a = infinity                                 */
                if (MANT_INFINITY(manta)) {
                   if (((*c)[A_STATUS] & A_QUIETNAN)==ZERO)
                      {
                      if (((*c)[A_STATUS] &
                         (A_MINFINITY | A_PINFINITY))==ZERO)
                         {
                         e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,4,E_TMSG,13,
                                E_TDBL,&a);
                         if (vza!=FALSE)
                            (*c)[A_STATUS] |= A_PINFINITY;
                         else
                            (*c)[A_STATUS] |= A_MINFINITY;
                         }
                      else if (NOT(vza) &&
                               ((*c)[A_STATUS] & A_PINFINITY))
                         {
                         e_trap(INV_OP+E_IEEE,6,E_TMSG,9,E_TDTP,c,E_TDBL,&a);
                         (*c)[A_STATUS] |= A_QUIETNAN;
                         (*c)[A_LOWNAN] = INV_OP;
                         }
                      else if (vza && ((*c)[A_STATUS] & A_MINFINITY))
                         {
                         e_trap(INV_OP+E_IEEE,6,E_TMSG,9,E_TDTP,c,E_TDBL,&a);
                         (*c)[A_STATUS] |= A_QUIETNAN;
                         (*c)[A_LOWNAN] = INV_OP;
                         }
                      }
                   E_TPOPP("d_rsub")
                   return;
                   }

                /* a = NaN                                      */
                if (SIGNALING(manta[0]))
                   {
                   e_trap(INV_OP+E_IEEE,4,E_TMSG,5,E_TDBL,&a);
                   (*c)[A_STATUS] |= A_QUIETNAN;
                   (*c)[A_LOWNAN] = INV_OP;
                   }
                else if (((*c)[A_STATUS] & A_QUIETNAN)==ZERO)
                   {
                   e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,4,E_TMSG,14,E_TDBL,&a);
                   (*c)[A_STATUS] |= A_QUIETNAN;
                   (*c)[A_LOWNAN] = ((a_btyp *)&a)[B_LPART];
                   }
                E_TPOPP("d_rsub")
                return;
                }

        if (((*c)[A_STATUS] &
             (A_QUIETNAN | A_PINFINITY | A_MINFINITY ))!=ZERO)
           {
           E_TPOPP("d_rsub")
           return;
           }

        /* allign mantissa to dotprecision                      */
        k = A_D_P-B_ASHR(expoa,LOG_B_LENGTH);
        if ( (i = (expoa & (B_LENGTH-1))-EXPO_SHIFT) < 0)
                {
                manta[2] = manta[1]<<(B_LENGTH+i);
                manta[1] = (manta[1]>>-i) |
                           (manta[0]<<(B_LENGTH+i));
                manta[0] = (manta[0]>>-i);
                }
        else if (i>0)
                {
                manta[0] =
                   (manta[0]<<i) | (manta[1]>>(B_LENGTH-i));
                manta[1] <<= i;
                }

        /* accu is empty                                        */
        if ((*c)[A_BEGIN]==ZERO) {
                (*c)[A_END] = (D_U_RATIO)+k;
                (*c)[A_BEGIN] = k;
                for (i=0;i<D_U_RATIO+1;i++)
                    (*c)[k+i] = manta[i];
                (*c)[A_SIGN] = NOT(vza);
                }

        /* subtract value from dotprecision with different sign */
        else if (vza!=(*c)[A_SIGN]) {
                if ((*c)[A_END]<(D_U_RATIO)+k)
                    (*c)[A_END] = (D_U_RATIO)+k;
                if ((*c)[A_BEGIN]>k) (*c)[A_BEGIN] = k;

                if (b_addm(D_U_RATIO+1,&((*c)[k]),manta)) {
                    b_addc(&((*c)[(-1)+k]));

                    /* carry extends used area of accu          */
                    if ((*c)[(*c)[A_BEGIN]-1]) (*c)[A_BEGIN]--;
                    }
                }

        /* subtract value from dotprecision with same sign      */
        else    {

                int m;
    
                /* consider leading zeros for denormalized numbers */
                for (m=0;manta[m]==0;m++) k++;

                /* compare dotprecision value with real value  */
                if ((c_ge_a=((*c)[A_BEGIN]<=k) ? TRUE : FALSE)!=FALSE)
                    {
                    if ((*c)[A_BEGIN]==k) {
                        for (i=m;i<D_U_RATIO+1;i++)
                            if ((*c)[k+i-m]>manta[i]) break;
                            else if ((*c)[k+i-m]<manta[i]) {
                                c_ge_a = FALSE;
                                break;
                                }
                        }
                    }
                if (c_ge_a) {
                    if (b_subm(D_U_RATIO+1-m,(*c)+k,manta+m))
                        b_subc((*c)+(-1)+k);
                    }
                else {
                    cout = 0;
                    for (i=(*c)[A_END];i>D_U_RATIO+k-m;i--)
                        b_subu(ZERO,(*c)[i],cout,
                               (*c)+i,&cout);

                    for (i=D_U_RATIO-m;i>=0;i--)
                        b_subu(manta[i+m],(*c)[k+i],cout,
                              (*c)+k+i,&cout);

                    (*c)[A_SIGN] = 1-(*c)[A_SIGN];
                    }

                if ((*c)[A_END]<(D_U_RATIO)+k-m)
                    (*c)[A_END] = (D_U_RATIO)+k-m;
                if ((*c)[A_BEGIN]>k) (*c)[A_BEGIN] = k;
                }

        /* eliminate leading zeros                  */
        /* required because of denormalized numbers */
        d = &((*c)[(*c)[A_BEGIN]]);
        while (*d++==ZERO)
                {
                if (++((*c)[A_BEGIN])>(*c)[A_END])
                    {
                    (*c)[A_BEGIN] = (*c)[A_END] =
                        (*c)[A_SIGN] = ZERO;
                    break;
                    }
                }

        /* eliminate trailing zeros     */
        if ((*c)[A_BEGIN])
            {
            d = (*c)+(*c)[A_END];
            while (*d--==ZERO) (*c)[A_END]--;
            }

        E_TPOPP("d_rsub")
        }





