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

/* CVS $Id: r_succ.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_succ.c                              */
/*                                                              */
/*      Entries         : a_real r_succ(a)                      */
/*                        a_real a;                             */
/*                                                              */
/*      Arguments       : a = real value                        */
/*                                                              */
/*      Function value  : Succeding normalized IEEE value.      */
/*                                                              */
/*      Description     : Succeding normalized IEEE value.      */
/*                        +infinity returned if overflow.       */
/*                        smallest positive normalized or       */
/*                        denormalized number returned if a==0. */
/*                                                              */
/*      Note            : infinity and NaN are returned         */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_bool f_ppdn;
#endif

#ifdef LINT_ARGS
local a_real r_succ(a_real a)
#else
local a_real r_succ(a)

a_real a;
#endif
        {
        a_intg expoa;
        a_bool vza;
        a_btyp manta[D_U_RATIO];

        E_TPUSH("r_succ")

        /* decompose IEEE value                                 */
        if (b_deko(a,&expoa,manta,&vza))
           {
           vza = FALSE;
           expoa = EXPO_MIN;
#ifdef NORMALIZE_ENABLE
           if (f_ppdn)
              manta[0] = HIDDEN_BIT;
           else
#endif
              manta[1] = LSB;
           }
        else if (expoa>EXPO_MAX)
           {
           E_TPOPP("r_succ")
           return(a);
           }
        else if (NOT(vza))
           {
           b_addc(manta+(D_U_RATIO-1));
           if (SHFT_MASK & *manta) {
#if C_P_3
              manta[1] =
                 (manta[1]>>1) | (manta[0]<<(B_LENGTH-1));
              manta[0] >>= 1;
#else
              b_shr1(manta,D_U_RATIO);
#endif
              expoa++;
              }
           }
        else
           {
           b_subc(manta+(D_U_RATIO-1));
           if ((HIDDEN_BIT & *manta)==ZERO)
              {
              if (--expoa<EXPO_MIN)
                 {
#ifdef NORMALIZE_ENABLE
                 if (f_ppdn)
                    {
#if C_P_3
                    manta[0] = manta[1] = ZERO;
#else
                    B_CLEAR(manta,(a_intg)D_U_RATIO);
#endif
                    vza = FALSE;
                    }
                 else
#endif
                    {
                    expoa = EXPO_MIN;
                    }
                 }
              else
                 {
#if C_P_3
                 manta[0] = (manta[0]<<1) |
                            (manta[1]>>(B_LENGTH-1));
                 manta[1] <<= 1;
#else
                 b_shl1(manta,(a_intg)D_U_RATIO);
#endif
                 manta[D_U_RATIO-1] |= LSB;
                 }
              }
           }

        /* compose result                                       */
        b_comp(&a,expoa,manta,vza);

        E_TPOPP("r_succ")
        return(a);
        }





