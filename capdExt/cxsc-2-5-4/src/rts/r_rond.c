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

/* CVS $Id: r_rond.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_rond.c                              */
/*                                                              */
/*      Entries         : a_intg r_rond(a)                      */
/*                        a_real a;                             */
/*                                                              */
/*      Arguments       : a = IEEE value                        */
/*                                                              */
/*      Description     : Round real value to next integer      */
/*                        number.                               */
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
local a_intg r_rond(a_real a)
#else
local a_intg r_rond(a)

a_real a;
#endif
        {
        a_intg expoa;
        a_intg res = 0;
        a_btyp manta[BSIZE];
        a_bool vza;

        E_TPUSH("r_rond")

        B_CLEAR(manta, BSIZE);

        if (b_deko(a,&expoa,manta,&vza))
           {
           }
        else if (expoa<-1)
           {
           }
        else
           {
           if (expoa>EXPO_MAX)
              {
              if (MANT_INFINITY(manta))
                 e_trap(INV_ARG,4,E_TMSG,13,E_TDBL+E_TEXT(7),&a);
              else if (SIGNALING(manta[0]))
                 e_trap(INV_OP+E_IEEE,4,E_TMSG,5,E_TDBL+E_TEXT(7),&a);
              else
                 e_trap(INV_ARG,4,E_TMSG,14,E_TDBL+E_TEXT(7),&a);
              }
           else
              {

              /* shift to 1 position before decimal point */
              b_shru(manta,BSIZE,(MANTL-1-1)-expoa);

              /* perform rounding by incrementing */
              b_addc(manta+(D_U_RATIO-1));

              /* shift mantissa to decimal point */
              manta[1] =
                 (manta[1]>>1) | (manta[0]<<(B_LENGTH-1));
              manta[0] >>= 1;

              if (expoa>=I_LENGTH-1)
                 {
                 if (vza && expoa==I_LENGTH-1 &&
                     manta[D_U_RATIO-1]<MSB)
                    e_trap(OVERFLOW,2,E_TDBL+E_TEXT(7),&a);
                 else
                    res = -manta[D_U_RATIO-1];
                 }
              else if (vza)
                 res = -manta[D_U_RATIO-1];
              else
                 res = manta[D_U_RATIO-1];
              }
           }

        E_TPOPP("r_rond")
        return(res);
        }





