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

/* CVS $Id: r_sign.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_sign.c                              */
/*                                                              */
/*      Entries         : a_intg r_sign(a)                      */
/*                        a_real a;                             */
/*                                                              */
/*      Arguments       : a = real value                        */
/*                                                              */
/*      Description     : Sign of normalized IEEE value.        */
/*                                                              */
/*      Function value  : -1 if value is negative               */
/*                         0 if value is zero or error occurred */
/*                         1 if value is positive               */
/*                       error if quiet or signaling NaN        */
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
local a_intg r_sign(a_real a)
#else
local a_intg r_sign(a)

a_real a;
#endif
        {
        a_intg expoa;
        a_bool vza;
        a_btyp manta[D_U_RATIO];
        register a_btyp *u;

        E_TPUSH("r_sign")

        /* decompose IEEE value                                 */

        u = (a_btyp *)&a;


        /* sign of real value   */
        vza = (u[B_HPART] & MSB) ? TRUE : FALSE;


        /* exponent of number   */
        if ((expoa = ((u[B_HPART] & EXPO_MASK)>>EXPO_SHIFT)-CHARAC)
            ==-CHARAC)
           {
		 /* zero exponent field*/

           /* high part of mantissa */
           manta[0] = (u[B_HPART] & MANT_HIGH) | HIDDEN_BIT;

           /* value is zero     */
           if (!((manta[0] &= NOT_HIDDEN_BIT) | u[B_LPART]))
		 {
              expoa=0;
		    goto _ready;
	 	 }
           /* denormalized number  */
           expoa = -CHARAC+1;
           }



        if (expoa>EXPO_MAX)
           {
           /* low part of mantissa */
           manta[1] = u[B_LPART];

           /* high part of mantissa */
           manta[0] = (u[B_HPART] & MANT_HIGH) | HIDDEN_BIT;

           /* check for infinity */
           if (MANT_INFINITY(manta))
              expoa = (vza) ? -1 : 1;
           else
              {
              if (SIGNALING(manta[0]))
                 e_trap(INV_OP+E_IEEE,4,E_TMSG,5,E_TDBL+E_TEXT(7),&a);
              else
                 e_trap(INV_OP+E_IEEE,4,E_TMSG,14,E_TDBL+E_TEXT(7),&a);
              expoa = 0;
              }
           }
        else
           expoa = (vza) ? -1 : 1;
_ready:
        E_TPOPP("r_sign")
        return(expoa);
        }






