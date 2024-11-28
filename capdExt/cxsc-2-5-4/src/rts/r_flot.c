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

/* CVS $Id: r_flot.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_flot.c                              */
/*                                                              */
/*      Entries         : a_real r_flot(i)                      */
/*                        a_intg i;                             */
/*                                                              */
/*      Arguments       :                                       */
/*                        i = integer value to be converted     */
/*                                                              */
/*      Description     : Convert integer value to real.        */
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
local a_real r_flot(a_intg i)
#else
local a_real r_flot(i)

a_intg i;
#endif
        {
        a_real m;
        a_btyp mant[D_U_RATIO];
        a_bool vz = FALSE;
        a_intg expo;

        E_TPUSH("r_flot")

        B_CLEAR(mant,D_U_RATIO)

        if (i==0) expo = -CHARAC;
        else
           {
           if (i<0)
              {
              vz = TRUE;
              mant[0] = (i==MININT) ? (a_btyp)i : -i;
              }
           else
              mant[0] = i;

           expo = EXPO_SHIFT;

           /* normalization of mantissa                         */
           if ( SHFT_MASK & mant[0] )
              {
              do
                 {
                 b_shr1(mant,D_U_RATIO);
                 expo++;
                 }
              while (SHFT_MASK & mant[0]);
              }
           else
              {
              while ( (HIDDEN_BIT & mant[0])==ZERO )
                 {
                 mant[0] <<= 1;
                 expo--;
                 }
              }
           }

        b_comp(&m,expo,mant,vz);

        E_TPOPP("r_flot")
        return(m);
        }





