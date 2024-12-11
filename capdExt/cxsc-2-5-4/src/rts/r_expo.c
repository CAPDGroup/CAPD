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

/* CVS $Id: r_expo.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_expo.c                              */
/*                                                              */
/*      Entries         : a_intg r_expo(m)                      */
/*                        a_real m;                             */
/*                                                              */
/*      Arguments       : m = mantissa                          */
/*                                                              */
/*      Description     : Base 2 exponent of mantissa with      */
/*                        0.5 <= |mantissa| < 1.0.              */
/*                                                              */
/*      Note            : -MAXINT is returned if r==0.0.        */
/*                         MAXINT is returned if r=+/-infinity  */
/*                        qNaN causes an INV_ARG exception      */
/*                        sNaN causes an INV_OP  exception      */
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
local a_intg r_expo(a_real m)
#else
local a_intg r_expo(m)

a_real m;
#endif
        {
        a_bool vz;
        a_intg expo;
        a_btyp mant[D_U_RATIO];

        E_TPUSH("r_expo")

        if (b_deko(m,&expo,mant,&vz))
           {
           expo = -MAXINT;
           }
        else if (expo>EXPO_MAX)
           {
           if (MANT_INFINITY(mant))
              expo = MAXINT;
           else if (SIGNALING(mant[0]))
              e_trap(INV_OP+E_IEEE,4,E_TMSG,5,E_TDBL+E_TEXT(7),&m);
           else
              e_trap(INV_ARG,2,E_TDBL+E_TEXT(7),&m);
           }
        else
           {
           expo++;
           if ((mant[0] & HIDDEN_BIT)==ZERO)
              {
              do {
                 b_shl1(mant,D_U_RATIO);
                 expo--;
                 }
              while ((mant[0] & HIDDEN_BIT)==ZERO);
              }
           }

        E_TPOPP("r_expo")
        return(expo);
        }





