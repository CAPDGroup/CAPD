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

/* CVS $Id: r_ne.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_ne.c                                */
/*                                                              */
/*      Entries         : a_bool r_ne(a,b)                      */
/*                        a_real a,b;                           */
/*                                                              */
/*      Arguments       : a = IEEE operand                      */
/*                        b = IEEE operand                      */
/*                                                              */
/*      Function value  : TRUE if values are not equal          */
/*                        FALSE if values are equal             */
/*                                                              */
/*      Description     : compare IEEE values                   */
/*                                                              */
/*      Note            : FALSE returned in case of error       */
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
local a_bool r_ne(a_real a,a_real b)
#else
local a_bool r_ne(a,b)

a_real a;
a_real b;
#endif
        {
        a_bool res,vza,vzb,zeroa,zerob;
        a_intg expoa,expob;
        a_btyp manta[D_U_RATIO];
        a_btyp mantb[D_U_RATIO];

        E_TPUSH("r_ne")

        zeroa = b_deko(a,&expoa,manta,&vza);
        zerob = b_deko(b,&expob,mantb,&vzb);

        if (expoa>=EXPO_MAX+1 || expob>=EXPO_MAX+1)
           {
           if (expoa==EXPO_MAX+1 && expob==EXPO_MAX+1)
              {
#if C_P_3
              if (manta[0]==mantb[0] && manta[1]==mantb[1])
#else
              if (b_bmcm(D_U_RATIO,manta,mantb)==0)
#endif
                 {
                 if (MANT_INFINITY(manta))
                    res = (vza!=vzb) ? TRUE : FALSE;
                 else res = FALSE;
                 }
              else res = TRUE;
              }
           else
              res = TRUE;
           }
        else if (zeroa) res = (zerob) ? FALSE : TRUE;
        else if (zerob) res = TRUE;
        else if (vza!=vzb) res = TRUE;
        else if (expoa!=expob) res = TRUE;
#if C_P_3
        else if ((manta[0]!=mantb[0]) ||
                 (manta[1]!=mantb[1])) res = TRUE;
#else
        else if (b_bmcm(D_U_RATIO,manta,mantb[1])!=0) res = TRUE;
#endif
        else res = FALSE;

        E_TPOPP("r_ne")
        return(res);
        }





