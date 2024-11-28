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

/* CVS $Id: r_lt.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_lt.c                                */
/*                                                              */
/*      Entries         : a_bool r_lt(a,b)                      */
/*                        a_real a,b;                           */
/*                                                              */
/*      Arguments       : a = IEEE operand                      */
/*                        b = IEEE operand                      */
/*                                                              */
/*      Function value  : TRUE if a<b                           */
/*                        FALSE if a>=b                         */
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
extern a_bool e_efio;
extern a_bool e_ofio;
#endif

#ifdef LINT_ARGS
local a_bool r_lt(a_real a,a_real b)
#else
local a_bool r_lt(a,b)

a_real a;
a_real b;
#endif
        {
        a_intg expoa,expob,i;
        a_bool res = FALSE,vza,vzb,zeroa,zerob;
        a_btyp manta[D_U_RATIO];
        a_btyp mantb[D_U_RATIO];

        E_TPUSH("r_lt")

        zeroa = b_deko(a,&expoa,manta,&vza);
        zerob = b_deko(b,&expob,mantb,&vzb);

        if (expoa>=EXPO_MAX+1 || expob>=EXPO_MAX+1)
           {
           if (expoa==EXPO_MAX+1 && MANT_INFINITY(manta))
              {
              if (expob==EXPO_MAX+1)
                 {
                 if (MANT_INFINITY(mantb))
                    res = (NOT(vzb) && vza!=vzb) ? TRUE : FALSE;
                 else
                    {
                    res = FALSE;
                    if (e_efio)
                       e_trap(INV_OP+E_IEEE,6,E_TMSG,5,
                                              E_TDBL+E_TEXT(1),&a,
                                              E_TDBL+E_TEXT(2),&b);
                    else
                       e_ofio = TRUE;
                    }
                 }
              else res = NOT(vzb);
              }
           else if (expob==EXPO_MAX+1 && MANT_INFINITY(mantb))
              {
              if (expoa!=EXPO_MAX+1) res = NOT(vzb);
              else
                 {
                 res = FALSE;
                 if (e_efio)
                    e_trap(INV_OP+E_IEEE,6,E_TMSG,5,
                                           E_TDBL+E_TEXT(1),&a,
                                           E_TDBL+E_TEXT(2),&b);
                 else
                    e_ofio = TRUE;
                 }
              }
           else
              {
              res = FALSE;
              if (e_efio)
                 e_trap(INV_OP+E_IEEE,6,E_TMSG,5,
                                        E_TDBL+E_TEXT(1),&a,
                                        E_TDBL+E_TEXT(2),&b);
              else
                 e_ofio = TRUE;
              }
           }
        else if (zeroa) res = (zerob) ? FALSE : NOT(vzb);
        else if (zerob) res = vza;
        else if (vza!=vzb) res = vza;
        else if (expoa>expob) res = vza;
        else if (expoa<expob) res = NOT(vza);
        else
           {
           for (i=0;i<D_U_RATIO;i++)
              if (manta[i]>mantb[i]) { res = vza; break; }
              else if (manta[i]<mantb[i]) { res = NOT(vza); break; }
           if (i==D_U_RATIO) res = FALSE;
           }

        E_TPOPP("r_lt")
        return(res);
        }





