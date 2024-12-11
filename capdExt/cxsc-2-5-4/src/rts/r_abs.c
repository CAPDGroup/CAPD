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

/* CVS $Id: r_abs.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_abs.c                               */
/*                                                              */
/*      Entries         : a_real r_abs(a)                       */
/*                        a_real a;                             */
/*                                                              */
/*      Arguments       : a = IEEE value                        */
/*                                                              */
/*      Description     : Return positive value of argument.    */
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
local a_real r_abs(a_real a)
#else
local a_real r_abs(a)

a_real a;
#endif
        {
        a_real res;
        a_intg expoa;
        a_bool vza;
        a_btyp manta[D_U_RATIO];

        E_SPUSH("r_abs")

        if (b_deko(a,&expoa,manta,&vza))
           b_comp(&res,expoa,manta,FALSE);
        else if (expoa>EXPO_MAX)
           {
           if (MANT_INFINITY(manta))
              b_comp(&res,expoa,manta,FALSE);
           else
              {
              if (SIGNALING(manta[0]))
                 e_trap(INV_OP+E_IEEE,2,E_TDBL+E_TEXT(7),&a);
              else
                 e_trap(INV_ARG,2,E_TDBL+E_TEXT(7),&a);
              }
           }
        else
           b_comp(&res,expoa,manta,FALSE);

        E_SPOPP("r_abs")
        return(res);
        }





