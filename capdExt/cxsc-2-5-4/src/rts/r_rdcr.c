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

/* CVS $Id: r_rdcr.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_rdcr.c                              */
/*                                                              */
/*      Entries         : void r_rdcr(device,s,rnd,i)           */
/*                        FILE *device;                         */
/*                        a_real *s;                            */
/*                        a_intg rnd;                           */
/*                        int i;                                */
/*                                                              */
/*      Arguments       : device - input device reading string  */
/*                        s - resultant IEEE variable           */
/*                        rnd - rounding mode                   */
/*                              -1 = rounding downwards         */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        i - window character                  */
/*                                                              */
/*      Description     : Read a character string from input    */
/*                        device and convert to real format     */
/*                                                              */
/*      Note            : rounding ignored                      */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern dotprecision b_cp__;
#endif

#ifdef LINT_ARGS
local void r_rdcr(FILE *device,a_real *s,a_intg rnd,int i)
#else
local void r_rdcr(device,s,rnd,i)

FILE *device;
a_real *s;
a_intg rnd;
int i;
#endif
        {
        a_intg dp,expo,length;
        a_bool sign;

        E_TPUSH("r_rdcr")

        switch (b_cscn(device,(char *)b_cp__,&expo,&dp,&length,&sign,i))
           {
           case  1:
           case  2:
           case  3:
           case  6:
           case  7:
           case  8: e_trap(I_O_ERROR,2,E_TMSG,58);
                    E_TPOPP("r_rdcr")
                    return;
           case  4:
           case  5: e_trap(I_O_BUFFER,2,E_TMSG,56);
                    E_TPOPP("r_rdcr")
                    return;
           case  9: e_trap(OVERFLOW,2,E_TMSG,56);
                    E_TPOPP("r_rdcr")
                    return;
           case 10: e_trap(UNDERFLOW,2,E_TMSG,56);
                    E_TPOPP("r_rdcr")
                    return;
           default: ;
           }

        ((char *)b_cp__)[length] = '\0';

        if (expo<DBL_MIN_10_EXP)
           e_trap(UNDERFLOW,2,E_TMSG,56);
        else if (expo>DBL_MAX_10_EXP-1)
           e_trap(OVERFLOW,2,E_TMSG,56);
        else if (sscanf((char *)b_cp__,"%le",s)==EOF)
           e_trap(I_O_ERROR,2,E_TMSG,20);

        E_TPOPP("r_rdcr")
        return;
        }





