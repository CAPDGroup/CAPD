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

/* CVS $Id: r_read.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_read.c                              */
/*                                                              */
/*      Entries         : void r_read(device,s,rnd,i)           */
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
/*                        device and convert to IEEE double     */
/*                        format                                */
/*                                                              */
/*                   remove e_trap() in case 5:                 */
/*                   b_scan argument changed                    */
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
local void r_read(FILE *device,a_real *s,a_intg rnd,int i)
#else
local void r_read(device,s,rnd,i)

FILE *device;
a_real *s;
a_intg rnd;
int i;
#endif
        {
        a_intg dp,expo,length,size;
        a_bool sign;
        a_btyp rc;
        char *buffer;

        E_TPUSH("r_read")

        buffer = (char *)b_cp__;
        size = BUFFERSIZE;

        switch (b_scan(device,&buffer,&size,&expo,&dp,&length,&sign,i))
           {
           case  1: e_trap(ALLOCATION,2,E_TMSG,56);
                    break;
           case  2:
           case  3:
           case  4: e_trap(I_O_ERROR,2,E_TMSG,58);
                    break;
           case  5: /* e_trap(NO_ERROR+E_EMSG+E_ECNT,2,E_TMSG,64); */
           default:
              if (expo>DBL_MAX_10_EXP+2)
                 expo = DBL_MAX_10_EXP+2;
              else if (expo<DBL_MIN_10_EXP-DBL_MANT_DIG-2)
                 expo = DBL_MIN_10_EXP-DBL_MANT_DIG-2;
              if ((rc = b_form((a_btyp *)buffer,&size,
                               expo,dp,length,sign,rnd,s))!=0)
                 e_trap(rc,6,E_TMSG,56,E_TDBL+E_TEXT(3),s,
                                       E_TDBL+E_TRES,s);
           }

        if (size!=BUFFERSIZE)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&buffer,(a_char *)buffer,(a_char *)"r_read");
#endif
           free(buffer);
           }

        E_TPOPP("r_read")
        }





