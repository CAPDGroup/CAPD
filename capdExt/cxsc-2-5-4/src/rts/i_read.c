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

/* CVS $Id: i_read.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_read.c                              */
/*                                                              */
/*      Entries         : void i_read(device,s,i)               */
/*                        FILE *device;                         */
/*                        a_intv *s;                            */
/*                        int i;                                */
/*                                                              */
/*      Arguments       : device - input device reading string  */
/*                        s - resultant interval                */
/*                        i - window character                  */
/*                                                              */
/*      Description     : Read a real value in form of a        */
/*                        character string from input device    */
/*                        convert to an interval                */
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
local void i_read(FILE *device,a_intv *s,int i)
#else
local void i_read(device,s,i)

FILE *device;
a_intv *s;
int i;
#endif
        {
        a_intg dp,expo,length,size;
        a_bool sign;
        a_btyp rc;
        char *buffer;

        E_TPUSH("i_read")

        buffer = (char *)b_cp__;
        size = BUFFERSIZE;

        switch (b_scan(device,&buffer,&size,&expo,&dp,&length,&sign,i))
           {
           case  1: e_trap(ALLOCATION,2,E_TMSG,66);
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
              if ((rc = b_ifrm((a_btyp *)buffer,expo,dp,length,sign,s))!=0)
                 e_trap(rc,6,E_TMSG,56,E_TDBL+E_TEXT(5),&s->INF,
                                       E_TDBL+E_TEXT(6),&s->SUP);
           }

        if (size!=BUFFERSIZE)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&buffer,(a_char *)buffer,(a_char *)"i_read");
#endif
           free(buffer);
           }

        E_TPOPP("i_read")
        }





