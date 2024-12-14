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

/* CVS $Id: r_conv.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_conv.c                              */
/*                                                              */
/*      Entries         : void r_conv(str,s,rnd,next)           */
/*                        a_char *str;                          */
/*                        a_real *s;                            */
/*                        a_intg rnd;                           */
/*                        a_char **next;                        */
/*                                                              */
/*      Arguments       : str - input string                    */
/*                        s - resultant IEEE variable           */
/*                        rnd - rounding mode                   */
/*                              -1 = rounding downwards         */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        next - address of last character      */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to IEEE double format                 */
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
local void r_conv(a_char *str,a_real *s,a_intg rnd,a_char **next)
#else
local void r_conv(str,s,rnd,next)

a_char *str;
a_real *s;
a_intg rnd;
a_char **next;
#endif
        {
        a_intg dp,expo,length,size;
        a_bool sign;
        char *buffer;
        a_btyp rc;

        E_TPUSH("r_conv")

        buffer = (char *)b_cp__;
        size = BUFFERSIZE;

        switch (b_chck(str,&buffer,&size,&expo,&dp,&length,&sign,next))
           {
           case  1: e_trap(ALLOCATION,2,E_TMSG,56);
                    break;
           case  2:
           case  3:
           case  4: e_trap(I_O_ERROR,4,E_TMSG,58,E_TSTR+E_TEXT(10),str);
                    break;
           case  5: e_trap(NO_ERROR+E_EMSG+E_ECNT,2,E_TMSG,64);
           default: if ((rc = b_form((a_btyp *)buffer,&size,
                                     expo,dp,length,sign,rnd,s))!=0)
                       e_trap(rc,2,E_TSTR+E_TEXT(10),str);
           }

        if (size!=BUFFERSIZE)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&buffer,(a_char *)buffer,(a_char *)"r_conv");
#endif
           free(buffer);
           }


        E_TPOPP("r_conv")
        return;
        }





