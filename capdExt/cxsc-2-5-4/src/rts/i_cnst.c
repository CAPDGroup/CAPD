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

/* CVS $Id: i_cnst.c,v 1.21 2014/01/30 17:24:08 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_cnst.c                              */
/*                                                              */
/*      Entries         : a_intv i_cnst(str)                    */
/*                        a_char *str;                          */
/*                                                              */
/*      Arguments       : str - input string                    */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to IEEE double format rounded         */
/*                        to interval.                          */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_zero;
extern dotprecision b_cp__;
#endif

#ifdef LINT_ARGS
local a_intv i_cnst(a_char *str)
#else
local a_intv i_cnst(str)

a_char *str;
#endif
        {
        a_intv s;
        a_char *next;
        a_intg dp,expo,length,size;
        a_bool sign;
        char *buffer;
        a_btyp rc;

        E_TPUSH("i_cnst")

        R_ASSIGN(s.INF,*r_zero);
        R_ASSIGN(s.SUP,*r_zero);

        buffer = (char *)b_cp__;
        size = BUFFERSIZE;

        switch (b_chck(str,&buffer,&size,&expo,&dp,&length,&sign,&next))
           {
           case  1: e_trap(ALLOCATION,2,E_TMSG,66);
                    E_TPOPP("i_cnst")
                    return(s);
           case  2:
           case  3:
           case  4: e_trap(I_O_ERROR,4,E_TMSG,58,E_TSTR+E_TEXT(10),str);
                    E_TPOPP("i_cnst")
                    return(s);
           case  5: e_trap(NO_ERROR+E_EMSG+E_ECNT,2,E_TSTR+E_TEXT(10),str);
           default: ;
           }

        if ((rc = b_form((a_btyp *)buffer,&size,
                         expo,dp,length,sign,(a_intg)-1,&s.INF))!=0)
           e_trap(rc,2,E_TSTR+E_TEXT(10),str);

        else if ((rc = b_form((a_btyp *)buffer,&size,
                         expo,dp,length,sign,(a_intg)1,&s.SUP))!=0)
           e_trap(rc,2,E_TSTR+E_TEXT(10),str);

        if (size!=BUFFERSIZE)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&buffer,(a_char *)buffer,(a_char *)"i_cnst");
#endif
           free(buffer);
           }

        E_TPOPP("i_cnst")
        return(s);
        }





