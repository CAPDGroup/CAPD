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

/* CVS $Id: r_cnst.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_cnst.c                              */
/*                                                              */
/*      Entries         : a_real r_cnst(str)                    */
/*                        char *str;                            */
/*                                                              */
/*      Arguments       : str - input string                    */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to IEEE double format rounded         */
/*                        to nearest.                           */
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
extern a_bool f_ppcc;
extern a_bool e_ofie;
#endif

#ifdef LINT_ARGS
local a_real r_cnst(char *str)
#else
local a_real r_cnst(str)

char *str;
#endif
        {
        a_real s;
        a_char *next;
        a_bool ie_flag;

        E_TPUSH("r_cnst")

        R_ASSIGN(s,*r_zero);

        ie_flag = e_ofie;
        e_ofie = FALSE;
        r_conv((a_char*) str,&s,(a_intg)0,&next);
        if (e_ofie)
           {
           if (f_ppcc)
              e_trap(NO_ERROR+E_EMSG+E_EARG+E_ECNT,6,E_TMSG,68,E_TMSG,59,
                     E_TSTR+E_TEXT(10),str);
           e_ofie = TRUE;
           }
        else if (ie_flag) e_ofie = TRUE;

        E_TPOPP("r_cnst")
        return(s);
        }





