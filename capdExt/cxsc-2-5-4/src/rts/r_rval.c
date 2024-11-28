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

/* CVS $Id: r_rval.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_rval.c                              */
/*                                                              */
/*      Entries         : a_real r_rval(s,rnd)                  */
/*                        s_trng s;                             */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : s   - input string                    */
/*                        rnd - rounding mode                   */
/*                              -1 = rounding downwards         */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
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
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_real r_rval(s_trng s,a_intg rnd)
#else
local a_real r_rval(s,rnd)

s_trng s;
a_intg rnd;
#endif
        {
        a_real res;
        a_char *next;
        char ch;

        E_TPUSH("r_rval")

        if (s.clen<=0)
           {
           R_ASSIGN(res,*r_zero);
           e_trap(I_O_ERROR,2,E_TMSG,58);
           }
        else if (s.clen>=BUFFERSIZE)
           {
           R_ASSIGN(res,*r_zero);
           e_trap(I_O_BUFFER,2,E_TMSG,56);
           }
        else
           {
           if (s.suba) s_asgn(&s,s);
           ch = s.ptr[s.clen];
           s.ptr[s.clen] = '\0';

           r_conv((a_char *)s.ptr,&res,rnd,&next);

           s.ptr[s.clen] = ch;
           }

        if (s.tmp) s_free(&s);

        E_TPOPP("r_rval")
        return(res);
        }





