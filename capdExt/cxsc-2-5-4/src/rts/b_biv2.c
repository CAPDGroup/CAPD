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

/* CVS $Id: b_biv2.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_biv2.c                              */
/*                                                              */
/*      Entries         : a_btyp b_biv2(func,a1,a2,rlb,rub)     */
/*                        int (*func)();                        */
/*                        a_real a1, a2, *rlb, *rub             */
/*                                                              */
/*      Arguments       : func = function to be evaluated for   */
/*                               given real arguments a1,a2     */
/*                        a1,a2= arguments of function          */
/*                        rlb  = lower bound of the result      */
/*                        rub  = upper bound of the result      */
/*                                                              */
/*      Function value  : function dependent return code.       */
/*                                                              */
/*      Description     : evaluate bounds for a function        */
/*                        of two real arguments                 */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_btyp b_maxl;
#endif

#ifdef LINT_ARGS
local a_btyp b_biv2(int (*func)(multiprecision,
                                  multiprecision,multiprecision),
                      a_real a1, a_real a2, a_real *rlb, a_real *rub)
#else
local a_btyp b_biv2(func, a1, a2, rlb, rub)

int (*func)();
a_real a1;
a_real a2;
a_real *rlb;
a_real *rub;
#endif
        {
        a_intg code;
        multiprecision l_a1,l_a2,l_res;
        a_btyp rc;
        a_btyp Oldb_maxl;
        a_intg down = -1, up = 1;

        l_init(&l_a1);
        l_init(&l_a2);
        l_init(&l_res);

        if (b_rtol(a1,&l_a1,down)) return(ALLOCATION);
        if (b_rtol(a2,&l_a2,up)) return(ALLOCATION);

        Oldb_maxl = b_maxl;
        b_maxl = 3;

        code = (*func)(l_a1,l_a2,l_res);

        b_maxl = Oldb_maxl;

        rc = b_ltor(l_res,rlb,down);
        rc += b_ltor(l_res,rub,up);

        l_free(&l_a1);
        l_free(&l_a2);
        l_free(&l_res);

        return((code)?(a_btyp)code:rc);
        }





