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

/* CVS $Id: b_inv1.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_inv1.c                              */
/*                                                              */
/*      Entries         : a_btyp b_inv1(func,arg,res,rnd)       */
/*                        int (*func)();                        */
/*                        a_real arg,*res;                      */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : func = function to be evaluated for   */
/*                               a given real argument          */
/*                        arg  = argument of function           */
/*                        res  = result of function             */
/*                        rnd  = rounding                       */
/*                               <0 : downwards                 */
/*                               =0 : nearest                   */
/*                               >0 : upwards                   */
/*                                                              */
/*      Function value  : function dependent return code.       */
/*                                                              */
/*      Description     : evaluate a function for a given       */
/*                        real argument.                        */
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
local a_btyp b_inv1(int (*func)(multiprecision,multiprecision),
                      a_real arg,a_real *res,a_intg rnd)
#else
local a_btyp b_inv1(func,arg,res,rnd)

int (*func)();
a_real arg;
a_real *res;
a_intg rnd;
#endif
        {
        a_intg code;
        multiprecision l_arg,l_res;
        a_btyp rc;
        a_btyp Oldb_maxl;

        l_init(&l_arg);
        l_init(&l_res);

        if (b_rtol(arg,&l_arg,(a_intg)0)) return(ALLOCATION);

        Oldb_maxl = b_maxl;
        b_maxl = 3;

        code = (*func)(l_arg,l_res);

        b_maxl = Oldb_maxl;

        rc = b_ltor(l_res,res,rnd);

        l_free(&l_arg);
        l_free(&l_res);

        return((code)?(a_btyp)code:rc);
        }





