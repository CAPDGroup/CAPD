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

/* CVS $Id: b_inv2.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_inv2.c                              */
/*                                                              */
/*      Entries         : a_btyp b_inv2(func,a1,a2,res,rnd)     */
/*                        int (*func)();                        */
/*                        a_real a1,a2,*res;                    */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : func = function to be evaluated for   */
/*                               given real arguments a1,a2     */
/*                        a1,a2= arguments of function          */
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
local a_btyp b_inv2(int (*func)(multiprecision,
                                  multiprecision,multiprecision),
                      a_real a1,a_real a2,a_real *res,a_intg rnd)
#else
local a_btyp b_inv2(func,a1,a2,res,rnd)

int (*func)();
a_real a1;
a_real a2;
a_real *res;
a_intg rnd;
#endif
        {
        a_intg code;
        multiprecision l_a1,l_a2,l_res;
        a_btyp rc;
        a_btyp Oldb_maxl;

        l_init(&l_a1);
        l_init(&l_a2);
        l_init(&l_res);

        if (b_rtol(a1,&l_a1,(a_intg)0)) return(ALLOCATION);
        if (b_rtol(a2,&l_a2,(a_intg)0)) return(ALLOCATION);

        Oldb_maxl = b_maxl;
        b_maxl = 3;

        code = (*func)(l_a1,l_a2,l_res);

        b_maxl = Oldb_maxl;

        rc = b_ltor(l_res,res,rnd);

        l_free(&l_a1);
        l_free(&l_a2);
        l_free(&l_res);

        return((code)?(a_btyp)code:rc);
        }





