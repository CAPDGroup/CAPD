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

/* CVS $Id: b_blgx.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_blgx.c                              */
/*                                                              */
/*      Entries         : a_btyp b_blgx(a,res)                  */
/*                        a_real a;                             */
/*                        a_real res;                           */
/*                                                              */
/*      Arguments       : a = real argument of function         */
/*                      : res = real result of function         */
/*                                                              */
/*      Description     : logarithm to base 10 (exact)          */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_ten_,*r_one_;
#endif

#ifdef LINT_ARGS
local a_btyp b_blgx( a_real a, a_real *res )
#else
local a_btyp b_blgx( a, res )
a_real a;
a_real *res;
#endif

{
        a_real          k;
        int             i;
        a_btyp          rc;

        E_TPUSH("b_blgx")

        rc= 0;
        R_ASSIGN(k,*r_one_);
        for (i=1; i<=23; i++) {
          R_ASSIGN(k,r_muld(k,*r_ten_));
          if (r_eq(k,a)) {
            R_ASSIGN(*res,r_flot((a_intg)i));
            rc= 1;
          }
        } /* end for */

        E_TPOPP("b_blgx")
        return(rc);  /*result not exactly representable */
        }





