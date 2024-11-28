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

/* CVS $Id: d_stan.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_stan.c                              */
/*                                                              */
/*      Entries         : a_real d_stan(a)                      */
/*                        d_otpr a;                             */
/*                                                              */
/*      Arguments       : a = dotprecision variable             */
/*                                                              */
/*      Description     : Round dotprecision variable to        */
/*                        nearest IEEE double value according   */
/*                        to IEEE rounding conventions          */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_sero;
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_real d_stan(d_otpr a)
#else
local a_real d_stan(a)

d_otpr a;
#endif
        {
        a_btyp rc,lang[BSIZE];
        a_intg expo;
        a_bool vz;
        a_real res;

        E_TPUSH("d_stan")

        if (a[A_STATUS] & A_QUIETNAN)
           {
           ((a_btyp *)&res)[B_HPART] = SET_QUIET(EXPO_MASK);
           ((a_btyp *)&res)[B_LPART] = a[A_LOWNAN];
           }
        else if (a[A_STATUS] & A_PINFINITY)
           {
           ((a_btyp *)&res)[B_HPART] = EXPO_MASK;
           ((a_btyp *)&res)[B_LPART] = ZERO;
           }
        else if (a[A_STATUS] & A_MINFINITY)
           {
           ((a_btyp *)&res)[B_HPART] = EXPO_MASK | MSB;
           ((a_btyp *)&res)[B_LPART] = ZERO;
           }

        /* get accu value in IEEE length                        */
        else if (b_geta( a, lang, &expo, &vz ))
           {
           R_ASSIGN(res,
              ((a[A_STATUS] & A_MZERO)==ZERO) ? *r_sero : *r_zero);
           }
        else
           {

           /* adjust denormalized number                        */
           if ( (rc = b_adj( lang, &expo )) ==NO_ERROR)
              rc = b_rndn( lang, &expo );
           else
              (void)b_rndn( lang, &expo );

           /* composition                                       */
           b_comp( &res, expo, lang, vz );

           if (rc) e_trap(rc+E_IEEE,2,E_TDTP,&a);
           }

        /* free temporary dotprecison variables                 */
        if (a[A_STATUS] & A_TEMPORARY) d_free(&a);

        E_TPOPP("d_stan")
        return(res);
        }





