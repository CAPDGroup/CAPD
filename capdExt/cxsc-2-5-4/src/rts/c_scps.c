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

/* CVS $Id: c_scps.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : c_scps.c                              */
/*                                                              */
/*      Entries         : a_cmpx c_scps(r,s,n,rnd)              */
/*                        a_cmpx r[];                           */
/*                        a_cmpx s[];                           */
/*                        a_intg n;                             */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : n   = size of arrays r and s          */
/*                        r   = complex vector                  */
/*                        s   = complex vector                  */
/*                        rnd = rounding mode                   */
/*                              +1 rounding upwards             */
/*                               0 round to nearest             */
/*                              -1 rounding downwards           */
/*                          addition of 4 does not clear accu   */
/*                                                              */
/*      Function value  : Rounded complex dotproduct            */
/*                                                              */
/*      Description     : Dotproduct for static complex arrays  */
/*                        with directed rounding                */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern dotprecision b_acrl;
extern dotprecision b_acil;
#endif

#ifdef LINT_ARGS
local a_cmpx c_scps(a_cmpx r[],a_cmpx s[],a_intg n,a_intg rnd)
#else
local a_cmpx c_scps(r,s,n,rnd)

a_cmpx r[];
a_cmpx s[];
a_intg n;
a_intg rnd;
#endif
        {
        a_intg i;
        a_cmpx res;

        E_TPUSH("c_scps")

        /* force linkage of module containing global variables only */
#if VAX_VMS_C
        b_accu();
#endif

        if (rnd<3)
           {
           d_clr(&b_acrl);
           d_clr(&b_acil);
           }

        for (i=0;i<n;i++)
           c_padd(&b_acrl,&b_acil,r[i],s[i]);

        if (rnd==0)     res = c_stan(b_acrl,b_acil);
        else if (rnd>0) res = c_stau(b_acrl,b_acil);
        else            res = c_stad(b_acrl,b_acil);

        E_TPOPP("c_scps")
        return(res);
        }





