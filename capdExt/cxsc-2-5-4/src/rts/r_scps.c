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

/* CVS $Id: r_scps.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_scps.c                              */
/*                                                              */
/*      Entries         : a_real r_scps(r,s,n,rnd)              */
/*                        a_real r[];                           */
/*                        a_real s[];                           */
/*                        a_intg n;                             */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : n   = size of arrays r and s          */
/*                        r   = real vector                     */
/*                        s   = real vector                     */
/*                        rnd = rounding mode                   */
/*                              +1 rounding upwards             */
/*                               0 round to nearest             */
/*                              -1 rounding downwards           */
/*                          addition of 4 does not clear accu   */
/*                                                              */
/*      Function value  : Rounded real dotproduct               */
/*                                                              */
/*      Description     : Dotproduct for static real arrays     */
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
#endif

#ifdef LINT_ARGS
local a_real r_scps(a_real r[],a_real s[],a_intg n,a_intg rnd)
#else
local a_real r_scps(r,s,n,rnd)

a_real r[];
a_real s[];
a_intg n;
a_intg rnd;
#endif
        {
        a_intg i;
        a_real res;

        E_TPUSH("r_scps")

        /* force linkage of module containing global variables only */
#if VAX_VMS_C
        b_accu();
#endif

        if (rnd<3)
           {
           d_clr(&b_acrl);
           }

        for (i=0;i<n;i++)
           d_padd(&b_acrl,r[i],s[i]);

        if (rnd==0)     res = d_stan(b_acrl);
        else if (rnd>0) res = d_stau(b_acrl);
        else            res = d_stad(b_acrl);

        E_TPOPP("r_scps")
        return(res);
        }





