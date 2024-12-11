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

/* CVS $Id: r_scpy.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_scpy.c                              */
/*                                                              */
/*      Entries         : a_real r_scpy(r,s,rnd)                */
/*                        y_desc *r;                            */
/*                        y_desc *s;                            */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : r   = descriptor of real array        */
/*                        s   = descriptor of real array        */
/*                        rnd = rounding mode                   */
/*                              +1 rounding upwards             */
/*                               0 round to nearest             */
/*                              -1 rounding downwards           */
/*                          addition of 4 does not clear accus  */
/*                                                              */
/*      Function value  : Rounded real dotproduct               */
/*                                                              */
/*      Description     : Dotproduct for dynamic real arrays    */
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
extern d_otpr b_acrl;
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_real r_scpy(y_desc *r,y_desc *s,a_intg rnd)
#else
local a_real r_scpy(r,s,rnd)

y_desc *r;
y_desc *s;
a_intg rnd;
#endif
        {
        a_intg i,k;
        a_real res;

        E_TPUSH("r_scpy")

        /* force linkage of module containing global variables only */
#if VAX_VMS_C
        b_accu();
#endif

        R_ASSIGN(res,*r_zero);

        if (r->numdim!=1 || s->numdim!=1)
           {
           i = r->numdim;
           k = s->numdim;
           e_trap(INV_ARG,6,E_TMSG,26,E_TINT+E_TEXT(11),&i,
                                      E_TINT+E_TEXT(11),&k);
           }
        else if (r->elnum!=s->elnum)
           {
           i = r->elnum;
           k = s->elnum;
           e_trap(INV_ARG,6,E_TMSG,27,E_TINT+E_TEXT(12),&i,
                                      E_TINT+E_TEXT(12),&k);
           }
        else
           {
           if (rnd<3) d_clr(&b_acrl);

           for (i=0;i<r->elnum;i++)
              /* possible pointer allignment problems              */
              d_padd(&b_acrl,
                     *(a_real *)y_inxn(r,(r->fd[0]).lbound+i),
                     *(a_real *)y_inxn(s,(s->fd[0]).lbound+i));

           if (rnd==0)     res = d_stan(b_acrl);
           else if (rnd>0) res = d_stau(b_acrl);
           else            res = d_stad(b_acrl);
           }

        E_TPOPP("r_scpy")
        return(res);
        }





