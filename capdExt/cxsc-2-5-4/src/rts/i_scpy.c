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

/* CVS $Id: i_scpy.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_scpy.c                              */
/*                                                              */
/*      Entries         : a_intv i_scpy(r,s,rnd)                */
/*                        y_desc *r;                            */
/*                        y_desc *s;                            */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : r   = dynamic interval vector         */
/*                        s   = dynamic interval vector         */
/*                        rnd = mode                            */
/*                              0  clears accus first           */
/*                              4  does not clear accus first   */
/*                                                              */
/*      Function value  : Real interval dotproduct              */
/*                                                              */
/*      Description     : Dotproduct for dynamic real interval  */
/*                        arrays.                               */
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
extern dotprecision b_acru;
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_intv i_scpy(y_desc *r,y_desc *s,a_intg rnd)
#else
local a_intv i_scpy(r,s,rnd)

y_desc *r;
y_desc *s;
a_intg rnd;
#endif
        {
        a_intg i,k;
        a_intv res;

        E_TPUSH("i_scpy")

        /* force linkage of module containing global variables only */
#if VAX_VMS_C
        b_accu();
#endif

        R_ASSIGN(res.INF,*r_zero);
        R_ASSIGN(res.SUP,*r_zero);

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
           if (rnd<3)
              {
              d_clr(&b_acrl);
              d_clr(&b_acru);
              }

           for (i=0;i<r->elnum;i++)
              /* possible pointer allignment problems              */
              i_padd(&b_acrl,&b_acru,
                  *(a_intv *)y_inxn(r,(r->fd[0]).lbound+i),
                  *(a_intv *)y_inxn(s,(s->fd[0]).lbound+i));

           res = i_ista(b_acrl,b_acru);
           }

        E_TPOPP("i_scpy")
        return(res);
        }





