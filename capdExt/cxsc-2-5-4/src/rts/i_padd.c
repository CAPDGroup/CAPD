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

/* CVS $Id: i_padd.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_padd.c                              */
/*                                                              */
/*      Entries         : void i_padd(cl,cu,a,b)                */
/*                        dotprecision *cl,*cu;                 */
/*                        a_intv a,b;                           */
/*                                                              */
/*      Arguments       : cl= dotprecision variable(lower bound)*/
/*                        cu= dotprecision variable(upper bound)*/
/*                        a = interval value                    */
/*                        b = interval value                    */
/*                                                              */
/*      Description     : Add product to dotprecision variable  */
/*                        c = c+a*b                             */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
/* --- replaced by r_pcmp() ---
extern a_real *r_zero;
*/
#endif

#ifdef LINT_ARGS
local void i_padd(dotprecision *cl,dotprecision *cu,a_intv a,a_intv b)
#else
local void i_padd(cl,cu,a,b)

dotprecision *cl;
dotprecision *cu;
a_intv a;
a_intv b;
#endif
        {
        E_TPUSH("i_padd")

        if (r_sign(a.INF)>=0)
           {
           if (r_sign(b.INF)>=0)
              {
              d_padd(cl,a.INF,b.INF);
              d_padd(cu,a.SUP,b.SUP);
              }
           else if (r_sign(b.SUP)<=0)
              {
              d_padd(cl,a.SUP,b.INF);
              d_padd(cu,a.INF,b.SUP);
              }
           else
              {
              d_padd(cl,a.SUP,b.INF);
              d_padd(cu,a.SUP,b.SUP);
              }
           }
        else if (r_sign(a.SUP)<=0)
           {
           if (r_sign(b.INF)>=0)
              {
              d_padd(cl,a.INF,b.SUP);
              d_padd(cu,a.SUP,b.INF);
              }
           else if (r_sign(b.SUP)<=0)
              {
              d_padd(cl,a.SUP,b.SUP);
              d_padd(cu,a.INF,b.INF);
              }
           else
              {
              d_padd(cl,a.INF,b.SUP);
              d_padd(cu,a.INF,b.INF);
              }
           }
        else
           {
           if (r_sign(b.INF)>=0)
              {
              d_padd(cl,a.INF,b.SUP);
              d_padd(cu,a.SUP,b.SUP);
              }
           else if (r_sign(b.SUP)<=0)
              {
              d_padd(cl,a.SUP,b.INF);
              d_padd(cu,a.INF,b.INF);
              }
           else
              {
/* --- replaced by r_pcmp() ---
              dotprecision q;
              d_init(&q);
              d_clr(&q);

              d_padd(&q,a.INF,b.SUP);
              d_psub(&q,a.SUP,b.INF);
              if (r_gt(d_stau(q),*r_zero))
*/
              if (r_pcmp(a.INF,b.SUP,a.SUP,b.INF)>0)
                 {
                 d_padd(cl,a.SUP,b.INF);
                 }
              else
                 {
                 d_padd(cl,a.INF,b.SUP);
                 }

/* --- replaced by r_pcmp() ---
              d_clr(&q);
              d_padd(&q,a.INF,b.INF);
              d_psub(&q,a.SUP,b.SUP);
              if (r_gt(d_stau(q),*r_zero))
*/
              if (r_pcmp(a.INF,b.INF,a.SUP,b.SUP)>0)
                 {
                 d_padd(cu,a.INF,b.INF);
                 }
              else
                 {
                 d_padd(cu,a.SUP,b.SUP);
                 }

/* --- replaced by r_pcmp() ---
              d_free(&q);
*/
              }
           }

        E_TPOPP("i_padd")
        return;
        }





