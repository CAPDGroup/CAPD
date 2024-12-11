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

/* CVS $Id: i_loga.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_loga.c                              */
/*                                                              */
/*      Entries         : a_intv i_loga(x,base)                 */
/*                        a_intv x;                             */
/*                        a_real base;                          */
/*                                                              */
/*      Arguments       : base = base                           */
/*                      : x = interval argument of function     */
/*                                                              */
/*      Description     : Interval logarithm to base base       */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_two_,*r_ten_,*r_one_;
#endif

#ifdef LINT_ARGS
local a_intv i_loga( a_intv x, a_real base)
#else
local a_intv i_loga( x, base )
a_intv x;
a_real base;
#endif

{
        a_btyp        rc;
        a_intv          res;
        a_real            dummy;

        E_SPUSH("i_loga")

        rc= 0;

        /* printf("x in i_loga:  %lf   %lf %lf\n", x.INF, x.SUP, base); */

        if (r_eq(base,*r_two_)) {
          res = i_log2(x);
        }
        else if (r_eq(base,*r_ten_)) {
          res = i_log(x);
        }
        /* --- invalid argument --- */
        else if (r_sign(x.INF)<=0) {
          rc = 1;
        }
        /* --- invalid base --- */
        else if ( r_sign(base)<=0 || r_eq(base,*r_one_) ) {
          rc = 1;
        }
        /* --- point argument --- */
        else if ( i_point(x) ) {
            /* printf("point argument!\n"); */
            rc =  i_inv2(Lloga,x.INF,base,&res.INF,&res.SUP);
        }

        /* --- interval argument --- */
        else if ((rc = i_iv(x))!=0) {
            /* printf("interval argument!\n"); */
            rc =  i_inv2(Lloga,x.INF,base,&res.INF,&dummy);
            rc += i_inv2(Lloga,x.SUP,base,&dummy,&res.SUP);
        }
        /* --- invalid argument --- */
        else rc = 1;

        /* --- error --- */
        if (rc) {
            e_trap(INV_ARG,5,E_TDBL+E_TEXT(17),&base,
                             E_TDBL+E_TEXT(5),&x.INF,
                             E_TDBL+E_TEXT(6),&x.SUP);
        }

        /* --- exit --- */
        E_SPOPP("i_loga")
        return(res);
        }





