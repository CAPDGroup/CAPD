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

/* CVS $Id: i_ep10.c,v 1.21 2014/01/30 17:24:08 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_ep10.c                10**e         */
/*                                                              */
/*      Entries         : a_intv i_ep10(e)                      */
/*                        a_intv e;                             */
/*                                                              */
/*      Arguments                                               */
/*                      : e = interval argument exponent        */
/*                                                              */
/*      Description     : Interval power function 10**e         */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_ten_;
#endif

#ifdef LINT_ARGS
local a_intv i_ep10( a_intv e )
#else
local a_intv i_ep10( e )
a_intv e;
#endif

{
        a_btyp        rc;
        a_real            dummy;
        a_intv          res;

        E_SPUSH("i_ep10")

        if (0==i_iv(e) || 0==i_iv(e)) rc = 1; /* error no interval */

        /* --- 2 types --- */
        else if ( i_point(e) ) {
            rc =  i_inv2(Lpower,*r_ten_,e.INF,&res.INF,&res.SUP);
        }
        else if ( i_iv(e) ) {
            rc =   i_inv2(Lpower,*r_ten_,e.INF,&res.INF,&dummy);
            rc +=  i_inv2(Lpower,*r_ten_,e.SUP,&dummy,&res.SUP);
        }
        else rc = 1;

        /* --- error --- */
        if (rc)
            e_trap(INV_ARG,4,E_TDBL+E_TEXT(5),&e.INF,
                             E_TDBL+E_TEXT(6),&e.SUP);

        /* --- exit --- */
        E_SPOPP("i_ep10")
        return(res);
}






