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

/* CVS $Id: a_ixch.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_ixch.c                              */
/*                                                              */
/*      Entries         : a_btyp a_ixch(i,l,u)                  */
/*                        a_intg i,l,u;                         */
/*                                                              */
/*      Arguments       : i = index                             */
/*                        l = lower bound of index range        */
/*                        u = upper bound of index range        */
/*                                                              */
/*      Description     : Index check                           */
/*                                                              */
/*      Function value  : Offset of index i from lower bound l. */
/*                                                              */
/*      Note            : -1 returned in case of error          */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local a_btyp a_ixch(a_intg i,a_intg l,a_intg u)
#else
local a_btyp a_ixch(i,l,u)

a_intg i;
a_intg l;
a_intg u;
#endif
        {
        a_btyp res;

        E_TPUSH("a_ixch")

        if (l<=i && i<=u)
           res = (a_btyp)(i-l);
        else
           {
           e_trap(INDEX_RANGE,6,E_TINT+E_TEXT(4),&i,
                                E_TINT+E_TEXT(5),&l,
                                E_TINT+E_TEXT(6),&u);
           res = (a_btyp) -1;
           }

        E_TPOPP("a_ixch")
        return(res);
        }





