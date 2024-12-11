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

/* CVS $Id: y_ixch.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_ixch.c                              */
/*                                                              */
/*      Entry           : a_btyp y_ixch(index,range)            */
/*                        a_intg index;                         */
/*                        y_bnds range;                         */
/*                                                              */
/*      Arguments       : index  - index of element             */
/*                        range  - range of dimension           */
/*                                                              */
/*      Return value    : offset to indexed element             */
/*                                                              */
/*      Description     : check for valid index and determine   */
/*                        offset to indexed element.            */
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
local a_btyp y_ixch(a_intg index,y_bnds range)
#else
local a_btyp y_ixch(index,range)

a_intg index;
y_bnds range;
#endif
        {
        E_TPUSH("y_ixch")

        if (index<range.lbound || index>range.ubound)
           {
           e_trap(INDEX_RANGE,6,
                  E_TINT+E_TEXT(4),&index,
                  E_TINT+E_TEXT(5),&range.lbound,
                  E_TINT+E_TEXT(6),&range.ubound);
           }

        E_TPOPP("y_ixch")
        return((index-range.lbound)*range.stride);
        }





