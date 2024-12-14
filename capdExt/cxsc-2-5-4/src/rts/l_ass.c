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

/* CVS $Id: l_ass.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_ass.c                               */
/*                                                              */
/*      Entries         : void l_ass(a,b)                       */
/*                        multiprecision *a;                    */
/*                        multiprecision b;                     */
/*                                                              */
/*      Arguments       : a = multiprecision variable           */
/*                        b = multiprecision value              */
/*                                                              */
/*      Description     : Assign multiprecision value.          */
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
local void l_ass(multiprecision *a,multiprecision b)
#else
local void l_ass(a,b)

multiprecision *a;
multiprecision b;
#endif
        {
        E_TPUSH("l_ass")

        if (b_bcpy(b,*a))
           e_trap(ALLOCATION,2,E_TMSG,65);

        if (b->f) l_free(&b);

        E_TPOPP("l_ass")
        return;
        }





