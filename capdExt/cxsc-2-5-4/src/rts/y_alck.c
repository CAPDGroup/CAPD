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

/* CVS $Id: y_alck.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_alck.c                              */
/*                                                              */
/*      Entry           : y_dscp y_alck(d)                      */
/*                        y_dscp d;                             */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                                                              */
/*      Return value    : d (input value) null if not allocated */
/*                                                              */
/*      Description     : check for allocation, engage exception*/
/*                        if array is not allocated.            */
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
local y_dscp y_alck(y_dscp d)
#else
local y_dscp y_alck(d)
y_dscp d;
#endif

{ 
        E_TPUSH("y_alck")
        if (!y_aled(d)) 
	   {
		 e_trap(ALLOCATION,2,E_EMSG,42);
		 d = (y_dscp) NULL;
        }
        E_TPOPP("y_alck")
        return d;
}





