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

/* CVS $Id: d_chs.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_chs.c                               */
/*                                                              */
/*      Entries         : dotprecision d_chs(a)                 */
/*                        dotprecision a;                       */
/*                                                              */
/*      Arguments       : a = dotprecision value                */
/*                                                              */
/*      Function value  : negative dotprecision value           */
/*                                                              */
/*      Description     : Change sign of dotprecision value     */
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
local dotprecision d_chs(dotprecision a)
#else
local dotprecision d_chs(a)

dotprecision a;
#endif
        {
        dotprecision c;

        d_init(&c);

        d_ass(&c,a);

        if (c[A_BEGIN]) c[A_SIGN] = 1-c[A_SIGN];

        c[A_STATUS] |= A_TEMPORARY;
        return(c);
        }





