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

/* CVS $Id: d_free.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_free.c                              */
/*                                                              */
/*      Entries         : void d_free(a)                        */
/*                        dotprecision *a;                      */
/*                                                              */
/*      Arguments       : a = dotprecision variable             */
/*                                                              */
/*      Description     : Free allocated dotprecision value     */
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
local void d_free(dotprecision *a)
#else
local void d_free(a)

dotprecision *a;
#endif
        {
        if ((*a)!=NO_DOTPRECISION)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)a,(a_char *)(*a),(a_char *)"d_free");
#endif
           B_FREE(*a)
           *a = NO_DOTPRECISION;
           }
#ifdef HEAP_CHECK
else
b_freh((a_char *)a,(a_char *)(*a),(a_char *)"d_free");
#endif

        return;
        }





