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

/* CVS $Id: y_free.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_free.c                              */
/*                                                              */
/*      Entry           : void y_free(d)                        */
/*                        y_dscp d;                             */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                                                              */
/*      Description     : free dynamic array                    */
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
local void y_free(y_dscp d)
#else
local void y_free(d)

y_dscp d;
#endif
        {
        E_TPUSH("y_free")

        if ((char *)((y_desc *)d)->array!=NULL && NOT(((y_desc *)d)->subarr))
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&((y_desc *)d)->array,
(a_char *)((y_desc *)d)->array,(a_char *)"y_free");
#endif
           B_FREE(((y_desc *)d)->array)
           ((y_desc *)d)->array = NULL;
           }

#ifdef HEAP_CHECK
else if (NOT(((y_desc *)d)->subarr))
b_freh((a_char *)&((y_desc *)d)->array,
(a_char *)((y_desc *)d)->array,(a_char *)"y_free");
#endif

        E_TPOPP("y_free")
        return;
        }





