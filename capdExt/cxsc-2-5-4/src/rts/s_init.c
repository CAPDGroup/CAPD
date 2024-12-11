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

/* CVS $Id: s_init.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_init.c                              */
/*                                                              */
/*      Entry           : void s_init(d,size)                   */
/*                        s_trng *d;                            */
/*                        size_t size;                          */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                        size - size of string                 */
/*                                                              */
/*      Description     : initialization of a dynamic array     */
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
local void s_init(s_trng *d,size_t size)
#else
local void s_init(d,size)

s_trng *d;
size_t size;
#endif
        {
        E_TPUSH("s_init")

        d->alen = d->clen = 0;
        d->tmp = d->suba = d->fix = FALSE;
        if (size<=0)
           {
           size = 0;
           d->ptr = NULL;
           }

        /* one more allocated for direct use of '\0'    */
        else if ((d->ptr = (char*) malloc(size+1))==NULL)
           e_trap(ALLOCATION,2,E_TMSG,54);
        else
           {
           d->fix = TRUE;
           d->alen = size;
#ifdef HEAP_CHECK
b_geth((a_char *)&d->ptr,(a_char *)d->ptr,(a_char *)"s_init");
#endif
           }

        E_TPOPP("s_init")
        return;
        }





