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

/* CVS $Id: d_vlcp.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_vlcp.c                              */
/*                                                              */
/*      Entry           : void d_vlcp(d)                        */
/*                        d_otpr *d;                            */
/*                                                              */
/*      Arguments       : d - dotprecision variable             */
/*                                                              */
/*      Description     : copy value-parameter if temporary-flag*/
/*                        is FALSE, otherwise set temporary-flag*/
/*                        FALSE                                 */
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
local void d_vlcp(d_otpr *d)
#else
local void d_vlcp(d)

d_otpr *d;
#endif
        {
        char *p;
        size_t n;

        E_TPUSH("d_vlcp")

        if (((*d)[A_STATUS] & A_TEMPORARY)==FALSE)
           {
           if ((p = (char *)malloc(n = A_LENGTH*sizeof(a_btyp)))==NULL)
              {
              e_trap(ALLOCATION,2,E_TMSG,40);

              E_TPOPP("d_vlcp")
              return;
              }
#ifdef HEAP_CHECK
b_freh((a_char *)d,(a_char *)p,(a_char *)"d_vlcp");
#endif
           C_COPY(p,(char *)(*d),n)
           *d = (d_otpr)p;
           }
        else
           (*d)[A_STATUS] &= ~A_TEMPORARY;

        E_TPOPP("d_vlcp")
        }





