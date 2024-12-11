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

/* CVS $Id: s_sub.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_sub.c                               */
/*                                                              */
/*      Entries         : a_VOID s_sub(res,s,t)                 */
/*                        s_etof s,t,res;                       */
/*                                                              */
/*      Arguments       : s = first set                         */
/*                        t = second set                        */
/*                        res = result set                      */
/*                                                              */
/*      Function value  : pointer to result set                 */
/*                                                              */
/*      Description     : res = difference of two sets s and t. */
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
local a_VOID s_sub(s_etof res,s_etof s,s_etof t)
#else
local a_VOID s_sub(res,s,t)

s_etof res;
s_etof s;
s_etof t;
#endif
        {
        size_t i;

        E_TPUSH("s_sub")

        for (i=0;i<s_SIZE;i++) res[i] = s[i] & (~t[i]);

        E_TPOPP("s_sub")
        return((a_VOID)&res[0]);
        }





