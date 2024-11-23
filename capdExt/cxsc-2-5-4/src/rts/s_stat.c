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

/* CVS $Id: s_stat.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_stat.c                              */
/*                                                              */
/*      Entry           : s_trng s_stat(s,n)                    */
/*                        a_char s[];                           */
/*                        a_intg n;                             */
/*                                                              */
/*      Arguments       : s - static string                     */
/*                        n - length of static string           */
/*                                                              */
/*      Function value  : dynamic string descriptor             */
/*                                                              */
/*      Description     : copy static string to dynamic string  */
/*                        (descriptor is initialized)           */
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
local s_trng s_stat(a_char s[],a_intg n)
#else
local s_trng s_stat(s,n)

a_char s[];
a_intg n;
#endif
        {
        s_trng res;

        /* tmp- and fix-flag must not be TRUE since static array   */
        /* is referenced by pointer of dynamic string              */
        res.fix = res.tmp = FALSE;
        res.suba = TRUE;
        res.ptr = (char *) &s[0];
        res.alen = res.clen = (size_t) n;

        return(res);
        }





