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

/* CVS $Id: s_suba.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_suba.c                              */
/*                                                              */
/*      Entries         : s_trng s_suba(s,pos,end)              */
/*                        s_trng s;                             */
/*                        a_intg pos,end;                       */
/*                                                              */
/*      Arguments       : s = string                            */
/*                        pos = position                        */
/*                        end = end position                    */
/*                                                              */
/*      Return value    : substring                             */
/*                                                              */
/*      Description     : Substring of string beginning at      */
/*                        position pos up to an including       */
/*                        position end                          */
/*                                                              */
/*                   memcpy call                                */
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
local s_trng s_suba(s_trng s,a_intg pos,a_intg end)
#else
local s_trng s_suba(s,pos,end)

s_trng s;
a_intg pos;
a_intg end;
#endif
        {
        s_trng res;
        size_t n;

        if (pos<=s.clen && end>=pos)
           {
           if (pos<1) pos = 1;
           n = (size_t) (((end<=s.clen) ? end : s.clen)-pos+1);
           s_init(&res,n);
           if (res.ptr!=NULL) memcpy((void*) res.ptr, (void*) (s.ptr+pos-1),n);
           res.clen = n;
           }
        else
           s_init(&res,(size_t)0);

        res.tmp = TRUE;

        if (s.tmp) s_free(&s);

#ifdef HEAP_CHECK
b_tmph((a_char *)&res.ptr);
#endif

        return(res);
        }





