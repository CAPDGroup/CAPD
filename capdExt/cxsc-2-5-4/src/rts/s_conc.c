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

/* CVS $Id: s_conc.c,v 1.21 2014/01/30 17:24:13 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_conc.c                              */
/*                                                              */
/*      Entry           : s_trng s_conc(d,s)                    */
/*                        s_trng d;                             */
/*                        s_trng s;                             */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic string      */
/*                        s - descriptor of dynamic string      */
/*                                                              */
/*      Function value  : dynamic string descriptor             */
/*                                                              */
/*      Description     : concatenate dynamic strings           */
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
local s_trng s_conc(s_trng d,s_trng s)
#else
local s_trng s_conc(d,s)

s_trng d;
s_trng s;
#endif
        {
        s_trng res;
        size_t i,n;

        E_TPUSH("s_conc")

        n = d.clen+s.clen;
        if (d.tmp && n<=d.alen)
           {
           res = d;
           memcpy(res.ptr+d.clen,s.ptr,s.clen);
           if (s.tmp) s_free(&s);
           }
        else if (s.tmp && n<=s.alen)
           {
           res = s;
           for (i=s.clen;i>0;i--)
              res.ptr[d.clen+(i-1)] = res.ptr[i-1];
           memcpy(res.ptr,d.ptr,d.clen);
           if (d.tmp) s_free(&d);
           }
        else
           {
           res.tmp = TRUE;
           if ((res.ptr = (char *)malloc(n+1))==NULL)
              {
              e_trap(ALLOCATION,2,E_TMSG,54);
              res.clen = res.alen = 0;
              res.suba = res.fix = FALSE;
              res.ptr = NULL;
              E_TPOPP("s_conc")
              return(res);
              }
#ifdef HEAP_CHECK
b_geth((a_char *)&res.ptr,(a_char *)res.ptr,(a_char *)"s_conc");
#endif
           res.alen = n;
           memcpy(res.ptr,d.ptr,d.clen);
           memcpy(res.ptr+d.clen,s.ptr,s.clen);
           if (d.tmp) s_free(&d);
           if (s.tmp) s_free(&s);
           }

        res.suba = res.fix = FALSE;
        res.clen = n;

#ifdef HEAP_CHECK
b_tmph((a_char *)&res.ptr);
#endif

        E_TPOPP("s_conc")
        return(res);
        }





