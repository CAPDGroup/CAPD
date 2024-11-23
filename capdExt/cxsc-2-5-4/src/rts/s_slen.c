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

/* CVS $Id: s_slen.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_slen.c                              */
/*                                                              */
/*      Entry           : void s_slen(d,n)                      */
/*                        s_trng *d;                            */
/*                        a_intg n;                             */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic string      */
/*                        n - length of string                  */
/*                                                              */
/*      Description     : Set current length.                   */
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
local void s_slen(s_trng *d,a_intg n)
#else
local void s_slen(d,n)

s_trng *d;
a_intg n;
#endif
        {
        a_intg lb,ub;
        char *ptr;

        E_TPUSH("s_slen")

        if (n<0 || (n>d->alen && d->fix))
           {
           lb = 0;
           ub = d->alen;
           e_trap(INDEX_RANGE,6,E_TINT+E_TEXT(15),&n,
                  E_TINT+E_TEXT(5),&lb,
                  E_TINT+E_TEXT(6),&ub);
           }
        else if (n<=d->alen)
           d->clen = (size_t) n;
        else
           {

           /* one more allocated for direct use of '\0' */
           if ((ptr = (char*) malloc((size_t)(n+1)))==NULL)
              e_trap(ALLOCATION,2,E_TMSG,54);
           else
              {
              if (d->alen>0)
                 {
                 memcpy(ptr,d->ptr,d->alen);
#ifdef HEAP_CHECK
b_freh((a_char *)&d->ptr,(a_char *)d->ptr,(a_char *)"s_slen");
#endif
                 B_FREE(d->ptr)
                 }
#ifdef HEAP_CHECK
b_geth((a_char *)&d->ptr,(a_char *)ptr,(a_char *)"s_slen");
#endif
              d->ptr = ptr;
              d->clen = d->alen = (size_t) n;
              }
           }

        E_TPOPP("s_slen")
        return;
        }





