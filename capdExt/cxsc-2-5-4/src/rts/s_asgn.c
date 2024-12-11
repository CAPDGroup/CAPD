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

/* CVS $Id: s_asgn.c,v 1.21 2014/01/30 17:24:13 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_asgn.c                              */
/*                                                              */
/*      Entry           : void s_asgn(d,s)                      */
/*                        s_trng *d;                            */
/*                        s_trng s;                             */
/*                                                              */
/*      Arguments       : d - descriptor of destination string  */
/*                        s - descriptor of source string       */
/*                                                              */
/*      Description     : Assign contents of a dynamic string to*/
/*                        another dynamic string                */
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
local void s_asgn(s_trng *d,s_trng s)
#else
local void s_asgn(d,s)

s_trng *d;
s_trng s;
#endif
        {
        char *ptr;

        E_TPUSH("s_asgn")

        /* destination is smaller than source */
        if (d->alen<s.clen || d->suba)
           {

           /* destination is dynamic */
           if (d->fix==FALSE)
              {
              s_free(d);
              d->alen = d->clen = s.clen;

              /* source is temporary */
              if (s.tmp==TRUE)
                 {
                 d->alen = s.alen;
#ifdef HEAP_CHECK
b_varh((a_char *)&d->ptr,(a_char *)s.ptr);
#endif
                 d->ptr = s.ptr;
                 s.ptr = NULL;
                 s.clen = s.alen = 0;
                 s.tmp = FALSE;
                 }

              /* new allocation required */
              else
                 {

                 /* one more allocated for direct use of '\0'   */
                 if ((ptr=(char*) malloc(s.clen+1))==NULL)
                    e_trap(ALLOCATION,2,E_TMSG,54);
                 else
                    {
#ifdef HEAP_CHECK
b_geth((a_char *)&d->ptr,(a_char *)ptr,(a_char *)"s_asgn");
#endif
                    memcpy(ptr,s.ptr,s.clen);
                    d->ptr = ptr;
                    if (d->suba)
                       {
                       d->suba = FALSE;
                       d->tmp = TRUE;
                       }
                    }
                 }
              }

           /* destination is not dynamic == chopping */
           else
              {
              memcpy(d->ptr,s.ptr,d->alen);
              d->clen = d->alen;
              if (s.tmp) s_free(&s);
              }
           }

        /* destination large enough to hold source */
        else
           {
           d->clen = s.clen;
           memcpy(d->ptr,s.ptr,s.clen);
           if (s.tmp) s_free(&s);
           }

        E_TPOPP("s_asgn")
        return;
        }





