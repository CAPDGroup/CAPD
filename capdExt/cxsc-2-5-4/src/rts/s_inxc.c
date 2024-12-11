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

/* CVS $Id: s_inxc.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_inxc.c                              */
/*                                                              */
/*      Entry           : a_char *s_inxc(s,i)                   */
/*                        s_trng s;                             */
/*                        a_intg i;                             */
/*                                                              */
/*      Arguments       : s - string                            */
/*                        i - index of string component         */
/*                                                              */
/*      Description     : Determine address of string component */
/*                        with index checking.                  */
/*                                                              */
/*      Note            : NULL returned in case of error        */
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
local a_char *s_inxc(s_trng s,a_intg i)
#else
local a_char *s_inxc(s,i)

s_trng s;
a_intg i;
#endif
        {
        a_intg l,u;
        a_char *res;
        char *ptr;

        E_TPUSH("s_inxc")

        if (i<=0 || (i>s.alen && s.fix==1))
           {
           l = 1;
           u = s.alen;
           res = NULL;
           e_trap(INDEX_RANGE,6,E_TINT+E_TEXT(4),&i,
                                E_TINT+E_TEXT(5),&l,
                                E_TINT+E_TEXT(6),&u);
           }
        else
           {
           if (i>s.alen)
              {

              /* one more allocated for direct use of '\0'      */
              if ((ptr = (char*) malloc((size_t)(i+1)))==NULL)
                 {
                 res = NULL;
                 e_trap(ALLOCATION,2,E_TMSG,54);
                 }
              else
                 {
                 if (s.alen>0)
                    {
                    memcpy(ptr,s.ptr,s.alen);
#ifdef HEAP_CHECK
b_freh((a_char *)&s.ptr,(a_char *)s.ptr,(a_char *)"s_inxc");
#endif
                    B_FREE(s.ptr)
                    }
#ifdef HEAP_CHECK
b_geth((a_char *)&s.ptr,(a_char *)ptr,(a_char *)"s_inxc");
#endif
                 s.ptr = ptr;
                 s.alen = (size_t) i;
                 res = (a_char *)(s.ptr+(i-1));
                 }
              }
           else
              res = (a_char *)(s.ptr+(i-1));
           }

        if (s.tmp) s_free(&s);

        E_TPOPP("s_inxc")
        return(res);
        }





