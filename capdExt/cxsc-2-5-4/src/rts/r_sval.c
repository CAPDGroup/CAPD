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

/* CVS $Id: r_sval.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_sval.c                              */
/*                                                              */
/*      Entries         : a_real r_sval(s,rnd,r)                */
/*                        s_trng s;                             */
/*                        a_intg rnd;                           */
/*                        s_trng *r;                            */
/*                                                              */
/*      Arguments       : s   - input string                    */
/*                        rnd - rounding mode                   */
/*                              -1 = rounding downwards         */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        r   - remainder string                */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to IEEE double format                 */
/*                                                              */
/*                   type of ch -> char                         */
/*                   b_cp__ removed                             */
/*                   declared length of remainder considered    */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_real r_sval(s_trng s,a_intg rnd,s_trng *r)
#else
local a_real r_sval(s,rnd,r)

s_trng s;
a_intg rnd;
s_trng *r;
#endif
        {
        a_real res;
        a_char *next;
        char ch;

        E_TPUSH("r_sval")

        if (s.clen<=0)
           {
           R_ASSIGN(res,*r_zero);
           e_trap(I_O_ERROR,2,E_TMSG,58);
           }
        else
           {
           if (s.suba) s_asgn(&s,s);
           ch = s.ptr[s.clen];
           s.ptr[s.clen] = '\0';

           r_conv((a_char *)s.ptr,&res,rnd,&next);

           if (r->alen<(r->clen = s.clen-((char *)next-s.ptr)))
              {
              if (r->fix==FALSE)
                 {
                 if (r->alen>0)
                    {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->ptr,(a_char *)r->ptr,(a_char *)"r_sval");
#endif
                    free(r->ptr);
                    r->alen = 0;
                    }

                 /* one more allocated for direct use of '\0'   */
                 if ((r->ptr=(char*) malloc(r->clen+1))==NULL)
                    {
                    e_trap(ALLOCATION,2,E_TMSG,54);
                    r->alen = r->clen = 0;
                    }
                 else
                    r->alen = r->clen;
#ifdef HEAP_CHECK
b_geth((a_char *)&r->ptr,(a_char *)r->ptr,(a_char *)"r_sval");
#endif
                 }
              else
                 r->clen = r->alen;
              }
           if (r->clen>0) memcpy(r->ptr,next,r->clen);
           s.ptr[s.clen] = ch;
           }

        if (s.tmp) s_free(&s);

        E_TPOPP("r_sval")
        return(res);
        }





