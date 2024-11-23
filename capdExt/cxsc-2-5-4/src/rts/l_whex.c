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

/* CVS $Id: l_whex.c,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_whex.c                              */
/*                                                              */
/*      Entry           : void l_whex(desc,r,mode)              */
/*                        f_text *desc;                         */
/*                        multiprecision r;                     */
/*                        a_char mode;                          */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        r      - multiprecision value         */
/*                        mode   - mode                         */
/*                                 'x' = hexadecimal (small)    */
/*                                 'X' = hexadecimal (capital)  */
/*                                                              */
/*      Description     : write multiprecision value in special */
/*                        format                                */
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
local void l_whex(f_text *desc,multiprecision r,a_char mode)
#else
local void l_whex(desc,r,mode)

f_text *desc;
multiprecision r;
a_char mode;
#endif
        {
        int i,k;
        a_char *p;

        E_TPUSH("l_whex")

        if (b_text(desc,FALSE))
           {
           if (mode!='x' && mode!='X')
              e_trap(I_O_ERROR,4,E_TMSG,51,E_TCHR,&mode);
           else
              {
              f_putc((a_char)'|',desc);
              p = (a_char *)&r;
#if INTEL
              for (i=(int)sizeof(r)-1;i>=0;i--)
#else
              for (i=0;i<(int)sizeof(r);i++)
#endif
                 f_bhex(desc,p[i],mode);
              f_putc((a_char)'-',desc);
              f_putc((a_char)'>',desc);
              f_putc((a_char)'z',desc);
              f_putc((a_char)'=',desc);
              f_putc((r->z) ? (a_char)'1' : (a_char)'0',desc);
              f_putc((a_char)' ',desc);
              f_putc((a_char)'s',desc);
              f_putc((a_char)'=',desc);
              f_putc((r->s) ? (a_char)'1' : (a_char)'0',desc);
              f_putc((a_char)' ',desc);
              f_putc((a_char)'r',desc);
              f_putc((a_char)'=',desc);
              f_putc((r->r) ? (a_char)'1' : (a_char)'0',desc);
              f_putc((a_char)' ',desc);
              f_putc((a_char)'f',desc);
              f_putc((a_char)'=',desc);
              f_putc((r->f) ? (a_char)'1' : (a_char)'0',desc);
              f_putc((a_char)' ',desc);
              p = (a_char *)&(r->e);
              f_putc((a_char)'e',desc);
              f_putc((a_char)'=',desc);
#if INTEL
              for (i=(int)sizeof(r->e)-1;i>=0;i--)
#else
              for (i=0;i<(int)sizeof(r->e);i++)
#endif
                 {
                 f_bhex(desc,p[i],mode);
                 }
              f_putc((a_char)' ',desc);
              p = (a_char *)&(r->l);
              f_putc((a_char)'l',desc);
              f_putc((a_char)'=',desc);
#if INTEL
              for (i=(int)sizeof(r->l)-1;i>=0;i--)
#else
              for (i=0;i<(int)sizeof(r->l);i++)
#endif
                 {
                 f_bhex(desc,p[i],mode);
                 }
              f_putc((a_char)' ',desc);
              p = (a_char *)&(r->m);
              f_putc((a_char)'m',desc);
              f_putc((a_char)'=',desc);
#if INTEL
              for (i=(int)sizeof(r->m)-1;i>=0;i--)
#else
              for (i=0;i<(int)sizeof(r->m);i++)
#endif
                 {
                 f_bhex(desc,p[i],mode);
                 }
              f_putc((a_char)'-',desc);
              f_putc((a_char)'>',desc);
              p = (a_char *)(r->m);
              for (i=0;i<r->l;i++)
                 {
#if INTEL
                 for (k=(int)sizeof(a_btyp)-1;k>=0;k--)
#else
                 for (k=0;k<(int)sizeof(a_btyp);k++)
#endif
                    f_bhex(desc,p[i*sizeof(a_btyp)+k],mode);
                 if (i!=r->l-1) f_putc((a_char)' ',desc);
                 }
              f_putc((a_char)'|',desc);
              }
           }

        if (r->f) l_free(&r);

        E_TPOPP("l_whex")
        return;
        }





