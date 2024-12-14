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

/* CVS $Id: l_succ.c,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_succ.c                              */
/*                                                              */
/*      Entries         : multiprecision l_succ(a)              */
/*                        multiprecision a;                     */
/*                                                              */
/*      Description     : Successor of multiprecision value     */
/*                        within screen of mantissa length Maxl.*/
/*                        If a=0, +2**(-32*Maxl) is returned.   */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_btyp b_maxl;
#endif

#ifdef LINT_ARGS
local multiprecision l_succ(multiprecision a)
#else
local multiprecision l_succ(a)

multiprecision a;
#endif
        {
        multiprecision res;

        E_TPUSH("l_succ")

        l_init(&res);

        if (res==NULL)
           e_trap(ALLOCATION,2,E_TMSG,65);
        else if (a->z)
           {
           res->z = FALSE;
           res->s = FALSE;
           res->e = -b_maxl;
           res->l = 1;
           if (b_ball(1,&res->m))
              {
              e_trap(ALLOCATION,2,E_TMSG,65);
              res->z = TRUE;
              res->l = 0;
              }
           else
              res->m[0] = LSB;
           }
        else
           {
           if (b_ball(b_maxl,&res->m))
              {
              e_trap(ALLOCATION,2,E_TMSG,65);
              res->z = TRUE;
              res->l = 0;
              E_TPOPP("l_succ")
              return(res);
              }
           res->z = FALSE;
           res->s = a->s;
           res->l = b_maxl;
           res->e = a->e;
           if (a->l<b_maxl)
              B_COPY(res->m,a->m,a->l)
           else
              B_COPY(res->m,a->m,b_maxl)

           if (NOT(res->s))
              {
              if (b_bcad(b_maxl,res->m))
                 {
                 if (res->e==MAXINT)
                    e_trap(OVERFLOW,2,E_TMLT+E_TEXT(7),&a);
                 else
                    {
                    res->e++;
                    res->m[0] = LSB;
                    }
                 }
              }
           else if (b_test(a->l-b_maxl,&a->m[b_maxl]))
              {
              (void)b_bcsu(b_maxl,res->m);
              if (res->m[0]==ZERO)
                 {
                 if (res->e==-MAXINT)
                    e_trap(UNDERFLOW,2,E_TMLT+E_TEXT(7),&a);
                 else
                    {
                    res->e--;
                    res->m[0] = MAX_BASETYPE;
                    }
                 }
              }
           }

        if (a->f) l_free(&a);

        E_TPOPP("l_succ")
        return(res);
        }





