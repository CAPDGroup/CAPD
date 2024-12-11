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

/* CVS $Id: l_exct.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */


/****************************************************************/
/*                                                              */
/*      Filename        : l_exct.c                              */
/*                                                              */
/*      Entries         : void l_exct(a,b,r,l)                  */
/*                        multiprecision *a;                    */
/*                        multiprecision b;                     */
/*                        a_intg *r;                            */
/*                        a_intg *l;                            */
/*                                                              */
/*      Arguments       : a = multiprecision variable           */
/*                        b = multiprecision value              */
/*                        r = rounding bits                     */
/*                        l = length of mantissa                */
/*                                                              */
/*      Description     : Assign exact multiprecision value     */
/*                        without considering Maxl and without  */
/*                        changing rounding flags or mantissa   */
/*                        length.                               */
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
local void l_exct(multiprecision *a,multiprecision b,a_intg *r,a_intg *l)
#else
local void l_exct(a,b,r,l)

multiprecision *a;
multiprecision b;
a_intg *r;
a_intg *l;
#endif
        {
        E_TPUSH("l_exct")

        *l = (b->z) ? -((a_intg)MAXINT) : (a_intg)b->l;
        *r = b->r;

        if (*a==b) { E_TPOPP("l_exct") return; }

        if ((*a)->l!=0)
           {
           (*a)->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&(*a)->m,(a_char *)(*a)->m,(a_char *)"l_exct");
#endif
           B_FREE((*a)->m)
           }

        if (((*a)->z = b->z)==0)
           {
           if (b_ball(b->l,&(*a)->m)) 
              e_trap(ALLOCATION,2,E_TMSG,65);
           else
              {
              (*a)->e = b->e;
              (*a)->s = b->s;
              (*a)->l = b->l;
              B_COPY((*a)->m,b->m,b->l)
              }
           }

        (*a)->r = b->r;
 
        if (b->f) l_free(&b);

        E_TPOPP("l_exct")
        return;
        }





