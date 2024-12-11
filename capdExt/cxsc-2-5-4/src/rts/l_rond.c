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

/* CVS $Id: l_rond.c,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_rond.c                              */
/*                                                              */
/*      Entries         : a_intg l_rond(a)                      */
/*                        multiprecision a;                     */
/*                                                              */
/*      Arguments       : a = multiprecision value              */
/*                                                              */
/*      Description     : Round INTERN value to next integer    */
/*                        number.                               */
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
local a_intg l_rond(multiprecision a)
#else
local a_intg l_rond(a)

multiprecision a;
#endif
        {
        a_intg res = 0;

        E_TPUSH("l_rond")

        if (a->z)
           {
           }
        else if (a->e>0 || (a->e==0 && (a->m[0] & MSB)))
           e_trap(OVERFLOW,4,E_TMSG,15,E_TMLT+E_TEXT(7),&a);
        else if (a->e<0)
           {
           if (a->e==-1 && (a->m[0] & MSB)) res = 1-2*a->s;
           }
        else
           {
           res = a->m[0];

           if (a->l>1 && (a->m[1] & MSB))
              {
              if (res>MAXINT)
                 {
                 e_trap(OVERFLOW,4,E_TMSG,15,E_TMLT+E_TEXT(7),&a);
                 res = 0;
                 }
              else
                 res++;
              }

           if (a->s)
              {
#if C_P_2
              res = -res;
#else
              res = ~res+1;
#endif
              }
           }

        if (a->f) l_free(&a);

        E_TPOPP("l_rond")
        return(res);
        }





