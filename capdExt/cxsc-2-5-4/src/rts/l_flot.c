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

/* CVS $Id: l_flot.c,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_flot.c                              */
/*                                                              */
/*      Entries         : multiprecision l_flot(n)              */
/*                        a_intg n;                             */
/*                                                              */
/*      Description     : Generate multiprecision value from    */
/*                        integer value.                        */
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
local multiprecision l_flot(a_intg n)
#else
local multiprecision l_flot(n)

a_intg n;
#endif
        {
        multiprecision res;

        E_TPUSH("l_flot")

        l_init(&res);

        if (res==NULL)
           e_trap(ALLOCATION,2,E_TMSG,65);
        else
           {
           res->f = 1;
           if (n==0)
              res->z = 1;
           else if (b_ball(1,&res->m))
              e_trap(ALLOCATION,2,E_TMSG,65);
           else
              {
              res->e = res->z = res->r = 0;
              res->l = 1;
              if (n>0)
                 {
                 res->m[0] = n;
                 res->s = 0;
                 }
              else
                 {
                 res->m[0] = (n==MININT) ? MSB : -n;
                 res->s = 1;
                 }
              }
           }

        E_TPOPP("l_flot")
        return(res);
        }





