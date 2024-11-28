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

/* CVS $Id: l_abs.c,v 1.21 2014/01/30 17:24:09 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_abs.c                               */
/*                                                              */
/*      Entries         : multiprecision l_abs(a)               */
/*                        multiprecision a;                     */
/*                                                              */
/*      Arguments       : a = mulitprecision value              */
/*                                                              */
/*      Description     : Return positive value of argument.    */
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
local multiprecision l_abs(multiprecision a)
#else
local multiprecision l_abs(a)

multiprecision a;
#endif
        {
        multiprecision res;

        E_TPUSH("l_abs")

        if (a->f && a->l<=b_maxl)
           {
           res = a;
           res->s = 0;
           }
        else
           {
           l_init(&res);
           if (res==NULL)
              e_trap(ALLOCATION,2,E_TMSG,65);
           else if (b_bcpy(a,res))
              e_trap(ALLOCATION,2,E_TMSG,65);
           else
              {
              res->s = 0;
              res->f = 1;
              }
           }

        E_TPOPP("l_abs")
        return(res);
        }





