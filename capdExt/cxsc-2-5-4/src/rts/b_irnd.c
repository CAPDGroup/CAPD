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

/* CVS $Id: b_irnd.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_irnd.c                              */
/*                                                              */
/*      Entries         : void b_irnd(value,lb,ub);             */
/*                        multiprecision value;                 */
/*                        multiprecision *lb;                   */
/*                        multiprecision *ub;                   */
/*                                                              */
/*      Arguments       : value - being rounded to an interval  */
/*                        lb    - lower bound of interval       */
/*                        ub    - upper bound of interval       */
/*                                                              */
/*      Description     : Determine multiprecision interval     */
/*                        that is formed by a multiprecision    */
/*                        value and its rounding information.   */
/*                                                              */
/*      Note            : zero returned if zero was input!      */
/*                        temporary value is freed!             */
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
local void b_irnd(multiprecision value,multiprecision *lb,multiprecision *ub)
#else
local void b_irnd(value,lb,ub)

multiprecision value;
multiprecision *lb;
multiprecision *ub;
#endif
        {
        int rc1,rc2;       /* !!! must be int !!! */
        multiprecision bb;

        E_TPUSH("b_irnd")

        /* determine bounds not considering the sign of the value */
        rc2 = b_brnd(value,*ub);
        if ((rc1 = b_bcpy(value,*lb))!=0 || rc2!=0)
           {
           if (rc1==ALLOC || rc2==ALLOC)
              e_trap(ALLOCATION,2,E_TMLT,&value);
           else
              e_trap(OVERFLOW,2,E_TMLT,&value);
           }

        /* interchange bounds if negative value */
        else if (value->z==0 && value->s==1)
           {
           bb = *lb;
           *lb = *ub;
           *ub = bb;
           }

        if (value->f) l_free(&value);

        E_TPOPP("b_irnd")
        return;
        }





