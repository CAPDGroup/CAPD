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

/* CVS $Id: a_div_.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_div_.c                              */
/*                                                              */
/*      Entries         : a_intg a_div_(i,j)                    */
/*                        a_intg i,j;                           */
/*                                                              */
/*      Arguments       : i = integer numerator                 */
/*                        j = integer denominator               */
/*                                                              */
/*      Description     : Integer division according to         */
/*                        PASCAL Standard:                      */
/*                           error if j is zero                 */
/*                 abs(i)-abs(j) < abs((i div j)*j) <= abs(j)   */
/*                 sign(i div j) = sign(i)*sign(j)              */
/*                                                              */
/*      Note            : 0 returned in case of error           */
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
local a_intg a_div_(a_intg i,a_intg j)
#else
local a_intg a_div_(i,j)

a_intg i;
a_intg j;
#endif
        {
        a_intg res;

        E_TPUSH("a_div_")

        if (j==0)
           {
           e_trap(DIV_BY_ZERO,4,E_TINT+E_TEXT(1),&i,
                                E_TINT+E_TEXT(2),&j);
           res = 0;
           }
        else
           {
           if (i>=0)
              {
              if (j>0)
                 res = i/j;
              else if (j==MININT)
                 res = 0;
              else
                 res = -(i/(-j));
              }
           else if (i==MININT)
              {
              if (j>0)
                 {
                 i += j;
                 res = -((-i)/j)-1;
                 }
              else if (j==MININT)
                 res = 1;
              else if (j==-1)
                 {
                 e_trap(OVERFLOW,4,E_TINT+E_TEXT(1),&i,
                                   E_TINT+E_TEXT(2),&j);
                 res = 0;
                 }
              else
                 {
                 i -= j;
                 res = (-i)/(-j)+1;
                 }
              }
           else
              {
              if (j>0)
                 res = -((-i)/j);
              else
                 res = (j==MININT) ? 0 : (-i)/(-j);
              }
           }

        E_TPOPP("a_div_")
        return(res);
        }





