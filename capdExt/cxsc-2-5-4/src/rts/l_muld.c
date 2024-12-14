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

/* CVS $Id: l_muld.c,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_muld.c                              */
/*                                                              */
/*      Entries         : multiprecision l_muld(i1,i2)          */
/*                        multiprecision i1,i2;                 */
/*                                                              */
/*      Description     : Muliplication of multi-precision      */
/*                        numbers with rounding downwards.      */
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
local multiprecision l_muld(multiprecision i1,multiprecision i2)
#else
local multiprecision l_muld(i1,i2)

multiprecision i1;
multiprecision i2;
#endif
        {
        int rc;
        multiprecision res1,res2;

        E_TPUSH("l_muld")

        l_init(&res1);

        if (res1==NULL)
                e_trap(ALLOCATION,2,E_TMSG,65);
        else
            {
            l_init(&res2);

            if (res2==NULL)
                e_trap(ALLOCATION,2,E_TMSG,65);
            else if ((rc = b_bmul(i1,i2,res1,res2))!=0)
                {
                switch(rc)
                   {
                   case OFLOW : e_trap(OVERFLOW,4,E_TMLT+E_TEXT(1),&i1,
                                                  E_TMLT+E_TEXT(2),&i2);
                                break;
                   case UFLOW : e_trap(UNDERFLOW,4,E_TMLT+E_TEXT(1),&i1,
                                                   E_TMLT+E_TEXT(2),&i2);
                                break;
                   case ALLOC : e_trap(ALLOCATION,2,E_TMSG,65);
                                break;
                   }
                }
            else
               {
               res1->f = 1;
               if (res1->s)
                  {
                  if ((rc = b_brnd(res1,res1))!=0)
                     {
                     switch(rc)
                        {
                   case OFLOW : e_trap(OVERFLOW,4,E_TMLT+E_TEXT(1),&i1,
                                                  E_TMLT+E_TEXT(2),&i2);
                                break;
                   case ALLOC : e_trap(ALLOCATION,2,E_TMSG,65);
                                break;
                        }
                     }
                  }
               }

            l_free(&res2);
            }

        if (i1->f) l_free(&i1);
        if (i2->f) l_free(&i2);

        E_TPOPP("l_muld")
        return(res1);
        }





