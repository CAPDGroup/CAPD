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

/* CVS $Id: b_dtol.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_dtol.c                              */
/*                                                              */
/*      Entries         : a_btyp b_dtol(d,i,rnd)                */
/*                        d_otpr d;                             */
/*                        multiprecision *i;                    */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : i    = long variable                  */
/*                        d    = dotprecision value             */
/*                        rnd  = rounding mode                  */
/*                                                              */
/*      Function value  : ALLOCATION - allocation error         */
/*                                                              */
/*      Description     : convert dotprecision to long.         */
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
local a_btyp b_dtol(d_otpr d,multiprecision *i,a_intg rnd)
#else
local a_btyp b_dtol(d,i,rnd)

d_otpr d;
multiprecision *i;
a_intg rnd;
#endif
        {
        a_btyp j,l;

        E_TPUSH("b_dtol")

        (*i)->r = 0;
        (*i)->f = 0;

        if (d[A_BEGIN]==ZERO)
           {
           (*i)->z = TRUE;
           }
        else
           {
           (*i)->z = FALSE;
           (*i)->s = (d[A_SIGN]) ? TRUE : FALSE;

           /* required length of long value mantissa = l+1         */
           if ((l = d[A_END]-d[A_BEGIN])>=b_maxl) 
              {
              l = b_maxl-1;
              (*i)->r = 1;

              /* remove trailing zeros */
              while (d[d[A_BEGIN]+l]==ZERO) l--;
              }

           /* allocate and copy mantissa */
           if (l+1!=(*i)->l)
              {
              if ((*i)->l)
                 {
                 (*i)->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&(*i)->m,(a_char *)(*i)->m,(a_char *)"b_dtol");
#endif
                 B_FREE((*i)->m)
                 }
              if (b_ball(l+1,&(*i)->m)) return(ALLOCATION);
              (*i)->l = l+1;
              }
           for (j=0;j<=l;j++) (*i)->m[j] = d[d[A_BEGIN]+j];

           (*i)->e = A_D_P-d[A_BEGIN];
           }

        return(NO_ERROR);
        }





