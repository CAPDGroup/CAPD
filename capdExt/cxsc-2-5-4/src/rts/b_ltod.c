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

/* CVS $Id: b_ltod.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_ltod.c                              */
/*                                                              */
/*      Entries         : a_btyp b_ltod(i,d,rnd)                */
/*                        multiprecision i;                     */
/*                        d_otpr *d;                            */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : i    = long value                     */
/*                        d    = dotprecision variable          */
/*                        rnd  = rounding                       */
/*                               <0 : downwards                 */
/*                               =0 : nearest                   */
/*                               >0 : upwards                   */
/*                                                              */
/*      Function value  : error code                            */
/*                                                              */
/*      Description     : convert long to dotprecision.         */
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
local a_btyp b_ltod(multiprecision i,d_otpr *d,a_intg rnd)
#else
local a_btyp b_ltod(i,d,rnd)

multiprecision i;
d_otpr *d;
a_intg rnd;
#endif
        {
        register a_btyp *rega,k;

        E_TPUSH("b_ltod")

        rega = *d;
        B_CLEAR(rega,A_LENGTH) 

        /* number is zero */
        if (i->z)
           {
           E_TPOPP("b_ltod")
           return(NO_ERROR);
           }

        /* sign and exponent of real value */
        rega[A_SIGN] = i->s;
        rega[A_BEGIN] = A_D_P-i->e;
        rega[A_END] = A_D_P-i->e+i->l-1;
        rega[A_STATUS] |= A_MZERO | A_PZERO;
 
        /* overflow occurred */
        if (rega[A_BEGIN]<A_START)
           {
           e_trap(OVERFLOW,2,E_TMSG,48); 
           rega[A_STATUS] |= (i->s) ? A_MINFINITY : A_PINFINITY;

           E_TPOPP("b_ltod")
           return(OVERFLOW);
           }
        /* underflow occurred */
        else if (rega[A_BEGIN]>=A_LENGTH)
           {
           e_trap(UNDERFLOW,0);
           if ((i->s && rnd<0) || (!i->s && rnd>0))
              {
              rega[A_BEGIN] = rega[A_END] = A_LENGTH-1;
              rega[A_LENGTH-1] = LSB;
              }
           else
              {
              rega[A_BEGIN] = rega[A_END] = ZERO;
              }
           }
        /* copy mantissa of multiprecision value */
        else
           {
           /* chop length of mantissa */
           if (rega[A_END]>=A_LENGTH)
              {
              /* rnd holds rounding information */
              rnd = ((i->s && rnd<0) || (!i->s && rnd>0));
              rega[A_END] = A_LENGTH-1;
              }
           else
              rnd = 0;
       
           for (k=rega[A_BEGIN];k<=rega[A_END];k++) 
              rega[k] = i->m[k-rega[A_BEGIN]]; 

           /* propagate rounding carry */
           if (rnd)
              {
              if (b_bcad(rega[A_END]-rega[A_BEGIN]+1,
                         &rega[rega[A_END]]))
                 {
                 rega[A_END] = ++rega[A_BEGIN];
                 /* overflow is not considered */
                 }
              }

           /* eliminate trailing zeros */
           while (rega[rega[A_END]]==ZERO) rega[A_END]--;
           }

        E_TPOPP("b_ltod")
        return(NO_ERROR);
        }





