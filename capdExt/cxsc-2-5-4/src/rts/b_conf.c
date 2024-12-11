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

/* CVS $Id: b_conf.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_conf.c                              */
/*                                                              */
/*      Entries         : void b_conf(dp,pd,a_b,a_e,d,bits)     */
/*                        a_intg dp;                            */
/*                        a_btyp *pd;                           */
/*                        a_intg *a_b;                          */
/*                        a_intg *a_e;                          */
/*                        d_otpr d;                             */
/*                        a_intg *bits;                         */
/*                                                              */
/*      Arguments       : dp      - length of pd                */
/*                        pd      - packed decimals             */
/*                        a_b     - first position in mantissa  */
/*                        a_e     - last  position in mantissa  */
/*                        d       - dotprecision value          */
/*                        bits    - number of required bits     */
/*                                                              */
/*      Description     : Convert fractional part of number     */
/*                                                              */
/*      Note            : No range check on array pd.           */
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
local void b_conf(a_intg dp,a_btyp *pd,a_intg *a_b,a_intg *a_e,
                  d_otpr d,a_intg *bits)
#else
local void b_conf(dp,pd,a_b,a_e,d,bits)

a_intg dp;
a_btyp *pd;
a_intg *a_b;
a_intg *a_e;
d_otpr d;
a_intg *bits;
#endif
        {
        a_bool new_v;
        a_intg block;
        a_btyp carry,*q,*p;

        /* initialize                                   */
        block = B_LENGTH;
        new_v = (*a_b==0) ? TRUE : FALSE;
        if (new_v) *a_b = A_D_P+1;
        *a_e = A_D_P+1;
        q = &d[A_D_P+1];

        /* convert by repeated multiplication           */
        while (*bits>-BIN_PACK)
           {
           if ( (block -= BIN_PACK) <0)
              {
              block = B_LENGTH-BIN_PACK;
              q++;
              (*a_e)++;
              }

           carry = ZERO;
           for (p=pd+dp-1;p>=pd;p--)
              {
              *p = ((*p)<<BIN_PACK)+carry;
              carry = (*p) / DEC_PACK_POWER;
              *p %= DEC_PACK_POWER;
              }
           *q |= carry<<block;

           if (new_v)
              {
              if (carry) new_v = FALSE;
              }
           else *bits -= BIN_PACK;
           }

        while (--dp>=0)
           {
           if (pd[dp])
               {
               d[*a_e] |= LSB;
               break;
               }
           }

        return;
        }





