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

/* CVS $Id: b_rtol.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_rtol.c                              */
/*                                                              */
/*      Entries         : a_btyp b_rtol(r,i,rnd)                */
/*                        a_real r;                             */
/*                        multiprecision *i;                    */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : i    = long value                     */
/*                        r    = real value                     */
/*                        rnd  = rounding mode                  */
/*                                                              */
/*      Function value  : ALLOCATION - allocation error         */
/*                                                              */
/*      Description     : convert real to long value.           */
/*                                                              */
/*      Note            : Conversion is always exact.           */
/*                        Infinity and NaNs not considered.     */
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
local a_btyp b_rtol(a_real r,multiprecision *i,a_intg rnd)
#else
local a_btyp b_rtol(r,i,rnd)

a_real r;
multiprecision *i;
a_intg rnd;
#endif
        {
        a_btyp mant[D_U_RATIO+1];
        a_intg expo,j,k,l,m;
        a_bool sign;

        /* decompose real value                                 */
        if (b_deko(r,&expo,mant,&sign)) return(NO_ERROR);
        mant[D_U_RATIO] = ZERO;

        /* allign mantissa to long variable                     */
        k = B_ASHR(expo,LOG_B_LENGTH);
        if ( (m = (expo & (B_LENGTH-1))-EXPO_SHIFT) < 0)
                {
#if C_P_3
                mant[2] = mant[1]<<(B_LENGTH+m);
                mant[1] = (mant[1]>>-m) |
                          (mant[0]<<(B_LENGTH+m));
                mant[0] >>= -m;
#else
                b_shru(mant,D_U_RATIO+1,-m);
#endif
                }
        else if (m>0)
                {
#if C_P_3
                mant[0] = (mant[0]<<m) |
                          (mant[1]>>(B_LENGTH-m));
                mant[1] <<= m;
#else
                b_shlu(mant,D_U_RATIO,m);
#endif
                }

        (*i)->z = 0;
        (*i)->s = sign;
        (*i)->r = 0;
        (*i)->f = 0;

        /* required length of long value mantissa               */
        for (m=0;m<=D_U_RATIO;m++)
            if (mant[m]) break;
        for (l=D_U_RATIO;l>=0;l--)
            if (mant[l]) break;

        /* allocate and copy mantissa                           */
        if (l-m+1!=(*i)->l)
           {
           if ((*i)->l)
              {
              (*i)->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&(*i)->m,(a_char *)(*i)->m,(a_char *)"b_rtol");
#endif
              B_FREE((*i)->m)
              }
           if (b_ball(l-m+1,&(*i)->m)) return(ALLOCATION);
           (*i)->l = l-m+1;
           }
        for (j=m;j<=l;j++) (*i)->m[j-m] = mant[j];

        (*i)->e = k-m;

        return(NO_ERROR);
        }





