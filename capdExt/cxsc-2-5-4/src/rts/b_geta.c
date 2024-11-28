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

/* CVS $Id: b_geta.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_geta.c                              */
/*                                                              */
/*      Entries         : a_bool b_geta(a,result,expo,vz)       */
/*                        dotprecision a;                       */
/*                        a_btyp *result;                       */
/*                        a_intg *expo;                         */
/*                        a_bool *vz;                           */
/*                                                              */
/*      Arguments       : a = dotprecision variable             */
/*                        result = a_btyp array[BSIZE]          */
/*                        expo = exponent of number             */
/*                        vz = sign of number                   */
/*                                                              */
/*      Function value  : TRUE if accu is empty                 */
/*                        FALSE otherwise                       */
/*                                                              */
/*      Description     : Get contents of accu ready for use    */
/*                        as IEEE number                        */
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
local a_bool b_geta(dotprecision a,a_btyp *result,a_intg *expo,
                     a_bool *vz)
#else
local a_bool b_geta(a,result,expo,vz)

dotprecision a;
a_btyp *result;
a_intg *expo;
a_bool *vz;
#endif
        {
        a_intg i,k;
        a_btyp q;

        /* accu is empty                                        */
        if (a[A_BEGIN]==ZERO) return(TRUE);

        /* exponent determined from begin of accu               */
        *expo = (A_D_P-a[A_BEGIN])*B_LENGTH+EXPO_SHIFT;

        /* copy accu to result                                  */
        k = a[A_END]-a[A_BEGIN];
        for (i=0;i<=k && i<BSIZE-1;i++)
           result[i] = a[a[A_BEGIN]+i];

        while (i<BSIZE) { result[i] = ZERO; i++; }

        /* normalization of mantissa */
        if (SHFT_MASK & *result)
           {
           b_shru(result,BSIZE,ZERO_BITS-1);
           (*expo) += ZERO_BITS-1;
           }
        q = *result;
        for(i=0;(HIDDEN_BIT & q)==ZERO;i++) q <<= 1;
        if (i>0)
           {
           b_shlu(result,BSIZE,i);
           (*expo) -= i;
           }

        /* set any insignificant bit in result if inexact       */
        if (k>=BSIZE-1) result[BSIZE-1] |= LSB;

        /* sign of accu                                         */
        *vz = (a_bool) a[A_SIGN];

        return(FALSE);
        }





