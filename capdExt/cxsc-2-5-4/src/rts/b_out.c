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

/* CVS $Id: b_out.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_out.c                               */
/*                                                              */
/*      Entries         : void b_out                            */
/*                        (mant,expo,digits,buffer,bdp,dexpo)   */
/*                        a_btyp *mant;                         */
/*                        a_intg *dexpo,digits,expo,*bdp;       */
/*                        char *buffer;                         */
/*                                                              */
/*      Arguments       : expo - exponent of IEEE value         */
/*                        mant - mantissa of IEEE value         */
/*                        digits - number of digits             */
/*                        buffer - output string                */
/*                        bdp - position of decimal point       */
/*                        dexpo - decimal exponent of first     */
/*                                non-zero digit                */
/*                                dexpo<0 : float format        */
/*                                dexpo>=0: fixed format        */
/*                                                              */
/*                                                              */
/*      Description     : Decimal representation determined from*/
/*                        exponent and mantissa of IEEE value   */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern dotprecision b_cm__;
#endif

#ifdef LINT_ARGS
local void b_out(a_btyp *mant,a_intg expo,a_intg digits,char *buffer,
                 a_intg *bdp,a_intg *dexpo)
#else
local void b_out(mant,expo,digits,buffer,bdp,dexpo)

a_btyp *mant;
a_intg expo;
a_intg digits;
char *buffer;
a_intg *bdp;
a_intg *dexpo;
#endif
        {
        a_intg i,k;

        /* allign mantissa to dotprecision              */
        k = B_ASHR(expo,LOG_B_LENGTH);
        if ( (i = (expo & (B_LENGTH-1))-EXPO_SHIFT) <0)
             {
#if C_P_3
             mant[2] = (mant[2]>>-i) | (mant[1]<<(B_LENGTH+i));
             mant[1] = (mant[1]>>-i) | (mant[0]<<(B_LENGTH+i));
             mant[0] = (mant[0]>>-i);
#else
             b_shru(mant,D_U_RATIO+1,-i);
#endif
             }
        else if (i>0)
             {
#if C_P_3
             mant[0] = (mant[0]<<i) | (mant[1]>>(B_LENGTH-i));
             mant[1] <<= i;
#else
             b_shlu(mant,D_U_RATIO,i);
#endif
             }

        b_cm__[A_BEGIN] = A_D_P-k;
        b_cm__[A_END] = (A_D_P+D_U_RATIO)-k;

   /* copy value to I/O-buffer                                  */
        for (i=D_U_RATIO;i>=0;i--)
             b_cm__[A_D_P-k+i] = mant[i];

   /* clear accu contents between number and decimal point      */
        for (i=(A_D_P+D_U_RATIO+1)-k;i<=A_D_P;i++)
            b_cm__[i] = ZERO;
        for (i=A_D_P+1;i<=(A_D_P-1)-k;i++)
            b_cm__[i] = ZERO;

/*                                                                   */
/* conversion of integer part                                        */
/*                                                                   */
        if (expo>=0) {
            b_outi(&digits,buffer,bdp,dexpo,b_cm__);
            }

/*                                                                   */
/* conversion of fraction part                                       */
/*                                                                   */
        if (digits>0) {
            b_outf(&digits,buffer,bdp,dexpo,b_cm__);
            }

        return;
        }





