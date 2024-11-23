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

/* CVS $Id: b_outm.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_outm.c                              */
/*                                                              */
/*      Entries         : void b_outm                           */
/*                      (mant,len,expo,digits,buffer,bdp,dexpo) */
/*                        a_intg len,*dexpo,digits,expo,*bdp;   */
/*                        a_btyp *mant;                         */
/*                        char *buffer;                         */
/*                                                              */
/*      Arguments       : mant   - mantissa of intern value     */
/*                        len    - length of mantissa           */
/*                        expo   - exponent of intern value     */
/*                        digits - number of digits             */
/*                        buffer - output string                */
/*                        bdp    - position of decimal point    */
/*                        dexpo  - decimal exponent of first    */
/*                                non-zero digit                */
/*                                dexpo<0 : float format        */
/*                                dexpo>=0: fixed format        */
/*                                                              */
/*      Description     : Decimal representation determined from*/
/*                        exponent and mantissa of multi-       */
/*                        precision value                       */
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
local void b_outm(a_btyp *mant,a_intg len,a_intg expo,a_intg digits,
                  char *buffer,a_intg *bdp,a_intg *dexpo)
#else
local void b_outm(mant,len,expo,digits,buffer,bdp,dexpo)

a_btyp *mant;
a_intg len;
a_intg expo;
a_intg digits;
char *buffer;
a_intg *bdp;
a_intg *dexpo;
#endif
        {
        a_intg i;

        E_TPUSH("b_outm")

        b_cm__[A_BEGIN] = A_D_P-expo;
        b_cm__[A_END] = (A_D_P+len-1)-expo;

        /* insufficient buffer size                             */
        if (b_cm__[A_END]>A_LENGTH || b_cm__[A_BEGIN]<A_START)
           e_trap(I_O_BUFFER,2,E_TMSG,39);

        /* copy value to I/O-buffer                             */
        for (i=len-1;i>=0;i--)
           b_cm__[A_D_P-expo+i] = mant[i];

        /* clear accu contents between number and decimal point */
        for (i=A_D_P+len-expo;i<=A_D_P;i++)
           b_cm__[i] = ZERO;
        for (i=A_D_P+1;i<=(A_D_P-1)-expo;i++)
           b_cm__[i] = ZERO;

        /* conversion of integer part    */
        if (expo>=0)
           {
           b_outi(&digits,buffer,bdp,dexpo,b_cm__);
           }

        /* conversion of fraction part                          */
        if (digits>0)
           {
           b_outf(&digits,buffer,bdp,dexpo,b_cm__);
           }

        E_TPOPP("b_outm")
        return;
        }





