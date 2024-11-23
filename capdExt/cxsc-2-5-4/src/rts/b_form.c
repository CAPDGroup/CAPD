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

/* CVS $Id: b_form.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_form.c                              */
/*                                                              */
/*      Entries         : a_btyp b_form                         */
/*                           (dc,size,expo,dp,length,sign,rnd,s)*/
/*                        a_btyp *dc;                           */
/*                        a_intg size;                          */
/*                        a_intg expo;                          */
/*                        a_intg dp;                            */
/*                        a_intg length;                        */
/*                        a_intg rnd;                           */
/*                        a_bool sign;                          */
/*                        a_real *s;                            */
/*                                                              */
/*      Arguments       : dc - decimal digits                   */
/*                        size - length of dc                   */
/*                        expo - converted exponent             */
/*                        dp - position of decimal point        */
/*                        length - number of dc                 */
/*                        sign - sign of number                 */
/*                        rnd - rounding mode                   */
/*                              -1 = rounding downwards         */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        s - resultant IEEE variable           */
/*                                                              */
/*      Function value  : error indicator                       */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to IEEE double format                 */
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
extern a_real *r_sero;
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_btyp b_form(a_btyp *dc,a_intg *size,a_intg expo,a_intg dp,
                    a_intg length,a_bool sign,a_intg rnd,a_real *s)
#else
local a_btyp b_form(dc,size,expo,dp,length,sign,rnd,s)

a_btyp *dc;
a_intg *size;
a_intg expo;
a_intg dp;
a_intg length;
a_bool sign;
a_intg rnd;
a_real *s;
#endif
        {
        a_intg k;
        a_intg bits;
        a_btyp rc,rc1,mant[BSIZE];
        a_intg p_dp,p_length;
        char *buffer;

        E_TPUSH("b_form")

        /* number is zero                               */
        if (dp==0 && length==0)
           {
           R_ASSIGN(*s,((sign) ? *r_sero : *r_zero));

           E_TPOPP("b_form")
           return(NO_ERROR);
           }

        /* initialize data                              */
        B_CLEAR(mant,BSIZE)
        bits = MANTL+2;
        buffer = (char *)dc;
/*
printf("b_form :   size=%d\n",*size);
printf("b_form :   expo=%d\n",expo);
printf("b_form :     dp=%d\n",dp);
printf("b_form : length=%d\n",length);
printf("b_form :   sign=%d\n",sign);
printf("b_form :    rnd=%d\n",rnd);
{int i;
 printf("b_form :     dc=");
 for (i=0;i<length;i++) printf("%2.2x",buffer[i]);
 printf("\n");
}
*/

        /* adjust and pack decimal digits */
        /* k = index of start of packed mantissa */
        if (b_adpp((a_btyp **)&buffer,size,
                   expo,dp,length,&k,&p_dp,&p_length))
           {
           E_TPOPP("b_form")
           return(ALLOCATION);
           }

        B_CLEAR(b_cm__,A_LENGTH)
        b_cm__[A_SIGN] = (sign) ? LSB : ZERO;

        /* convert integer part                         */
        if (p_dp>k) b_coni(p_dp-k,&dc[k],
                           (a_intg *)&b_cm__[A_BEGIN],
                           (a_intg *)&b_cm__[A_END],b_cm__,&bits);

        /* convert fraction part                        */
        if (p_length>p_dp)
           {
           if (bits>0) b_conf(p_length-p_dp,&dc[p_dp],
                              (a_intg *)&b_cm__[A_BEGIN],
                              (a_intg *)&b_cm__[A_END],b_cm__,&bits);
           else b_cm__[b_cm__[A_END]] |= LSB;
           }

        /* eliminate trailing zeros                                  */
        while (b_cm__[b_cm__[A_END]]==ZERO) b_cm__[A_END]--;

        /* eliminate leading zeros                                   */
        while (b_cm__[b_cm__[A_BEGIN]]==ZERO) b_cm__[A_BEGIN]++;

/*
{int i;
 printf("b_form :   akku=");
 for (i=b_cm__[A_BEGIN];i<=b_cm__[A_END];i++)
    printf("%08.8lx ",b_cm__[i]);
 printf("\n");
}
*/

        /* get accu                                                  */
        if (b_geta(b_cm__,mant,&expo,&sign))
           {
           R_ASSIGN(*s,((sign) ? *r_sero : *r_zero));

           E_TPOPP("b_form")
           return(NO_ERROR);
           }

/*
printf("b_form :   mant=%08.8lx %08.8lx %08.8lx %08.8lx %08.8lx\n",
mant[0],mant[1],mant[2],mant[3],mant[4]);
*/

        /* adjust denormalized number                                */
        rc = b_adj(mant,&expo);

/*
printf("b_form :   mant=%08.8lx %08.8lx %08.8lx %08.8lx %08.8lx\n",
mant[0],mant[1],mant[2],mant[3],mant[4]);
*/

        /* round                                                     */
        if (rnd==0) rc1 = b_rndn(mant,&expo);
        else if (rnd<0) rc1 = b_rndd(mant,&expo,sign);
        else rc1 = b_rndu(mant,&expo,sign);

/*
printf("b_form :   mant=%08.8lx %08.8lx\n",mant[0],mant[1]);
printf("b_form :   expo=%d\n",expo);
*/

        if (rc==NO_ERROR) rc = rc1;

        /* composition                                               */
        b_comp(s,expo,mant,sign);

        E_TPOPP("b_form")
        return(rc);
        }





