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

/* CVS $Id: b_ifrm.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_ifrm.c                              */
/*                                                              */
/*      Entries         : a_btyp b_ifrm                         */
/*                           (dc,expo,dp,length,sign,s)         */
/*                        a_btyp *dc;                           */
/*                        a_intg expo,dp,length;                */
/*                        a_bool sign;                          */
/*                        a_intv *s;                            */
/*                                                              */
/*      Arguments       : dc - decimal digits                   */
/*                        expo - converted exponent             */
/*                        dp - position of decimal point        */
/*                        length - number of dc                 */
/*                        sign - sign of number                 */
/*                        s - resultant interval                */
/*                                                              */
/*      Function value  : error indicator                       */
/*                                                              */
/*      Description     : Convert a character string of a real  */
/*                        number to the smallest enclosing      */
/*                        interval.                             */
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
local a_btyp b_ifrm(a_btyp *dc,a_intg expo,a_intg dp,a_intg length,
                      a_bool sign,a_intv *s)
#else
local a_btyp b_ifrm(dc,expo,dp,length,sign,s)

a_btyp *dc;
a_intg expo;
a_intg dp;
a_intg length;
a_bool sign;
a_intv *s;
#endif
        {
        a_intg k;
        a_intg bits,size;
        a_btyp rc,rc1,mant[BSIZE],mantb[BSIZE];
        a_intg expob,p_dp,p_length;
        char *buffer;

        E_TPUSH("b_ifrm")

        /* number is zero       */
        if (dp==0 && length==0) {
           if (sign)
              {
              R_ASSIGN(s->INF,*r_sero);
              R_ASSIGN(s->SUP,*r_sero);
              }
           else
              {
              R_ASSIGN(s->INF,*r_zero);
              R_ASSIGN(s->SUP,*r_zero);
              }

           E_TPOPP("b_ifrm")
           return(NO_ERROR);
           }

        /* initialize data      */
        B_CLEAR(b_cm__,A_LENGTH)
        B_CLEAR(mant,BSIZE)
        b_cm__[A_SIGN] = (sign) ? LSB : ZERO;
        bits = MANTL+2;
        buffer = (char *)dc;
        size = BUFFERSIZE;

        /* adjust and pack decimal digits */
        /* k = index of start of packed mantissa */
        if (b_adpp((a_btyp **)&buffer,&size,
                   expo,dp,length,&k,&p_dp,&p_length))
           {
           E_TPOPP("b_form")
           return(ALLOCATION);
           }

        /* convert integer part         */
        if (p_dp>k) b_coni(p_dp-k,&dc[k],
                       (a_intg *)&b_cm__[A_BEGIN],
                       (a_intg *)&b_cm__[A_END],b_cm__,&bits);

        /* convert fraction part        */
        if (p_length>p_dp)
           {
           if (bits>0) b_conf(p_length-p_dp,&dc[p_dp],
                              (a_intg *)&b_cm__[A_BEGIN],
                              (a_intg *)&b_cm__[A_END],
                              b_cm__,&bits);
           else b_cm__[b_cm__[A_END]] |= LSB;
           }

        /* eliminate trailing zeros     */
        while (b_cm__[b_cm__[A_END]]==ZERO) b_cm__[A_END]--;

        /* eliminate leading zeros      */
        while (b_cm__[b_cm__[A_BEGIN]]==ZERO) b_cm__[A_BEGIN]++;

        /* get accu     */
        if (b_geta(b_cm__,mant,&expo,&sign))
           {
           if (sign)
              {
              R_ASSIGN(s->INF,*r_sero);
              R_ASSIGN(s->SUP,*r_sero);
              }
           else
              {
              R_ASSIGN(s->INF,*r_zero);
              R_ASSIGN(s->SUP,*r_zero);
              }
           rc = NO_ERROR;
           }
        else
           {

           /* adjust denormalized number   */
           rc = b_adj(mant,&expo);

           /* round        */
           expob = expo;
           memcpy(mantb,mant,BSIZE*sizeof(a_btyp));
           rc1 = b_rndd(mant,&expo,sign);
           if (rc==NO_ERROR) rc = rc1;
           rc1 = b_rndu(mantb,&expob,sign);
           if (rc==NO_ERROR) rc = rc1;
           if (rc==INEXACT) rc = NO_ERROR;

           /* composition  */
           b_comp(&(s->INF),expo,mant,sign);
           b_comp(&(s->SUP),expob,mantb,sign);
           }

        E_TPOPP("b_ifrm")
        return(rc);
        }





