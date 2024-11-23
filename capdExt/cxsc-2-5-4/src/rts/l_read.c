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

/* CVS $Id: l_read.c,v 1.22 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_read.c                              */
/*                                                              */
/*      Entries         : void l_read(device,s,rnd,ii)          */
/*                        FILE *device;                         */
/*                        multiprecision *s;                    */
/*                        a_intg rnd;                           */
/*                        a_intg ii;                            */
/*                                                              */
/*      Arguments       : device - input device reading string  */
/*                        s - resultant multiprecision value    */
/*                        rnd - rounding mode                   */
/*                              -1 = rounding downwards         */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        ii - first input character            */
/*                                                              */
/*      Description     : Read a character string from input    */
/*                        device and convert to multiprecision  */
/*                                                              */
/*      Note            : Length of input string is restricted  */
/*                        by BUFFERSIZE.                        */
/*                                                              */
/*                   b_adpp() used                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern d_otpr b_cp__;
extern a_btyp b_maxl;
#endif

#ifdef LINT_ARGS
local void l_read(FILE *device,multiprecision *s,a_intg rnd,a_intg ii)
#else
local void l_read(device,s,rnd,ii)

FILE *device;
multiprecision *s;
a_intg rnd;
a_intg ii;
#endif
        {
        a_intg bits,dp,expo,length;
        a_intg p_dp,p_length,res_length;
        a_intg a_b,a_e,i,k,size,offset;
        a_bool sign;
        a_btyp *res;
        char *buffer;

        E_TPUSH("l_read")

        buffer = (char *)b_cp__;
        size = BUFFERSIZE;

        /* scan digits from device and store decimal value in buffer */
        switch (b_scan(device,&buffer,&size,&expo,&dp,&length,&sign,(int)ii))
           {
           case  1: e_trap(ALLOCATION,2,E_TMSG,56);
                    E_TPOPP("l_read")
                    return;
           case  2:
           case  3:
           case  4: e_trap(I_O_ERROR,2,E_TMSG,58);
                    E_TPOPP("l_read")
                    return;
           case  5: e_trap(NO_ERROR+E_EMSG+E_ECNT,2,E_TMSG,64);
           default: ;
           }

        /* number is zero                                       */
        if (dp==0 && length==0)
           {
           (*s)->z = TRUE;
           (*s)->r = 0;
           E_TPOPP("l_read")
           return;
           }

        /* initialize data                                      */
        (*s)->s = sign;
        (*s)->z = FALSE;

        /* adjust and pack decimal digits */
        /* k = index of start of packed mantissa */
        if (b_adpp((a_btyp **)&buffer,&size,
                   expo,dp,length,&k,&p_dp,&p_length))
           {
           e_trap(ALLOCATION,2,E_TMSG,56);
           E_TPOPP("l_read")
           return;
           }

/*
printf("l_read : size=%d p_start=%d p_dp=%d p_length=%d\n",
size,k,p_dp,p_length);
printf("         buffer=");
for (i=k;i<p_length;i++) printf("%08.8lx",((a_btyp *)buffer)[i]);
printf("\n");
*/
        /* minimum number of required bits */
        bits = B_LENGTH*b_maxl+2;
/*
printf("l_read : bits=%d\n",bits);
*/
        /* allocate output buffer */
        offset = (p_dp-k>A_D_P) ? p_dp-k-A_D_P : 0;
        i = offset+A_D_P+((p_length<A_LENGTH) ? A_LENGTH : p_length);
        if (b_maxl+2>i) i = b_maxl+2;
        if ((res = (a_btyp *)malloc((size_t)(i*sizeof(a_btyp))))==NULL)
           {
           e_trap(ALLOCATION,2,E_TMSG,56);
           E_TPOPP("l_read")
           return;
           }

#ifdef HEAP_CHECK
b_geth((a_char *)&res,(a_char *)res,(a_char *)"l_read");
#endif

        /* clear output buffer */
        B_CLEAR(res,i)
        a_e = a_b = 0;

        /* convert integer part if there is an integer part */
        if (p_dp>k) b_coni(p_dp-k,&((a_btyp *)buffer)[k],
                           &a_b,&a_e,&res[offset],&bits);
/*
printf("l_read : a_b=%d a_e=%d res=",a_b,a_e);
for (i=a_b;i<=a_e;i++) printf("%08.8lx",res[i]);
printf("\n");

printf("l_read : bits=%d\n",bits);
*/
        /* convert fraction part if any and there are still bits required */
        if (p_length>p_dp)
           {
           if (bits>0)
              b_conf(p_length-p_dp,&((a_btyp *)buffer)[p_dp],
                     &a_b,&a_e,&res[offset],&bits);
           else res[a_e] |= LSB;
           }
/*
printf("l_read : A_BEGIN=%d A_END=%d res=",a_b,a_e);
for (i=a_b;i<=a_e;i++) printf("%08.8lx",res[i]);
printf("\n");
*/
        /* eliminate trailing zeros */
        while (res[a_e]==ZERO) a_e--;

        /* eliminate leading zeros */
        while (res[a_b]==ZERO) a_b++;

        (*s)->e = A_D_P-a_b;

        res_length = a_e-a_b+1;
        length = (res_length<b_maxl) ? res_length : b_maxl;

        if (length!=(*s)->l)
           {
           if ((*s)->l)
              {
              (*s)->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&(*s)->m,(a_char *)(*s)->m,(a_char *)"l_read");
#endif
              B_FREE((*s)->m)
              }
           if (b_ball(length,&(*s)->m))
              {
              e_trap(ALLOCATION,2,E_TMSG,56);

              E_TPOPP("l_read")
              return;
              }
           (*s)->l = length;
           }

        /* copy mantissa digits */
        for (i=0;i<length;i++)
           (*s)->m[i] = res[i+A_D_P-(*s)->e];

        /* determine rounded result */
        (*s)->r = NOT(b_test(res_length-b_maxl,&res[a_b+b_maxl]));
        if ((((rnd<0 && sign) || (rnd>0 && NOT(sign))) && (*s)->r) ||
            (rnd==0 && res_length>b_maxl && (res[a_b+b_maxl] & MSB)))
           {
           if (b_bcad(b_maxl,(*s)->m))
              {
              if ((*s)->e==MAXINT)
                 e_trap(OVERFLOW,2,E_TMSG,56);
              (*s)->e++;
              }
           }

        /* free output buffer */
#ifdef HEAP_CHECK
b_freh((a_char *)&res,(a_char *)res,(a_char *)"l_read");
#endif
        free((char *)res);

        if (size!=BUFFERSIZE)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&buffer,(a_char *)buffer,(a_char *)"l_read");
#endif
           free(buffer);
           }

        E_TPOPP("l_read")
        return;
        }





