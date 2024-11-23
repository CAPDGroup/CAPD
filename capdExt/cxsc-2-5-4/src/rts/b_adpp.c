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

/* CVS $Id: b_adpp.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_adpp.c                              */
/*                                                              */
/*      Entries         : a_btyp b_adpp                         */
/*                           (dc,size,expo,dp,length,           */
/*                            p_start,p_dp,p_length)            */
/*                        a_btyp **dc;                          */
/*                        a_intg *size;                         */
/*                        a_intg expo;                          */
/*                        a_intg dp;                            */
/*                        a_intg length;                        */
/*                        a_intg *p_start;                      */
/*                        a_intg *p_dp;                         */
/*                        a_intg *p_length;                     */
/*                                                              */
/*      Arguments       : dc     - decimal digits               */
/*                        size   - size of dc                   */
/*                        expo   - converted exponent           */
/*                        dp     - position of decimal point    */
/*                        length - number of decimal digits     */
/*                        p_start - start of packed mantissa    */
/*                        p_dp   - packed decimal point         */
/*                        p_length - length of packed decimals  */
/*                                                              */
/*      Function value  : error indicator                       */
/*                                                              */
/*      Description     : Adjust decimal point and pack decimal */
/*                        digits.                               */
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
local a_btyp b_adpp(a_btyp **dc,a_intg *size,a_intg expo,a_intg dp,
                    a_intg length,a_intg *p_start,
                    a_intg *p_dp,a_intg *p_length)
#else
local a_btyp b_adpp(dc,size,expo,dp,length,p_start,p_dp,p_length)

a_btyp **dc;
a_intg *size;
a_intg expo;
a_intg dp;
a_intg length;
a_intg *p_start;
a_intg *p_dp;
a_intg *p_length;
#endif
        {
        a_intg i,j,k,m;
        char *buffer,*ptr;

        E_TPUSH("b_adpp")

        *p_start = *p_dp = *p_length = 0;

        /* number is non-zero   */
        if (dp!=0 || length!=0)
           {
           buffer = (char *)(*dc);

           /* move decimal point left to first digit    */
           if (0<(i = -dp-expo))
              {
              if (i+length>*size)
                 {
                 if ((ptr = (char*) malloc((size_t) (i+length)))==NULL)
                    return(ALLOCATION);

                 memcpy(ptr,buffer,(size_t)length);

                 /* if *size==BUFFERSIZE then static object assumed */
                 if (*size!=BUFFERSIZE)
                    {
#ifdef HEAP_CHECK
b_freh((a_char *)dc,(a_char *)buffer,(a_char *)"b_adpp");
#endif
                    free(buffer);
                    }
#ifdef HEAP_CHECK
b_geth((a_char *)dc,(a_char *)ptr,(a_char *)"b_adpp");
#endif
                 buffer = ptr;
                 *dc = (a_btyp *)ptr;
                 *size = i+length;
                 }
              for (k=length-1;k>=0;k--) buffer[k+i] = buffer[k];
              for (k=0;k<i;k++) buffer[k] = 0x00;
              length += i;
              dp = 0;
              }

           /* move decimal point right to last digit    */
           else if (-i>length)
              {
              if (-i>*size)
                 {
                 if ((ptr = (char*) malloc((size_t) -i))==NULL)
                    return(ALLOCATION);

                 memcpy(ptr,buffer,(size_t) *size);

                 /* if *size==BUFFERSIZE then static object assumed */
                 if (*size!=BUFFERSIZE)
                    {
#ifdef HEAP_CHECK
b_freh((a_char *)dc,(a_char *)buffer,(a_char *)"b_adpp");
#endif
                    free(buffer);
                    }
#ifdef HEAP_CHECK
b_geth((a_char *)dc,(a_char *)ptr,(a_char *)"b_adpp");
#endif
                 buffer = ptr;
                 *dc = (a_btyp *)ptr;
                 *size = -i;
                 }
              while (length<-i) buffer[length++] = 0x00;
              dp = length;
              }

           /* move decimal point within mantissa digits */
           else
              dp += expo;

           /* p_length = maximum number of a_btyp for packed number */
#if CRAY_UNIX_C
           k = (dp) ? 2*dp/sizeof(a_btyp)+2 : 1;
           k += ((length!=dp) ? 2*(length-dp)/sizeof(a_btyp)+2 : 1);
#else
           k = (dp) ? dp/sizeof(a_btyp)+2 : 1;
           k += ((length!=dp) ? (length-dp)/sizeof(a_btyp)+2 : 1);
#endif

           /* check for sufficient buffer space */
           if (k*sizeof(a_btyp)>*size)
              {
              if ((ptr = (char*) malloc((size_t) k*sizeof(a_btyp)))==NULL)
                 return(ALLOCATION);

              memcpy(ptr,buffer,(size_t) length);

              /* if *size==BUFFERSIZE then static object assumed */
              if (*size!=BUFFERSIZE)
                 {
#ifdef HEAP_CHECK
b_freh((a_char *)dc,(a_char *)buffer,(a_char *)"b_adpp");
#endif
                 free(buffer);
                 }
#ifdef HEAP_CHECK
b_geth((a_char *)dc,(a_char *)ptr,(a_char *)"b_adpp");
#endif
              buffer = ptr;
              *dc = (a_btyp *)ptr;
              *size = k*sizeof(a_btyp);
              }
           *p_length = k;

           /* pack digits of fraction part */
           if ((i = j = (length-dp) % DEC_PACK)!=0)
              {
              (*dc)[--k] = buffer[length-j];
              while (--i>0)
                 (*dc)[k] = (*dc)[k]*10+buffer[length-i];
              i = j;
              while (i++<DEC_PACK) (*dc)[k] *= 10;
              }
           for (i=length-j;i>dp;i -= DEC_PACK)
              {
              (*dc)[--k] = buffer[i-DEC_PACK];
              for (m=DEC_PACK-1;m>0;m--)
                 (*dc)[k] = (*dc)[k]*10+buffer[i-m];
              }

           /* position of packed decimal point */
           *p_dp = k;

           /* pack digits of integer part */
           j = dp % DEC_PACK;
           for (i=dp;i>j;i -= DEC_PACK)
              {
              (*dc)[--k] = buffer[i-DEC_PACK];
              for (m=DEC_PACK-1;m>0;m--)
                 (*dc)[k] = (*dc)[k]*10+buffer[i-m];
              }
           if (j)
              {
              (*dc)[--k] = buffer[0];
              for (i=1;i<j;i++)
                 (*dc)[k] = (*dc)[k]*10+buffer[i];
              }

           *p_start = k;
           }

        E_TPOPP("b_adpp")
        return(NO_ERROR);
        }





