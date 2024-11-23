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

/* CVS $Id: b_chck.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_chck.c                              */
/*                                                              */
/*      Entries         : int b_chck                            */
/*                          (strng,buffer,size,expo,dp,         */
/*                           length,sign,next)                  */
/*                        a_char *strng;                        */
/*                        char **buffer;                        */
/*                        a_intg *size;                         */
/*                        a_intg *expo;                         */
/*                        a_intg *dp;                           */
/*                        a_intg *length;                       */
/*                        a_bool *sign;                         */
/*                        a_char **next;                        */
/*                                                              */
/*      Arguments       : strng  - input string                 */
/*                        buffer - buffer to hold string        */
/*                        size   - size of buffer               */
/*                        expo   - exponent                     */
/*                        dp     - digits before decimal point  */
/*                        length - number of digits read        */
/*                        sign   - sign of number               */
/*                        next   - address of last character    */
/*                                                              */
/*      Description     : Scan input string for real number     */
/*                        with exponent range restricted to     */
/*                        -9999 <= expo <= 9999.                */
/*                                                              */
/*      Function value  : 0 = real number scanned               */
/*                        1 = allocation error                  */
/*                        2 = no digit in mantissa              */
/*                        3 = no digit in fraction              */
/*                        4 = no digit in exponent              */
/*                        5 = restricted exponent range         */
/*                                                              */
/*      Note            : Variable 'buffer' points to a         */
/*                        d_otpr field (b_cp__) which cannot    */
/*                        be freed. 'buffer' must be alligned   */
/*                        according to type a_btyp.             */
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
local int b_chck(a_char *strng,char **buffer,a_intg *size,
                 a_intg *expo,a_intg *dp,
                 a_intg *length,a_bool *sign,a_char **next)
#else
local int b_chck(strng,buffer,size,expo,dp,length,sign,next)

a_char *strng;
char **buffer;
a_intg *size;
a_intg *expo;
a_intg *dp;
a_intg *length;
a_bool *sign;
a_char **next;
#endif
        {
        a_intg esign;
        char *ptr;

        *sign = FALSE;
        *length = *expo = *dp = 0;

        /* skip leading blanks and control codes        */
        while (*strng==' ' || *strng==EOLN) strng++;

        /* sign of number                               */
        if (*strng=='-' || *strng=='+')
           {
           *sign = (*strng=='-') ? TRUE : FALSE;
           strng++;
           }
        else
           *sign = FALSE;

        /* first character of mantissa must be a digit */
        if (!isdigit(*strng))
           {
           *next = strng;
           return(2);
           }

        /* skip leading zeros                           */
        while (*strng=='0') strng++;

        /* scan digits of integer part                  */
        while (isdigit(*strng))
           {
           if (*length>=*size)
              {
              if ((ptr = (char*) malloc((size_t)(*size+BUFFERSIZE)))==NULL)
                 {
                 *next = strng;
                 return(1);
                 }

              memcpy(ptr,*buffer,(size_t) *size);

              if (*size!=BUFFERSIZE)
                 {
#ifdef HEAP_CHECK
b_freh((a_char *)buffer,(a_char *)*buffer,(a_char *)"b_chck");
#endif
                 free(*buffer);
                 }
#ifdef HEAP_CHECK
b_geth((a_char *)buffer,(a_char *)ptr,(a_char *)"b_chck");
#endif
              *buffer = ptr;
              *size += BUFFERSIZE;
              }
           (*buffer)[*length] = *strng-'0';
           (*length)++;
           strng++;
           }

        /* save position of decimal point */
        *dp = *length;

        /* skip decimal point                           */
        if (*strng=='.')
           {
           strng++;

           /* first character of fraction must be a digit */
           if (!isdigit(*strng))
              {
              *next = strng;
              return(3);
              }

           /* scan digits of fractional part            */
           while (isdigit(*strng))
              {
              if (*length>=*size)
                 {
                 if ((ptr = (char*) malloc((size_t)(*size+BUFFERSIZE)))==NULL)
                    {
                    *next = strng;
                    return(1);
                    }

                 memcpy(ptr,*buffer,(size_t) *size);

                 if (*size!=BUFFERSIZE)
                    {
#ifdef HEAP_CHECK
b_freh((a_char *)buffer,(a_char *)*buffer,(a_char *)"b_chck");
#endif
                    free(*buffer);
                    }
#ifdef HEAP_CHECK
b_geth((a_char *)buffer,(a_char *)ptr,(a_char *)"b_chck");
#endif
                 *buffer = ptr;
                 *size += BUFFERSIZE;
                 }
              (*buffer)[*length] = *strng-'0';
              (*length)++;
              strng++;
              }

           /* skip trailing blanks of fraction in buffer */
           while (*length>*dp && (*buffer)[*length-1]=='\0') (*length)--;
           }

        /* skip letter of exponential part              */
        if (*strng=='e' || *strng=='E' )
           {
           strng++;

           /* sign of exponent                          */
           if (*strng=='-' || *strng=='+')
              {
              esign = (*strng=='-') ? -1 : 1;
              strng++;
              }
           else esign = 1;

           /* no digits in exponent                     */
           if (!isdigit(*strng))
              {
              *next = strng;
              return(4);
              }

           /* get value of exponent                     */
           while (isdigit(*strng) && *expo<1000)
              {
              *expo = *expo*10+*strng-'0';
              strng++;
              }
           *expo *= esign;
           }

        *next = strng;

        return(isdigit(*strng) ? 5 : 0);
        }





