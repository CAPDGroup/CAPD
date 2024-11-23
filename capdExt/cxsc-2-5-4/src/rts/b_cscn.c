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

/* CVS $Id: b_cscn.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_cscn.c                              */
/*                                                              */
/*      Entries         : int b_cscn                            */
/*                        (device,buffer,expo,dp,length,sign,i) */
/*                        FILE *device;                         */
/*                        a_intg *expo,*dp,*length;             */
/*                        a_bool *sign;                         */
/*                        char *buffer;                         */
/*                        int i;                                */
/*                                                              */
/*      Arguments       : device - input device reading string  */
/*                        buffer - buffer to hold string        */
/*                        expo - exponent                       */
/*                        dp - digits before decimal point      */
/*                        length - number of digits read        */
/*                        sign - sign of number                 */
/*                        i - first character scanned           */
/*                                                              */
/*      Description     : Scan input device for real number     */
/*                                                              */
/*      Function value  : error code                            */
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
local int b_cscn(FILE *device,char *buffer,a_intg *expo,a_intg *dp,
                 a_intg *length,a_bool *sign,int i)
#else
local int b_cscn(device,buffer,expo,dp,length,sign,i)

FILE *device;
char *buffer;
a_intg *expo;
a_intg *dp;
a_intg *length;
a_bool *sign;
int i;
#endif
        {
        a_intg digits;
        int c;  /* !!! must be int !!! */

        *sign = FALSE;
        *length = *expo = *dp = digits = 0;
        c = i;

/*                                                                   */
/* skip leading blanks and control codes                             */
/*                                                                   */
        while (c==' ' || c==EOLN || c==TAB)
           if ( (c = getc(device)) ==EOF)
              {
              ungetc(c,device);
              return(1);
              }
/*                                                                   */
/* first character is not a valid character                          */
/*                                                                   */
        if ( (c<'0' || c>'9') && c!='-' && c!='+' )
           {
           ungetc(c,device);
           return(2);
           }
/*                                                                   */
/* sign of number                                                    */
/*                                                                   */
        if (c=='-')
           {
           if (*length>=BUFFERSIZE) return(4);
           buffer[(*length)++] = c;
           c = getc(device);
           *sign = TRUE;
           }
        else if (c=='+') c = getc(device);
/*                                                                   */
/* end of input (no digits found)                                    */
/*                                                                   */
        if (c==EOF)
           {
           ungetc(c,device);
           return(3);
           }
/*                                                                   */
/* scan digits of integer part                                       */
/*                                                                   */
        while (c>='0' && c<='9')
           {
           if (*length>=BUFFERSIZE) return(4);
           buffer[*length] = c;
           (*length)++;
           digits++;
           c = getc(device);
           }
        *dp = *length;
/*                                                                   */
/* no digits in mantissa                                             */
/*                                                                   */
        if (digits==0)
           {
           return(6);
           }
/*                                                                   */
/* read decimal point                                                */
/*                                                                   */
        if (c=='.')
           {
           if (*length>=BUFFERSIZE) return(4);
           buffer[(*length)++] = c;
           c = getc(device);
/*                                                                   */
/* no digits in fraction                                             */
/*                                                                   */
           if ( c<'0' || c>'9')
              {
              ungetc(c,device);
              return(6);
              }
/*                                                                   */
/* scan digits of fractional part                                    */
/*                                                                   */
           do
              {
              if (*length>=BUFFERSIZE) return(5);
              buffer[*length] = c;
              (*length)++;
              digits++;
              c = getc(device);
              }
           while ( c>='0' && c<='9');
/*                                                                   */
/* skip trailing blanks in buffer                                    */
/*                                                                   */
           while (*length>*dp && buffer[*length-1]=='\0') (*length)--;
           }

         digits = 1;
/*                                                                   */
/* skip letter of exponential part                                   */
/*                                                                   */
        if (c=='e' || c=='E' )
           {
           if (*length>=BUFFERSIZE) return(5);
           buffer[(*length)++] = c;
           if ( (c = getc(device)) ==EOF)
              {
              ungetc(c,device);
              return(7);
              }
/*                                                                   */
/* sign of exponent                                                  */
/*                                                                   */
           if (c=='-')
              {
              if (*length>=BUFFERSIZE) return(5);
              buffer[(*length)++] = c;
              digits = -1;
              c = getc(device);
              }
           else if (c=='+') c = getc(device);
/*                                                                   */
/* no digits in exponent                                             */
/*                                                                   */
           if (c<'0' || c>'9')
              {
              ungetc(c,device);
              return(8);
              }
/*                                                                   */
/* get value of exponent                                             */
/*                                                                   */
           while (c>='0' && c<='9')
              {
              if (digits>0)
                 {
                 if (*expo>DBL_MAX_10_EXP) return(9);
                 }
              else
                 {
                 if (*expo>-(DBL_MIN_10_EXP)) return(10);
                 }
              *expo = *expo*10+c-'0';
              if (*length>=BUFFERSIZE) return(5);
              buffer[(*length)++] = c;
              c = getc(device);
              }
           *expo *= digits;
           }

        ungetc(c,device);

        return(0);
        }





