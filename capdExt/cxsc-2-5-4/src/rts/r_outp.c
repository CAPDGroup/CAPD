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

/* CVS $Id: r_outp.c,v 1.22 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_outp.c                              */
/*                                                              */
/*      Entries         : void r_outp                           */
/*                         (buffer,s,Totalwidth,                */
/*                          FracDigits,rnd,length)              */
/*                        char *buffer;                         */
/*                        a_real s;                             */
/*                        a_intg rnd,Totalwidth,FracDigits;     */
/*                        a_intg *length;                       */
/*                                                              */
/*      Arguments       : buffer - output buffer holding string */
/*                        s - IEEE value                        */
/*                        Totalwidth - total length of string   */
/*                        FracDigits - number of fraction digits*/
/*                        rnd - rounding mode                   */
/*                              -1 = round downwards            */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                        length - size of buffer string        */
/*                                                              */
/*      Description     : Convert an IEEE double format number  */
/*                        to a character string.                */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern char *o_text[];
#ifdef SIGNED_ENABLE
extern a_bool f_ppsz;
#endif
#endif

#ifdef LINT_ARGS
local void r_outp(char *buffer,a_real s,a_intg TotalWidth,a_intg FracDigits,
                  a_intg rnd,a_intg *length)
#else
local void r_outp(buffer,s,TotalWidth,FracDigits,rnd,length)

char *buffer;
a_real s;
a_intg TotalWidth;
a_intg FracDigits;
a_intg rnd;
a_intg *length;
#endif
        {
        a_intg ActWidth,DecPlaces,expo,IntDigits,MinNumChars;
        a_intg dexpo,digits,blanks,bdp;
        a_intg k;
        a_bool vz,zero;
        a_btyp mant[BSIZE];

        E_TPUSH("r_outp")

        blanks = *length = 0;

        zero = b_deko(s,&expo,mant,&vz);
/*                                                                   */
/* infinity or NaN                                                   */
/*                                                                   */
        if (expo>EXPO_MAX)
           {
           *length = TotalWidth;

           /* Infinity                                  */
           if (MANT_INFINITY(mant))
              {
              for (k=0;k<TotalWidth-strlen(o_text[36])-1;k++)
                 buffer[k] = ' ';
              buffer[k++] = (vz) ? MINUS_SIGN : PLUS_SIGN;
              memcpy(buffer+k,o_text[36],(size_t)(TotalWidth-k));
              }

           /* NaN                                       */
           else if (SIGNALING(mant[0]))
              {
              for (k=0;k<TotalWidth-strlen(o_text[37]);k++) 
                 buffer[k]=' ';
              memcpy(buffer+k,o_text[37],(size_t)(TotalWidth-k));
              }
           else
              {
              for (k=0;k<TotalWidth-strlen(o_text[38]);k++) 
                 buffer[k]=' ';
              memcpy(buffer+k,o_text[38],(size_t)(TotalWidth-k));
              }

           E_TPOPP("r_outp")
           return;
           }

        for (k=D_U_RATIO;k<BSIZE;k++) mant[k] = ZERO;

        if (vz) rnd = -rnd;

        if (expo>800) bdp = (BUFFERSIZE-8)-FracDigits;
        else if (expo<-800) bdp = 8;
        else bdp = B_D_P;
/*                                                                   */
/* restrict TotalWidth to BUFFERSIZE-8                               */
/*                                                                   */
        if (TotalWidth>BUFFERSIZE-8) TotalWidth = BUFFERSIZE-8;
/*                                                                   */
/* floating-point representation                                     */
/*                                                                   */
        if (FracDigits==0)
           {
/*                                                                   */
/* ActWidth                                                          */
/*                                                                   */
           ActWidth = (TotalWidth>=ExpDigits+6) ? TotalWidth : ExpDigits+6;
           DecPlaces = ActWidth-ExpDigits-5;
/*                                                                   */
/* number is zero                                                    */
/*                                                                   */
           if (zero)
              {
              *buffer++ =
#ifdef SIGNED_ENABLE
                          (vz && f_ppsz ) ? MINUS_SIGN :
#endif
                                                  ' ';
              *buffer++ = '0';
              *buffer++ = XSC_DECIMAL_POINT;
              for (k=0;k<DecPlaces;k++) *buffer++ = '0';
              *buffer++ = EXPONENT_E;
              *buffer++ = PLUS_SIGN;
              for (k=0;k<ExpDigits;k++) *buffer++ = '0';

              *length = 5+DecPlaces+ExpDigits;

              E_TPOPP("r_outp")
              return;
              }
/*                                                                   */
/* determine output string                                           */
/*                                                                   */
           digits = DecPlaces+(1+2);
           dexpo = -1;
           b_out(mant,expo,digits,buffer,&bdp,&dexpo);

           if (dexpo>0 && dexpo>DecPlaces+2) digits = dexpo+1;

           b_rnd(rnd,buffer,digits,DecPlaces+1,&bdp,&dexpo);

           buffer[bdp-dexpo-1] = buffer[bdp-dexpo];
           buffer[bdp-dexpo] = XSC_DECIMAL_POINT;
           buffer[bdp-dexpo+DecPlaces+1] = EXPONENT_E;
           buffer[bdp-dexpo+DecPlaces+2] =
              (dexpo<0) ? MINUS_SIGN : PLUS_SIGN;
           expo = (dexpo<0) ? -dexpo : dexpo;
           for (k=ExpDigits;k>0;k--) {
              buffer[bdp-dexpo+DecPlaces+2+k] = (char) (expo%10+'0');
              expo /= 10;
              }
           *length = 3+DecPlaces+2+ExpDigits;

           buffer[bdp-2-dexpo] = (vz) ? MINUS_SIGN : ' ';

           for (k=0;k<*length;k++)
              buffer[k] = buffer[bdp-2-dexpo+k];
           }
        else
           {
/*                                                                   */
/* fixed-point representation                                        */
/*                                                                   */

/*                                                                   */
/* number is zero                                                    */
/*                                                                   */
           if (zero)
              {
              *length = TotalWidth;
              for (k=TotalWidth-FracDigits-3;k>0;k--) *buffer++ = ' ';
              if (k==0)
                 {
                 *buffer++ =
#ifdef SIGNED_ENABLE
                          (vz && f_ppsz) ? MINUS_SIGN :
#endif
                                                 ' ';
                 }
#ifdef SIGNED_ENABLE
              else if (vz && f_ppsz)
                 {
                 *buffer++ = MINUS_SIGN;
                 (*length)++;
                 }
#endif
              *buffer++ = '0';
              *buffer++ = XSC_DECIMAL_POINT;
              for (k=0;k<FracDigits;k++) *buffer++ = '0';
              E_TPOPP("r_outp")
              return;
              }

           /* estimate number of decimal digits */
           if (expo>=0) {
               IntDigits = ((expo+1)*61)/200+1;
               }
           else {
               IntDigits = 0;
               dexpo = 0;
               }

           /* fill fractional part with zeros                   */
           for (k=0;k<=FracDigits+2;k++) buffer[bdp+k] = '0';

           digits = IntDigits+FracDigits+2;
           b_out(mant,expo,digits,buffer,&bdp,&dexpo);

           /* correct setting of IntDigits */
           if (expo>=0) {
               IntDigits = dexpo+1;
               }
           else {
               IntDigits = 1;
               digits++;
               }

           b_rnd(rnd,buffer,digits,IntDigits+FracDigits,&bdp,&dexpo);

           /* correct setting of IntDigits after rounding */
           IntDigits = dexpo+1;

           /* value is zero after rounding;                     */
           /* sign is changed if f_ppsz is FALSE                */
           if (vz
#ifdef SIGNED_ENABLE
                  && NOT(f_ppsz)
#endif
                                 ) {
              for (k=(bdp+1)-IntDigits;k<(bdp+1)+FracDigits;k++)
                 if (buffer[k]!='0') break;
              if (k==(bdp+1)+FracDigits) vz = FALSE;
              }

           MinNumChars = IntDigits+FracDigits+1;
           if (vz) MinNumChars++;

           if ( (blanks = TotalWidth-MinNumChars)<=0 )
              blanks = 0;
           else
              for (k=0;k<blanks;k++) buffer[k] = ' ';

           /* generate decimal point                            */
           for (k=bdp-IntDigits;k<bdp;k++) buffer[k] = buffer[k+1];
           buffer[bdp] = XSC_DECIMAL_POINT;

           *length = MinNumChars;
/*                                                                   */
/* sign of number                                                    */
/*                                                                   */
           buffer[bdp-2-dexpo] = (vz) ? MINUS_SIGN : ' ';

           for (k=0;k<*length;k++)
              buffer[blanks+k] = buffer[bdp-1-vz-dexpo+k];
           *length += blanks;
           }

        E_TPOPP("r_outp")
        return;
        }





