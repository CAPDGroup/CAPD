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

/* CVS $Id: l_writ.c,v 1.21 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : l_writ.c                              */
/*                                                              */
/*      Entries         : void l_writ(desc,s,width,frac,rnd)    */
/*                        f_text *desc;                         */
/*                        INTERN *s;                            */
/*                        a_intg width;                         */
/*                        a_intg frac;                          */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : desc   - output device                */
/*                        s      - multiprecision value         */
/*                        width  - field width                  */
/*                        frac   - number of fraction digits    */
/*                        rnd    - rounding mode                */
/*                              -1 = round downwards            */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                                                              */
/*      Description     : Convert a multiprecision value        */
/*                        to a character string and write to    */
/*                        output device                         */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_btyp b_maxl;
#endif
#ifdef LINT_ARGS
local void l_writ(f_text *desc,multiprecision s,a_intg width,
                  a_intg frac,a_intg rnd)
#else
local void l_writ(desc,s,width,frac,rnd)

f_text *desc;
multiprecision s;
a_intg width;
a_intg frac;
a_intg rnd;
#endif
        {
        char *istr = NULL;
        char *fstr = NULL;
        char sticky;
        char round;
        char last;
        a_btyp *tmp;
        int i;
        int ilen = 0;
        int flen = 0;
        int decimal_exponent = 0;
        size_t isize;
        size_t fsize;
        size_t tmp_size;
        a_btyp idigits;
        a_btyp fdigits;
        a_btyp iminlen;
        a_btyp fminlen;
        a_btyp msb;
        a_btyp lsb;
        a_btyp rem;
        a_btyp digits;

#define DEFAULT_POW10 1000

        E_TPUSH("l_writ")

        /*   Number of mantissa digits required for integer conversion   */
        iminlen = (s->e >= 0) ? (a_btyp)s->e+1 : 0;

        /*   Number of mantissa digits required for fraction conversion   */
        fminlen = ((s->e>=0 && s->l>s->e+1) || s->e<0) ? s->l-s->e-1 : 0;

        /*   Number of mantissa digits stored for integer part   */
        idigits = (s->e >= 0 && (a_btyp)s->e+1 >= s->l) ? s->l : iminlen;

        /*   Number of mantissa digits stored for fraction part   */
        fdigits = s->l-idigits;

        /*   Size of result string for integer part.    */
        /*   One more for additional digit in rounding. */
        isize = 10*iminlen + 1;

        /*   Allocation of buffer for integer result   */
        if ((istr = (char *)malloc(isize))==NULL)
           {
           }

        /*   Conversion of integer part */
        if (s->z || s->e<0)
           {
           istr[0] = '0';
           ilen = 1;
           }
        else
           {

           /*   Size of temporary buffer for conversion of integer part.  */
           /*   Mantissa digits of integer part are splitted in 2 parts.  */
           iminlen *= 2;
           tmp_size = iminlen*sizeof(a_btyp);

           /*   Allocation of temporary buffer    */
           if ((tmp = (a_btyp *)malloc(tmp_size))==NULL)
              {
              }

           /*   Copy and split integer mantissa to temporary buffer   */
           for (i=0; i<idigits; i++)
             {
             tmp[2*i] = s->m[i]>>16;
             tmp[2*i+1] = s->m[i] & 0x0000ffff;
             }

           /*   Fill trailing buffer positions with zero   */
           for (i=2*idigits; i<iminlen; i++)
              {
              tmp[i] = 0x00000000;
              }

           /*   Conversion of mantissa digits of integer part     */
           /*   is done by repeated division.                     */
           /*   Resulting output characters are stored in reverse */
           /*   in order in result string istr.                   */

           /*   Index of temporary digit that holds most sigificant bit   */
           msb = 0;

           /*   As long as there are digits perform conversion   */
           ilen = 0;
           while (msb<iminlen)
              {

              /*   Initialize remainder   */
              rem = 0;

              /*   Divide all digits in temporary buffer   */
              for (i=msb; i<iminlen; i++)
                 {
                 rem = (tmp[i] += rem<<16) % 10000;
                 tmp[i] /= 10000;
                 }

              /*   Determine result characters   */
              for (i=0; i<4 && ilen<isize; i++)
                 {
                 istr[ilen++] = (rem % 10) + '0';
                 rem /= 10;
                 }

              /*   Adjust index of most significant bit   */
              if (tmp[msb]==0x00000000L) msb++;
              }
           
           /*    Ignore leading zeros in result string   */
           while (istr[ilen-1]=='0') ilen--;

           /*    Free temporary buffer    */
           free((char *)tmp);
           }

/*---------------------------------------------------------------------------*/

        /*   Set default field width.                                 */
        /*   The value "digits" holds the number of required digits   */
        /*   including 2 digits for exact rounding.                   */
        /*   Decimal exponent set to non-negative value if integer    */
        /*   part is non-zero.                                        */
        if (width<0) width = -width;
        if (frac==0)
           {
           decimal_exponent = (ilen==1 && istr[0]=='0' && !s->z) ? -1 : ilen-1;
           if (width<=0) width = ExpDigits + 4 + (963*b_maxl/100+1);
		 else if (width<9) width = 9;
           digits = width-7+2;
           if (s->e<0 && !s->z) 
              {
              digits -= 10*s->e;
              }
           }
        else
           {
           if (width<ilen+frac+1) width = ilen+frac+1;
           digits = ilen+frac+2;
           }

/*---------------------------------------------------------------------------*/

        /*   Size of result string for fraction part */
        fsize = (ilen >= digits) ? 0 : digits-ilen;

        /*   Fraction digits are required   */
        if (fsize>0) 
           {

           /*   Allocation of buffer for fraction result   */
           if ((fstr = (char *)malloc(fsize))==NULL)
              {
              }

           /*   Conversion of fraction part if there exist fraction digits.   */
           if (fdigits>0)
              {

              /*   Size of temporary buffer for conversion of fraction part. */
              /*   Mantissa digits of fraction part are splitted in 2 parts. */
              fminlen = 2*fminlen;
              tmp_size = fminlen*sizeof(a_btyp);

              /*   Allocation of temporary buffer    */
              if ((tmp = (a_btyp *)malloc(tmp_size))==NULL)
                 {
                 }

              /*   Initialize leading buffer digits with zero   */
              for (i= -2*s->e-3; i>=0; i--)
                 {
                 tmp[i] = 0x00000000;
                 }

              /*   Copy and split integer mantissa to temporary buffer   */
              for (i=0; i<fdigits; i++)
                 {
                 tmp[2*(i-fdigits)+fminlen] = s->m[idigits+i]>>16;
                 tmp[2*(i-fdigits)+1+fminlen] = s->m[idigits+i] & 0x0000ffff;
                 }

              /*   Conversion of mantissa digits of fraction part    */
              /*   is done by repeated multiplication.               */
              /*   Resulting output characters are stored in correct */
              /*   order in result string fstr.                      */

              /*   Index of temporary digit that holds least sigificant bit   */
              lsb = fminlen-1;

              /*   As long as digits are required perform conversion   */
              flen = 0;
              while (flen<fsize)
                 {

                 /*   Initialize remainder   */
                 rem = 0;

                 /*   Multiply all digits in temporary buffer   */
                 for (i=lsb; i>=0; i--)
                    {
                    tmp[i] = tmp[i]*10000+rem;
                    rem = tmp[i]>>16;
                    tmp[i] &= 0x0000ffff;
                    }

                 /*   Determine result characters   */
                 fstr[flen++] = (rem / 1000) + '0';
                 rem %= 1000;
                 if (flen<fsize) 
                    fstr[flen++] = (rem / 100) + '0';
                 else
                    fstr[flen-1] |= rem / 100;
                 rem %= 100;
                 if (flen<fsize) 
                    fstr[flen++] = (rem / 10) + '0';
                 else
                    fstr[flen-1] |= rem / 10;
                 rem %= 10;
                 if (flen<fsize) 
                    fstr[flen++] = rem + '0';
                 else
                    fstr[flen-1] |= rem;

                 /*   Adjust index of least significant bit   */
                 if (tmp[lsb]==0x00000000) lsb--;
                 }
           
              /*   Check for non-zero fraction digits   */
              for (i=lsb; i>=0; i--)
                 {
                 if (tmp[i])
                    {
                    fstr[flen-1] |= 1;
                    break;
                    }
                 }

              /*    Free temporary buffer    */
              free((char *)tmp);
              }

           /*   Set all fraction digits to zero.   */
           else
              {
              for (flen=0; flen<fsize; flen++)
                 {
                 fstr[flen] = '0';
                 }
              }

           }

        /*    Fraction digits are not required   */
        else
           {

           /*   Set sticky bit information in integer result   */
           /*   if there exist fraction digits.                */
           if (fdigits>0)
              {
              istr[ilen-1] |= 1;
              }

           }

        /*   Determine decimal exponent    */
        if (frac==0 && decimal_exponent<0)
           {
           for (i=0; i<flen; i++)
              {
              if (fstr[i]!='0')
                 {
                 decimal_exponent = -i-1;
                 break;
                 }
              }
           digits = -decimal_exponent+width-7+2;
           }

/*---------------------------------------------------------------------------*/

        /*    Rounding of result string   */
        if (ilen>=digits)
           {

           /*    Gather sticky bit information of integer digits    */
           sticky = istr[ilen-digits];
           round = istr[ilen-digits+1];
           last = istr[ilen-digits+2];
           for (i=0; i<ilen-digits && sticky=='0'; i++)
              {
              sticky |= istr[i];
              }

           /*   Check for increment required   */
           if ((rnd==0 && (round>'5' ||
                           (round=='5' && 
                            (sticky!='0' ||
                             (sticky=='0' && 
                              (last & 1)
                          ))))
               ) ||
               (((rnd>0 && s->s==0) || (rnd<0 && s->s!=0)) && 
                (round | sticky)!='0'
              ))
              {
              for (i=ilen-digits+2; i<ilen; i++)
                 {
                 if (istr[i]=='9')
                    {
                    istr[i] = '0';
                    }
                 else
                    {
                    istr[i]++;
                    break;
                    }
                 }
              if (i==ilen)
                 {
                 istr[i] = '1';
                 ilen++;
                 decimal_exponent++;
                 }
              }
           }
        else
           {

           /*   Gather sticky bit information of fraction digits   */
           sticky = fstr[digits-ilen-1];
           round = (digits>=ilen+2) ? fstr[digits-ilen-2]:istr[ilen+1-digits];
           last  = (digits>=ilen+3) ? fstr[digits-ilen-3]:istr[ilen+2-digits];
           for (i=digits-ilen; i<flen && sticky=='0'; i++)
              {
              sticky |= fstr[i];
              }

           /*   Check for increment required   */
           if ((rnd==0 && (round>'5' ||
                           (round=='5' && 
                            (sticky!='0' ||
                             (sticky=='0' && 
                              (last & 1)
                          ))))
               ) ||
               (((rnd>0 && s->s==0) || (rnd<0 && s->s!=0)) && 
                (sticky | round)!='0'
              ))
              {
              for (i=digits-ilen-3; i>=0; i--)
                 {
                 if (fstr[i]=='9')
                    {
                    fstr[i] = '0';
                    }
                 else
                    {
                    fstr[i]++;
                    break;
                    }
                 }
              if (i<0)
                 {
                 for (i=0; i<ilen; i++)
                    {
                    if (istr[i]=='9')
                       {
                       istr[i] = '0';
                       }
                    else
                       {
                       istr[i]++;
                       break;
                       }
                    }
                 if (i==ilen)
                    {
                    istr[i] = '1';
                    ilen++;
                    decimal_exponent++;
                    }
                 }
              }
           }

/*---------------------------------------------------------------------------*/

        /*    Generating floating-point number string     */
        if (frac==0)
           {
           f_putc((s->s) ? '-' : ' ',desc);
           if (decimal_exponent>=0)
              f_putc(istr[decimal_exponent],desc);
           else
              {
              f_putc(fstr[-decimal_exponent-1],desc);

              /*    Adjust number of digits for leading zeros are not used   */
              digits += decimal_exponent;
              }
           f_putc('.',desc);
           for (i=1; i<digits-2; i++)
              {
              if (decimal_exponent-i>=0)
                 f_putc(istr[decimal_exponent-i],desc);
              else
                 f_putc(fstr[i-decimal_exponent-1],desc);
              }
           f_putc('E',desc);
           if (decimal_exponent>=0)
              {
              f_putc('+',desc);
              }
           else
              {
              decimal_exponent = -decimal_exponent;
              f_putc('-',desc);
              }
           for (i=DEFAULT_POW10; i<=decimal_exponent; i *= 10)
              {
              /*   Do nothing   */
              }
           do
              {
              i /= 10;
              f_putc(decimal_exponent/i+'0',desc);
              decimal_exponent %= i;
              }
           while (i>1);
           }

        /*    Generating fixed-point number string     */
        else
           {
           for (i=ilen+frac+3; i<=width; i++)
              {
              f_putc(' ',desc);
              }
           if (s->s) f_putc('-',desc);
           else if (ilen+frac+2<=width) f_putc(' ',desc);
           for (i=ilen-1; i>=0; i--)
              f_putc(istr[i],desc);
           f_putc('.',desc);
           for (i=0; i<frac; i++)
              f_putc(fstr[i],desc);
           }

/*---------------------------------------------------------------------------*/

        /*    Free buffer for fraction result    */
        if (fsize>0) free((char *)fstr);

        /*   free buffer for integer result    */
        free((char *)istr);

        E_TPOPP("l_writ")
        return;
        }





