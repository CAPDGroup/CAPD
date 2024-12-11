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

/* CVS $Id: b_outi.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_outi.c                              */
/*                                                              */
/*      Entries         : void b_outi(digits,buffer,bdp,dexpo,c)*/
/*                        char *buffer;                         */
/*                        a_intg *digits,*bdp,*dexpo;           */
/*                        dotprecision c;                       */
/*                                                              */
/*      Arguments       : digits  - number of required digits   */
/*                        buffer  - output buffer               */
/*                        bdp - position of decimal point       */
/*                        dexpo   - exponent of first non-zero  */
/*                                  digit                       */
/*                        c       - accu holding IEEE value     */
/*                                                              */
/*      Description     : Convert integer part of IEEE value    */
/*                        to a character string                 */
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
local void b_outi(a_intg *digits,char *buffer,a_intg *bdp,a_intg *dexpo,
                  dotprecision c)
#else
local void b_outi(digits,buffer,bdp,dexpo,c)

a_intg *digits;
char *buffer;
a_intg *bdp;
a_intg *dexpo;
dotprecision c;
#endif
        {
        a_intg i,mod;
        a_btyp *s,*p,h,hh;
        char *q;

        /* initialize                                   */
        q = buffer+*bdp;
        s = &c[c[A_BEGIN]];

        /* convert by repeated division                 */
        while (s<=&c[A_D_P]) {
           mod = 0;
           for (p=s;p<=&c[A_D_P];p++)
               {
               h = GETHIGH(*p) | MOVEHIGH(mod);
               hh = GETLOW(*p) | MOVEHIGH(h%B2D_POWER);
               mod = hh%B2D_POWER;
               *p = MOVEHIGH(h/B2D_POWER) | (hh/B2D_POWER);
               }

           /* store decimal digits                              */
           for (i=0;i<B2D_LOG10-1;i++) {
               *q-- = (char) (mod%10+'0');
               mod /= 10;
               }
           *q-- = (char) (mod+'0');

           /* reduce length of binary mantissa                  */
           if (*s==ZERO) s++;
           }

        /* ignore leading zeros */
        while (*(++q)=='0') { /* empty */ }

        /* decimal exponent of first non-zero digit     */
        *dexpo = (buffer+*bdp)-q;

        if (*digits>*dexpo+1) {
           *digits -= (*dexpo+1);
           }
        else
           {                                                 /*(1)*/
           *digits = 0;
           if (!b_test((a_intg)c[A_END]-A_D_P,c+(A_D_P+1))) {/*(1)*/
                 buffer[*bdp] = '1';                         /*(1)*/
                 }                                           /*(1)*/
           }                                                 /*(1)*/

        return;
        }





