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

/* CVS $Id: b_outf.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_outf.c                              */
/*                                                              */
/*      Entries         : void b_outf(digits,buffer,bdp,dexpo,c)*/
/*                        char *buffer;                         */
/*                        dotprecision c;                       */
/*                        a_intg *digits,*bdp,*dexpo;           */
/*                                                              */
/*      Arguments       : digits  - number of required digits   */
/*                        buffer  - output buffer               */
/*                        bdp - position of decimal point       */
/*                        dexpo   - exponent of first non-zero  */
/*                                  digit                       */
/*                        c       - accu holding IEEE value     */
/*                                                              */
/*      Description     : Convert fractional part of an IEEE    */
/*                        value to a character string           */
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
local void b_outf(a_intg *digits,char *buffer,a_intg *bdp,a_intg *dexpo,
                  dotprecision c)
#else
local void b_outf(digits,buffer,bdp,dexpo,c)

a_intg *digits;
char *buffer;
a_intg *bdp;
a_intg *dexpo;
dotprecision c;
#endif
        {
        a_intg i,carry,cont,dg;
        a_btyp *p,*s,h,hh;
        char *q;

        /* initialize                                   */
        dg = *digits;
        cont = *dexpo;
        q = buffer+(*bdp+1);
        s = &c[c[A_END]];
        while (s>&c[A_D_P] && *s==ZERO) s--;

        /* convert by repeated multiplication           */
        while (dg>0 && s>&c[A_D_P])
           {
           carry = 0;
           for (p=s;p>&c[A_D_P];p--)
                {
                hh = GETLOW(*p)*B2D_POWER+carry;
                h = GETHIGH(*p)*B2D_POWER+GETHIGH(hh);
                carry = GETHIGH(h);
                *p = MOVEHIGH(h) | GETLOW(hh);
                }

           /* store decimal digits                              */
           for (i=B2D_LOG10-1;i>0;i--) {
               q[i] = (char) (carry%10+'0');
               carry /= 10;
               }
           *q = (char) (carry+'0');

           /* adjust dexpo and digits                           */
           for (i=0;i<B2D_LOG10;i++) {
               if (cont<0) {
                   if (*q!='0') {
                       cont = 0;
                       *dexpo = (buffer+*bdp)-q;
                       dg--;
                       }
                   }
               else dg--;
               q++;
               }

           if (*s==ZERO) s--;
           }

        /* pad with zeros if fraction part is zero remainder */
        if (dg>0) {
           while (dg--) {
              *q++ = '0';
              }
           }

        /* set least significant digit to 1 if non-zero remainder */
        else {
           if (s!=&c[A_D_P]) q[dg-1] = '1';                 /*(1)*/
           }

   return;
   }





