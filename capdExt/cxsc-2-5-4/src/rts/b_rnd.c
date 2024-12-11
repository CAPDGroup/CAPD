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

/* CVS $Id: b_rnd.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_rnd.c                               */
/*                                                              */
/*      Entries         : void b_rnd                            */
/*                         (rnd,buffer,digits,pos,bdp,dexpo)    */
/*                        a_intg *dexpo,rnd,digits,pos,*bdp;    */
/*                        char *buffer;                         */
/*                                                              */
/*      Arguments       : rnd - rounding mode for unsigned      */
/*                        buffer - mantissa                     */
/*                        digits - number of mantissa digits    */
/*                        pos - position where to round         */
/*                        bdp - position of decimal point       */
/*                        dexpo - exponent of result            */
/*                                                              */
/*      Description     : Round decimal representation of number*/
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
local void b_rnd(a_intg rnd,char *buffer,a_intg digits,a_intg pos,
                 a_intg *bdp,a_intg *dexpo)
#else
local void b_rnd(rnd,buffer,digits,pos,bdp,dexpo)

a_intg rnd;
char *buffer;
a_intg digits;
a_intg pos;
a_intg *bdp;
a_intg *dexpo;
#endif
        {
        char *s,*q;

        s = buffer+(*bdp-*dexpo);

        if (rnd>0) {
                for (q=s+pos;q<s+digits;q++)
                   if (*q!='0') break;
                if (q==s+digits) return;
                }
        else if (rnd==0) {
                if (s[pos]<'5') return;
                if (s[pos]=='5')
                   {
                   for (q=s+pos+1;q<s+digits;q++)
                      if (*q!='0') break;
                   if (q>=s+digits && (s[pos-1] & 1)==0) return;
                   }
                }
        else return;

        /* advance to next decimal number                       */
        for (q=s+pos-1;q>=s;q--)
                if (*q=='9') *q = '0';
                else {
                    (*q)++;
                    break;
                    }

        /* carry occurs                                         */
        if (q<s) {
                (*dexpo)++;
                *q = '1';
                }

        return;
        }





