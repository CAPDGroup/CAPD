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

/* CVS $Id: r_writ.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_writ.c                              */
/*                                                              */
/*      Entries         : void r_writ                           */
/*                         (device,s,TotalWidth,FracDigits,rnd) */
/*                        FILE *device;                         */
/*                        a_real s;                             */
/*                        a_intg rnd,TotalWidth,FracDigits;     */
/*                                                              */
/*      Arguments       : device - output device writing string */
/*                        s - IEEE value                        */
/*                        Totalwidth - total length of string   */
/*                        FracDigits - number of fraction digits*/
/*                        rnd - rounding mode                   */
/*                              -1 = round downwards            */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                                                              */
/*      Description     : Convert an IEEE double format number  */
/*                        to a character string and write to    */
/*                        output device                         */
/*                                                              */
/*                   TotalWidth-2<FracDigits considered.        */
/*                   No special consideration of TotalWidth=0.  */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern dotprecision b_cp__;
#endif

#ifdef LINT_ARGS
local void r_writ(FILE *device,a_real s,a_intg TotalWidth,a_intg FracDigits,
                  a_intg rnd)
#else
local void r_writ(device,s,TotalWidth,FracDigits,rnd)

FILE *device;
a_real s;
a_intg TotalWidth;
a_intg FracDigits;
a_intg rnd;
#endif
        {
        a_intg length;
        int k,sign;     /* !!! must be int !!! */
        char *buffer;

        E_TPUSH("r_writ")

        if (TotalWidth<0)
           {
           sign = 1;
           TotalWidth = -TotalWidth;
           }
        else
           sign = 0;

        if (FracDigits<0) FracDigits = -FracDigits;

/* Set default TotalWidth according to data format                   */
        if (FracDigits==0)
           {
/* --- No special consideration of case TotalWidth=0 ---
           if (TotalWidth==0)
              TotalWidth = ExpDigits+4+(301*MANTL/1000+1);
           else
  ---                                       --- */
           if (TotalWidth<9)
              TotalWidth = 9;
           }
        else if (TotalWidth-2<FracDigits)
           {
           TotalWidth = 2+FracDigits;
           }

        buffer = (char *)b_cp__;

/* generate output string                                            */
        r_outp(buffer,s,TotalWidth,FracDigits,rnd,&length);

/*                                                                   */
/* output of string                                                  */
/*                                                                   */
        if (sign)
           {
           for (k=0;buffer[k]==' ';k++) { /* empty */ }
           if (k>0 && buffer[k]!='-') k--;
           sign = k;
           while (k<(int)length)
              fputc(buffer[k++],device);
           for (k=0;k<sign;k++)
              fputc(' ',device);
           }
        else
           {
           for (k=0;k<(int)length;k++)
              fputc(buffer[k],device);
           }

        E_TPOPP("r_writ")
        return;
        }





