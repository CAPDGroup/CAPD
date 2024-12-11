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

/* CVS $Id: s_int_.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_int_.c                              */
/*                                                              */
/*      Entries         : s_trng s_int_(s,Totalwidth)           */
/*                        a_intg s,Totalwidth;                  */
/*                                                              */
/*      Arguments       : s - integer value                     */
/*                        Totalwidth - total length of string   */
/*                                                              */
/*      Description     : Convert an integer value              */
/*                        to a character string.                */
/*                                                              */
/*                   and size_t at length calculation           */
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
local s_trng s_int_(a_intg s,a_intg TotalWidth)
#else
local s_trng s_int_(s,TotalWidth)

a_intg s;
a_intg TotalWidth;
#endif
        {
        a_intg i,IntDigits;
        a_bool sign;
        char *buffer;
        a_btyp b;
        s_trng str;

        buffer = (char *)b_cp__;

        IntDigits = 0;
        if (s<0)
           {
           if (s==MININT) b = MSB;
           else b = -s;
           sign = TRUE;
           }
        else
           {
           b = s;
           sign = FALSE;
           }

        do
           {
           buffer[IntDigits++] = (char)(b%10+'0');
           b /= 10;
           }
        while (b>0);

        i = 0;
        if (TotalWidth>=IntDigits+1)
           {
           s_init(&str,(size_t)TotalWidth);
           if (str.ptr!=NULL)
              {
              for (;i<TotalWidth-IntDigits-1;i++) str.ptr[i] = ' ';
              str.ptr[i++] = ((sign) ? '-' : ' ');
              str.clen = (size_t) TotalWidth;
              }
           }
        else if (sign)
           {
           s_init(&str,(size_t)(IntDigits+1));
           if (str.ptr!=NULL)
              {
              str.ptr[0] = '-';
              i = 1;
              str.clen = (size_t) (IntDigits+1);
              }
           }
        else
           {
           s_init(&str,(size_t)IntDigits);
           if (str.ptr!=NULL) str.clen = (size_t) IntDigits;
           }

        if (str.ptr!=NULL)
           while (IntDigits-->0) str.ptr[i++] = buffer[IntDigits];
        str.tmp = TRUE;

#ifdef HEAP_CHECK
b_tmph((a_char *)&str.ptr);
#endif

        return(str);
        }





