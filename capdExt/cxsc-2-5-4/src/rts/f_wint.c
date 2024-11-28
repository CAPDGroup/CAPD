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

/* CVS $Id: f_wint.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_wint.c                              */
/*                                                              */
/*      Entry           : void f_wint(device,i,TotalWidth)      */
/*                        f_text *device;                       */
/*                        a_intg i,TotalWidth;                  */
/*                                                              */
/*      Arguments       : device - file descriptor              */
/*                        i - integer value to be written       */
/*                        TotalWidth - number of characters     */
/*                                                              */
/*      Description     : perform PASCAL integer write.         */
/*                                                              */
/*      Note            : BS PASCAL Standard 6.9.3.2            */
/*                                                              */
/*                   left adjust output if width is negative    */
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
local void f_wint(f_text *device,a_intg i,a_intg TotalWidth)
#else
local void f_wint(device,i,TotalWidth)

f_text *device;
a_intg i;
a_intg TotalWidth;
#endif
        {
        a_intg IntDigits;
        a_bool sign;
        char *buffer;
        a_btyp b;

        buffer = (char *)b_cp__;

        IntDigits = 0;
        if (i<0)
           {
           if (i==MININT) b = MSB;
           else b = -i;
           sign = TRUE;
           }
        else
           {
           b = i;
           sign = FALSE;
           }

        do
           {
           buffer[IntDigits++] = (char) (b%10+'0');
           b /= 10;
           }
        while (b>0);

        if (TotalWidth>0)
           {
           if (TotalWidth>=IntDigits+1)
              {
              for (i=0;i<TotalWidth-IntDigits-1;i++)
                 f_putc((a_char)' ',device);
              f_putc((a_char)((sign) ? MINUS_SIGN : ' '),device);
              }
           else
              if (sign) f_putc((a_char)MINUS_SIGN,device);

           while (--IntDigits>=0)
              f_putc((a_char)buffer[IntDigits],device);
           }
        else
           {
           if ((TotalWidth = -TotalWidth-IntDigits-1)>=0)
              f_putc((a_char)((sign) ? MINUS_SIGN : ' '),device);
           else if (sign)
              f_putc((a_char)MINUS_SIGN,device);
           for (i = --IntDigits;i>=0;i--)
              f_putc((a_char)buffer[i],device);
           while (++i<TotalWidth)
              f_putc((a_char)' ',device);
           }

        return;
        }





