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

/* CVS $Id: s_real.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_real.c                              */
/*                                                              */
/*      Entries         : s_trng s_real(s,Totalwidth,           */
/*                                    FracDigits,rnd)           */
/*                        a_real s;                             */
/*                        a_intg rnd,Totalwidth,FracDigits;     */
/*                                                              */
/*      Arguments       : s - IEEE value                        */
/*                        Totalwidth - total length of string   */
/*                        FracDigits - number of fraction digits*/
/*                        rnd - rounding mode                   */
/*                              -1 = round downwards            */
/*                               0 = round to nearest           */
/*                               1 = round upwards              */
/*                                                              */
/*      Description     : Convert an IEEE double format number  */
/*                        to a character string.                */
/*                                                              */
/*                   replace strncpy by memcpy                  */
/*                   set default TotalWidth                     */
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
local s_trng s_real(a_real s,a_intg TotalWidth,a_intg FracDigits,a_intg rnd)
#else
local s_trng s_real(s,TotalWidth,FracDigits,rnd)

a_real s;
a_intg TotalWidth;
a_intg FracDigits;
a_intg rnd;
#endif
        {
        s_trng str;
        a_intg length;
        char *buffer;

        E_TPUSH("s_real")

        buffer = (char *)&b_cp__[0];

        if (TotalWidth<=0)
           TotalWidth = ExpDigits+4+3*MANTL/10;

        r_outp(buffer,s,TotalWidth,FracDigits,rnd,&length);

        s_init(&str,(size_t)length);
        str.clen = str.alen;
        if (length>0 && str.ptr!=NULL)
           memcpy(str.ptr,buffer,(size_t)length);
        str.tmp = TRUE;

#ifdef HEAP_CHECK
b_tmph((a_char *)&str.ptr);
#endif

        E_TPOPP("s_real")
        return(str);
        }





