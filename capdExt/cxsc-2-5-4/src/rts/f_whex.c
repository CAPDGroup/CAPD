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

/* CVS $Id: f_whex.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_whex.c                              */
/*                                                              */
/*      Entry           : void f_whex(desc,r,mode)              */
/*                        f_text *desc;                         */
/*                        a_real r;                             */
/*                        a_char mode;                          */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        r      - real value                   */
/*                        mode   - read mode                    */
/*                                 'x' = hexadecimal (small)    */
/*                                 'X' = hexadecimal (capital)  */
/*                                                              */
/*      Description     : write real value in special format    */
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
local void f_whex(f_text *desc,a_real r,a_char mode)
#else
local void f_whex(desc,r,mode)

f_text *desc;
a_real r;
a_char mode;
#endif
        {
        size_t i;
        a_char *p;

        E_TPUSH("f_whex")

        if (b_text(desc,FALSE))
           {
           if (mode!='x' && mode!='X')
              e_trap(I_O_ERROR,4,E_TMSG,51,E_TCHR,&mode);
           else
              {
#if INTEL
              p = ((a_char *)&r)+(sizeof(r)-1);
              for (i=0;i<sizeof(r);i++) f_bhex(desc,*p--,mode);
#else
              p = (a_char *)&r;
              for (i=0;i<sizeof(r);i++) f_bhex(desc,*p++,mode);
#endif
              }
           }

        E_TPOPP("f_whex")
        return;
        }





