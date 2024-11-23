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

/* CVS $Id: f_bhex.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_bhex.c                              */
/*                                                              */
/*      Entry           : void f_bhex(desc,r,mode)              */
/*                        f_text *desc;                         */
/*                        a_char r;                             */
/*                        a_char mode;                          */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        ch     - character value              */
/*                        mode   - mode value for letter        */
/*                                 generation                   */
/*                                 'x' = 'abcdef' is used       */
/*                                 'X' = 'ABCDEF' is used       */
/*                                                              */
/*      Description     : write character value in hexadecimal  */
/*                        format                                */
/*                                                              */
/*      Note            : No checking for valid descriptor done */
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

static char lc[] = "0123456789abcdef";
static char uc[] = "0123456789ABCDEF";

#ifdef LINT_ARGS
local void f_bhex(f_text *desc,a_char r,a_char mode)
#else
local void f_bhex(desc,r,mode)

f_text *desc;
a_char r;
a_char mode;
#endif
        {
        a_char c;

        E_TPUSH("f_bhex")

        c = r>>4;
        f_putc((mode=='x') ? lc[c] : uc[c],desc);
        c = r & 0x0f;
        f_putc((mode=='x') ? lc[c] : uc[c],desc);

        E_TPOPP("f_bhex")
        return;
        }





