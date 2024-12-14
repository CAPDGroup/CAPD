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

/* CVS $Id: f_get_.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_get_.c                              */
/*                                                              */
/*      Entry           : void f_get_(desc)                     */
/*                        f_text *desc;                         */
/*                                                              */
/*      Arguments       : desc   - descriptor of device         */
/*                                                              */
/*      Description     : PASCAL get function.                  */
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
local void f_get_(f_text *desc)
#else
local void f_get_(desc)

f_text *desc;
#endif
        {
        E_TPUSH("f_get_")

        if (desc->fp==NULL || desc->asgd==FALSE)
           e_trap(I_O_ERROR,4,E_TMSG,17,E_TSTR+E_TEXT(8),desc->name);
        else if (desc->infl==FALSE)
           e_trap(I_O_ERROR,4,E_TMSG,18,E_TSTR+E_TEXT(8),desc->name);
        else if (desc->eof==TRUE)
           e_trap(I_O_ERROR,4,E_TMSG,20,E_TSTR+E_TEXT(8),desc->name);
        else
           f_getc(desc);

        E_TPOPP("f_get_")
        }





