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

/* CVS $Id: f_wri1.c,v 1.21 2014/01/30 17:24:08 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_wri1.c                              */
/*                                                              */
/*      Entry           : void f_wri1(desc,i)                   */
/*                        f_text *desc;                         */
/*                        a_intg i;                             */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        i      - integer value                */
/*                                                              */
/*      Description     : perform PASCAL write(integer).        */
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
local void f_wri1(f_text *desc,a_intg i)
#else
local void f_wri1(desc,i)

f_text *desc;
a_intg i;
#endif
        {
        E_TPUSH("f_wri1")

        if (b_text(desc,FALSE))
           f_wint(desc,i,0);

        E_TPOPP("f_wri1")
        }





