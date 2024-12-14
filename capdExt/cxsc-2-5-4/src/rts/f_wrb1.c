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

/* CVS $Id: f_wrb1.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_wrb1.c                              */
/*                                                              */
/*      Entry           : void f_wrb1(desc,b)                   */
/*                        f_text *desc;                         */
/*                        a_bool b;                             */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        b      - logical value                */
/*                                                              */
/*      Description     : perform PASCAL WRITE(BOOLEAN).        */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local void f_wrb1(f_text *desc,a_bool b)
#else
local void f_wrb1(desc,b)

f_text *desc;
a_bool b;
#endif
        {
        E_TPUSH("f_wrb1")

        if (b_text(desc,FALSE))
           {
           if (b)
              f_wrc1(desc,(a_char *)o_text[35],
                          strlen(o_text[35]));
           else
              f_wrc1(desc,(a_char *)o_text[34],
                          strlen(o_text[34]));
           }

        E_TPOPP("f_wrb1")
        return;
        }





