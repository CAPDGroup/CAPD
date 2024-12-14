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

/* CVS $Id: f_wrb2.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_wrb2.c                              */
/*                                                              */
/*      Entry           : void f_wrb2(desc,b,w)                 */
/*                        f_text *desc;                         */
/*                        a_bool b;                             */
/*                        a_intg w;                             */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        b      - logical value                */
/*                        w      - width of output field        */
/*                                                              */
/*      Description     : perform PASCAL WRITE(BOOLEAN:W).      */
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
local void f_wrb2(f_text *desc,a_bool b,a_intg w)
#else
local void f_wrb2(desc,b,w)

f_text *desc;
a_bool b;
a_intg w;
#endif
        {
        int p;

        E_TPUSH("f_wrb2")

        if (b_text(desc,FALSE))
           {
           p = 0;
           if (b)
              {
              if ((w<0) && (w>-strlen(o_text[35])))
                 p = (int)(strlen(o_text[35])+w);
              f_wrc2(desc,(a_char *)o_text[35],
                          strlen(o_text[35])-p,w);
              }
           else
              {
              if ((w<0) && (w>-strlen(o_text[34])))
                 p = (int)(strlen(o_text[34])+w);
              f_wrc2(desc,(a_char *)o_text[34],
                          strlen(o_text[34])-p,w);
              }
           }

        E_TPOPP("f_wrb2")
        return;
        }





