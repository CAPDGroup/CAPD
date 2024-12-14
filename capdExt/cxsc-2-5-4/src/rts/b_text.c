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

/* CVS $Id: b_text.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_text.c                              */
/*                                                              */
/*      Entry           : a_bool b_text(desc,in)                */
/*                        f_text *desc;                         */
/*                        a_bool in;                            */
/*                                                              */
/*      Arguments       : desc   - file descriptor              */
/*                        in     - test for input file          */
/*                                                              */
/*      Description     : Display message if not text file      */
/*                        and not opened for input (in==TRUE)   */
/*                        or not opened for output (in==FALSE). */
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
local a_bool b_text(f_text *desc,a_bool in)
#else
local a_bool b_text(desc,in)

f_text *desc;
a_bool in;
#endif
        {
        int msg = 0;

        E_TPUSH("b_text")

        if (desc->asgd==FALSE || desc->fp==NULL) msg = 17;
#ifdef CHECK_TEXT_FILE
        else if (desc->text==FALSE) msg = 19;
#endif
        else if (in)
           {
           if (desc->infl==FALSE) msg = 18;
           if (desc->eof==TRUE) msg = 20;
           }
        else
           {
           if (desc->outf==FALSE) msg = 34;
           if (desc->err==TRUE) msg = 35;
           }

        if (msg) e_trap(I_O_ERROR,4,E_TMSG,msg,E_TSTR+E_TEXT(8),desc->name);

        E_TPOPP("b_text")
        return((msg) ? FALSE : TRUE);
        }





