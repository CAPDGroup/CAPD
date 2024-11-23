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

/* CVS $Id: e_xarg.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_xarg.c                              */
/*                                                              */
/*      Entry           : void e_xarg(code,e_argc,e_argv)       */
/*                        a_btyp code;                          */
/*                        va_list e_argv;                       */
/*                                                              */
/*      Arguments       : code - error code                     */
/*                        e_argc - number of arguments          */
/*                        e_argv - reference to arguments       */
/*                                                              */
/*      Description     : Invalid argument.                     */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern f_text f_errr;
#endif

#ifdef LINT_ARGS
local void e_xarg(a_btyp code,int e_argc,va_list e_argv)
#else
local void e_xarg(code,e_argc,e_argv)

a_btyp code;
int e_argc;
va_list e_argv;
#endif
        {
        a_bool print;

        if (!(code & E_EQIE))
           {
           if ((print = ((code & E_EMSG) ? TRUE : FALSE))==TRUE)
              e_tmsg(12);

           if (code & E_EARG) e_tprt(e_argc,e_argv);
           else if (code & E_EMSG) e_tmrt(e_argc,e_argv,print);

           if (code & E_ETBC) e_back(f_errr.fp);
           else if (print) e_bmsg(f_errr.fp);
           }

        if (code & E_EXIT)
           {
           e_tmsg(25);
           exit(EXIT_VALUE);
           }

        return;
        }





