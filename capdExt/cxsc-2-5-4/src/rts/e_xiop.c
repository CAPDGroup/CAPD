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

/* CVS $Id: e_xiop.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_xiop.c                              */
/*                                                              */
/*      Entry           : void e_xiop(code,e_argc,e_argv)       */
/*                        a_btyp code;                          */
/*                        va_list e_argv;                       */
/*                                                              */
/*      Arguments       : code - error code                     */
/*                        e_argc - number of arguments          */
/*                        e_argv - reference to arguments       */
/*                                                              */
/*      Description     : Invalid operation.                    */
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
extern int e_rtyp;
extern a_VOID e_rptr;
extern a_bool e_efio;
extern a_bool e_ofio;
#endif

#ifdef LINT_ARGS
local void e_xiop(a_btyp code,int e_argc,va_list e_argv)
#else
local void e_xiop(code,e_argc,e_argv)

a_btyp code;
int e_argc;
va_list e_argv;
#endif
        {
        a_bool print;

        if (!(code & E_EQIE))
           {
           if (code & E_IEEE)
              {
              if ((print = e_efio ||
                           ((code & E_EXIT) ? TRUE : FALSE))==TRUE)
                 e_tmsg(3);
              e_ofio = TRUE;
              }
           else if ((print = ((code & E_EMSG) ? TRUE : FALSE))==TRUE)
              e_tmsg(3);

           if (code & E_EARG) e_tprt(e_argc,e_argv);
           else if (code & E_EMSG) e_tmrt(e_argc,e_argv,print);

           if (code & E_ETBC) e_back(f_errr.fp);
           else if (print) e_bmsg(f_errr.fp);
           }
        else if (code & E_IEEE) e_ofio = TRUE;

        if (e_rtyp==E_TDBL)
           {
           ((a_btyp *)e_rptr)[B_HPART] = SET_QUIET(EXPO_MASK);
           ((a_btyp *)e_rptr)[B_LPART] = INV_OP;
           }
        else if (e_rtyp==E_TLNG)
           {
           int i;
           a_btyp x = INV_OP;
           for (i=0; i<12; i++)
              ((char *)e_rptr)[i] = 0x00;
           ((char *)e_rptr)[0] = 0x3f;
           ((char *)e_rptr)[1] = 0xf0;
           for (i=0; i<4; i++)
              ((char *)e_rptr)[i+4] = ((char *)&x)[i];
           }

        if (code & E_EXIT)
           {
           e_tmsg(25);
           exit(EXIT_VALUE);
           }

        return;
        }





