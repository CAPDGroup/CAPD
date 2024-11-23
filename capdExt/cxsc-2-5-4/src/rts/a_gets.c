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

/* CVS $Id: a_gets.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_gets.c                              */
/*                                                              */
/*      Entries         : a_btyp a_gets(rc)                     */
/*                        a_btyp *rc;                           */
/*                                                              */
/*      Arguments       : rc = error condition code             */
/*                                                              */
/*      Description     : Get hardware status and determine     */
/*                        error condition code.                 */
/*                                                              */
/*      Note            : Version for IBM_RT_C                  */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_bool e_efdz;
extern a_bool e_efie;
extern a_bool e_efio;
extern a_bool e_efof;
extern a_bool e_efuf;
extern a_bool e_ofdz;
extern a_bool e_ofie;
extern a_bool e_ofio;
extern a_bool e_ofof;
extern a_bool e_ofuf;
#endif

#ifdef IEEE_HARDWARE
#if IBM_RT_C

#undef MININT
#include <sys/FP.h>
#undef MININT

#ifdef LINT_ARGS
local a_btyp a_gets(a_btyp *rc)
#else
local a_btyp a_gets(rc)

a_btyp *rc;
#endif
        {
        FP_STATUS status;

        /* get status */
        *((unsigned *)&status) =
           (*(unsigned (*)())_fpfpf[(int)FP_getst])();

        /* test status flag for exception occurrence */
        if (status.xcp_flag==0)
           {

           /* test status flag for inexact result */
           if (status.ir_flag==0) return(*rc = NO_ERROR);
           if (e_efie) return(*rc = INEXACT);
           e_ofie = TRUE;
           }

        /* test status flag for overflow occurrence */
        else if (status.of_flag)
           {
           if (e_efof) return(*rc = OVERFLOW);
           e_ofof = TRUE;
           }

        /* test status flag for underflow occurrence */
        else if (status.uf_flag)
           {
           if (e_efuf) return(*rc = UNDERFLOW);
           e_ofuf = TRUE;
           }

        /* test status flag for invalid-operation occurrence */
        else if (status.io_flag)
           {
           if (e_efio) return(*rc = INV_OP);
           e_ofio = TRUE;
           }

        /* test status flag for division-by-zero occurrence */
        else if (status.dz_flag)
           {
           if (e_efdz) return(*rc = DIV_BY_ZERO);
           e_ofdz = TRUE;
           }

        return(*rc = NO_ERROR);
        }
#else
#ifdef LINT_ARGS
local a_btyp a_gets(a_btyp *rc)
#else
local a_btyp a_gets(rc)

a_btyp *rc;
#endif
        {
        return(NO_ERROR);
        }
#endif
#else
#ifdef LINT_ARGS
local a_btyp a_gets(a_btyp *rc)
#else
local a_btyp a_gets(rc)

a_btyp *rc;
#endif
        {
        return(NO_ERROR);
        }
#endif





