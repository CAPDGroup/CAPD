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

/* CVS $Id: a_sets.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_sets.c                              */
/*                                                              */
/*      Entries         : void a_sets(rounding)                 */
/*                        int rounding;                         */
/*                                                              */
/*      Arguments       : rounding = rounding mode              */
/*                                                              */
/*      Description     : Set rounding mode in hardware status. */
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
/*
extern a_bool e_efdz;
extern a_bool e_efie;
extern a_bool e_efio;
extern a_bool e_efof;
extern a_bool e_efuf;
*/
#endif

#ifdef IEEE_HARDWARE
#if IBM_RT_C

#undef MININT
#include <sys/FP.h>
#undef MININT

#ifdef LINT_ARGS
void a_sets(int rounding)
#else
void a_sets(rounding)

int rounding;
#endif
        {
        FP_STATUS status;

        /* get status from processor */
        *((unsigned *)&status) =
           (*(unsigned (*)())_fpfpf[(int)FP_getst])();

        /* set/reset trap enable flags if kill flag is set */
/*
        *((unsigned *)&status) |= 0x80000000L;
        if (e_efio) *((unsigned *)&status) |= 0x10000000L;
        else        *((unsigned *)&status) &= 0xefffffffL;
        if (e_efdz) *((unsigned *)&status) |= 0x04000000L;
        else        *((unsigned *)&status) &= 0xfbffffffL;
        if (e_efof) *((unsigned *)&status) |= 0x01000000L;
        else        *((unsigned *)&status) &= 0xfeffffffL;
        if (e_efuf) *((unsigned *)&status) |= 0x00400000L;
        else        *((unsigned *)&status) &= 0xffbfffffL;
        if (e_efie) *((unsigned *)&status) |= 0x00000020L;
        else        *((unsigned *)&status) &= 0xffffffdfL;
*/

        /* reset exception occurrence flags and rounding */
        /* reset kill flag */
/*
        *((unsigned *)&status) &= 0x957ffe3fL;
        *((unsigned *)&status) &= 0x7fffffffL;
*/
        *((unsigned *)&status) &= 0x157ffe3fL;

        /* set rounding mode */
        status.rnd_mode = rounding;

        /* set status in processor */
        (*(unsigned (*)())_fpfpf[(int)FP_setst])(status);
        return;
        }
#else
#ifdef LINT_ARGS
void a_sets(int rounding)
#else
void a_sets(rounding)

int rounding;
#endif
        {
        }
#endif
#else
#ifdef LINT_ARGS
void a_sets(int rounding)
#else
void a_sets(rounding)

int rounding;
#endif
        {
        }
#endif





