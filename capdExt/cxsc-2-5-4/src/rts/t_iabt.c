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

/* CVS $Id: t_iabt.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/*
**   ieee_abort[ri][12]:
**     werden aus den PXSC-aufrufbaren Funktionen aus aufgerufen,
**     falls ein Laufzeitfehler erzeugt werden soll.
**     Die Aufrufe stehen in makehead.h.
**
**     pxsccode: setzt die ieee-eigenen Fehlercodes in die 
**               PXSC-Fehlercodes fuer e_trap um.
*/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#ifdef LINT_ARGS
a_btyp pxsccode(int rc)
#else
a_btyp pxsccode(rc)
int rc;
#endif /* LINT_ARGS */
{
    switch (rc)
    {
        case DOMAIN:    return(INV_ARG);
        case SING:      return(INV_ARG);
        case OVER_FLOW: return(OVERFLOW);   /* +E_IEEE ? */
        case UNDER_FLOW: return(UNDERFLOW);
        case TLOSS:     return(INEXACT);
        case PLOSS:     return(INEXACT);
        case ExcNoI:    return(INV_ARG);
        case ExcInvalid: return(INV_ARG);
        case ExcISing:  return(INV_ARG);
        case ExcDBZ:    return(DIV_BY_ZERO);
        case ExcDIZ:    return(DIV_BY_ZERO);
        default:        return(INV_ARG);
    }
}


#ifdef LINT_ARGS
void ieee_abortr1(int rc, void *arg)
#else
void ieee_abortr1(rc, arg)
int rc;
int *arg;       /* IBM PC/RT cc kennt keine void * */
#endif /* LINT_ARGS */
{
    a_btyp ecode;
    ecode = pxsccode(rc);
    e_trap(ecode, 2, E_TDBL, arg);
}

#ifdef LINT_ARGS
void ieee_abortr2(int rc, void *arg1, void *arg2)
#else
void ieee_abortr2(rc, arg1, arg2)
int rc;
int *arg1, *arg2;
#endif /* LINT_ARGS */
{
    a_btyp ecode;
    ecode = pxsccode(rc);
    e_trap(ecode, 4, E_TDBL, arg1, E_TDBL, arg2);
}

#ifdef LINT_ARGS
void ieee_aborti1(int rc, void *ai)
#else
void ieee_aborti1(rc, ai)
int rc;
int *ai;
#endif /* LINT_ARGS */
{
    a_btyp ecode;
    ecode = pxsccode(rc);
    e_trap(ecode, 4,
           E_TDBL, &(((a_intv *)ai)->INF), E_TDBL, &(((a_intv *)ai)->SUP));
}

#ifdef LINT_ARGS
void ieee_aborti2(int rc, void *ai, void *bi)
#else
void ieee_aborti2(rc, ai, bi)
int rc;
int *ai, *bi;
#endif /* LINT_ARGS */
{
    a_btyp ecode;
    ecode = pxsccode(rc);
    e_trap(ecode, 8,
           E_TDBL, &(((a_intv *)ai)->INF), E_TDBL, &(((a_intv *)ai)->SUP),
           E_TDBL, &(((a_intv *)bi)->INF), E_TDBL, &(((a_intv *)bi)->SUP));
}





