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

/* CVS $Id: e_ienv.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_ienv.c                              */
/*                                                              */
/*      Entry           : void e_ienv(action,code,mode)         */
/*                        a_intg action;                        */
/*                        a_intg code;                          */
/*                        a_bool mode;                          */
/*                                                              */
/*      Arguments       : code - error code                     */
/*                        action - action to be performed       */
/*                        mode - set or reset                   */
/*                                                              */
/*      Description     : Interface routine for PASCAL.         */
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
local void e_ienv(a_intg action,a_intg code,a_bool mode)
#else
local void e_ienv(action,code,mode)

a_intg action;
a_intg code;
a_bool mode;
#endif
        {
        a_btyp e_action;

        e_action = E_CHNG+E_ACTIVE+E_DEFAULT+1;
        if (((a_btyp)action & (E_CONT+E_BACK))==ZERO) return;
        if (((a_btyp)action & E_CONT) && mode) e_action += E_CONT;
        if (((a_btyp)action & E_BACK) && mode) e_action += E_BACK;

        *(a_btyp *)&code = (a_btyp)code & E_MASK;

        if ((a_btyp)code==INV_OP)
           e_actn(e_action,INV_OP,NO_TRAP);
        else if ((a_btyp)code==DIV_BY_ZERO)
           e_actn(e_action,DIV_BY_ZERO,NO_TRAP);
        else if ((a_btyp)code==OVERFLOW)
           e_actn(e_action,OVERFLOW,NO_TRAP);
        else if ((a_btyp)code==UNDERFLOW)
           e_actn(e_action,UNDERFLOW,NO_TRAP);
        else if ((a_btyp)code==INEXACT)
           e_actn(e_action,INEXACT,NO_TRAP);

        return;
        }





