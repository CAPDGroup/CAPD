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

/* CVS $Id: b_subm.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_subm.c                              */
/*                                                              */
/*      Entry           : a_bool b_subm(k,a,b)                  */
/*                        a_intg k;                             */
/*                        a_btyp *a,*b;                         */
/*                                                              */
/*      Arguments       : a = first operand of subtraction      */
/*                        b = second operand of subtraction     */
/*                        k = size of operands                  */
/*                                                              */
/*      Function value  : TRUE if borrow remains                */
/*                                                              */
/*      Description     : Subtraction of a_btyp arrays          */
/*                        a = a-b                               */
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
local a_bool b_subm(a_intg k,a_btyp *a,a_btyp *b)
#else
local a_bool b_subm(k,a,b)

a_intg k;
a_btyp *a;
a_btyp *b;
#endif
        {
        a_intg borrow = 0;
        a_btyp c;

        while (--k>=0) {
#if C_P_2
                c = a[k]-b[k]-borrow;
                borrow = (a[k]<b[k] ||
                         (a[k]==b[k] && borrow)) ? 1 : 0;
                a[k] = c;
#else
                b_subu(a[k],b[k],borrow,a+k,&borrow);
#endif
                }

        return((borrow) ? TRUE : FALSE);
        }





