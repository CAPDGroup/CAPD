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

/* CVS $Id: s_ixch.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_ixch.c                              */
/*                                                              */
/*      Entry           : a_intg s_ixch(i,length)               */
/*                        a_intg i;                             */
/*                        size_t length                         */
/*                                                              */
/*      Arguments       : i - index to component of string      */
/*                        length - length of string             */
/*                                                              */
/*      Description     : Determine offset to component of      */
/*                        string with index checking.           */
/*                                                              */
/*      Note            : -1 is returned in case of error       */
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
local a_intg s_ixch(a_intg i,size_t length)
#else
local a_intg s_ixch(i,length)

a_intg i;
size_t length;
#endif
        {
        a_intg res,l,u;

        E_TPUSH("s_ixch")

        if (i>0 && i<=length) res = i-1;
        else
           {
           l = 1;
           u = length;
           e_trap(INDEX_RANGE,6,E_TINT+E_TEXT(4),&i,
                                E_TINT+E_TEXT(5),&l,
                                E_TINT+E_TEXT(6),&u);
           res = -1;
           }

        E_TPOPP("s_ixch")
        return(res);
        }





