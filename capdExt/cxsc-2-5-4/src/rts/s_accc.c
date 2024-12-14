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

/* CVS $Id: s_accc.c,v 1.21 2014/01/30 17:24:13 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_accc.c                              */
/*                                                              */
/*      Entries         : a_VOID s_accc(s,t,u)                  */
/*                        a_char s[],t[];                       */
/*                        a_char u;                             */
/*                                                              */
/*      Arguments       : s = concatenated string               */
/*                        t = string                            */
/*                        u = character                         */
/*                                                              */
/*      Return value    : concatenated string                   */
/*                                                              */
/*      Description     : Concatenate string with character.    */
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
local a_VOID s_accc(a_char s[],a_char t[],a_char u)
#else
local a_VOID s_accc(s,t,u)

a_char s[];
a_char t[];
a_char u;
#endif
        {
        a_char *p,*q;

        E_TPUSH("s_accc")

        p = &s[0];
        q = &t[0];

        while ((*p++ = *q++)!=0) { /* empty */ }
        p[-1] = u;
        *p = '\0';

        E_TPOPP("s_accc")
        return((a_VOID)&s[0]);
        }





