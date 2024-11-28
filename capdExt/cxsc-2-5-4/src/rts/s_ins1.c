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

/* CVS $Id: s_ins1.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_ins1.c                              */
/*                                                              */
/*      Entries         : a_VOID s_ins1(s,e)                    */
/*                        s_etof s;                             */
/*                        a_intg e;                             */
/*                                                              */
/*      Arguments       : s = set                               */
/*                        e = element                           */
/*                                                              */
/*      Function value  : pointer to updated set                */
/*                                                              */
/*      Description     : Add element e to set s.               */
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
local a_VOID s_ins1(s_etof s,a_intg e)
#else
local a_VOID s_ins1(s,e)

s_etof s;
a_intg e;
#endif
        {
        E_TPUSH("s_ins1")

        if (e<0 || e>s_SIZE*BITS_PER_CHAR-1)
           e_trap(INDEX_RANGE,2,E_TINT+E_TEXT(14),&e);
        else
           s[e/BITS_PER_CHAR] |=(a_char)(
              (((unsigned char)0x80)>>(e%BITS_PER_CHAR)));

        E_TPOPP("s_ins1")
        return((a_VOID)&s[0]);
        }





