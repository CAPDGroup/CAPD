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

/* CVS $Id: s_rhex.c,v 1.22 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_rhex.c                              */
/*                                                              */
/*      Entry           : a_real s_rhex(s_trng)                 */
/*                        s_trng s;                             */
/*                                                              */
/*      Arguments       : s      - input string                 */
/*                                                              */
/*      Description     : convert hex real string to real       */
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
local a_real s_rhex(s_trng s)
#else
local a_real s_rhex(s)
 
s_trng s;
#endif
        {
        int i;
        a_real r ;
        a_btyp b = ZERO;
	   int reallen= 2*sizeof(a_real); /* number of nibbles a word */
 
        E_TPUSH("s_rhex")

        /* check string lenght */
        if (reallen != s.alen)
	   {
                 e_trap(INV_ARG    ,6,E_TMSG,56,E_TSTR,&s.ptr[0],
				    E_TINT,&s.alen);
        }
        for (i=0;i<reallen;i++)
        {
                    b <<= 4;
                    if (isdigit((int)s.ptr[i]))
                       b += s.ptr[i]-'0';
                    else if (isalpha((int)s.ptr[i]))
                       b += toupper(s.ptr[i])-('A'-10);
                    else
                       {
                       e_trap(I_O_ERROR,4,E_TMSG,52,
                              E_TCHR+E_TEXT(10),&s.ptr[i]);
                       }
 
                    if (i==7)
                       {
                       ((a_btyp *)&r)[B_HPART] = b;
                       b = ZERO;
                       }
                    else if (i==15) ((a_btyp *)&r)[B_LPART] = b;
           }
 
        E_TPOPP("s_rhex")
        return(r);
        }





