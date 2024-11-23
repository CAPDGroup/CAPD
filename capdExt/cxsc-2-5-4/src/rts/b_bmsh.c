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

/* CVS $Id: b_bmsh.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bmsh.c          Processor : C                   */
/*                                                                      */
/* Shift unsigned array to the left by several positions.               */
/* Bits shifted out are lost. Shift value : 0<n<B_LENGTH                */
/*                                                                      */
/* Function value : int     0 -                                         */
/*                                                                      */
/* Function references : none                                           */
/*                                                                      */
/************************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local int b_bmsh(a_intg l,a_btyp *s,a_intg n)
#else
local int b_bmsh(l,s,n)

a_intg l;       /* length of unsigned array                             */
a_btyp *s;      /* pointer to unsigned array                            */
a_intg n;       /* number of bits to be shifted                         */
#endif
        {
/*----------------------------------------------------------------------*/
        while (--l>0)
           {
           *s = ((*s)<<n) | (*(s+1)>>(B_LENGTH-n));
           s++;
           }
/*----------------------------------------------------------------------*/
        *s <<= n;
/*----------------------------------------------------------------------*/
        return(0);
        }





