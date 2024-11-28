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

/* CVS $Id: b_bmsp.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bmsp.c          Processor : C                   */
/*                                                                      */
/* Subtract product of internal mantissa with unsigned                  */
/* from the result unsigned array.                                      */
/* The internal and the unsigned number must be nonzero.                */
/*                                                                      */
/* Function value : int       - borrow value                            */
/*                                                                      */
/* Function references : b_busp - subtract product of unsigned numbers  */
/*                       b_bcsu - propagate borrow                      */
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
local int b_bmsp(a_intg n,a_btyp *i,a_btyp u,a_btyp *r)
#else
local int b_bmsp(n,i,u,r)

a_intg n;               /* length of internal mantissa                  */
a_btyp *i;            /* pointer to intern mantissa of length n       */
a_btyp u;             /* unsigned factor                              */
a_btyp *r;            /* pointer to result mantissa of length n+1     */
#endif
        {
        a_intg borrow = 0; /* borrow value                              */
/*----------------------------------------------------------------------*/
        while (--n>=0)
        if (*(i+n))
        if (b_busp(*(i+n),u,r+n))
        if (b_bcsu(n,r)) borrow = 1;
                        /* subtract product of unsigned factors         */
/*----------------------------------------------------------------------*/
        return((int)borrow);
        }





