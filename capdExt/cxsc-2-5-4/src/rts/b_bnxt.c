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

/* CVS $Id: b_bnxt.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bnxt.c          Processor : C                   */
/*                                                                      */
/* Determine next number away from zero.                                */
/* If zero was input zero is returned.                                  */
/*                                                                      */
/* Function value : int     0 - next number returned                    */
/*                      OFLOW - exponent overflow                       */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Function references :                                                */
/*                       b_banx - next number in absolute value         */
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
local int b_bnxt(multiprecision i,multiprecision r)
#else
local int b_bnxt(i,r)

multiprecision i;               /* pointer to intern variable           */
multiprecision r;               /* pointer to intern result             */
#endif
        {
/*----------------------------------------------------------------------*/
        if (i->z)
           {
           r->r = 0;
           r->z = 1;
           return(0);
           }            /* input is zero   zero returned                */
/*----------------------------------------------------------------------*/
        r->r = 0;
/*----------------------------------------------------------------------*/
        r->z = 0;       /* result is non-zero                           */
/*----------------------------------------------------------------------*/
        r->s = i->s;    /* copy sign of number                          */
/*----------------------------------------------------------------------*/
        return(b_banx(i,r));
                        /* next absolute number                         */
        }





