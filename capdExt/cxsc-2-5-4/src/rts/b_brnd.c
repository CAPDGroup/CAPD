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

/* CVS $Id: b_brnd.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_brnd.c          Processor : C                   */
/*                                                                      */
/* Determine intern floating point representation which is greater      */
/* or equal to input number + rounding ulps. Zero is returned if        */
/* zero was input.                                                      */
/*                                                                      */
/* Function value : int     0 - next number returned                    */
/*                      OFLOW - exponent overflow                       */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Function references : b_bcpy - copy intern numbers                   */
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
local int b_brnd(multiprecision i,multiprecision r)
#else
local int b_brnd(i,r)

multiprecision i;       /* pointer to intern variable                   */
multiprecision r;       /* pointer to intern result                     */
#endif
        {
        int rc;            /* !!! must be int !!! */
        a_intg k;          /* loop variable                                */
/*----------------------------------------------------------------------*/
        if (i->r==0)
           {
           if ((rc = b_bcpy(i,r))!=0) return(rc);
           if (r->r==0) return(rc);
           return(b_banx(r,r));
           }
                        /* rounding indicator not set                   */
/*----------------------------------------------------------------------*/
        if (i->z)
           {
           r->r = 0;
           r->z = 1;
           return(0);
           }            /* input is zero   zero returned                */
/*----------------------------------------------------------------------*/
        r->z = 0;       /* result is non-zero                           */
/*----------------------------------------------------------------------*/
        r->s = i->s;    /* copy sign of number                          */
/*----------------------------------------------------------------------*/
        if ((rc = b_banx(i,r))!=0) return(rc);
/*----------------------------------------------------------------------*/
        for (k=1;k<i->r;k++) if ((rc = b_banx(r,r))!=0) return(rc);
/*----------------------------------------------------------------------*/
        r->r = 0;
/*----------------------------------------------------------------------*/
        return(rc);     /* next absolute number                         */
        }





