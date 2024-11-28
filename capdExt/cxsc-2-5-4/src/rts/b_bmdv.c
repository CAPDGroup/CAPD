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

/* CVS $Id: b_bmdv.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bmdv.c          Processor : C                   */
/*                                                                      */
/* Divide 3 digit number by 2 digit number to get a 1 digit result.     */
/*                                                                      */
/* Note : First digit of dividend < first digit of divisor.             */
/*                                                                      */
/* Entry type          : static                                         */
/* Function type       : int                                            */
/* Function value      : 0 - division done                              */
/* Error codes         : -                                              */
/* Include files       : defs.h                                         */
/* Function references : -                                              */
/* Globals referenced  : -                                              */
/* Rounding            : -                                              */
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
local int b_bmdv(a_btyp *i1,a_btyp *i2,a_btyp *r)
#else
local int b_bmdv(i1,i2,r)

a_btyp *i1;                           /* pointer to dividend          */
a_btyp *i2;                           /* pointer to divisor           */
a_btyp *r;                            /* pointer to result digit      */
#endif
        {
        a_intg k;                       /* loop variable                */
        a_btyp s;                     /* indicates subtraction        */
        a_btyp d,e,f;                 /* digits of dividend           */
/*----------------------------------------------------------------------*/
        *r = 0;                         /* initialize result            */
/*----------------------------------------------------------------------*/
        if ((d = *i1)==ZERO) return(0); /* result is zero               */
/*----------------------------------------------------------------------*/
        e = *(i1+1);                    /* copy dividend i1 => d,e,f    */
        f = *(i1+2);
/*----------------------------------------------------------------------*/
        if ((d==*i2) && e==*(i2+1))
           {
           *r = MAX_BASETYPE;
           return(0);
           }
/*----------------------------------------------------------------------*/
        for (k=0;k<B_LENGTH;k++)
           {
           *r <<= 1;                    /* shift result 1 bit left      */
           s = d & MSB;                 /* MSB of dividend              */
                                        /* shift dividend 1 bit left    */
           d = (d<<1) | (e>>(B_LENGTH-1));
           e = (e<<1) | (f>>(B_LENGTH-1));
           f <<= 1;
           if (s==ZERO)                 /* MSB not set, divid.>=divis.  */
              if (d>*i2 || (d==*i2 && e>=*(i2+1))) s = 1;
           if (s)
              {
              (*r)++;                   /* increment result             */
              if (e<*(i2+1))
                 {
                 if (d) d--; else d = MAX_BASETYPE;
                 }
              d -= *i2;
              e -= *(i2+1);
              }
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





