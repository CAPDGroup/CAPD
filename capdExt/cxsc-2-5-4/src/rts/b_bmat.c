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

/* CVS $Id: b_bmat.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bmat.c          Processor : C                   */
/*                                                                      */
/* Carry from addition of mantissas                                     */
/*                                                                      */
/* Function value :   int     - value of carry                          */
/*                                                                      */
/* Functione references : none                                          */
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
local int b_bmat(a_intg n,a_btyp *i1,a_btyp *i2,a_intg carry,a_intg *r)
#else
local int b_bmat(n,i1,i2,carry,r)

a_intg n;               /* length of mantissas added                    */
a_btyp *i1;
a_btyp *i2;           /* pointer to intern mantissas                  */
a_intg carry;           /* carry value at lsb position                  */
a_intg *r;              /* non-zero mantissa                            */
#endif
        {
        a_btyp rn;    /* sum of digits                                */
/*----------------------------------------------------------------------*/
        i1 += n;
        i2 += n;
        rn = 0;
        while (--n>=0)
        if (~*(--i1)>*(--i2))
           {
           rn |= *i1+*i2+carry;
           carry = 0;
           }
        else if (~*i1<*i2)
           {
           if (*i1 & MSB)
              if (*i2 & MSB)
                 rn |= (*i1 ^ MSB) + (*i2 ^ MSB) + carry;
              else
                 rn |= ((*i1 ^ MSB) + *i2 + carry) ^ MSB;
           else
              rn |= (*i1 + (*i2 ^ MSB) + carry) ^ MSB;
           carry = 1;
           }
        else
           rn |= ((carry) ? ZERO : MAX_BASETYPE);
/*----------------------------------------------------------------------*/
        *r = (rn) ? 1 : 0;
/*----------------------------------------------------------------------*/
        return((int)carry);
        }





