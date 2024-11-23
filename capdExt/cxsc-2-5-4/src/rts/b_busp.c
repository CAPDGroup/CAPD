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

/* CVS $Id: b_busp.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_busp.c          Processor : C                   */
/*                                                                      */
/* Subtract product of unsigned factors from unsigned array.            */
/* Factors must be non-zero.                                            */
/*                                                                      */
/* Function value : int       - borrow value                            */
/*                                                                      */
/* Function references : b_bms1 - subtract 1 digit number from 1 digit  */
/*                       b_bms2 - subtract 1 digit number from 2 digit  */
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
local int b_busp(a_btyp i,a_btyp u,a_btyp *r)
#else
local int b_busp(i,u,r)

a_btyp i;
a_btyp u;             /* unsigned factors                             */
a_btyp *r;            /* pointer to result mantissa of length 2       */
#endif
        {
        a_btyp ih;    /* upper half of i                              */
        a_btyp uh;    /* upper half of u                              */
        a_btyp p;     /* intermediate product                         */
        a_intg borrow;     /* borrow value                                 */
/*----------------------------------------------------------------------*/
        ih = i >> (B_LENGTH/2);
        i &= LOW_MASK;
        uh = u >> (B_LENGTH/2);
        u &= LOW_MASK;
                        /* split unsigned                               */
/*----------------------------------------------------------------------*/
        borrow = 0;
/*----------------------------------------------------------------------*/
        if (u)
           {
           if (i) borrow += b_bms2(i*u,r);
           if (ih)
              {
              p = ih*u;
              borrow += b_bms2(p<<(B_LENGTH/2),r)+
                        b_bms1(p>>(B_LENGTH/2),r);
              }
           }
/*----------------------------------------------------------------------*/
        if (uh)
           {
           if (i)
              {
              p = i*uh;
              borrow += b_bms2(p<<(B_LENGTH/2),r)+
                        b_bms1(p>>(B_LENGTH/2),r);
              }
           if (ih) borrow += b_bms1(ih*uh,r);
           }
/*----------------------------------------------------------------------*/
        return((int)borrow);
        }





