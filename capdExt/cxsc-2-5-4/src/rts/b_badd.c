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

/* CVS $Id: b_badd.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_badd.c          Processor : C                   */
/*                                                                      */
/* Add numbers in intern representation                                 */
/*                                                                      */
/* Function value : int     0 - intern numbers added                    */
/*                      OFLOW - exponent overflow                       */
/*                      UFLOW - exponent underflow                      */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Function references : b_bacm - compare absolute values of intern     */
/*                       b_bcpy - copy intern numbers                   */
/*                       b_baad - add absolute intern values            */
/*                       b_basu - subtract absolute intern values       */
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
local int b_badd(multiprecision i1,multiprecision i2,multiprecision r)
#else
local int b_badd(i1,i2,r)

multiprecision i1;
multiprecision i2;      /* pointer to intern variables                  */
multiprecision r;       /* pointer to intern result                     */
#endif
        {
        a_intg cv;         /* compare value                                */
/*----------------------------------------------------------------------*/
        if (i1->z) return(b_bcpy(i2,r));
                        /* first operand is zero                        */
/*----------------------------------------------------------------------*/
        if (i2->z) return(b_bcpy(i1,r));
                        /* second operand is zero                       */
/*----------------------------------------------------------------------*/
        cv = b_bacm(i1,i2);
                        /* compare absolute values of intern numbers    */
/*----------------------------------------------------------------------*/
        if (i1->s==i2->s)
           {
           r->s = i1->s;
           if (cv>=0)
              return(b_baad(i1,i2,r));
           else
              return(b_baad(i2,i1,r));
           }            /* addition of equaly signed numbers            */
/*----------------------------------------------------------------------*/
        if (cv<0)
           {
           r->s = i2->s;
           return(b_basu(i2,i1,r));
           }            /* subtract first from second absolute value    */
/*----------------------------------------------------------------------*/
        if (cv)
           {
           r->s = i1->s;
           return(b_basu(i1,i2,r));
           }            /* subtract second from first absolute value    */
/*----------------------------------------------------------------------*/
        if (r->l) { r->l = 0; B_FREE(r->m) }
        r->z = 1;
        r->r = 0;       /* subtraction of equal absolute values         */
/*----------------------------------------------------------------------*/
        return(0);
        }





