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

/* CVS $Id: b_bdvn.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bdvn.c          Processor : C                   */
/*                                                                      */
/* Divide numbers in internal format                                    */
/*                                                                      */
/* Function value : int     0 - division executed                       */
/*                      OFLOW - overflow                                */
/*                      UFLOW - underflow                               */
/*                      ZEROD - divide by zero                          */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Rounding : towards zero                                              */
/*                                                                      */
/*                   cp,Maxl renamed                    */
/*                   a_btyp *cp -> dotprecision cp      */
/********************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_btyp b_maxl;
extern dotprecision b_cp__;
#endif

#ifdef LINT_ARGS
local int b_bdvn(multiprecision i1,a_btyp n,multiprecision r)
#else
local int b_bdvn(i1,n,r)

multiprecision i1;      /* pointer to intern variables                  */
a_btyp n;             /* integer divisor                              */
multiprecision r;       /* pointer to intern result                     */
#endif
        {
        a_intg k;          /* loop variable                                */
        a_intg cv;         /* compare value                                */
        a_intg m;          /* length of intermediate array                 */
        a_btyp *q;    /* quotient                                     */
        a_btyp *d;    /* dividend                                     */
/*----------------------------------------------------------------------*/
        if (n==0) return(ZEROD);
                        /* division by zero                             */
/*----------------------------------------------------------------------*/
        if (i1->z) { r->z = 1; r->r = 0; return(0); }
                        /* result is zero                               */
/*----------------------------------------------------------------------*/
        r->z = 0;       /* result is non-zero                           */
/*----------------------------------------------------------------------*/
        r->s = i1->s;
                        /* sign of quotient                             */
/*----------------------------------------------------------------------*/
        r->e = i1->e;
                        /* determine exponent of first part of product  */
/*----------------------------------------------------------------------*/
        if (*i1->m<n)
           {
           if (r->e<MININT+1) return(UFLOW);
           r->e--;
           cv = 0;
           }
        else cv = 1;
                        /* adjust result exponent                       */
/*----------------------------------------------------------------------*/
        if ( (m = b_maxl+i1->l+1)+2>A_LENGTH)
           {
           if (b_ball(m+2,&d)) return(ALLOC);
           }
        else
           {
           d = b_cp__;
           B_CLEAR(d,m+2)
           }
/*----------------------------------------------------------------------*/
        if (b_ball(b_maxl,&q)) return(ALLOC);
/*----------------------------------------------------------------------*/
        for (k=0;k<i1->l;k++) d[cv+k] = i1->m[k];
                        /* copy dividend                                */
/*----------------------------------------------------------------------*/
        for (k=0;(n & MSB)==0;k++) n <<= 1;
        if (k) b_bmsh(i1->l+1,d,k);
                        /* shift dividend and divisor                   */
/*----------------------------------------------------------------------*/
        *(d+m) = n;
                        /* copy divisor                                 */
/*----------------------------------------------------------------------*/
        for (k=0;k<b_maxl;k++)
        if (*(d+k))
           {
           (void)b_bmdv(d+k,d+m,q+k);
           (void)b_busp(n,*(q+k),d+k);
           }
        else if (*(d+k+1)>=n)
           {
           *(d+k+1) -= n;
           *(q+k) = 1;
           }
/*----------------------------------------------------------------------*/
        r->r = b_bmts(m-b_maxl,d+b_maxl);
/*----------------------------------------------------------------------*/
        if (r->l)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_bdvn");
#endif
           B_FREE(r->m)
           }
        r->l = b_maxl;
        r->m = q;
/*----------------------------------------------------------------------*/
        if (d!=b_cp__)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&d,(a_char *)d,(a_char *)"b_bdvn");
#endif
           B_FREE(d)
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





