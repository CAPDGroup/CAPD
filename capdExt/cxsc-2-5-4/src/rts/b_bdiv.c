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

/* CVS $Id: b_bdiv.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/* Descriptive Name : b_bdiv.c          Processor : C           */
/*                                                              */
/* Divide numbers in internal format                            */
/*                                                              */
/* Function value : int     0 - division executed               */
/*                      OFLOW - overflow                        */
/*                      UFLOW - underflow                       */
/*                      ZEROD - divide by zero                  */
/*                      ALLOC - allocation error                */
/*                                                              */
/* Rounding : towards zero                                      */
/*                                                              */
/****************************************************************/

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
local int b_bdiv(multiprecision i1,multiprecision i2,multiprecision r)
#else
local int b_bdiv(i1,i2,r)

multiprecision i1;
multiprecision i2;      /* pointer to intern variables          */
multiprecision r;       /* pointer to intern result             */
#endif
        {
        a_intg k;          /* loop variable                     */
        a_intg cv;         /* compare value                     */
        a_intg l;          /* minimum of lengths of input numbers */
        a_intg m;          /* length of intermediate array      */
        a_btyp f;     /* factor                                 */
        a_btyp *q;    /* quotient                               */
        a_btyp *d;    /* dividend                               */
/*----------------------------------------------------------------------*/
        if (i2->z) return(ZEROD);
                        /* division by zero                             */
/*----------------------------------------------------------------------*/
        if ((r->z = i1->z)!=0) { r->r = 0; return(0); }
                        /* result is zero                               */
/*----------------------------------------------------------------------*/
        r->s = (i1->s==i2->s) ? 0 : 1;
                        /* sign of quotient                             */
/*----------------------------------------------------------------------*/
        if (i2->e>=0)
           {
           if (MININT+i2->e>i1->e) return(UFLOW);
           }
        else
           if (MAXINT+i2->e<i1->e) return(OFLOW);
                        /* underflow/overflow checking                  */
/*----------------------------------------------------------------------*/
        r->e = i1->e-i2->e;
                        /* determine exponent of quotient               */
/*----------------------------------------------------------------------*/
        l = (i1->l<i2->l)?i1->l:i2->l;
        if ((cv = b_bmcm(l,i1->m,i2->m))==0)
           cv = -b_bmts(i2->l-l,i2->m+l);
/*----------------------------------------------------------------------*/
        if (cv>=0) cv = 1;
        else
           {
           if (r->e<MININT+1) return(UFLOW);
           r->e--;
           cv = 0;
           }            /* adjust result exponent                       */
/*----------------------------------------------------------------------*/
                        /* reserved space for division operation        */
                        /* 1     : numerator larger than denominator    */
                        /* i1->l : first mantissa                       */
                        /* i2->l : length of second mantissa            */
                        /* 2     : force minimum length 3               */
                        /* b_maxl-1: at most b_maxl-1 subtractions      */
                        /*           required                           */
                        /* i2->l : second mantissa                      */
                        /* 1     : force minimum length 2               */
        if ( (m = b_maxl+i1->l+i2->l+2)+i2->l+1>A_LENGTH)
           {
           if (b_ball(m+i2->l+1,&d)) return(ALLOC);
           }
        else
           {
           d = b_cp__;
           B_CLEAR(d,m+i2->l+1)
           }
/*----------------------------------------------------------------------*/
        if (b_ball(b_maxl,&q)) return(ALLOC);
/*----------------------------------------------------------------------*/
        for (k=0;k<i1->l;k++) d[cv+k] = i1->m[k];
        for (k=0;k<i2->l;k++) d[m+k] = i2->m[k];
                        /* copy dividend and divisor                    */
/*----------------------------------------------------------------------*/
        f = *(d+m);
        for (k=0;(f & MSB)==ZERO;k++) f <<= 1;
        if (k)
           {
           b_bmsh(i1->l+1,d,k);
           b_bmsh(i2->l,d+m,k);
           }
                        /* shift dividend and divisor                   */
/*----------------------------------------------------------------------*/
        for (k=0;k<b_maxl;k++)
        if (*(d+k))
           {
           (void)b_bmdv(d+k,d+m,q+k);
           if (b_bmsp(i2->l,d+m,*(q+k),d+k))
              {
              (*(q+k))--;
              (void)b_addm(i2->l,d+k+1,d+m);
              *(d+k) = 0;
              }
           }
        else if (b_bmcm(i2->l,d+k+1,d+m)>=0)
           {
           (void)b_subm(i2->l,d+k+1,d+m);
           *(q+k) = 1;
           }
/*----------------------------------------------------------------------*/
        r->r = b_bmts(m-b_maxl,d+b_maxl);
/*----------------------------------------------------------------------*/
        if (q[(m = b_maxl-1)]==ZERO)
           {
           while (q[m-1]==ZERO) m--;
                                /* remove trailing blanks               */

           if (b_badj(m,r)) return(ALLOC);
           while (m-->0) r->m[m] = q[m];
                                /* copy mantissa                        */

#ifdef HEAP_CHECK
b_freh((a_char *)&q,(a_char *)q,(a_char *)"b_bdiv");
#endif
           B_FREE(q)
           }
        else
           {
           if (r->l)
              {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_bdiv");
#endif
              B_FREE(r->m)
              }
           r->l = b_maxl;
           r->m = q;
           }
/*----------------------------------------------------------------------*/
        if (d!=b_cp__)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&d,(a_char *)d,(a_char *)"b_bdiv");
#endif
           B_FREE(d)
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





