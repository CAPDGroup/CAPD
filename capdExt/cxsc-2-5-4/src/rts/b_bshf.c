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

/* CVS $Id: b_bshf.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/* Descriptive Name : b_bshf.c          Processor : C           */
/*                                                              */
/* Multiply intern number by power of two                       */
/*                                                              */
/* Function value : int     0 - intern numbers added            */
/*                      OFLOW - exponent overflow               */
/*                      UFLOW - exponent underflow              */
/*                      ALLOC - allocation error                */
/*                                                              */
/* Globals referenced : Maxl                                    */
/*                                                              */
/* Rounding : towards zero                                      */
/*                                                              */
/*                   reduce length if trailing zero mantissa    */
/*                   explicit compare                           */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_btyp b_maxl;
#endif

#ifdef LINT_ARGS
local int b_bshf(a_intg n,multiprecision i,multiprecision r)
#else
local int b_bshf(n,i,r)

a_intg n;               /* number of shifts                     */
multiprecision i;       /* pointer to intern variable           */
multiprecision r;       /* pointer to intern result             */
#endif
        {
        a_intg k;       /* loop variable                        */
        a_intg m;       /* length of array to be allocated      */
        a_btyp *s;      /* allocated array                      */
/*--------------------------------------------------------------*/
        r->r = 0;       /* initially no rounding                */
/*--------------------------------------------------------------*/
        if ((r->z = i->z)!=0) return(0);
                        /* number is zero                       */
/*--------------------------------------------------------------*/
        r->s = i->s;    /* copy sign of number                  */
/*--------------------------------------------------------------*/
        k = B_ASHR(n,LOG_B_LENGTH);
                /* number of unsigned to be shifted (virtually) */
        n -= k*B_LENGTH;
/*--------------------------------------------------------------*/
        if (k>0)
           {
           if (MAXINT-k<i->e) return(OFLOW);
           }
        else if (k<0)
           {
           if (MININT-k>i->e) return(UFLOW);
           }            /* underflow/overflow checking          */
/*--------------------------------------------------------------*/
        r->e = i->e+k;  /* determine resultant exponent         */
/*--------------------------------------------------------------*/
        m = (i->l<b_maxl) ? i->l : b_maxl;
/*--------------------------------------------------------------*/
        if (n>0)
           {
           if (i->m[0]>>(B_LENGTH-n))
              {
              if (r->e>MAXINT-1) return(OFLOW);
              r->e++;
              n -= B_LENGTH;
              }
           }
        else if (n<0)
           if ((i->m[0]>>(-n))==0)
              {
              if (r->e<MININT+1) return(UFLOW);
              r->e--;
              n += B_LENGTH;
              }         /* consider bits in leading digit       */
/*--------------------------------------------------------------*/
        if (n<0)
           {
           if (m<b_maxl && i->m[m-1]<<(B_LENGTH+n))
              m++;
           else if (m>=b_maxl) {
              if (i->m[m-1]<<(B_LENGTH+n)) r->r = 1;
              else {
                 r->r = 0;
                 for (k=m;k<i->l;k++)
                 if (i->m[k]) {
                    r->r = 1;
                    break;
                    }
                 }
              }
           }
                /* some bits drop out of result mantissa        */
/*--------------------------------------------------------------*/
        if (b_ball(m,&s)) return(ALLOC);
                        /* allocate unsigned array              */
/*--------------------------------------------------------------*/
        if (n<0)
           {
           n = -n;
           if (m<=i->l)
              {
              if (m!=1)
                 *(s+m-1) = (i->m[m-1]>>n) |
                            (i->m[m-2]<<(B_LENGTH-n));
              }
           else
              *(s+m-1) = i->m[m-2]<<(B_LENGTH-n);
           for (k=m-2;k>0;k--)
              *(s+k) = (i->m[k]>>n) |
                       (i->m[k-1]<<(B_LENGTH-n));
           *s = i->m[0]>>n;
           }
        else if (n>0)
           {
           for (k=0;k<m-1;k++)
              *(s+k) = (i->m[k]<<n) | (i->m[k+1]>>(B_LENGTH-n));
           *(s+m-1) = i->m[m-1]<<n;
           if (i->l>m) *(s+m-1) |= (i->m[m]>>(B_LENGTH-n));
           }
        else
           {
           for (k=0;k<m;k++)
              *(s+k) = i->m[k];
           }
                        /* shift digits                         */
/*--------------------------------------------------------------*/
        if (r->l)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_bshf");
#endif
           B_FREE(r->m)
           }
        r->l = m;
        r->m = s;
                        /* update pointer of result             */
/*--------------------------------------------------------------*/
        while (s[r->l-1]==ZERO) r->l--;
                        /* remove trailing zeros by reducing r->l */
/*--------------------------------------------------------------*/
        return(0);
        }





