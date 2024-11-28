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

/* CVS $Id: b_baad.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_baad.c          Processor : C                   */
/*                                                                      */
/* Addition of absolute values of intern numbers                        */
/* abs(first number) > abs(second number) > 0                           */
/*                                                                      */
/* Function value : int     0 - absolute value of intern numbers added  */
/*                      OFLOW - exponent overflow                       */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Rounding : towards zero                                              */
/*                                                                      */
/************************************************************************/

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
local int b_baad(multiprecision i1,multiprecision i2,multiprecision r)
#else
local int b_baad(i1,i2,r)

multiprecision i1;
multiprecision i2;      /* pointer to intern variables          */
multiprecision r;       /* pointer to result variable           */
#endif
        {
        a_intg k;       /* loop variable                         */
        a_intg l;       /* length of significant second mantissa */
        a_intg m;       /* length of resultant mantissa          */
        a_intg cv;      /* exponent difference                   */
        a_intg carry;   /* carry from choped addition            */
        a_intg rr;      /* non-zero mantissa                     */
        a_btyp *s;      /* pointer to result mantissa            */
/*---------------------------------------------------------------*/
        r->z = 0;       /* result is non-zero                    */
/*---------------------------------------------------------------*/
        cv = i1->e-i2->e;
                        /* exponent difference                   */
/*---------------------------------------------------------------*/
        l = cv+i2->l;   /* intermediate value                    */
/*---------------------------------------------------------------*/
        if (i1->l>=b_maxl) m = b_maxl;
        else if (l>=b_maxl) m = b_maxl;
        else if (i1->l>l) m = i1->l+1;
        else m = l+1;
                        /* resultant length of mantissa          */
/*---------------------------------------------------------------*/
        if (b_ball(m,&s)) return(ALLOC);
                        /* allocate and clear result mantissa    */
/*---------------------------------------------------------------*/
        B_COPY(s,i1->m,(i1->l>m)?m:i1->l)
                        /* copy mantissa of first operand to result     */
/*----------------------------------------------------------------------*/
        r->e = i1->e;   /* copy exponent of i1 to r                     */
/*----------------------------------------------------------------------*/
        if (cv>=m)      /* exponents differ by at least m               */
                        /* mantissas do not overlap                     */
           if (i1->l<=cv)
              r->r = 1; /* second mantissa is non-zero                  */
           else         /* mantissas overlap                            */
              {
              if (b_bmat((i1->l<l)?
                         i1->l-cv:
                         i2->l,&i1->m[cv],i2->m,0,&rr))
                 if (b_bcat(cv-m,&i1->m[m]))
                    {
                    if (b_bcad(m,s))
                       {
                       for (k=m-1;k>0;k--) *(s+k) = 0;
                       *s = 1;
                       if (r->e==MAXINT) return(OFLOW);
                       r->e++;
                       }
                    if (!rr)
                    rr = (i1->l<l) ? b_bmts(l-i1->l,&i2->m[i1->l])
                                   : b_bmts(i1->l-l,&i1->m[l]);
                    }
                 else rr = 1;
              else rr = 1;
              r->r = (unsigned int)rr;
              }
        else
           {
                        /* test choped digits of operands       */
           if (i1->l<=m)
              {
              r->r = b_bmts(l-m,&i2->m[m-cv]);
              carry = 0;
              }
           else
              {
                        /* determine carry from overlapping part */
              if (i1->l<l) {
                 carry = b_bmat(i1->l-m,&i1->m[m],&i2->m[m],0,&rr);
                 if ((r->r = (unsigned int)rr)==0)
                    r->r = b_bmts(l-i1->l,&i2->m[i1->l]);
                 }
              else {
                 carry = b_bmat(l-m,&i1->m[m],&i2->m[m],0,&rr);
                 if ((r->r = (unsigned int)rr)==0)
                    r->r = b_bmts(i1->l-l,&i1->m[l]);
                 }
              if (carry)
                 carry = b_bcad(m,s);
              }

           if (b_addm((l>m)?m-cv:i2->l,s+cv,i2->m))
           if (b_bcad(cv,s) || carry)
              {
              if (*(s+m-1)) r->r = 1;
              for (k=m-1;k>0;k--) *(s+k) = *(s+k-1);
              *s = 1;
              if (r->e==MAXINT) return(OFLOW);
              r->e++;
              }
           }
/*--------------------------------------------------------------*/
        if (r->l)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_baad");
#endif
           B_FREE(r->m)
           }
        r->l = m;
        r->m = s;
/*--------------------------------------------------------------*/
        while (s[r->l-1]==ZERO) r->l--;
                        /* remove trailing zeros by reducing r->l */
/*--------------------------------------------------------------*/
        return(0);
        }





