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

/* CVS $Id: b_basu.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_basu.c          Processor : C                   */
/*                                                                      */
/* Subtraction of absolute values of intern numbers                     */
/* abs(first number) > abs(second number) > 0                           */
/*                                                                      */
/* Function value : int     0 - absolute value of intern numbers added  */
/*                      ALLOC - allocation error                        */
/*                      UFLOW - exponent underflow                      */
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
extern dotprecision b_cp__;
#endif

#ifdef LINT_ARGS
local int b_basu(multiprecision i1,multiprecision i2,multiprecision r)
#else
local int b_basu(i1,i2,r)

multiprecision i1;
multiprecision i2;      /* pointer to intern variables                  */
multiprecision r;       /* pointer to result variable                   */
#endif
        {
        a_intg i;          /* loop variable                                */
        a_intg k;          /* loop variable                                */
        a_intg m;          /* length of resultant mantissa                 */
        a_intg cv;         /* exponent difference                          */
        a_btyp *s;    /* pointer to intermediate result               */
/*----------------------------------------------------------------------*/
        cv = i1->e-i2->e;
/*----------------------------------------------------------------------*/
        r->z = 0;       /* result is non-zero                           */
        r->e = i1->e;   /* result exponent                              */
/*----------------------------------------------------------------------*/
        if ( (m = (b_maxl>i1->l) ? b_maxl+1 : i1->l+1) <=cv)
           {
                        /* mantissas do not overlap                     */
           if (m>A_LENGTH)
              {
              if (b_ball(m,&s)) return(ALLOC);
              }
           else
              {
              s = b_cp__;
              B_CLEAR(s,m)
              }         /* intermediate result                          */
           B_COPY(s,i1->m,i1->l)
                        /* copy mantissa                                */
           b_bcsu(m,s);
                        /* subtract least significant bit               */
           k = 0;
           if (*s==0)
              {
              if (r->e<MININT+1) return(UFLOW);
              r->e--;
              k = 1;
              }
           r->r = 1;
           b_badj(b_maxl,r);              /* adjust result                */
           for (i=0;i<b_maxl;i++) r->m[i] = s[k+i];
                                        /* copy result mantissa         */
           if (s!=b_cp__)
              {
#ifdef HEAP_CHECK
b_freh((a_char *)&s,(a_char *)s,(a_char *)"b_basu");
#endif
              B_FREE(s)
              }
           return(0);
           }
/*----------------------------------------------------------------------*/
        m = (cv+i2->l<i1->l) ? i1->l : cv+i2->l;
                        /* length of intermediate result                */
/*----------------------------------------------------------------------*/
        if (m>A_LENGTH)
           {
           if (b_ball(m,&s)) return(ALLOC);
           }
        else
           {
           s = b_cp__;
           B_CLEAR(s,m)
           }            /* intermediate result                          */
/*----------------------------------------------------------------------*/
        B_COPY(s,i1->m,i1->l)
                        /* copy mantissa                                */
/*----------------------------------------------------------------------*/
        if (b_subm(i2->l,s+cv,i2->m)) b_bcsu(cv,s);
                        /* subtract mantissas                           */
/*----------------------------------------------------------------------*/
        k = 0;
        while (*(s+k)==0) k++;
        if (k)
           {
           if (r->e<MININT+k) return(UFLOW);
           r->e -= k;
           m -= k;
           }            /* adjust exponent                              */
/*----------------------------------------------------------------------*/
        if (m>b_maxl)
           {
           r->r = b_bmts(m-b_maxl,s+k+b_maxl);
           m = b_maxl;
           }
        else
           r->r = 0;    /* test for rounded result                      */
/*----------------------------------------------------------------------*/
        while (s[k+m-1]==ZERO) m--;
                        /* eliminate trailing zeros from mantissa       */
/*----------------------------------------------------------------------*/
        b_badj(m,r);    /* adjust result                                */
/*----------------------------------------------------------------------*/
        for (i=0;i<m;i++) r->m[i] = s[k+i];
                        /* copy result                                  */
/*----------------------------------------------------------------------*/
        if (s!=b_cp__)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&s,(a_char *)s,(a_char *)"b_basu");
#endif
           B_FREE(s)
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





