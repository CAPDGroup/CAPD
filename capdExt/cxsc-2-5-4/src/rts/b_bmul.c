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

/* CVS $Id: b_bmul.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bmul.c          Processor : C                   */
/*                                                                      */
/* Multiply numbers in internal format                                  */
/*                                                                      */
/* Function value : int     0 - multiplication done                     */
/*                      OFLOW - overflow                                */
/*                      UFLOW - underflow                               */
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
extern dotprecision b_cp__;
#endif

#ifdef LINT_ARGS
local int b_bmul(multiprecision i1,multiprecision i2,multiprecision r1,
                 multiprecision r2)
#else
local int b_bmul(i1,i2,r1,r2)

multiprecision i1;
multiprecision i2;   /* pointer to intern variables                  */
multiprecision r1;
multiprecision r2;   /* pointer to intern result                     */
#endif
        {
        a_intg i;          /* current length of mantissa part              */
        a_intg m;          /* length of exact product                      */
        a_intg k;          /* loop variable                                */
        a_intg cv;         /* count value                                  */
        a_btyp *s;    /* array holding complete product of mantissas  */
/*----------------------------------------------------------------------*/
        r1->r = r2->r = 0;
                        /* no rounding expected                         */
/*----------------------------------------------------------------------*/
        if (i1->z || i2->z) { r1->z = r2->z = 1; return(0); }
                        /* result is zero                               */
/*----------------------------------------------------------------------*/
        r1->z = r2->z = 0;
                        /* result is non-zero                           */
/*----------------------------------------------------------------------*/
        r1->s = r2->s = (i1->s==i2->s) ? 0 : 1;
                        /* sign of product                              */
/*----------------------------------------------------------------------*/
        if (i1->e>=0)
           {
           if ((MAXINT-i1->e)-1<i2->e) return(OFLOW);
           }
        else
           {
           if (MININT-i1->e>i2->e) return(UFLOW);
           }
                        /* check underflow/overflow                     */
/*----------------------------------------------------------------------*/
        r2->e = r1->e = i1->e+i2->e+1;
                        /* determine exponent of first part of product  */
/*----------------------------------------------------------------------*/
        m = i1->l+i2->l;/* length of exact product mantissa             */
/*----------------------------------------------------------------------*/
        if (m>A_LENGTH)
           {
           if (b_ball(m,&s)) return(ALLOC);
           }
        else
           {
           s = &b_cp__[0];
           B_CLEAR(s,m)
           }
/*----------------------------------------------------------------------*/
        for (k=i2->l-1;k>=0;k--)
        if (i2->m[k])
           {
           for (i=i1->l-1;i>=0;i--)
              {
              if (i1->m[i])
                 b_muad(i1->m[i],i2->m[k],s+k+i+1);
              }
           }
                        /* multiply and add products                    */
/*----------------------------------------------------------------------*/
        cv = 0;         /* digit count                                  */
/*----------------------------------------------------------------------*/
        if (!*s)
           {
           if (r1->e<MININT+1) return(UFLOW);
           r1->e--;
           cv = 1;
           }            /* adjust leading digit to be non-zero          */
/*----------------------------------------------------------------------*/
        i = (m-cv<b_maxl)?m-cv:b_maxl;
/*----------------------------------------------------------------------*/
        if (b_badj(i,r1)) return(ALLOC);
/*----------------------------------------------------------------------*/
        for (k=0;k<i;k++) r1->m[k] = s[cv+k];
        cv += i;
                        /* copy first result mantissa                   */
/*----------------------------------------------------------------------*/
        while (cv<m) if (s[cv]) break; else cv++;
                        /* count leading zero mantissas of second part  */
/*----------------------------------------------------------------------*/
        if (cv==m)
           r2->z = 1;   /* second part of mantissa is zero              */
        else
           {
           r1->r = 1;   /* first mantissa does not hold exact product   */
           if (r2->e<MININT+cv) return(UFLOW);
           r2->e -= cv;
           i = (m-cv<b_maxl)?m-cv:b_maxl;
           if (b_badj(i,r2)) return(ALLOC);
           for (k=0;k<i;cv++,k++) r2->m[k] = *(s+cv);
           r2->r = b_bmts(m-cv,s+cv);
           }            /* copy second result mantissa                  */
/*----------------------------------------------------------------------*/
        if (m>A_LENGTH)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&s,(a_char *)s,(a_char *)"b_bmul");
#endif
           B_FREE(s)
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





