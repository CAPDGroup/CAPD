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

/* CVS $Id: b_bmun.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bmun.c          Processor : C                   */
/*                                                                      */
/* Multiply number in internal format with unsigned integer             */
/*                                                                      */
/* Function value : int     0 - multiplication done                     */
/*                      OFLOW - overflow                                */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Rounding : towards zero                                              */
/*                                                                      */
/*                   a_btyp *cp -> dotprecision cp                      */
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
local int b_bmun(multiprecision i,a_btyp n,multiprecision r)
#else
local int b_bmun(i,n,r)

multiprecision i;       /* pointer to intern variables                  */
a_btyp n;             /* unsigned integer factor                      */
multiprecision r;       /* pointer to intern result                     */
#endif
        {
        a_intg ii;         /* current length of mantissa part              */
        a_intg m;          /* length of exact product                      */
        a_intg k;          /* loop variable                                */
        a_intg cv;         /* count value                                  */
        a_btyp *s;    /* array holding complete product of mantissas  */
/*----------------------------------------------------------------------*/
        r->r = 0;
                        /* no rounding expected                         */
/*----------------------------------------------------------------------*/
        if (i->z || n==0) { r->z = 1; return(0); }
                        /* result is zero                               */
/*----------------------------------------------------------------------*/
        r->z = 0;
                        /* result is non-zero                           */
/*----------------------------------------------------------------------*/
        r->s = i->s;
                        /* sign of product                              */
/*----------------------------------------------------------------------*/
        if (i->e>=0)
        if (MAXINT-i->e<1) return(OFLOW);
                        /* check overflow                               */
/*----------------------------------------------------------------------*/
        r->e = i->e+1;
                        /* determine exponent of first part of product  */
/*----------------------------------------------------------------------*/
        m = i->l+1;     /* length of exact product                      */
/*----------------------------------------------------------------------*/
        if (m>A_LENGTH)
           {
           if (b_ball(m,&s)) return(ALLOC);
           }
        else
           {
           s = b_cp__;
           B_CLEAR(s,m)
           }
/*----------------------------------------------------------------------*/
        for (k=i->l-1;k>=0;k--)
           b_muad(i->m[k],n,s+k+1);
                        /* multiply and add product                     */
/*----------------------------------------------------------------------*/
        cv = 0;         /* digit count                                  */
/*----------------------------------------------------------------------*/
        if (!*s)
           {
           r->e--;
           cv = 1;
           }            /* adjust leading digit to be non-zero          */
/*----------------------------------------------------------------------*/
        ii = (m-cv<b_maxl) ? m-cv : b_maxl;
/*----------------------------------------------------------------------*/
        b_badj(ii,r);
/*----------------------------------------------------------------------*/
        for (k=0;k<ii;k++) r->m[k] = s[cv+k];
                        /* copy first result mantissa                   */
/*----------------------------------------------------------------------*/
        r->r = b_bmts(m-ii-cv,s+ii+cv);
                        /* result rounded                               */
/*----------------------------------------------------------------------*/
        if (m>A_LENGTH)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&s,(a_char *)s,(a_char *)"b_bmun");
#endif
           B_FREE(s)
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





