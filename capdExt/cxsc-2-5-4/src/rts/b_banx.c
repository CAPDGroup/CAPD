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

/* CVS $Id: b_banx.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_banx.x          Processor : C                   */
/*                                                                      */
/* Determine next number of absolute value of non-zero number           */
/*                                                                      */
/* Function value : int     0 - next number returned                    */
/*                      OFLOW - exponent overflow                       */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Rounding : none                                                      */
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
local int b_banx(multiprecision i,multiprecision r)
#else
local int b_banx(i,r)

multiprecision i;       /* pointer to intern variable                   */
multiprecision r;       /* pointer to intern result                     */
#endif
        {
        a_intg k;
        a_btyp *s;    /* pointer to result mantissa                   */
/*----------------------------------------------------------------------*/
        if (i==r)                       /* input = output               */
           {
           if (i->l!=b_maxl)            /* input differs from b_maxl    */
              {
                                        /* allocate result mantissa     */
              if (b_ball(b_maxl,&s)) return(ALLOC);
              if (b_maxl<i->l)
                 B_COPY(s,i->m,b_maxl)
              else
                 {
                 for (k=b_maxl-1;k>=i->l;k--) s[k] = ZERO;
                 B_COPY(s,i->m,i->l)
                 }
                                        /* free old mantissa            */
#ifdef HEAP_CHECK
b_freh((a_char *)&i->m,(a_char *)i->m,(a_char *)"b_banx");
#endif
              B_FREE(i->m)
              r->m = s;
              r->l = b_maxl;
              }
                                        /* determine next number        */
           if (b_bcad(b_maxl,r->m))
              {
              if (r->e>MAXINT-1) return(OFLOW);
              else
                 {
                 r->e++;
                 r->m[0] = LSB;
                 }
              }
           }
        else                            /* input differs from output    */
           {
           if (r->l==b_maxl)
              B_CLEAR(r->m,b_maxl)
           else
              {
              if (r->l)
                 {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_banx");
#endif
                 B_FREE(r->m)
                 }
              if (b_ball(b_maxl,&r->m)) return(ALLOC);
              r->l = b_maxl;
              }
           r->e = i->e;
           B_COPY(r->m,i->m, ((b_maxl<i->l) ? b_maxl : i->l) )
           if (b_bcad(b_maxl,r->m))
              {
              if (r->e>MAXINT-1) return(OFLOW);
              else { r->e++; *r->m = 1; }
              }
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





