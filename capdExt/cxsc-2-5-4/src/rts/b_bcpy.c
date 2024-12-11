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

/* CVS $Id: b_bcpy.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bcpy.c          Processor : C                   */
/*                                                                      */
/* Copy number in intern representation                                 */
/*                                                                      */
/* Function value : int     0 - data copied                             */
/*                      ALLOC - allocation error                        */
/*                                                                      */
/* Function references :                                                */
/*                       b_ball - allocate unsigned array               */
/*                       b_bmts - test unsigned arrat to be non-zero    */
/*                                                                      */
/* Rounding : none                                                      */
/*                                                                      */
/********************************************************/

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
local int b_bcpy(multiprecision i,multiprecision r)
#else
local int b_bcpy(i,r)

multiprecision i;
multiprecision r;     /* pointer to intern variables                  */
#endif
        {
        a_btyp *s;    /* pointer to result mantissa                   */
        a_intg m;          /* length of result mantissa                    */
/*----------------------------------------------------------------------*/
        if ((r->z = i->z)!=0)
           {
           r->r = 0;
           if (r->l)
              {
              r->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_bcpy");
#endif
              B_FREE(r->m)
              }
           return(0);
           }
                        /* number to be copied is zero                  */
/*----------------------------------------------------------------------*/
        m = (i->l>b_maxl) ? b_maxl : i->l;
/*----------------------------------------------------------------------*/
        if (i!=r)
           {
           r->e = i->e;
           r->s = i->s;
           }
        else if (i->l==m) { r->r = 0; return(0); }
/*----------------------------------------------------------------------*/
        r->r = b_bmts(i->l-m,i->m+m);
/*----------------------------------------------------------------------*/
        if (b_ball(m,&s)) return(ALLOC);
/*----------------------------------------------------------------------*/
        B_COPY(s,i->m,m)
/*----------------------------------------------------------------------*/
        if (r->l)
           {
#ifdef HEAP_CHECK
b_freh((a_char *)&r->m,(a_char *)r->m,(a_char *)"b_bcpy");
#endif
           B_FREE(r->m)
           }
        r->l = m;
        r->m = s;
/*----------------------------------------------------------------------*/
        return(0);
        }





