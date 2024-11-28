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

/* CVS $Id: b_badj.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_badj.c          Processor : C                   */
/*                                                                      */
/* Adjust length of mantissa of an intern variable                      */
/*                                                                      */
/* Function value : int     0 - intern variable allocated               */
/*                      ALLOC - allocation error                        */
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
local int b_badj(a_btyp n,multiprecision i)
#else
local int b_badj(n,i)

a_btyp n;               /* length of mantissa to be allocated           */
multiprecision i;       /* pointer to intern variable                   */
#endif
        {
/*----------------------------------------------------------------------*/
        if (n<=0)
           {
           if (i->l)
              {
              i->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&i->m,(a_char *)i->m,(a_char *)"b_badj");
#endif
              B_FREE(i->m);
              }
           }
        else if (i->l!=n)
           {
           if (i->l)
              {
              i->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&i->m,(a_char *)i->m,(a_char *)"b_badj");
#endif
              B_FREE(i->m);
              }
           if (b_ball(n,&i->m)) return(ALLOC);
           i->l = n;
           }
/*----------------------------------------------------------------------*/
        return(0);
        }





