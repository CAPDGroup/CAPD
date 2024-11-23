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

/* CVS $Id: b_ball.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_ball.c          Processor : C                   */
/*                                                                      */
/* Allocation of a_btyp array and clear allocated array.                */
/*                                                                      */
/* Function value : int     0 - a_btyp array allocated                  */
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
local int b_ball(a_btyp n,a_btyp **i)
#else
local int b_ball(n,i)

a_btyp n;             /* length of mantissa to be allocated             */
a_btyp **i;           /* address of pointer to a_btyp array             */
#endif
        {
        /* possible allignment problem in case of incorrect calloc()    */
        if ((*i = (a_btyp *)calloc((size_t)n,sizeof(a_btyp)))==NULL)
           return(ALLOC);
#ifdef HEAP_CHECK
b_geth((a_char *)i,(a_char *)(*i),(a_char *)"b_ball");
#endif
        return(0);
        }





