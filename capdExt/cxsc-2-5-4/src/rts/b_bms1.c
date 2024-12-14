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

/* CVS $Id: b_bms1.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bms1.c          Processor : C                   */
/*                                                                      */
/* Subtract 1 digit number from 1 digit number                          */
/*                                                                      */
/* Function value : int       - borrow value                            */
/*                                                                      */
/* Function references : none                                           */
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
local int b_bms1(a_btyp i,a_btyp *r)
#else
local int b_bms1(i,r)

a_btyp i;             /* 1 digit number                               */
a_btyp *r;            /* 1 digit number                               */
#endif
        {
/*----------------------------------------------------------------------*/
        if (i<=*r)
           {
           *r -= i;
           return(0);
           }
#if C_P_2
        *r -= i;
#else
        *r = ( (*r & i & MSB) ? *r-(i^MSB) : (*r^MSB)-i ) ^ MSB;
#endif
        return(1);
        }





