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

/* CVS $Id: b_bcmp.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bcmp.c          Processor : C                   */
/*                                                                      */
/* Compare numbers in intern representation                             */
/*                                                                      */
/* Function value : int     0 - numbers are equal                       */
/*                         -1 - first number less than second number    */
/*                          1 - first number greater than second number */
/*                                                                      */
/* Function references : b_bacm - compare positive numbers              */
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
local int b_bcmp(multiprecision i1,multiprecision i2)
#else
local int b_bcmp(i1,i2)

multiprecision i1;
multiprecision i2;   /* pointer to intern variables                  */
#endif
        {
        if (i1->z)
           {
           if (i2->z) return(0);
           else return((int)(2*i2->s-1));
           }
        if (i2->z || i1->s!=i2->s) return((int)(1-2*i1->s));
        return((int)(1-2*i1->s)*b_bacm(i1,i2));
        }





