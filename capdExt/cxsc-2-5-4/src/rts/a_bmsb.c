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

/* CVS $Id: a_bmsb.c,v 1.21 2014/01/30 17:24:01 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_bmsb.c                              */
/*                                                              */
/*      Entries         : a_intg a_bmsb(a)                      */
/*                        a_intg a;                             */
/*                                                              */
/*      Arguments       : a = integer bit operand               */
/*                                                              */
/*      Description     : Bit position of MSB of an integer     */
/*                                                              */
/*      Function value  : position of most significant bit      */
/*                        of an integer with position 0 is      */
/*                        right-most bit of integer.            */
/*                        -1 is returned if no bits are set.    */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local a_intg a_bmsb(a_intg a)
#else
local a_intg a_bmsb(a)

a_intg a;
#endif
        {
        a_intg res;

        E_TPUSH("a_bmsb")

        if (a<0) res = 31;
        else
           {
           res = -1;
           while (a)
              {
              a >>= 1;
              res++;
              }
           }

        E_TPOPP("a_bmsb")
        return(res);
        }





