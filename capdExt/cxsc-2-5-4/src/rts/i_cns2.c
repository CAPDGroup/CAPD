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

/* CVS $Id: i_cns2.c,v 1.21 2014/01/30 17:24:08 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : i_cns2.c                              */
/*                                                              */
/*      Entries         : a_intv i_cns2(stl,str)                */
/*                        a_char *stl;                          */
/*                        a_char *str;                          */
/*                                                              */
/*      Arguments       : stl - input string for lower bound    */
/*                        str - input string for upper bound    */
/*                                                              */
/*      Description     : Convert a character string            */
/*                        to IEEE double format rounded         */
/*                        to interval.                          */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local a_intv i_cns2(a_char *stl,a_char *str)
#else
local a_intv i_cns2(stl,str)

a_char *stl;
a_char *str;
#endif
        {
        a_intv s;
        a_char *next;

        E_TPUSH("i_cns2")

        R_ASSIGN(s.INF,*r_zero);
        R_ASSIGN(s.SUP,*r_zero);

        r_conv(stl,&s.INF,(a_intg)-1,&next);
        r_conv(str,&s.SUP,(a_intg)1,&next);

        E_TPOPP("i_cns2")
        return(s);
        }





