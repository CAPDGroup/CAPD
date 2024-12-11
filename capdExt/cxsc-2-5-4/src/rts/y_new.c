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

/* CVS $Id: y_new.c,v 1.22 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_new.c                               */
/*                                                              */
/*      Entry           : void   y_new (d,e_args)               */
/*                        y_dscp d;                             */
/*                        e_list                                */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                        e_args - variable argument list       */
/*                                                              */
/*      Return value    : void                                  */
/*                                                              */
/*      Description     : allocate dynamic array, bounds are=b= */
/*                        trailing.                             */
/*                        check allocation of d, free if any,   */
/*                        check bounds and                      */
/*                        calls y_ init                         */
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
#if C_P_7
local void   y_new (y_dscp d,...)
#else
local void   y_new (d,e_args)
y_dscp d;
e_list
#endif
#else
local void   y_new (d,e_args)

y_dscp d;
e_list
#endif

{
/*        a_btyp v = 0; */
        a_intg i;
        va_list e_argv;


        /* check if field is already allocated */
        y_free(d);

        e_open(d);

        /* assign bounds for every dimension */
        for (i=0;i<((y_dscp)d)->numdim;i++)
	   {
           ((y_dscp)d)->fd[i].lbound= e_ref(a_intg);
           ((y_dscp)d)->fd[i].ubound= e_ref(a_intg);

           /* check lbound less equal ubound */
           if (((y_dscp)d)->fd[i].lbound > ((y_dscp)d)->fd[i].ubound)
              e_trap(INDEX_RANGE,6,
                  E_TINT+E_TEXT(4),&i,
                  E_TINT+E_TEXT(5),&((y_dscp)d)->fd[i].lbound,
                  E_TINT+E_TEXT(6),&((y_dscp)d)->fd[i].ubound);
        }


        y_init(d,d->numdim,d->elsize);

        e_shut; /* end of dynamic arg list */
        E_TPOPP("y_new")
        return;
}





