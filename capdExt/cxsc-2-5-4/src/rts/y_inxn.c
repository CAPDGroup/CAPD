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

/* CVS $Id: y_inxn.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_inxn.c                              */
/*                                                              */
/*      Entry           : a_VOID y_inxn(d,e_args)               */
/*                        y_dscp d;                             */
/*                        e_list                                */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                        e_args - variable argument list       */
/*                                                              */
/*      Return value    : char pointer to indexed element       */
/*                                                              */
/*      Description     : determine pointer to indexed element  */
/*                        without index checking.               */
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
local a_VOID y_inxn(y_dscp d,...)
#else
local a_VOID y_inxn(d,e_args)

y_dscp d;
e_list
#endif
#else
local a_VOID y_inxn(d,e_args)

y_dscp d;
e_list
#endif

        {
        a_btyp v = 0;
        a_intg i;
        va_list e_argv;

        e_open(d);

        for (i=0;i<((y_desc *)d)->numdim;i++)
           /* possible pointer allignment problems              */
           v += (e_ref(a_intg)-
((y_desc *)d)->fd[i].lbound)*((y_desc *)d)->fd[i].stride;

        e_shut;

        return((char *)(((y_desc *)d)->array)+v*((y_desc *)d)->elsize);
        }





