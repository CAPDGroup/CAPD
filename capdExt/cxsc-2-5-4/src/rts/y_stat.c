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

/* CVS $Id: y_stat.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_stat.c                              */
/*                                                              */
/*      Entry           : a_VOID y_stat(d,s,z,dim,e_args)       */
/*                        y_dscp d;                             */
/*                        a_VOID s;                             */
/*                        size_t z;                             */
/*                        a_byte dim;                           */
/*                        e_list                                */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                        s - static array                      */
/*                        z - size of static array              */
/*                        dim - dimension of static array       */
/*                        ... - list of a_intg pairs lb,ub      */
/*                                                              */
/*      Function value  : address of descriptor of dynamic array*/
/*                                                              */
/*      Description     : copy static array to dynamic array    */
/*                        (descriptor is initialized)           */
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
local a_VOID y_stat(y_dscp d,a_VOID s,size_t z,a_byte dim,...)
#else
local a_VOID y_stat(d,s,z,dim,e_args)

y_dscp d;
a_VOID s;
size_t z;
a_byte dim;
e_list
#endif
#else
local a_VOID y_stat(d,s,z,dim,e_args)

y_dscp d;
a_VOID s;
size_t z;
a_byte dim;
e_list
#endif
        {
        int k;        /* !!! must be int !!! */
        va_list e_argv;

        E_TPUSH("y_stat")

        /* destroy-flag must not be TRUE since address of       */
        /* static array is assigned to array-component of       */
        /* dynamic array descriptor                             */
        ((y_desc *)d)->destroy = ((y_desc *)d)->subarr = FALSE;

        ((y_desc *)d)->numdim = dim;
        ((y_desc *)d)->elsize = z;

        /* read lower and upper bounds from argument list       */
        e_open(dim);
        ((y_desc *)d)->fd[0].lbound = e_ref(a_intg);
        ((y_desc *)d)->fd[0].ubound = e_ref(a_intg);
        for (k=1;k<dim;k++)
           {
           ((y_desc *)d)->fd[k].lbound = e_ref(a_intg);
           ((y_desc *)d)->fd[k].ubound = e_ref(a_intg);
           }
        e_shut;

        /* determine stride and number of elements      */
        k = dim-1;
        ((y_desc *)d)->elnum = ((y_desc *)d)->fd[k].stride = 1;
        for (;k>0;k--)
           {
           ((y_desc *)d)->elnum *= (size_t) ((((y_desc *)d)->fd[k].ubound-
                                    ((y_desc *)d)->fd[k].lbound+1));
           ((y_desc *)d)->fd[k-1].stride = ((y_desc *)d)->elnum;
           }
        ((y_desc *)d)->elnum *= (size_t) ((((y_desc *)d)->fd[0].ubound-
                                 ((y_desc *)d)->fd[0].lbound+1));

        /* assign static array to array component       */
        ((y_desc *)d)->array = (char *)s;

        E_TPOPP("y_stat")
        return((char *)d);
        }





