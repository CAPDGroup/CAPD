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

/* CVS $Id: y_inid.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_inid.c                              */
/*                                                              */
/*      Entry           : void y_inid(d,dim,elsize)             */
/*                        y_dscp d;                             */
/*                        a_byte dim;                           */
/*                        size_t elsize;                        */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                        dim - dimension of dynamic array      */
/*                        elsize - size of an element in bytes  */
/*                                                              */
/*      Description     : initialization of a dynamic array     */
/*                        same as y _init but no allocation     */
/*                        'array' should have been allocated    */
/*                        elnum =0, strid init 1,0,1,0 ....     */
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
local void y_inid(y_dscp d,a_byte dim,size_t elsize)
#else
local void y_inid(d,dim,elsize)

y_dscp d;
a_byte dim;
size_t elsize;
#endif
{
        int i;  /* !!! must be int !!! */

        E_TPUSH("y_inid")

        ((y_dscp)d)->subarr = ((y_dscp)d)->destroy = FALSE;
        ((y_dscp)d)->elsize = elsize;
        ((y_dscp)d)->numdim = dim;
        ((y_dscp)d)->elnum  = 0;
        ((y_dscp)d)->fd[dim-1].stride = 1;
        ((y_dscp)d)->array  = NULL;
        for (i=dim-2;i>=0;i--)
        {
           ((y_dscp)d)->fd[i].stride = (((y_dscp)d)->fd[i+1].ubound-
                                        ((y_dscp)d)->fd[i+1].lbound+1)*
                                        ((y_dscp)d)->fd[i+1].stride;
        }
        E_TPOPP("y_inid")
        return;
}





