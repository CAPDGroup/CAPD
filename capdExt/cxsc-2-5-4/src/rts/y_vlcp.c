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

/* CVS $Id: y_vlcp.c,v 1.21 2014/01/30 17:24:18 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_vlcp.c                              */
/*                                                              */
/*      Entry           : void y_vlcp(d)                        */
/*                        y_dscp d;                             */
/*                                                              */
/*      Arguments       : d - descriptor of dynamic array       */
/*                                                              */
/*      Description     : copy value-parameter if destroy-flag  */
/*                        is FALSE, otherwise set destroy-flag  */
/*                        FALSE                                 */
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
local void y_vlcp(y_dscp d)
#else
local void y_vlcp(d)

y_dscp d;
#endif
        {
        char *p;
        size_t n;
        a_intg i,k,*help,prod;

        E_TPUSH("y_vlcp")

        if (((y_desc *)d)->destroy==FALSE)
           {
           if ((p = (char *)malloc(n =
                    ((y_desc *)d)->elsize*((y_desc *)d)->elnum))==NULL)
              {
              e_trap(ALLOCATION,2,E_TMSG,42);
              }
           else
              {
#ifdef HEAP_CHECK
b_geth((a_char *)&((y_desc *)d)->array,(a_char *)p,(a_char *)"y_vlcp");
#endif
              /* check if d is allocated */
		    if (!y_aled(d) )
		    {
                    e_trap(ALLOCATION,2,E_TMSG,42);
                    E_TPOPP("y_vlcp")
                    return;
		    }
              if (((y_desc *)d)->subarr==FALSE)
                 C_COPY(p,(char *)((y_desc *)d)->array,n)

              /* copy individual components if source is subarray  */
              else
                 {

                 /* Allocate index array help.  */
                 if ((help =
            (a_intg *)malloc(((y_desc *)d)->numdim*sizeof(a_intg)))==NULL)
                    {
                    e_trap(ALLOCATION,2,E_TMSG,42);
                    E_TPOPP("y_vlcp")
                    return;
                    }
#ifdef HEAP_CHECK
b_geth((a_char *)&help,(a_char *)help,(a_char *)"y_vlcp");
#endif
                 for (i=0;i<((y_desc *)d)->numdim;i++)
                    help[i] = ((y_desc *)d)->fd[i].lbound;

                 for (k=0;k<((y_desc *)d)->elnum;k++)
                    {
                    prod = 0;
                    for (i=0;i<((y_desc *)d)->numdim;i++)
                       prod += (help[i]-((y_desc *)d)->fd[i].lbound)*
                               ((y_desc *)d)->fd[i].stride;

                    C_COPY((char *)p+k*((y_desc *)d)->elsize,
                           (char *)((y_desc *)d)->array+
prod*((y_desc *)d)->elsize,((y_desc *)d)->elsize)

                    for (i=((y_desc *)d)->numdim-1;i>=0;i--)
                       if (++help[i]>((y_desc *)d)->fd[i].ubound)
                          help[i] = ((y_desc *)d)->fd[i].lbound;
                       else break;
                    }
#ifdef HEAP_CHECK
b_freh((a_char *)&help,(a_char *)help,(a_char *)"y_vlcp");
#endif
                 B_FREE(help)

                 /* determine strides for destination   */
                 i = ((y_desc *)d)->numdim-1;
                 ((y_desc *)d)->fd[i--].stride = 1;
                 for (;i>=0;i--)
                    {
                    ((y_desc *)d)->fd[i].stride = (size_t) (
                                      (((y_desc *)d)->fd[i+1].ubound-
                                       ((y_desc *)d)->fd[i+1].lbound+1)*
                                      ((y_desc *)d)->fd[i+1].stride);
                    }

                 ((y_desc *)d)->subarr = FALSE;
                 }

              /* if (y_aled(d) ) y_free(d);*/
              ((y_desc *)d)->array = p;
              }
           }
        else
            ((y_desc *)d)->destroy = FALSE;

        E_TPOPP("y_vlcp")
        }





