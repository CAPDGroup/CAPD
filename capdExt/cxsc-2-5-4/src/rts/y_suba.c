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

/* CVS $Id: y_suba.c,v 1.21 2014/01/30 17:24:18 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_suba.c                              */
/*                                                              */
/*      Entry           : a_VOID y_suba(m,s,mode,e_args)        */
/*                        a_VOID m;                             */
/*                        a_VOID s;                             */
/*                        a_char *mode;                         */
/*                        e_list                                */
/*                                                              */
/*      Arguments       : m - descriptor of main array          */
/*                        s - descriptor of subarray            */
/*                        mode - subarray mode                  */
/*                           "0" - index range of main array    */
/*                           "1" - constant index               */
/*                           "2" - index range from arg. list   */
/*                           "3" - lbound from arg. list        */
/*                           "4" - ubound from arg. list        */
/*                                                              */
/*                        e_args - variable argument list       */
/*                                                              */
/*      Note            : subarrays are never temporary         */
/*                                                              */
/*      Description     : generation of a subarray descriptor   */
/*                                                              */
/*                   cast to size_t at elnum asignment          */
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
local a_VOID y_suba(a_VOID m,a_VOID s,a_char *mode,...)
#else
local a_VOID y_suba(m,s,mode,e_args)
a_VOID m;
a_VOID s;
a_char *mode;
e_list
#endif
#else
local a_VOID y_suba(m,s,mode,e_args)
a_VOID m;
a_VOID s;
a_char *mode;
e_list
#endif

{
        a_byte dim,sdim,i,j;
        a_byte length = (a_byte) strlen((char*) mode);
        a_btyp offset;
        va_list e_argv;

        E_TPUSH("y_suba")

        e_open(mode);

        /* dimension of main array */
        sdim = dim = ((y_desc *)m)->numdim;

        /* computation of the new subarray descriptor out of    */
        /* the main array descriptor                            */

        /* j = actual index position of subarray */
        offset = j = i = 0;

        /* scan characters of mode string */
        for (;i<length;i++)
           {

           /* get bounds and stride of index range for subarray index */
           switch (mode[i])
              {

              /* the whole index range of the main array        */
              /* is taken for the subarray                      */
              case '0': ((y_desc *)s)->fd[j].lbound =
                           ((y_desc *)m)->fd[i].lbound;
                        ((y_desc *)s)->fd[j].ubound =
                           ((y_desc *)m)->fd[i].ubound;
                        ((y_desc *)s)->fd[j].stride =
                           ((y_desc *)m)->fd[i].stride;
                        j++;
                        break;

              /* no subarray dimension in case of a single      */
              /* constant index value                           */
              /* possible pointer allignment problems           */
              case '1': 
				  { a_intg index = e_ref(a_intg);

				    if (((y_dscp)m)->fd[i].lbound > index ||
					   ((y_dscp)m)->fd[i].ubound < index)
                        {
                          e_trap(INDEX_RANGE,8,E_TMSG,67,
                                 E_TINT+E_TEXT(11),&index,
                                 E_TINT+E_TEXT(5),&((y_desc *)m)->fd[i].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)m)->fd[i].ubound);
                          e_shut;
                          E_TPOPP("y_suba")
                          return((a_VOID) s);
			         }	

				    offset += (index-((y_desc *)m)->fd[i].lbound)*
                                  ((y_desc *)m)->fd[i].stride;
                        sdim--;
                        break;
                      }
              /* both indexranges are taken from the arg_list   */
              /* possible pointer allignment problems           */
              case '2': ((y_desc *)s)->fd[j].lbound = e_ref(a_intg);

              /* possible pointer allignment problems           */
                        ((y_desc *)s)->fd[j].ubound = e_ref(a_intg);
				    if (((y_dscp)s)->fd[j].lbound < 
					   ((y_dscp)m)->fd[i].lbound ||
				        ((y_dscp)s)->fd[j].ubound > 
					   ((y_dscp)m)->fd[i].ubound )
                        {
					 a_intg k;
					 k=j+1;
                          e_trap(INDEX_RANGE,12,E_TMSG,67,
                                 E_TINT+E_TEXT(11),&k,
                                 E_TINT+E_TEXT(5),&((y_desc *)s)->fd[j].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)s)->fd[j].ubound,
                                 E_TINT+E_TEXT(5),&((y_desc *)m)->fd[i].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)m)->fd[i].ubound);
                          e_shut;
                          E_TPOPP("y_suba")
                          return((a_VOID) s);
			         }	

                        ((y_desc *)s)->fd[j].stride =
                           ((y_desc *)m)->fd[i].stride;
                        offset += (((y_desc *)s)->fd[i].lbound-
                                   ((y_desc *)m)->fd[j].lbound)*
                                  ((y_desc *)m)->fd[i].stride;
                        j++;
                        break;

              /* the lower bound is taken from the arg_list     */
              /* possible pointer allignment problems           */
              case '3': ((y_desc *)s)->fd[j].lbound = e_ref(a_intg);
                        ((y_desc *)s)->fd[j].ubound =
                           ((y_desc *)m)->fd[i].ubound;

				    if (((y_dscp)s)->fd[j].lbound < 
					   ((y_dscp)m)->fd[i].lbound)
                        {
					 a_intg k;
					 k=j+1;
                          e_trap(INDEX_RANGE,12,E_TMSG,67,
                                 E_TINT+E_TEXT(11),&k,
                                 E_TINT+E_TEXT(5),&((y_desc *)s)->fd[j].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)s)->fd[j].ubound,
                                 E_TINT+E_TEXT(5),&((y_desc *)m)->fd[i].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)m)->fd[i].ubound);
                          e_shut;
                          E_TPOPP("y_suba")
                          return((a_VOID) s);
			         }	

                        ((y_desc *)s)->fd[j].stride =
                           ((y_desc *)m)->fd[i].stride;
                        offset += (((y_desc *)s)->fd[i].lbound-
                                   ((y_desc *)m)->fd[j].lbound)*
                                  ((y_desc *)m)->fd[i].stride;
                        j++;
                        break;

              /* the upper bound is taken from the arg_list     */
              /* possible pointer allignment problems           */
              case '4': ((y_desc *)s)->fd[j].ubound = e_ref(a_intg);
                        ((y_desc *)s)->fd[j].lbound =
                           ((y_desc *)m)->fd[i].lbound;

				    if ( ((y_dscp)s)->fd[j].ubound > 
					   ((y_dscp)m)->fd[i].ubound )
                        {
					 a_intg k;
					 k=j+1;
                          e_trap(INDEX_RANGE,12,E_TMSG,67,
                                 E_TINT+E_TEXT(11),&k,
                                 E_TINT+E_TEXT(5),&((y_desc *)s)->fd[j].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)s)->fd[j].ubound,
                                 E_TINT+E_TEXT(5),&((y_desc *)m)->fd[i].lbound,
                                 E_TINT+E_TEXT(6),&((y_desc *)m)->fd[i].ubound);
                          e_shut;
                          E_TPOPP("y_suba")
                          return((a_VOID) s);
			         }	
                        ((y_desc *)s)->fd[j].stride =
                           ((y_desc *)m)->fd[i].stride;
                        j++;
                        break;

              default : ;
              }
           }

        e_shut;

        /* The higher indexranges which are not in the arg_list */
        /* are taken from the mainarray.                        */
        for (i=length;i<dim;i++)
           {
           ((y_desc *)s)->fd[j].lbound = ((y_desc *)m)->fd[i].lbound;
           ((y_desc *)s)->fd[j].ubound = ((y_desc *)m)->fd[i].ubound;
           ((y_desc *)s)->fd[j].stride = ((y_desc *)m)->fd[i].stride;
           j++;
           }

        ((y_desc *)s)->subarr = TRUE;
        ((y_desc *)s)->destroy = FALSE;
        ((y_desc *)s)->elsize = ((y_desc *)m)->elsize;
        ((y_desc *)s)->numdim = sdim;

	   /* assign sub(m) to s if is m allocated */
	   if (!y_aled(m) ) ((y_desc *)s)->array = NULL;
        else ((y_desc *)s)->array = (char *)((y_desc *)m)->array+
                                    ((y_desc *)m)->elsize*offset;

        ((y_desc *)s)->elnum = 1;
        for (j=0;j<sdim;j++)
           ((y_desc *)s)->elnum *= (size_t) (((y_desc *)s)->fd[j].ubound-
                                   ((y_desc *)s)->fd[j].lbound+1);

        E_TPOPP("y_suba")
        return((a_VOID) s);
        }





