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

/* CVS $Id: y_asgn.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : y_asgn.c                              */
/*                                                              */
/*      Entry           : void y_asgn(d,s)                      */
/*                        y_dscp d;                             */
/*                        y_dscp s;                             */
/*                                                              */
/*      Arguments       : d - descriptor of destination array   */
/*                        s - descriptor of source array        */
/*                                                              */
/*      Description     : Assign contents of a dynamic array to */
/*                        another dynamic array.                */
/*                           f_ppib==FALSE: check index bounds  */
/*                           f_ppib==TRUE : copy index bounds   */
/*                                                              */
/*      Note            : elsize and numdim of s and d must     */
/*                        coincide.                             */
/*                                                              */
/*                   beautyfied text for full dynamic           */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#ifdef COPYIB_ENABLE
extern a_bool f_ppib;
#endif
#endif

#ifdef LINT_ARGS
local void y_asgn(y_dscp d,y_dscp s)
#else
local void y_asgn(d,s)

y_dscp d;
y_dscp s;
#endif
 {
 a_intg i,k;
 a_intg *help;
 a_btyp od = 0;
 a_btyp os = 0;
 a_btyp prod;
 char *pos,*pod;
 a_bool reverse = FALSE;

 E_TPUSH("y_asgn")

 /* checking right side for allocation,
  * error if not
  */
  if (!y_aled(s)) 
  {
       e_trap(ALLOCATION,2,E_TMSG,42);
       E_TPOPP("y_asgn")
       return;
  }

 /* checking left side for allocation , 
  * allocate , 
  * initialize if not 
  */
  if (!y_aled(d)) 
  {
     a_bool suba,destr;
	suba  = d->subarr;
	destr = d->destroy;

	/* left induct bounds to right side */
     for (i=0;i<((y_desc *)s)->numdim;i++)
     {
          ((y_desc *)d)->fd[i].lbound = ((y_desc *)s)->fd[i].lbound;
          ((y_desc *)d)->fd[i].ubound = ((y_desc *)s)->fd[i].ubound;
     }
	y_init(d, s->numdim, s->elsize); /* new allocate and save flags */
	d->destroy = destr;
	d->subarr  = suba ;
  }


 /* Copy matching positions only if destination is a     */
 /* subarray. Other positions are unaltered.             */
 if (((y_desc *)d)->subarr)
    {

    /* Allocate index array help.                        */
    if ((help =
(a_intg *)malloc(((y_desc *)s)->numdim*sizeof(a_intg)))==NULL)
       {
       e_trap(ALLOCATION,2,E_TMSG,42);
       E_TPOPP("y_asgn")
       return;
       }
#ifdef HEAP_CHECK
b_geth((a_char *)&help,(a_char *)help,(a_char *)"y_asgn");
#endif

    for (i=((y_desc *)s)->numdim-1;i>=0;i--)
       {

       /* Return if widths of index ranges are different */
       if (((y_desc *)s)->fd[i].ubound- ((y_desc *)s)->fd[i].lbound!=
           ((y_desc *)d)->fd[i].ubound- ((y_desc *)d)->fd[i].lbound)
          {
          k = i+1;
          e_trap(INDEX_RANGE,12,E_TMSG,67,
                 E_TINT+E_TEXT(11),&k,
                 E_TINT+E_TEXT(5),&((y_desc *)d)->fd[i].lbound,
                 E_TINT+E_TEXT(6),&((y_desc *)d)->fd[i].ubound,
                 E_TINT+E_TEXT(5),&((y_desc *)s)->fd[i].lbound,
                 E_TINT+E_TEXT(6),&((y_desc *)s)->fd[i].ubound);
#ifdef HEAP_CHECK
b_freh((a_char *)&help,(a_char *)help,(a_char *)"y_asgn");
#endif
          B_FREE(help)
          E_TPOPP("y_asgn")
          return;
          }

       /* help[i] holds actual i-th index.               */
       help[i] = ((y_desc *)s)->fd[i].lbound;
       }

    /* copy subarray to subarray                         */
    /* check for overlapping                             */
    if (((y_desc *)s)->subarr)
       {
       prod = ((y_desc *)s)->elsize*
          (((y_desc *)s)->fd[((y_desc *)s)->numdim-1].ubound-
           ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].lbound);

       /* source stride smaller than destination stride  */
       if (((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride<
           ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride)
          {
          pos = (char *)((y_desc *)s)->array;
          pod = (char *)((y_desc *)d)->array;
          while (pos>pod)
             {
             pos += ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride*
                    ((y_desc *)s)->elsize;
             pod += ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride*
                    ((y_desc *)s)->elsize;
             }

          if (pos<pod &&
              pod<=(char *)((y_desc *)s)->array+
                   ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride*prod)
             {
             do {
                pos += ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride*
			        ((y_desc *)s)->elsize;
                }
             while(pos<pod);
             if (pod==pos)
                {
                reverse = TRUE;
/* subarrays overlap */
                }
             }
          }

       /* source stride larger than destination stride   */
       else if (((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride>
                ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride)
          {
          pos = (char *)((y_desc *)s)->array+
((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride*prod;
          pod = (char *)((y_desc *)d)->array+
                ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride*prod;
          while (pos>pod)
             {
             pos -=
((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride*((y_desc *)s)->elsize;
             pod -=
((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride*((y_desc *)s)->elsize;
             }
          if (pos<pod && (char *)((y_desc *)d)->array<=pos)
             {
             do {
                pod -=
((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride*((y_desc *)d)->elsize;
                }
             while(pos<pod);
             if (pos==pod)
                {
                reverse = TRUE;
/* subarrays overlap */
                }
             }
          }/* if else stride */
       }/* if s subarray */

    if (reverse)
       {
       do
          {
          os += (((y_desc *)s)->fd[((y_desc *)s)->numdim-1].ubound-
                 ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].lbound)*
                 ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride;
          od += (((y_desc *)d)->fd[((y_desc *)s)->numdim-1].ubound-
                 ((y_desc *)d)->fd[((y_desc *)s)->numdim-1].lbound)*
                 ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride;

          /* copy element from source to destination     */
          for (i=((y_desc *)s)->fd[((y_desc *)s)->numdim-1].lbound;
               i<=((y_desc *)s)->fd[((y_desc *)s)->numdim-1].ubound;i++)
             {
             C_COPY((char *)((y_desc *)d)->array+
od*((y_desc *)s)->elsize,
                    (char *)((y_desc *)s)->array+
os*((y_desc *)s)->elsize,((y_desc *)s)->elsize)
             os -= ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride;
             od -= ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride;
             }

          if (((y_desc *)s)->numdim>1)
             {
             os += ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride+
                   ((y_desc *)s)->fd[((y_desc *)s)->numdim-2].stride;
             od += ((y_desc *)d)->fd[((y_desc *)s)->numdim-1].stride+
                   ((y_desc *)d)->fd[((y_desc *)s)->numdim-2].stride;
             help[((y_desc *)s)->numdim-2]++;

             /* advance to next index of source element     */
             for (i=(int)((y_desc *)s)->numdim-2;i>0;i--)
             if (help[i]>((y_desc *)s)->fd[i].ubound)
                {
                os -= ((y_desc *)s)->fd[i].stride*
				  (((y_desc *)s)->fd[i].ubound-
                       ((y_desc *)s)->fd[i].lbound+1)
                     -((y_desc *)s)->fd[i-1].stride;
                od -= ((y_desc *)d)->fd[i].stride*
                      (((y_desc *)d)->fd[i].ubound-
                       ((y_desc *)d)->fd[i].lbound+1)
                     -((y_desc *)d)->fd[i-1].stride;
                help[i] = ((y_desc *)s)->fd[i].lbound;
                help[i-1]++;
                }
             else break;
             }
          else
             break;
          }
       while (help[0]<=((y_desc *)s)->fd[0].ubound);
       }
    else
       {
       do
          {

          /* copy element from source to destination     */
          for (i=((y_desc *)s)->fd[((y_desc *)s)->numdim-1].lbound;
               i<=((y_desc *)s)->fd[((y_desc *)s)->numdim-1].ubound;i++)
             {
             C_COPY((char *)((y_desc *)d)->array+ od*((y_desc *)s)->elsize,
                    (char *)((y_desc *)s)->array+ os*((y_desc *)s)->elsize,
				((y_desc *)s)->elsize)
             os += ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride;
             od += ((y_desc *)d)->fd[((y_desc *)d)->numdim-1].stride;
             }

          if (((y_desc *)s)->numdim>1)
             {
             os -= ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].stride*
                   (((y_desc *)s)->fd[((y_desc *)s)->numdim-1].ubound-
                    ((y_desc *)s)->fd[((y_desc *)s)->numdim-1].lbound+1)
                  -((y_desc *)s)->fd[((y_desc *)s)->numdim-2].stride;
             od -= ((y_desc *)d)->fd[((y_desc *)s)->numdim-1].stride*
                   (((y_desc *)d)->fd[((y_desc *)s)->numdim-1].ubound-
                    ((y_desc *)d)->fd[((y_desc *)s)->numdim-1].lbound+1)
                  -((y_desc *)d)->fd[((y_desc *)s)->numdim-2].stride;
             help[((y_desc *)s)->numdim-2]++;

             /* advance to next index of source element     */
             for (i=(int)((y_desc *)s)->numdim-2;i>0;i--)
             if (help[i]>((y_desc *)s)->fd[i].ubound)
                {
                os -= ((y_desc *)s)->fd[i].stride* 
				  (((y_desc *)s)->fd[i].ubound-
                       ((y_desc *)s)->fd[i].lbound+1)-
				  ((y_desc *)s)->fd[i-1].stride;
                od -= ((y_desc *)d)->fd[i].stride*
                      (((y_desc *)d)->fd[i].ubound-
                       ((y_desc *)d)->fd[i].lbound+1)-
                      ((y_desc *)d)->fd[i-1].stride;
                help[i] = ((y_desc *)s)->fd[i].lbound;
                help[i-1]++;
                }
             else break;
             }
          else
             break;
          }
       while (help[0]<=((y_desc *)s)->fd[0].ubound);
       }

#ifdef HEAP_CHECK
b_freh((a_char *)&help,(a_char *)help,(a_char *)"y_asgn");
#endif
    B_FREE(help)
    }
 else
    {

    /* check index ranges if f_ppib==FALSE       */
#ifdef COPYIB_ENABLE
    if (f_ppib==FALSE)
#endif
       {
       for (i=((y_desc *)s)->numdim-1;i>=0;i--)
          {
          if (((y_desc *)d)->fd[i].ubound- ((y_desc *)d)->fd[i].lbound !=
              ((y_desc *)s)->fd[i].ubound- ((y_desc *)s)->fd[i].lbound)
             {
             k = i+1;
             e_trap(INDEX_RANGE,12,E_TMSG,67,
                    E_TINT+E_TEXT(11),&k,
                    E_TINT+E_TEXT(5),&((y_desc *)d)->fd[i].lbound,
                    E_TINT+E_TEXT(6),&((y_desc *)d)->fd[i].ubound,
                    E_TINT+E_TEXT(5),&((y_desc *)s)->fd[i].lbound,
                    E_TINT+E_TEXT(6),&((y_desc *)s)->fd[i].ubound);
             E_TPOPP("y_asgn")
             return;
             }
          }
       }

#ifdef COPYIB_ENABLE
    /* copy index bounds                         */
    else /* if (f_ppib) */
       {
       for (i=0;i<((y_desc *)s)->numdim;i++)
          {
          ((y_desc *)d)->fd[i].lbound = ((y_desc *)s)->fd[i].lbound;
          ((y_desc *)d)->fd[i].ubound = ((y_desc *)s)->fd[i].ubound;
          }
       }
#endif

    if (((y_desc *)s)->destroy==TRUE)
       {
       y_free(d);
#ifdef HEAP_CHECK
b_freh((a_char *)NULL,(a_char *)((y_desc *)s)->array,(a_char *)"y_asgn");
b_geth((a_char *)&((y_desc *)d)->array,
(a_char *)((y_desc *)s)->array,(a_char *)"y_asgn");
#endif
       ((y_desc *)d)->array = ((y_desc *)s)->array;
       ((y_desc *)s)->destroy = FALSE;
       ((y_desc *)s)->array = NULL;
       }
    else
       {
       if ((prod = ((y_desc *)s)->elnum*((y_desc *)s)->elsize)!=
                   ((y_desc *)d)->elnum*((y_desc *)d)->elsize)
          {
          y_free(d);
/*
          if (((a_char *)((y_desc *)d)->array = 
               (a_char *)malloc((size_t)prod))==NULL)
*/
          if ((((y_desc *)d)->array = (a_VOID) malloc((size_t)prod))==NULL)
             {
             e_trap(ALLOCATION,2,E_TMSG,42);
             E_TPOPP("y_asgn")
             return;
             }
#ifdef HEAP_CHECK
b_geth((a_char *)&((y_desc *)d)->array,
(a_char *)((y_desc *)d)->array,(a_char *)"y_asgn");
#endif
          }

       /* copy individual components if source is subarray  */
       if (((y_desc *)s)->subarr)
          {
          /* Allocate index array help.                        */
          if ((help =
(a_intg *)malloc(((y_desc *)s)->numdim*sizeof(a_intg)))==NULL)
             {
             e_trap(ALLOCATION,2,E_TMSG,42);
             E_TPOPP("y_asgn")
             return;
             }
#ifdef HEAP_CHECK
b_geth((a_char *)&help,(a_char *)help,(a_char *)"y_asgn");
#endif
          for (i=0;i<((y_desc *)s)->numdim;i++)
             help[i] = ((y_desc *)s)->fd[i].lbound;

          for (k=0;k<((y_desc *)s)->elnum;k++)
             {
             prod = 0;
             for (i=0;i<((y_desc *)s)->numdim;i++)
                prod += (help[i]-((y_desc *)s)->fd[i].lbound)*
                        ((y_desc *)s)->fd[i].stride;

             C_COPY((char *)((y_desc *)d)->array+ k*((y_desc *)s)->elsize,
                    (char *)((y_desc *)s)->array+
				prod*((y_desc *)s)->elsize,((y_desc *)s)->elsize)

             for (i=((y_desc *)s)->numdim-1;i>=0;i--)
                if (++help[i]>((y_desc *)s)->fd[i].ubound)
                   help[i] = ((y_desc *)s)->fd[i].lbound;
                else break;
             }
#ifdef HEAP_CHECK
b_freh((a_char *)&help,(a_char *)help,(a_char *)"y_asgn");
#endif
          B_FREE(help)
          }
       else
          C_COPY((char *)((y_desc *)d)->array,
			  (char *)((y_desc *)s)->array,prod)
       }

    /* determine strides for destination         */
    i = ((y_desc *)s)->numdim-1;
    ((y_desc *)d)->fd[i--].stride = 1;
    for (;i>=0;i--)
       {
       ((y_desc *)d)->fd[i].stride = (size_t) ((((y_desc *)d)->fd[i+1].ubound-
                                               ((y_desc *)d)->fd[i+1].lbound+1)*
                                               ((y_desc *)d)->fd[i+1].stride);
       }

    ((y_desc *)d)->elnum = ((y_desc *)s)->elnum;
    }

 if (((y_desc *)s)->destroy) y_free(s);

 E_TPOPP("y_asgn")
 return;
 }





