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

/* CVS $Id: f_rds1.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rds1.c                              */
/*                                                              */
/*      Entry           : void f_rds1(desc,s)                   */
/*                        f_text *desc;                         */
/*                        s_trng *s;                            */
/*                                                              */
/*      Arguments       : desc   - device descriptor            */
/*                        s      - string                       */
/*                                                              */
/*      Description     : perform PASCAL read(string).          */
/*                                                              */
/*                   one more allocated for direct use of '\0'  */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern dotprecision b_cm__;
#endif

#ifdef LINT_ARGS
local void f_rds1(f_text *desc,s_trng *s)
#else
local void f_rds1(desc,s)

f_text *desc;
s_trng *s;
#endif
        {
        size_t alen,i,k;
        char *buffer,*ptr,*dummy;

        E_TPUSH("f_rds1")

        if (b_text(desc,TRUE))
           {
           buffer = (char *)b_cm__;
           alen = (size_t)((s->fix==TRUE) ? s->alen : MAXINT);
           k = i = 0;
           ptr = NULL;

           while ((desc->eoln==FALSE) && (i<alen))
              {
              buffer[k] = desc->win.ch[0];
              i++;
              if (++k==BUFFERSIZE)
                 {
                 k = 0;
                 if ((dummy = (char*) malloc(i))==NULL)
                    {
                    e_trap(I_O_BUFFER,2,E_TMSG,55);
                    E_TPOPP("f_rds1")
                    return;
                    }
                 else
                    {
                    if (ptr)
                       {
                       memcpy(dummy,ptr,i-BUFFERSIZE);
#ifdef HEAP_CHECK
b_freh((a_char *)&ptr,(a_char *)ptr,(a_char *)"f_rds1");
#endif
                       free(ptr);
                       }
#ifdef HEAP_CHECK
b_geth((a_char *)&ptr,(a_char *)dummy,(a_char *)"f_rds1");
#endif
                    memcpy(dummy+(i-BUFFERSIZE),buffer,BUFFERSIZE);
                    ptr = dummy;
                    }
                 }
              f_getc(desc);
              }

           if (s->alen<i)
              {
              if (s->alen>0)
                 {
#ifdef HEAP_CHECK
b_freh((a_char *)&s->ptr,(a_char *)s->ptr,(a_char *)"f_rds1");
#endif
                 free((char *)s->ptr);
                 }

              /* one more for direct use of '\0'*/
              if ((s->ptr = (char*) malloc(i+1))==NULL)
                 {
                 e_trap(ALLOCATION,2,E_TMSG,54);
                 s->clen = s->alen = 0;
                 E_TPOPP("f_rds1")
                 return;
                 }
              else
                 {
#ifdef HEAP_CHECK
b_geth((a_char *)&s->ptr,(a_char *)s->ptr,(a_char *)"f_rds1");
#endif
                 s->alen = i;
                 }
              }

           s->clen = i;
           if (ptr)
              {
              memcpy(s->ptr,ptr,i-k);
#ifdef HEAP_CHECK
b_freh((a_char *)&ptr,(a_char *)ptr,(a_char *)"f_rds1");
#endif
              free(ptr);
              }
           if (k>0) memcpy(s->ptr+(i-k),buffer,k);
           }

        E_TPOPP("f_rds1")
        }





