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

/* CVS $Id: f_args.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_args.c                              */
/*                                                              */
/*      Entry           : void f_args(d)                        */
/*                        s_trng *d;                            */
/*                                                              */
/*      Arguments       : d - string for command line argument  */
/*                                                              */
/*      Description     : Copy command line argument string to  */
/*                        string variable.                      */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern int f_orgc;
extern char **f_orgv;
#endif

#ifdef LINT_ARGS
local void f_args(s_trng *d)
#else
local void f_args(d)

s_trng *d;
#endif
        {
        size_t arglen;
        static int index = 0;
 
        E_TPUSH("f_args")

        if (index>=f_orgc)
           {
           d->clen = 0;
      
           E_TPOPP("f_args")
           return; 
           }

        /* argument is empty string */
        if ((arglen = strlen(f_orgv[index]))==0)
           {
           d->clen = 0;
           }

        /* destination is smaller than source */
        else if (d->alen<arglen || d->suba)
           {

           /* destination is dynamic */
           if (d->fix==FALSE)
              {
              s_free(d);
              d->alen = d->clen = arglen;

              /* one more allocated for direct use of '\0'   */
              if ((d->ptr = (char*) malloc(arglen+1))==NULL)
                 e_trap(ALLOCATION,2,E_TMSG,54);
              else
                 {
#ifdef HEAP_CHECK
b_geth((a_char *)&d->ptr,(a_char *)d->ptr,(a_char *)"s_asgn");
#endif
                 strcpy(d->ptr,f_orgv[index]);
                 if (d->suba)
                    {
                    d->suba = FALSE;
                    d->tmp = TRUE;
                    }
                 }
              }

           /* destination is not dynamic == chopping */
           else
              {
              memcpy(d->ptr,f_orgv[index],d->alen);
              d->clen = d->alen;
              }
           }

        /* destination large enough to hold source */
        else
           {
           d->clen = arglen;
           strcpy(d->ptr,f_orgv[index]);
           }

        index++;

        E_TPOPP("f_args")
        return;
        }





