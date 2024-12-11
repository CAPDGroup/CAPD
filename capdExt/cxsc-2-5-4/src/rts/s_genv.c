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

/* CVS $Id: s_genv.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_genv.c                              */
/*                                                              */
/*      Entries         : void s_genv(name,value,exists)        */
/*                        s_trng name;                          */
/*                        s_trng *value;                        */
/*                        a_bool *exists;                       */
/*                                                              */
/*      Arguments       : name = name of environment variable   */
/*                        value = value of environment variable */
/*                        exists = TRUE if variable exists      */
/*                                                              */
/*      Description     : Determine value of environment        */
/*                        variable using call to getenv().      */
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
local void s_genv(s_trng name,s_trng *value,a_bool *exists)
#else
local void s_genv(name,value,exists)

s_trng name;
s_trng *value;
a_bool *exists;
#endif
        {
        char *p,ch;
        size_t len;

        E_TPUSH("s_genv")

        if (name.clen==0)
           {
           e_trap(I_O_ERROR,2,E_TMSG,61);
           *exists = FALSE;
           value->clen = 0;
           }
        else
           {
           if (name.suba) s_asgn(&name,name);
           ch = name.ptr[name.clen];
           name.ptr[name.clen] = '\0';
           if ((p = getenv(name.ptr))==NULL)
              {
              *exists = FALSE;
              value->clen = 0;
              }
           else
              {
              *exists = TRUE;
              if ((len = strlen(p))>value->alen)
                 {
                 if (value->fix)
                    len = value->alen;
                 else
                    {
                    if (value->alen>0)
                       {
#ifdef HEAP_CHECK
b_freh((a_char *)&value->ptr,(a_char *)value->ptr,(a_char *)"s_genv");
#endif
                       free(value->ptr);
                       }
                    if ((value->ptr = (char*) malloc(len+1))==NULL)
                       {
                       e_trap(ALLOCATION,2,E_TMSG,54);
                       len = 0;
                       }
                    else
                       {
#ifdef HEAP_CHECK
b_geth((a_char *)&value->ptr,(a_char *)value->ptr,(a_char *)"s_genv");
#endif
                       value->alen = len;
                       }
                    }
                 }
              if (len>0) memcpy(value->ptr,p,len);
              value->clen = len;
              }
           name.ptr[name.clen] = ch;
           }

        if (name.tmp) s_free(&name);

        E_TPOPP("s_genv")
        return;
        }





