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

/* CVS $Id: a_syst.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_syst.c                              */
/*                                                              */
/*      Entries         : a_intg a_syst(s)                      */
/*                        s_trng s;                             */
/*                                                              */
/*      Arguments       : command string                        */
/*                                                              */
/*      Description     : execute command in command string     */
/*                                                              */
/*                   correct return call in case of IBM_AT_MS_C */
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

/*------------------------------------------------------*/
/*  msc_system : system call that gives the correct     */
/*               return value of the given command.     */
/*------------------------------------------------------*/

#if IBM_AT_MS_C

#include <process.h>

static int msc_system (const char *s)

        {
        char *par[10];
        int   i = 0;
        char *p = (char *)b_cm__;

        (void)strcpy(p,s);

        while (*p==' ') p++;                    /* skip leading blanks */
        while (*p!='\0')
           {
           par[i++] = p;
           while (*p!=' ' && *p!='\0') p++;     /* search next blank */
           if (*p==' ')
              {
              *p = '\0';                        /* store null char */
              p++;
              while (*p==' ') p++;              /* skip leading blanks */
              }
           }
        par [i] = NULL ;

        return  spawnvp (P_WAIT, par[0],par ) ;
        }
#endif



#ifdef LINT_ARGS
local a_intg a_syst(s_trng s)
#else
local a_intg a_syst(s)

s_trng s;
#endif
        {
        a_intg n;

        E_TPUSH("a_syst")

        if (s.clen>0)
           {
           if (s.suba) s_asgn(&s,s);
           s.ptr[s.clen] = '\0';
#if IBM_AT_MS_C
           n = msc_system(s.ptr);
#else
           n = system(s.ptr);
#endif
           }
        else
#if IBM_AT_MS_C
           n = msc_system(NULL);
#else
           n = system(NULL);
#endif

        if (s.tmp) s_free(&s);

        E_TPOPP("a_syst")
        return(n);
        }







