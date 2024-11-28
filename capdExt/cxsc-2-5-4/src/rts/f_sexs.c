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

/* CVS $Id: f_sexs.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_sexs.c                              */
/*                                                              */
/*      Entry           : a_bool f_sexs(name)                   */
/*                        s_trng name;                          */
/*                                                              */
/*      Arguments       : name   - filename                     */
/*                                                              */
/*      Function value  : TRUE if file exists                   */
/*                                                              */
/*      Description     : Test for existing file                */
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
local a_bool f_sexs(s_trng name)
#else
local a_bool f_sexs(name)

s_trng name;
#endif

        {
        FILE *fp;
        a_bool res;
        char ch;

        E_TPUSH("f_sexs")

        res = FALSE;

        if (name.clen==0)
           {
           res = TRUE;
           }
        else if (name.ptr==NULL || name.ptr[0]=='\0')
           {
           e_trap(I_O_ERROR,2,E_TMSG,45);
           }
        else if (name.clen>=BUFFERSIZE-1)
           {
           e_trap(I_O_BUFFER,4,E_TMSG,29,E_TMSG,30);
           }
        else
           {
           if (name.suba) s_asgn(&name,name);
           ch = name.ptr[name.clen];
           name.ptr[name.clen] = '\0';
           if ((fp=fopen(name.ptr,"r"))!=NULL)
              {
              fclose(fp);
              res = TRUE;
              }
           name.ptr[name.clen] = ch;
           }

        if (name.tmp) s_free(&name);

        E_TPOPP("f_sexs")
        return(res);
        }





