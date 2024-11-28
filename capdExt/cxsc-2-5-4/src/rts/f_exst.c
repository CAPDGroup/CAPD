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

/* CVS $Id: f_exst.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_exst.c                              */
/*                                                              */
/*      Entry           : a_bool f_exst(name)                   */
/*                        a_char name[];                        */
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
local a_bool f_exst(a_char name[])
#else
local a_bool f_exst(name)

a_char name[];
#endif

        {
        FILE *fp;

        E_TPUSH("f_exst")

        if (name==NULL)
           {
           e_trap(I_O_ERROR,2,E_TMSG,62);
           E_TPOPP("f_exst")
           return(FALSE);
           }

        if (name[0]=='\0')
           {
           E_TPOPP("f_exst")
           return(TRUE);
           }

        if ((fp=fopen((char *)name,"r"))==NULL)
           {
           E_TPOPP("f_exst")
           return(FALSE);
           }

        fclose(fp);

        E_TPOPP("f_exst")
        return(TRUE);
        }





