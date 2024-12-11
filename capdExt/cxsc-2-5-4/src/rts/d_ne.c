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

/* CVS $Id: d_ne.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_ne.c                                */
/*                                                              */
/*      Entries         : a_bool d_ne(a,b)                      */
/*                        dotprecision a,b;                     */
/*                                                              */
/*      Arguments       : a = dotprecision value                */
/*                        b = dotprecision value                */
/*                                                              */
/*      Function value  : TRUE if a!=b                          */
/*                        FALSE if a==b                         */
/*                                                              */
/*      Description     : Compare dotprecision values           */
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
local a_bool d_ne(dotprecision a,dotprecision b)
#else
local a_bool d_ne(a,b)

dotprecision a;
dotprecision b;
#endif
        {
        int i;  /* !!! must be int !!! */
        a_bool res = FALSE;

        if (a[A_BEGIN]==ZERO)
           {
           if (b[A_BEGIN]==ZERO) res = FALSE;
           else res = TRUE;
           }
        else if (b[A_BEGIN]==ZERO) res = TRUE;
        else
           {
           if (a[A_SIGN]!=b[A_SIGN]) res = TRUE;
           else if (a[A_BEGIN]!=b[A_BEGIN]) res = TRUE;
           else
              {
              for (i = (int)a[A_BEGIN];
                   i<=a[A_END] && i<=b[A_END];i++)
                 if (a[i]!=b[i]) { res = TRUE; break; }

              if (i>a[A_END])
                 {
                 if (a[A_END]!=b[A_END]) res = TRUE;
                 else res = FALSE;
                 }
              }
           }

        if (a[A_STATUS] & A_TEMPORARY) d_free(&a);
        if (b[A_STATUS] & A_TEMPORARY) d_free(&b);

        return(res);
        }





