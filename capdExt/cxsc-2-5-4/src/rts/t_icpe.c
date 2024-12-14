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

/* CVS $Id: t_icpe.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | icopyee                                                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int icopyee (const IExtReal *arg, IExtReal *res)
#else
int icopyee (arg, res)
const IExtReal   *arg;
      IExtReal   *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int icopyee (IExtReal *arg, IExtReal *res)
#else
int icopyee (arg, res)
IExtReal   *arg;
IExtReal   *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int   ret;

   ret = copyee(&(arg->u), &(res->u));
   ret = copyee(&(arg->l), &(res->l));

   return ret;
} /* icopyee() */





