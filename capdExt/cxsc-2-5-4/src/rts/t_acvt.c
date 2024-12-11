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

/* CVS $Id: t_acvt.c,v 1.21 2014/01/30 17:24:14 cxsc Exp $ */


#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | acosviatan                                                   |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int acosviatan(const ExtReal *arg, ExtReal *res)
#else
int acosviatan(arg, res)
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int acosviatan(ExtReal *arg, ExtReal *res)
#else
int acosviatan(arg, res)
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal rp1;
   ExtReal rm1;
   ExtReal r;
   ExtReal rsqrt;
   ExtReal rdiv;
   int   ret;

   ret =  subee(&One, arg, &rm1);
   ret += addee(&One, arg, &rp1);
   ret += mulee(&rp1, &rm1, &r);
   ret += sqrtee(&r, &rsqrt);

   if(ret += divee(&rsqrt, arg, &rdiv)) return ret;

   ret = _s_atan(&rdiv, res);

   return ret;
} /* acosviatan() */





