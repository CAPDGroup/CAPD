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

/* CVS $Id: t_dnrm.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_dnrm.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

#ifdef LINT_ARGS
static int ret_zero(int sgn, DReal *res);
#else
static int ret_zero();
#endif
/*--------------------------------------------------------------*/
/* normalize DReal                                              */
/*--------------------------------------------------------------*/

#ifdef ANSI_C
#ifdef LINT_ARGS
int normd(const DReal *arg, DReal *res)
#else
int normd(arg, res)
const DReal *arg;
      DReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int normd(DReal *arg, DReal *res)
#else
int normd(arg, res)
DReal *arg;
DReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DReal d;
   DExp  e;
   Digit highdigit;

   if(arg->s == 0) return ret_zero(0, res);
   copydd(arg, &d);                     /* falls *res == *arg */

   if(d.m.digit[DMantLen]==0) {
      if(1==dmadjust(&d.m, DMantLen, &(res->m), &e))
         return ret_zero(d.s, res);
   }
   else {
      highdigit = d.m.digit[DMantLen];
      for(e=0; highdigit; highdigit >>=1, e++);
      dmshift(e, &d.m, &(res->m));
   }
   res->s = d.s;
   res->e = d.e + e;

   return NoErr;
} /* normd() */

/*--------------------------------------------------------------*/
#ifdef LINT_ARGS
static int ret_zero(int sgn, DReal *res)
#else
static int ret_zero(sgn, res)
int   sgn;
DReal *res;
#endif
{
   initd(res);
   if(sgn!=NEG)
      return ExcPZero;
   else
      return ExcMZero;
} /* ret_zero */





