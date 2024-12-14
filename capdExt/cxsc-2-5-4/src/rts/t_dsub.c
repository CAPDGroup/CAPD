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

/* CVS $Id: t_dsub.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_dsub.c                              */
/*                                                              */
/*                   dummy return 0                             */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

/*--------------------------------------------------------------*/
/* sub DReal arg to DReal res                                   */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int subdd(const DReal *arg1, const DReal *arg2, DReal *res)
#else
int subdd(arg1, arg2, res)
const DReal *arg1;
const DReal *arg2;
      DReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int subdd(DReal *arg1, DReal *arg2, DReal *res)
#else
int subdd(arg1, arg2, res)
DReal *arg1;
DReal *arg2;
DReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{

   /* --- ein arg == 0 --- */
   if(arg1->s==0) {
      copydd(arg2, res);
      res->s = -res->s;
      return 0;
   }

   if(arg2->s==0) {
      copydd(arg1, res);
      return 0;
   }

   /* --- ein arg < 0 --- */
   if(arg1->s!=arg2->s) {
      copydd(arg2, res);
      res->s = -res->s;
      return dadd(arg1, res, res);
   }

   /* --- beide gleiche sign --- */
   switch(cmpabsdd(arg1, arg2)) {
   case 1:
      return dsub(arg1, arg2, res);
   case 0:
      return initd(res);
   case -1:
      dsub(arg2, arg1, res);
      res->s *= (-1);
      return 0;
   }

   return 0;

/* return (1==cmpdd(arg1, arg2)*arg1->s ?
          dsub(arg1, arg2, res) :
          dsub(arg2, arg1, res), res->s = -res->s); */

   /*NOTREACHED*/
} /* subdd */





