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

/* CVS $Id: t_dcmp.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_dcmp.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

#ifdef ANSI_C
#ifdef LINT_ARGS
static int mantcmp(const DMant *arg1, const DMant *arg2);
#else
static int mantcmp();
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int mantcmp(DMant *arg1, DMant *arg2);
#else
static int mantcmp();
#endif /* LINT_ARGS */
#endif /* ANSI_C */

/*--------------------------------------------------------------*/
/* cmpdd return 1 if arg1>arg2, 0 if ==, -1 if <                */
/*--------------------------------------------------------------*/
#ifdef LINT_ARGS
int cmpdd(DReal *arg1, DReal *arg2)
#else
int cmpdd(arg1, arg2)
DReal   *arg1;
DReal   *arg2;
#endif /* LINT_ARGS */
{
   Digit highdigit, h1, h2;
   DMant m;
   register DExp  e1;
   register DExp  e2;

   /* --- compare sign --- */
   if(arg1->s!=arg2->s)
      return (int)(arg1->s > arg2->s ? 1 : -1);

   if(arg1->s==0)
      return 0;

   /* --- compare exponent --- */
   highdigit = h1 = arg1->m.digit[DMantLen];
   for(e1=arg1->e; highdigit; highdigit >>=1, e1++);

   highdigit = h2 = arg2->m.digit[DMantLen];
   for(e2=arg2->e; highdigit; highdigit >>=1, e2++);

   if (e1!=e2)
      return (int)(e1>e2 ? arg1->s : -arg1->s);

   /* --- compare mant --- */
   if(h1) {
      dmshift(h1, &(arg1->m), &m);
      memcpy(&(arg1->m), &m, sizeof(DMant));
      arg1->e = e1;
   }
   if(h2) {
      dmshift(h2, &(arg2->m), &m);
      memcpy(&(arg2->m), &m, sizeof(DMant));
      arg2->e = e2;
   }
   return mantcmp(&(arg1->m), &(arg2->m));

} /* cmpdd() */

/*--------------------------------------------------------------*/
/* cmpabsdd return 1 if abs(arg1) > abs(arg2), 0 if ==, -1 if < */
/*--------------------------------------------------------------*/
#ifdef LINT_ARGS
int cmpabsdd(DReal *arg1, DReal *arg2)
#else
int cmpabsdd(arg1, arg2)
DReal   *arg1;
DReal   *arg2;
#endif /* LINT_ARGS */
{
   Digit highdigit, h1, h2;
   DMant m;
   register DExp  e1;
   register DExp  e2;

   /* --- compare sign --- */
   if(arg1->s==0 && arg2->s==0)
      return 0;
   if(arg1->s==0)
      return -1;
   if(arg2->s==0)
      return 1;

   /* --- compare exponent --- */
   highdigit = h1 = arg1->m.digit[DMantLen];
   for(e1=arg1->e; highdigit; highdigit >>=1, e1++);

   highdigit = h2 = arg2->m.digit[DMantLen];
   for(e2=arg2->e; highdigit; highdigit >>=1, e2++);

   if (e1!=e2)
      return (int)(e1>e2 ? 1 : -1);

   /* --- compare mant --- */
   if(h1) {
      dmshift(h1, &(arg1->m), &m);
      memcpy(&(arg1->m), &m, sizeof(DMant));
      arg1->e = e1;
   }
   if(h2) {
      dmshift(h2, &(arg2->m), &m);
      memcpy(&(arg2->m), &m, sizeof(DMant));
      arg2->e = e2;
   }
   return mantcmp(&(arg1->m), &(arg2->m));

} /* cmpabsdd() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int mantcmp(const DMant *arg1, const DMant *arg2)
#else
static int mantcmp(arg1, arg2)
const DMant *arg1;
const DMant *arg2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int mantcmp(DMant *arg1, DMant *arg2)
#else
static int mantcmp(arg1, arg2)
DMant *arg1;
DMant *arg2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   register int i;

   for(i=(int)sizeof(DMant);
       --i>=0 && arg1->digit[i]==arg2->digit[i];);
   if(i<0) return 0;
   return (arg1->digit[i]>arg2->digit[i]? 1 : -1);
} /* mantcmp() */





