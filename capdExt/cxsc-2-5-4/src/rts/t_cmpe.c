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

/* CVS $Id: t_cmpe.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */


#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#ifdef LINT_ARGS
#ifdef ANSI_C
static int mantcmp(const Digit arg1[], const Digit arg2[]);
#else
static int mantcmp(Digit arg1[], Digit arg2[]);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int mantcmp();
#endif /* LINT_ARGS */
/*--------------------------------------------------------------*
 | cmpee                                                        |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int cmpee (const ExtReal *arg1, const ExtReal *arg2)
#else
int cmpee (arg1, arg2)
const ExtReal  *arg1;
const ExtReal  *arg2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int cmpee (ExtReal *arg1, ExtReal *arg2)
#else
int cmpee (arg1, arg2)
ExtReal  *arg1;
ExtReal  *arg2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   unsigned int   s1, s2;
   int   e1, e2;
   int   m;

   /* --- compare sign --- */
   s1 = arg1->s.exp&ExtSignMask;
   s2 = arg2->s.exp&ExtSignMask;
   if(s1!=s2)
      if(arg1->s.DIGIT(7)!=0&&arg2->s.DIGIT(7)!=0)
         return (int)(s1 > s2 ? -1 : 1);

   /* --- compare exponent --- */
   e1 = (arg1->s.exp&ExtExpMask);
   e2 = (arg2->s.exp&ExtExpMask);
   if (e1!=e2)
      return (int)(e1>e2 ? (int)(s1==0?1:-1) : (int)(s2==0?-1:1));

   /* --- compare mant --- */           /* 130391 cb */
/* m = mantcmp(&(arg1->s.digit[0]), &(arg2->s.digit[0])); */
   m = mantcmp(arg1->s.digit, arg2->s.digit);
   return (int)(s1==0?m:-m);
} /* cmpee() */

/*--------------------------------------------------------------*
 | cmpabsee                                                     |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int cmpabsee (const ExtReal *arg1, const ExtReal *arg2)
#else
int cmpabsee (arg1, arg2)
const ExtReal  *arg1;
const ExtReal  *arg2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int cmpabsee (ExtReal *arg1, ExtReal *arg2)
#else
int cmpabsee (arg1, arg2)
ExtReal  *arg1;
ExtReal  *arg2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   ExtReal a1, a2;

   /* --- abs --- */
   absee(arg1, &a1);
   absee(arg2, &a2);

   /* --- cmpabs --- */
   return cmpee(&a1, &a2);
} /* cmpabsee() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int mantcmp(const Digit arg1[], const Digit arg2[])
#else
static int mantcmp(arg1, arg2)
const Digit arg1[];
const Digit arg2[];
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int mantcmp(Digit arg1[], Digit arg2[])
#else
static int mantcmp(arg1, arg2)
Digit arg1[];
Digit arg2[];
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   register int i;

#if LSBFIRST
   for (i=ExtMantLen; --i>=0 && arg1[i]==arg2[i]; ) ;
   if (i<0) return(0);
   return (arg1[i]>arg2[i]? 1 : -1);
#else
   for (i=0; i<ExtMantLen && arg1[i]==arg2[i]; i++) ;
   if (i==ExtMantLen) return(0);
   return (arg1[i]>arg2[i] ? 1 : -1);
#endif
} /* mantcmp() */





