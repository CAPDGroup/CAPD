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

/* CVS $Id: t_das_.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

/*--------------------------------------------------------------*/
/* add two digits, return res and carry                         */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int addigit(Digit *d1, const Digit d2)
#else
static int addigit(d1, d2)
      Digit   *d1;
const Digit    d2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int addigit(Digit *d1, Digit d2)
#else
static int addigit(d1, d2)
Digit   *d1;
Digit    d2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int carry;

   if ((Digit)0xff - *d1 < d2) carry = 1;
   else carry = 0;

   *d1 += d2;

   return carry;
} /* addigit() */

/*--------------------------------------------------------------*/
/* add two DMants, Voraussetzung: highdigit.highbit = 0         */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int dmadd(DMant *arg1, const DMant *arg2)
#else
static int dmadd(arg1, arg2)
      DMant  *arg1;
const DMant  *arg2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int dmadd(DMant *arg1, DMant *arg2)
#else
static int dmadd(arg1, arg2)
DMant  *arg1;
DMant  *arg2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   register int i,j;

   for (i=0; i<sizeof(DMant); i++)
      if (addigit(&(arg1->digit[i]), arg2->digit[i])) {
         j = i;
         while (addigit(&(arg1->digit[++j]), 1));
      }

   return 0;
} /* dmadd() */

/*--------------------------------------------------------------*/
/* add two DReals, Voraussetzung: highdigit.highbit = 0         */
/* beide arg gleiches Sign                                      */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dadd(const DReal *arg1, const DReal *arg2, DReal *res)
#else
int dadd(arg1, arg2, res)
const DReal  *arg1;
const DReal  *arg2;
      DReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dadd(DReal *arg1, DReal *arg2, DReal *res)
#else
int dadd(arg1, arg2, res)
DReal  *arg1;
DReal  *arg2;
DReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DMant m;
   DExp e;

   e = arg1->e - arg2->e;
   if(e>0) {
      dmshift(e, &(arg2->m), &m);
      dmadd(&m, &(arg1->m));
      res->e = arg1->e;
   }
   else
   if(e<0) {
      dmshift(-e, &(arg1->m), &m);
      dmadd(&m, &(arg2->m));
      res->e = arg2->e;
   }
   else {
      memcpy(&m, &(arg1->m), sizeof(DMant));
      dmadd(&m, &(arg2->m));
      res->e = arg1->e;
   }

   /* --- */
   if(m.digit[DMantLen]&0x80) {
      dmshift(8, &m, &(res->m));
      e = 8;
   }
   else {
      memcpy(&(res->m), &m, sizeof(DMant));
      e=0;
   }
   res->s = arg1->s;
   res->e += e;

   return 0;
} /* dadd() */

/*--------------------------------------------------------------*/
/* sub two digits, return res and borrow                        */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int subdigit(Digit *d1, const Digit d2)
#else
static int subdigit(d1, d2)
      Digit *d1;
const Digit  d2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int subdigit(Digit *d1, Digit d2)
#else
static int subdigit(d1, d2)
Digit *d1;
Digit  d2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int borrow;

   if (*d1 < d2) borrow = 1;
   else borrow = 0;

   *d1 -= d2;

   return borrow;
} /* subdigit() */

/*--------------------------------------------------------------*/
/* sub two DMants, Voraussetzung: highdigit.highbit = 0         */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int dmsub(DMant *arg1, const DMant *arg2)
#else
static int dmsub(arg1, arg2)
      DMant  *arg1;
const DMant  *arg2;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int dmsub(DMant *arg1, DMant *arg2)
#else
static int dmsub(arg1, arg2)
DMant  *arg1;
DMant  *arg2;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   register int i,j;

   for (i=0; i<sizeof(DMant); i++)
      if (subdigit(&(arg1->digit[i]), arg2->digit[i])) {
         j = i;
         while (subdigit(&(arg1->digit[++j]), 1));
      }

   return 0;
} /* dmsub() */

/*--------------------------------------------------------------*/
/* sub two DReals, Voraussetzung: highdigit.highbit = 0         */
/* arg1>arg2                                                    */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dsub(const DReal *arg1, const DReal *arg2, DReal *res)
#else
int dsub(arg1, arg2, res)
const DReal  *arg1;
const DReal  *arg2;
      DReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dsub(DReal *arg1, DReal *arg2, DReal *res)
#else
int dsub(arg1, arg2, res)
DReal  *arg1;
DReal  *arg2;
DReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DMant m1;
   DMant m2;
   Digit highdigit;
   DExp e;

   e = arg1->e - arg2->e;
   if(e>0) {
      memcpy(&m1, &(arg1->m), sizeof(DMant));
      dmshift(e, &(arg2->m), &m2);
      dmsub(&m1, &m2);
      res->e = arg1->e;
   }
   else
   if(e<0) {
      dmshift(-e, &(arg1->m), &m1);
      dmsub(&m1, &(arg2->m));
      res->e = arg2->e;
   }
   else {
      memcpy(&m1, &(arg1->m), sizeof(DMant));
      dmsub(&m1, &(arg2->m));
      res->e = arg1->e;
   }

   /* --- */
   if(m1.digit[DMantLen]==0) {
      dmadjust(&m1, DMantLen, &(res->m), &e);
   }
   else {
      highdigit = m1.digit[DMantLen];
      for(e=0; highdigit; highdigit >>=1, e++);
      dmshift(e, &m1, &(res->m));
   }

   res->e += e;
   res->s = arg1->s;

   return 0;
} /* dsub() */





