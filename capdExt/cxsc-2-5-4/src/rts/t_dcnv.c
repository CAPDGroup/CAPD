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

/* CVS $Id: t_dcnv.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_dcnv.c                              */
/*                                                              */
/*                   remove unused variable i,highdigit         */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

/*--------------------------------------------------------------*/
/* convert DMant to ExtReal                                     */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dreal_to_extreal(const DReal *arg, ExtReal *res)
#else
int dreal_to_extreal(arg, res)
const DReal   *arg;
      ExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dreal_to_extreal(DReal *arg, ExtReal *res)
#else
int dreal_to_extreal(arg, res)
DReal   *arg;
ExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DReal d;
   register int i;
   int   rnd;
   int   ret;

   if(ExcPZero==(ret=normd(arg, &d))||ExcMZero==ret) {
      copyee(&Zero, res);
      return 0;
   }

#ifndef __CB_SUN4_OS4_HARDWARE
#if LSBFIRST                    /* 110391 cb */
   memcpy(res, &(d.m.digit[DMantLen-ExtMantLen]), ExtMantLen);
#else
   {
        unsigned char *p;
        int i;
        p = (unsigned char *)&(d.m.digit[DMantLen-ExtMantLen]);
        for (i=0; i<ExtMantLen; i++)
                res->s.DIGIT(i) = *(p++);
   }
#endif
#else
    {
	unsigned char *p, *q;
	int i;
	p = (unsigned char *)&(d.m.digit[DMantLen-1]);
	q = ((unsigned char *)res)+2;
	for (i=0; i<ExtMantLenTenbyte-1; i++, p--)
	    *q++ = (((*p)<<1)&0x0fe) + (((*(p-1))>>7)&0x001);
	*q++ = (((*p)<<1)&0x0fe);
	for (i=0; i<ExtMantLen-ExtMantLenTenbyte; i++)
	    *q++ = 0;
    }
#endif
   res->s.exp = ExtExpBias + d.e;

   if(d.s == -1) {
      res->s.exp |= ExtSignMask;
      if(DOWN==(rnd=getrndmode())) {
         for(i=(int)sizeof(DMant)-ExtMantLen-2;
             i>=0 && d.m.digit[i--]==0;);
         if(i>=0)
            SUBULP(res);
      }
      else
      if(NEAR==rnd)
         if(0x80 < d.m.digit[sizeof(DMant)-ExtMantLen-2]) {
            setrndmode(DOWN);
            SUBULP(res);
            setrndmode(NEAR);
         }
   }
   else { /* d->s == 1 */
      if(UP==(rnd=getrndmode())) {
         for(i=(int)sizeof(DMant)-ExtMantLen-2;
             i>=0 && d.m.digit[i--]==0;);
         if(i>=0)
            ADDULP(res);
      }
      else
      if(NEAR==rnd)
         if(0x80 < d.m.digit[sizeof(DMant)-ExtMantLen-2]) {
            setrndmode(UP);
            ADDULP(res);
            setrndmode(NEAR);
         }
   } /* else positiv */

   return 0;
} /* dreal_to_extreal() */

/*--------------------------------------------------------------*/
/* convert DMant to ExtReal.u, ExtReal.l                        */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dreal_to_2extreal(const DReal *d, ExtReal *resu, ExtReal *resl)
#else
int dreal_to_2extreal(d, resu, resl)
const DReal   *d;
      ExtReal *resu;
      ExtReal *resl;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dreal_to_2extreal(DReal *d, ExtReal *resu, ExtReal *resl)
#else
int dreal_to_2extreal(d, resu, resl)
DReal   *d;
ExtReal *resu;
ExtReal *resl;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   DMant mu;
   DMant ml;
   DExp  e;
   Digit highdigit;
   register int i;

   if(d->m.digit[DMantLen]==0) {
      if(1==dmadjust(&(d->m), DMantLen, &mu, &e)) {
         copyee(&Zero, resu);
         copyee(&Zero, resl);
         return 0;
      }
   }
   else {
      highdigit = d->m.digit[DMantLen];
      for(e=0; highdigit; highdigit >>=1, e++);
      dmshift(e, &(d->m), &mu);
   }

#ifndef __CB_SUN4_OS4_HARDWARE
#if LSBFIRST                    /* 110391 cb */
   memcpy(resu, &(mu.digit[DMantLen-ExtMantLen]), ExtMantLen);
#else
   {
        unsigned char *p;
        int i;
        p = (unsigned char *)&(mu.digit[DMantLen-ExtMantLen]);
        for (i=0; i<ExtMantLen; i++)
                resu->s.DIGIT(i) = *(p++);
   }
#endif
#else
    {
	unsigned char *p, *q;
	int i;
	p = (unsigned char *)&(mu.digit[DMantLen-1]);
	q = ((unsigned char *)resu)+2;
	for (i=0; i<ExtMantLenTenbyte-1; i++, p--)
	    *q++ = (((*p)<<1)&0x0fe) + (((*(p-1))>>7)&0x001);
	*q++ = (((*p)<<1)&0x0fe);
	for (i=0; i<ExtMantLen-ExtMantLenTenbyte; i++)
	    *q++ = 0;
    }
#endif
   resu->s.exp = ExtExpBias + d->e + e;
   if(d->s == -1) resu->s.exp |= ExtSignMask;

   /* --- resl --- */
   if(1==dmadjust(&mu, DMantLen-ExtMantLen, &ml, &e)) {
      copyee(&Zero, resl);
      return 0;
   }
#ifndef __CB_SUN4_OS4_HARDWARE
#if LSBFIRST                    /* 110391 cb */
   memcpy(resl, &(ml.digit[DMantLen-ExtMantLen]), ExtMantLen);
#else
   {
        unsigned char *p;
        int i;
        p = (unsigned char *)&(ml.digit[DMantLen-ExtMantLen]);
        for (i=0; i<ExtMantLen; i++)
                resl->s.DIGIT(i) = *(p++);
   }
#endif
#else
    {
	unsigned char *p, *q;
	int i;
	p = (unsigned char *)&(ml.digit[DMantLen-1]);
	q = ((unsigned char *)resl)+2;
	for (i=0; i<ExtMantLenTenbyte-1; i++, p--)
	    *q++ = (((*p)<<1)&0x0fe) + (((*(p-1))>>7)&0x001);
	*q++ = (((*p)<<1)&0x0fe);
	for (i=0; i<ExtMantLen-ExtMantLenTenbyte; i++)
	    *q++ = 0;
    }
#endif
   resl->s.exp = resu->s.exp + e - 64;
   if(d->s == -1) {
      if(DOWN==getrndmode()) {
         for(i=(int)sizeof(DMant)-ExtMantLen-2;
             i>=sizeof(DMant)-2*ExtMantLen-1 && ml.digit[i--]==0;);
         if(i>=sizeof(DMant)-2*ExtMantLen-1)
            SUBULP(resl);
      }
   }
   else
      if(UP==getrndmode()) {
         for(i=(int)sizeof(DMant)-ExtMantLen-2;
             i>=sizeof(DMant)-2*ExtMantLen-1 && ml.digit[i--]==0;);
         if(i>=sizeof(DMant)-2*ExtMantLen-1)
            ADDULP(resl);
      }

   return 0;
} /* dreal_to_2extreal() */

/*--------------------------------------------------------------*/
/* convert ExtReal to DMant                                     */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int extreal_to_dreal(const ExtReal *arg, DReal *d)
#else
int extreal_to_dreal(arg, d)
const ExtReal *arg;
      DReal   *d;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int extreal_to_dreal(ExtReal *arg, DReal *d)
#else
int extreal_to_dreal(arg, d)
ExtReal *arg;
DReal   *d;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{

   initd(d);
   if(0==cmpee(&Zero, arg)) return NoErr;

#ifndef __CB_SUN4_OS4_HARDWARE
#if LSBFIRST                    /* 110391 cb */
   memcpy(&(d->m.digit[DMantLen-ExtMantLen]), arg, ExtMantLen);
#else
   {
        unsigned char *p;
        int i;
        p = (unsigned char *)&(d->m.digit[DMantLen-ExtMantLen]);
        for (i=0; i<ExtMantLen; i++)
                *(p++) = arg->s.DIGIT(i);
   }
#endif
#else
    {
	unsigned char *p, *q;
	int i;
	p = (unsigned char *)&(d->m.digit[DMantLen-1]);
	q = ((unsigned char *)arg)+2;
	*p-- = 0x080 + (((*q)>>1)&0x07f);
	for (i=0; i<ExtMantLenTenbyte-1; i++, q++)
	    *p-- = (((*q)<<7)&0x080) + (((*(q+1))>>1)&0x07f);
    }
#endif
   xtrexpe(arg, &(d->e));
   if(arg->s.exp & ExtSignMask) d->s = -1; else d->s = 1;

   return NoErr;
} /* dreal_to_extreal() */





