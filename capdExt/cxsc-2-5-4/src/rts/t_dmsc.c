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

/* CVS $Id: t_dmsc.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/**************************************************************/
/* Name       : t_dmsc.c                                      */
/* Description: init, split, shift and adjust double mantissa */
/* Date       : 1992-11-26                                    */
/* 1992-11-26 : flipped conditions in calculate shift   =b=   */
/*              dmadjust                                      */
/*              add description tex                           */
/* 1992-07-22 : changed new to new_d      =b=                 */
/* 1992-07-21 : add sun4_gnu_c sun4_cpp_c =b=                 */
/* 1992-07-08 : =b= removed emshift                           */
/**************************************************************/
#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ddev.h"
#else
#include "t_ddev.h"
#endif

/* aus stdlib.h, da in c++ nicht enthalten */
#define max(a,b)    (((a) > (b)) ? (a) : (b))

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

/*--------------------------------------------------------------*/
/* initd                                                        */
/*--------------------------------------------------------------*/
#ifdef LINT_ARGS
int initd(DReal *d)
#else
int initd(d)
DReal *d;
#endif
{
   memset(&(d->m), 0x00, sizeof(DMant));
   d->e = 0;
   d->s = 0;

   return 0;
} /* initd() */

/*--------------------------------------------------------------*/
/* split mant von arg in mitte und lege in resh, bzw. resl ab   */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int msplitee(const ExtReal *arg, ExtReal *resh, ExtReal *resl)
#else
int msplitee(arg, resh, resl)
const ExtReal *arg;
      ExtReal *resh;
      ExtReal *resl;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int msplitee(ExtReal *arg, ExtReal *resh, ExtReal *resl)
#else
int msplitee(arg, resh, resl)
ExtReal *arg;
ExtReal *resh;
ExtReal *resl;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   copyee(arg, resh);
#ifndef __CB_SUN4_OS4_HARDWARE
#if LSBFIRST                    /* 110391 cb */
#if SHORTABTYP
   *((a_btyp *)resh) = 0;
#else
   *((long *)resh) = 0;
#endif
#else
   *((short int *)(((char *)resh)+6)) = 0;
   *((short int *)(((char *)resh)+8)) = 0;
#endif
#else
#if SHORTABTYP
   *(((a_btyp *)resh)+1) &= ((a_btyp)0x0fffe0000);
   *(((a_btyp *)resh)+2) = 0;
   *(((a_btyp *)resh)+3) = 0;
#else
   *(((long *)resh)+1) &= 0x0fffe0000;
   *(((long *)resh)+2) = 0;
   *(((long *)resh)+3) = 0;
#endif
#endif
   subee(arg, resh, resl);

   return 0;
} /* msplitee() */

/*----------------------------------------------------------------------*/
/* shift DMant m um fuehrende 0-Bits nach rechts und lege in mr ab      */
/*----------------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dmadjust(const DMant *arg, const int len, DMant *res, DExp *shiftr)
#else
int dmadjust(arg, len, res, shiftr)
const DMant *arg;
const int len;        /* adjust arg->digit[0] .. arg->digit[len] */
DMant *res;
DExp *shiftr;         /* Ergebnis, immer < 0 */
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dmadjust(DMant *arg, int len, DMant *res, DExp *shiftr)
#else
int dmadjust(arg, len, res, shiftr)
DMant *arg;
int len;           /* adjust arg->digit[0] .. arg->digit[len] */
DMant *res;
DExp *shiftr;      /* Ergebnis, immer < 0 */
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
#ifdef ANSI_C
   const
#endif
   Digit *ptr;
   Digit mask;
   Digit new_d;
   Digit old;
   DExp digits;
   DExp bits;
   int imax;
   register int i,j;

   /* calculate shifts */
   for (ptr = &(arg->digit[len-1]), digits=0;
                       digits < len && *ptr == 0 ; ptr --, digits++);
   if (digits == len) return 1; /* Zero */
   for (mask = 0x80, bits=0; (*ptr & mask) == 0; mask >>= 1, bits++);

   /* clear mant */
   memset(res, 0x00, sizeof(DMant));

   /* --- calc shiftr, leave in case no adjust --- */
   if(0==(*shiftr = (DExp) -(digits * BitsPerDigit + bits))) {
      memcpy(&(res->digit[DMantLen-len]), arg, len);
      return 0;
   }

   /* --- copy arg nach res, start high end --- */
   imax = len-digits;
   old = 0;
   for (i=0, j=digits+DMantLen-len; i<imax; i++, j++) {
      new_d =  arg->digit[i];
      res->digit[j] = old | (new_d << bits);
      old = new_d >> (BitsPerDigit-bits);
   }

   return 0;
} /* dmadjust() */

/*--------------------------------------------------------------*/
/* shift arg um e Stellen nach low und lege in DMant res ab     */
/* arg mindestens ein digit != 0x00                             */
/* high digit immer 0 fuer carry                                */
/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int dmshift(const DExp sh, const DMant *arg, DMant *res)
#else
int dmshift(sh, arg, res)
const DExp sh; /* Anzahl shifts, immer >= 0 (SHift) */
const DMant *arg;
DMant *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int dmshift(DExp sh, DMant *arg, DMant *res)
#else
int dmshift(sh, arg, res)
DExp sh; /* Anzahl shifts, immer >= 0 (SHift) */
DMant *arg;
DMant *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   Digit new_d;
   Digit old;
   DExp shdigits; /* Anzahl ganze Digits fuer shift */
   DExp shbits;   /* Anzahl bits fuer shift */
   int imin;
   register int i, j;

   /* clear mant */
   memset(res, 0x00, sizeof(DMant));

   if(sh>=sizeof(DMant)*BitsPerDigit)
      return 0;

   /* calculate shifts */
   shdigits = sh / BitsPerDigit;
   shbits = sh % BitsPerDigit;

   for(i=0; 0x00==arg->digit[i++];); /* Letztes Digit != 0 */
   j = (int)sizeof(DMant) - --i;          /* Anzahl Digits != 0 */

   if(shbits == 0) {
      if(0>(imin=(int)sizeof(DMant)-shdigits-j)) {
         i-=imin;
         j+=imin;
         imin=0;
      }
      memcpy(&(res->digit[imin]), &(arg->digit[i]), j);
   }
   else {
      /* --- copy arg nach m, start high end --- */
      if(0>(imin=(int)sizeof(DMant)-shdigits-j))
         imin = i-imin+1;
      else
         imin = (int)max(imin, i);
      old = 0;
      for (i=DMantLen, j=DMantLen-shdigits; i>=imin; ) {
         new_d = arg->digit[i--];
         res->digit[j--] = old | (new_d >> shbits);
         old = new_d << (BitsPerDigit-shbits);
      }
      if (j>=0) res->digit[j] = old;
   }

   return 0;
} /* dmshift() */

#if 0  /* =b= obsolete */
/*------------------------------------------------------------------*/
/* shift mant von arg um e Stellen nach low und lege in DMant m ab  */
/* high digit immer 0 fuer carry                                    */
/*------------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int emshift(const DExp sh, const ExtReal *arg, DMant *m)
#else
int emshift(sh, arg, m)
const DExp sh; /* Anzahl shifts, immer >= 0 (SHift) */
const ExtReal *arg;
DMant *m;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int emshift(DExp sh, ExtReal *arg, DMant *m)
#else
int emshift(sh, arg, m)
DExp sh; /* Anzahl shifts, immer >= 0 (SHift) */
ExtReal *arg;
DMant *m;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   Digit new_d;
   Digit old;
   DExp shdigits; /* Anzahl ganze Digits fuer shift */
   DExp shbits;   /* Anzahl bits fuer shift */
   int imin;
   register int i, j;

   /* clear mant */
   memset(&m, 0x00, sizeof(m));

   /* calculate shifts */
   shdigits = DMantLen - ExtMantLenTenbyte - sh / BitsPerDigit;
   shbits = sh % BitsPerDigit;

   /* --- copy arg nach m, start high end --- */
   imin = shdigits ? shdigits : 0;
   old = 0;
   for (i=ExtMantLenTenbyte-1, j=i+shdigits; i>imin; ) {
      new_d = arg->s.DIGIT(i--);
      m->digit[j--] = old | (new_d >> shbits);
      old = new_d << (BitsPerDigit-shbits);
   }
   if (shdigits>=0) m->digit[shdigits] = old;

   return 0;
} /* emshift() */
#endif /* obsolete */





