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

/* CVS $Id: t_cnst.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_cnst.c                              */
/*                                                              */
/****************************************************************/

/* -------------------------------------------------------------*/
/* Develop Konstanten fuer beide Versionen                      */
/* -------------------------------------------------------------*/
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#include "/u/p88c/runtime/tbyte/t_cnst.h"
#else
#include "o_defs.h"
#include "t_cnst.h"
#endif

#if SUN4_CPP_C
/*  switch off ansi_c flag caused by  initialization difficulties */
#undef ANSI_C
#endif

#ifdef ANSI_C

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#else  /* NOT ANSI_C */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_mach.h"
#else
#include "t_mach.h"
#endif

#define const

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

#ifdef __CB_SUN4_OS4_HARDWARE	/* fuer SUN4 ist LSBFIRST=false */
#  define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
	{  /* Exponent+VZ */	y0, y1,	\
	   /* 64bit Mant  */    (((y2<<1)&0x0fe)+((y3>>7)&0x001)),	\
				(((y3<<1)&0x0fe)+((y4>>7)&0x001)),	\
				(((y4<<1)&0x0fe)+((y5>>7)&0x001)),	\
				(((y5<<1)&0x0fe)+((y6>>7)&0x001)),	\
				(((y6<<1)&0x0fe)+((y7>>7)&0x001)),	\
				(((y7<<1)&0x0fe)+((y8>>7)&0x001)),	\
				(((y8<<1)&0x0fe)+((y9>>7)&0x001)),	\
				((y9<<1)&0x0fe),			\
	    /* Rest 0     */	0, 0, 0, 0, 0, 0  }
#else
#if LSBFIRST                    /* aus ieeedev.h */
#  define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
                 {y9,y8,y7,y6,y5,y4,y3,y2,y1,y0}
#else
#  define EXTREAL(y0,y1,y2,y3,y4,y5,y6,y7,y8,y9) \
                 {y0,y1,y2,y3,y4,y5,y6,y7,y8,y9}
#endif /* LSBFIRST */
#endif

typedef char ExtReal[t_size];

extern char arg_check;
#define On 1

#endif /* ANSI_C */

/* -------------------------------------------------------------*/
/* Globale Groessen Argument Pruefung                           */
/* -------------------------------------------------------------*/
char arg_check = On;

/* -------------------------------------------------------------*/
/* ExtReal                                                      */
/* -------------------------------------------------------------*/
/* const ExtReal Epsl    = EPS_L; */

/* const ExtReal IntMax =       */            /* (int)0x7fff = INT_MAX */
/* EXTREAL(0x40, 0x0D,
         0xFF, 0xFE, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00); */     /* +3.276700000000000E+004 */
/* const ExtReal IntMin =          */         /* (int)0xffff */
/* EXTREAL(0xC0, 0x0E,
         0x80, 0x00, 0x00, 0x00,
         0x00, 0x00, 0x00, 0x00); */     /* -3.276800000000000E+004 */

const ExtReal IntMax =                  /* 32-Bit-Integer !!! */
EXTREAL(0x40, 0x1d,
        0xff, 0xff, 0xff, 0xfe,
        0x00, 0x00, 0x00, 0x00);           /* +2.147483647000000E+09 */
const ExtReal IntMin =
EXTREAL(0xc0, 0x1e,
        0x80, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00);           /* -2.147483648000000E+09 */
const ExtReal LongRealMax =        /* = (double)DBL_MAX in float.h */
EXTREAL(0x43, 0xFE,
        0xFF, 0xFF, 0xFF, 0xFF,
        0xFF, 0xFF, 0xF8, 0x00);   /* +1.797693134862316E+308 */
const ExtReal LongRealDenormMin =  /* kleinste pos. denormalisierte Zahl */
EXTREAL(0x3b, 0xcd,
        0x80, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00);   /* +4.940656458412465E-324 */
const ExtReal LongRealMin =        /* = (double)DBL_MIN in float.h */
EXTREAL(0x3C, 0x01,
        0x80, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00);      /* +2.225073858507201E-308 */

/*  Konstanten fuer 80-bit Ergebnis */
/* const ExtReal MaxArgExp =
EXTREAL(0x40, 0x0C,
        0xB1, 0x72, 0x17, 0xF7,
        0xD1, 0xCF, 0x79, 0xAB); */ /*+1.135652340629414E+004*/
/* const ExtReal MinArgExp =
EXTREAL(0xC0, 0x0C,
        0xB1, 0x6C, 0x8C, 0x67,
        0x12, 0x10, 0xEB, 0x30); */ /*-1.135513711193302E+004*/

/*  Konstanten fuer 64-bit Ergebnis */
const ExtReal MaxArgExp =
EXTREAL(0x40, 0x08,
        0xB1, 0x72, 0x17, 0xF7,
        0xD1, 0xCF, 0x79, 0xAA); /* +7.097827128933840E+002 */
const ExtReal MinArgExp =
EXTREAL(0xC0, 0x08,
        0xBA, 0x1C, 0x2A, 0x23,
        0x6B, 0x8E, 0x1B, 0x1C); /* -7.444400719213812E+002 */

const ExtReal MaxArgHyp =
EXTREAL(0x40, 0x08,
        0xB1, 0x9E, 0x74, 0x7D,
        0xCF, 0xC3, 0xED, 0x87); /* +7.104758600739439E+002 */
const ExtReal MinArgHyp =
EXTREAL(0xC0, 0x08,
        0xB1, 0x9E, 0x74, 0x7D,
        0xCF, 0xC3, 0xED, 0x87); /* -7.104758600739439E+002 */

const ExtReal MaxArgLnp1 =               /* sqrt(2)/2               */
EXTREAL(0x3F, 0xFD,
        0x80, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00);      /* +2.500000000000000E-001 */
const ExtReal MaxArgSqrtm1 =
EXTREAL(0x3F, 0xFA,
        0x80, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00);      /* +3.125000000000000E-002 */

#if VAX_VMS_C
#ifdef LINT_ARGS
local void t_cnst(void)
#else
local void t_cnst()
#endif
        {
        }
#endif





