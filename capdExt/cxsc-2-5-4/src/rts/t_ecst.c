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

/* CVS $Id: t_ecst.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename: t_ecst.c                                      */
/*                                                              */
/****************************************************************/

/* -------------------------------------------------------------*/
/* Header Konstanten fuer EmulationsVersion                     */
/* -------------------------------------------------------------*/
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#include "/u/p88c/runtime/tbyte/t_cnst.h"
#else
#include "o_defs.h"
#include "t_cnst.h"
#endif

/* #include "esin.h"       */
/*---------------------------------------------------------------*/
/* Header Emulation Sinus                                        */
/*---------------------------------------------------------------*/

/* --- Auswahl zwischen einfacher und doppelter Mantisse ----- */
/* ACHTUNG ! bei Aenderung EpsSin in File econst.c verfolgen ! */
#define HPrec   /* falls       definiert: eps(sin) = 2.725e-19 */
                /* falls nicht definiert: eps(sin) = 3.993e-19 */

#define N 11                      /* (0..N) N = Grad Approx Polynom */
#define ApproxKoeff {        /* Koeffizienten Approximationspolynom */ \
           EXTREAL(0x3F, 0xFF,               \
                   0x80, 0x00, 0x00, 0x00,   \
                   0x00, 0x00, 0x00, 0x00),  /* +1.000000000000000E+000 */ \
           EXTREAL(0xBF, 0xFC,               \
                   0xAA, 0xAA, 0xAA, 0xAA,   \
                   0xAA, 0xAA, 0xAA, 0xAB),  /* -1.666666666666667E-001 */ \
           EXTREAL(0x3F, 0xF8,               \
                   0x88, 0x88, 0x88, 0x88,   \
                   0x88, 0x88, 0x88, 0x89),  /* +8.333333333333333E-003 */ \
           EXTREAL(0xBF, 0xF2,               \
                   0xD0, 0x0D, 0x00, 0xD0,   \
                   0x0D, 0x00, 0xD0, 0x0D),  /* -1.984126984126984E-004 */ \
           EXTREAL(0x3F, 0xEC,               \
                   0xB8, 0xEF, 0x1D, 0x2A,   \
                   0xB6, 0x39, 0x9C, 0x7D),  /* +2.755731922398589E-006 */ \
           EXTREAL(0xBF, 0xE5,               \
                   0xD7, 0x32, 0x2B, 0x3F,   \
                   0xAA, 0x27, 0x1C, 0x7F),  /* -2.505210838544172E-008 */ \
           EXTREAL(0x3F, 0xDE,               \
                   0xB0, 0x92, 0x30, 0x9D,   \
                   0x43, 0x68, 0x4B, 0xE5),  /* +1.605904383682161E-010 */ \
           EXTREAL(0xBF, 0xD6,               \
                   0xD7, 0x3F, 0x9F, 0x39,   \
                   0x9D, 0xC0, 0xF8, 0x8F),  /* -7.647163731819816E-013 */ \
           EXTREAL(0x3F, 0xCE,               \
                   0xCA, 0x96, 0x3B, 0x81,   \
                   0x85, 0x6A, 0x53, 0x59),  /* +2.811457254345521E-015 */ \
           EXTREAL(0xBF, 0xC6,               \
                   0x97, 0xA4, 0xDA, 0x34,   \
                   0x0A, 0x0A, 0xB9, 0x26),  /* -8.220635246624329E-018 */ \
           EXTREAL(0x3F, 0xBD,               \
                   0xB8, 0xDC, 0x77, 0xB6,   \
                   0xE7, 0xAB, 0x8C, 0x5F),  /* +1.957294106339126E-020 */ \
           EXTREAL(0xBF, 0xB4,               \
                   0xBB, 0x0D, 0xA0, 0x98,   \
                   0xB1, 0xC0, 0xCE, 0xCC)   /* -3.868170170630684E-023 */ \
} /* --- Ende Koeffizienten Approximationspolynom --- */

#if SUN4_CPP_C
/* switch off caused by initialization difficulties */
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

#endif /* ANSI_C */

#ifdef HPrec
const ExtReal EpsSin  = EPS_SIN_HPrec;
const ExtReal EpsCos  = EPS_SIN_HPrec;
#else
const ExtReal EpsSin  = EPS_SIN;
const ExtReal EpsCos  = EPS_SIN;
#endif

const ExtReal EpsTan  = EPS_TAN;
const ExtReal EpsCot  = EPS_TAN;
const ExtReal EpsASin = EPS_ASIN;
const ExtReal EpsACos = EPS_ACOS;
const ExtReal EpsATan = EPS_ATAN; /* mit EPS_CO keine Fehler gefunden! */
const ExtReal EpsACot = EPS_CO;
const ExtReal EpsSinh = EPS_SINH;
const ExtReal EpsCosh = EPS_COSH;
const ExtReal EpsTanh = EPS_TANH;
const ExtReal EpsCoth = EPS_COTH;
const ExtReal EpsExp  = EPS_EXP;
const ExtReal EpsLn   = EPS_LN;
const ExtReal EpsLnRel= EPS_LN_REL;
const ExtReal EpsLnAbs= EPS_LN_ABS;
const ExtReal EpsSqrt = EPS_SQRT;

const ExtReal EpsASinh = EPS_ASINH;
const ExtReal EpsACosh = EPS_ACOSH;
const ExtReal EpsATanh = EPS_ATANH;
const ExtReal EpsACoth = EPS_ACOTH;

const ExtReal EpsPow   = EPS_POW;

#if VAX_VMS_C
#ifdef LINT_ARGS
local void t_ecst(void)
#else
local void t_ecst()
#endif
        {
        }
#endif





