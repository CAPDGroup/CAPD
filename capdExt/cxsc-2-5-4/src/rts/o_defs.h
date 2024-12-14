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

/* CVS $Id: o_defs.h,v 1.22 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_defs.h                              */
/*                                                              */
/*      Description     : definitions included by all runtime   */
/*                        source files                          */
/*                                                              */
/****************************************************************/

#ifndef _O_DEFS_LOADED 
#define _O_DEFS_LOADED 

#ifdef AIX
#include "/u/p88c/runtime/o_slct.h"
#else
#include "o_slct.h"
#endif

#if IBM_RT_C+IBM_RS6000_C
#else
#undef AIX
#endif

#include <ctype.h>
#if SUN4_GNU_C
int toupper(int c);
#endif

#if GNU_C+ZORTECH_C
#include <varargs.h>
/* #define _VA_LIST_ */
#endif
#include <stdio.h>
#if SUN4_CPP_C
#include <sysent.h>
#define remove(path)    unlink(path)
#endif

#if IBM_RT_C+IBM_RS6000_C+HP_9000_C+CONVEX1_UNIX_C
#define remove(path)    unlink(path)
#else
#if SUN4_GNU_C
int fclose(FILE *stream);
int fflush(FILE *stream);
int fprintf(FILE *stream,const char *format,...);
int fgetc(FILE *stream);
int fputc(int c,FILE *stream);
int ungetc(int c, FILE *stream);
int sscanf(const char *s,const char *format,...);
int remove(const char *filename);
int _filbuf(FILE *stream);
#endif
#endif

#if IBM_RT_C+SUN4_OS4_C+CONVEX1_UNIX_C /* IBM_RS6000_C weg 140396 cb */
typedef unsigned long int size_t;
#else
#if IBM_370_C
#include <stdefs.h>
#else
#if GNU_C+ZORTECH_C+IBM_EMX_C
#else
#include <stddef.h>
#endif
#endif
#endif

#if IBM_RT_C+HP_9000_C+CONVEX1_UNIX_C /* IBM_RS6000_C weg 140396 cb */
#define remove(path)    unlink(path)
char *calloc();
char *malloc();
char *getenv();
int system();
int abs();
#else
#include <stdlib.h>
#if SUN4_GNU_C
int system(const char *string);
#endif
#endif

/* IBM_RS6000_C weg 140396 cb */
#if IBM_370_C+IBM_RT_C+HP_9000_C+SUN4_OS4_C+CONVEX1_UNIX_C+\
    GNU_C+ZORTECH_C+SUN4_GNU_C+SUN4_CPP_C+IBM_EMX_C
#define DBL_MAX_10_EXP          308
#define DBL_MIN_10_EXP        (-307)
#define DBL_MANT_DIG             53
#else
#include <float.h>
#endif

#define DUMMY_ARGS IBM_PS2_C+IBM_RT_C+IBM_RS6000_C+HP_9000_C+SUN4_OS4_C+CONVEX1_UNIX_C
#if DUMMY_ARGS+DEC_ULTRIX_C+VAX_VMS_C
#include <varargs.h>
#else
#include <stdarg.h>
#endif

#if CONVEX1_UNIX_C
#include <strings.h>
char *memchr();
char *memcpy();
char *memset();
#else
#include <string.h>
#endif
#if IBM_RT_C+IBM_RS6000_C+HP_9000_C+SUN4_OS4_C+SUN4_GNU_C+DEC_ULTRIX_C+\
    SUN4_CPP_C
#include <memory.h>
#endif

/* #if IBM_RT_C+IBM_RS6000_C+IBM_PS2_C+SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C */
#if IBM_RT_C+IBM_RS6000_C+IBM_PS2_C+SUN4_OS4_C+SUN4_CPP_C
#define CLOCKS_PER_SEC  (1000000L)
extern long clock();
#else

#if ZORTECH_C+IBM_EMX_C
#include <time.h>
#else
#if GNU_C
#define CLOCKS_PER_SEC  (0)
#else
#if HP_9000_C
/* #include <limits.h> */
#define CLOCKS_PER_SEC  (100)
#else
#include <time.h>
#if ATARI_TURBO_C+T800_HELIOS_C+IBM_370_C+IBM_AT_TURBO_C
#define CLOCKS_PER_SEC  CLK_TCK
#endif
#endif
#endif
#endif
#endif

#ifdef AIX
#include "/u/p88c/runtime/p88rts.h"
#include "/u/p88c/runtime/o_spec.h"
#include "/u/p88c/runtime/o_syst.h"
#include "/u/p88c/runtime/o_type.h"

#include "/u/p88c/runtime/int/a_defs.h"
#include "/u/p88c/runtime/base/b_defs.h"
#include "/u/p88c/runtime/dot/d_defs.h"
#include "/u/p88c/runtime/error/e_defs.h"
#include "/u/p88c/runtime/multi/l_defs.h"
#include "/u/p88c/runtime/real/r_defs.h"

#include "/u/p88c/runtime/o_name.h"

#include "/u/p88c/runtime/o_fcth.h"
#ifdef DEC_ARITH
#include "/u/p88c/runtime/g_fcth.h"
#endif
#else
#include "p88rts.h"
#include "o_spec.h"
#include "o_syst.h"
#include "o_type.h"

#include "a_defs.h"
#include "b_defs.h"
#include "d_defs.h"
#include "e_defs.h"
#include "l_defs.h"
#include "r_defs.h"

#include "o_name.h"

#include "o_fcth.h"
#ifdef DEC_ARITH
#include "g_fcth.h"
#endif
#endif

/* -------- 64 bit accumulator -------------------------------- */
/* redefines r_defs and d_defs macros suitable for 64 bit       */
/* by appending _64 to all constants' names                     */
#ifdef IS_64_BIT
#include "b_64bt.h"
#endif

/* -------- Test constant ------------------------------------- */

/* Heap checking is on if HEAP_CHECK is defined*/
/*#define HEAP_CHECK */

#endif /* _O_DEFS_LOADED */





