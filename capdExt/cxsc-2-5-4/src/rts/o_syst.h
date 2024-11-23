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

/* CVS $Id: o_syst.h,v 1.23 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_syst.h                              */
/*                                                              */
/****************************************************************/

#include "cxscconf.h"

#ifdef SUN4_FORTE
#undef SUN4_GNU_C
#define SUN4_GNU_C 1  
#endif

/* -------- Operating System Selection ------------------------ */

/* This option specifies the way data is stored in memory       */
/* cells.                                                       */
/* 0 - most significant byte at high address                    */
/* 1 - most significant byte at low address (Intel)             */
#define DUMMY_INTEL IBM_AT_MS_C+IBM_AT_TURBO_C+IBM_PS2_C+T800_HELIOS_C\
                    +IBM_AT_BOR_C+IBM_OS2_ICC_C+IBM_EMX_C+IBM_LINUX_C
#if DUMMY_INTEL+GNU_C+ZORTECH_C+VAX_VMS_C+VAX_UNIX_C+DEC_ULTRIX_C+DEC_ALPHA_C+\
    CHAM_32_C+CHAM_64_C+IBM_WATCOM_C
#define INTEL           1
#else
#define INTEL           0
#endif

/* -------- Operating System Properties ----------------------- */

/*    Special memory set operation available which initializes  */
/*    a certain amount of memory.                               */
/*    1 - char *memset(char *s,int c,size_t n)                  */
#define O_P_1           1

/*    Special memory copy operation available which copies      */
/*    a certain amount of memory.                               */
/*    1 - char *memcpy(char *s1,char *s2,size_t n)              */
#define O_P_2           1

/*    Special memory compare operation available which compares */
/*    a certain amount of memory.                               */
#define O_P_3           0

/*    Data for Intel processors is stored in reverse order, i.e.*/
/*    high storage addresses hold upper parts of data.          */
/*    This option is 1 if INTEL is selected.                    */
#define O_P_4           INTEL

/* -------- Operating System Dependent Macros ----------------- */

/*    Clear all bits of a 'a_btyp' array.                       */
/*    Actual argument s must not be an expression.              */
#if O_P_1
#define B_CLEAR(s,n)    { \
        (void)memset((char *)(s),(int)'\0',(size_t)((n)*sizeof(a_btyp)));}
#else
#define B_CLEAR(s,n) {size_t _; _=n; while(--_>=0) (s)[_]=ZERO;}
#endif

/*    Copy all bits of a 'a_btyp'/'char' array from left to right.     */
/*    Actual arguments d and s must not be expressions.                */
#if O_P_2
#define B_COPY(d,s,n)   { \
        (void)memcpy((char *)(d),(char *)(s),(size_t)((n)*sizeof(a_btyp)));}
#define C_COPY(d,s,n)   { \
        (void)memcpy((char *)(d),(char *)(s),(size_t)(n)); }
#else
#define B_COPY(d,s,n) {size_t _; _=n; while(--_>=0) (d)[_]=(s)[_];}
#define C_COPY(d,s,n) {size_t _; _=n; while(--_>=0) (d)[_]=(s)[_];}
#endif

/*    Compare all bits of 'a_btyp' arrays from left to right.   */
/*    Actual arguments d and s must not be expressions.         */
#if O_P_3
#else
#endif

/*    Specify indices of most significant and least significant */
/*    part of a a_btyp array which is equivalenced with a real  */
/*    number (D_U_RATIO = sizeof(a_real)/sizeof(a_btyp))        */
#if O_P_4
#define B_HPART         (D_U_RATIO-1)
#define B_LPART         0
#define B_TYPINI(a,b)   { b,a }       
#else
#define B_HPART         0
#define B_LPART         (D_U_RATIO-1)
#if CRAY_UNIX_C
#define B_TYPINI(a,b)   { (a<<B_LENGTH)+b }       
#else
#define B_TYPINI(a,b)   { a,b }       
#endif
#endif





