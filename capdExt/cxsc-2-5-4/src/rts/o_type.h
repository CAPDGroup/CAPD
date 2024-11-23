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

/* CVS $Id: o_type.h,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_type.h                              */
/*                                                              */
/*                   changed new to new_d in a_pptr             */
/****************************************************************/

typedef d_otpr dotprecision;

typedef struct bent {
                int line;
                char *function;
                char *filename;
                struct bent *pred;
                struct bent *succ;
               } bentry;

typedef struct aent {
                a_btyp   code;                  /* condition code group */
                unsigned int active : 1;        /* output generated     */
                unsigned int back : 1;          /* traceback generated  */
                unsigned int cont : 1;          /* execution continued  */
                unsigned int def : 1;           /* default action       */
                unsigned int stat : 1;          /* statically defined   */
                struct aent *succ;              /* linkage              */
                                                /* exception handler    */
#ifdef LINT_ARGS
#ifndef __cplusplus
                void (*action)(a_btyp,int,va_list);
#endif
#else
                void (*action)();
#endif
               } aentry;

typedef struct { int msgid; char *text; } e_mtyp;

struct intern
       {
       unsigned int z : 1; /* z    ==  0  number is non-zero           */
                           /* z    ==  1  identifies number zero       */
       unsigned int s : 1; /* s    ==  0 identifies positive numbers   */
                           /* s    ==  1 identifies negative numbers   */
       unsigned int r : 2; /* number of ulps for outer bound           */
       unsigned int f : 1; /* f    ==  0 permanent value               */
                           /* f    ==  1 temporary value               */
       a_intg e;           /* exponent of most significant a_btyp      */
                           /* mantissa digit m[0]                      */
                           /* the exponent of m[i] is e-i              */
       a_btyp l;           /* number of mantissa digits                */
       a_btyp *m;          /* mantissa digits                          */
       };

typedef struct intern *multiprecision;

typedef struct a_pptr
      {
      a_VOID new_d,  /* New Display. NULL if already installed. */
             old;    /* Old Display. */
      struct a_pptr * prev;
      } a_pptr;

#ifdef IEEE_HARDWARE
#if SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C
/* #define __CB_SUN4_OS4_HARDWARE */
#endif
#endif

#ifdef __CB_SUN4_OS4_HARDWARE
#define t_size          16
#define ExtMantLen      14
#define ExtMantLenTenbyte 8	/* fuer t_dcnv, t_dmsc, t_rtrg */
#else
#define t_size          10
#define ExtMantLen       8
#define ExtMantLenTenbyte 8
#endif

typedef unsigned short int  ExtExp;
typedef unsigned char Digit;

typedef union {
   unsigned char c[t_size];
#if SUN4_CPP_C
   struct {
#else
   struct s{
#endif
#if INTEL
      Digit  digit[ExtMantLen];
      ExtExp exp;
#else
      ExtExp exp;
      Digit  digit[ExtMantLen];
#endif
   } s;
} TByte;

typedef TByte tenbyte;

/*
typedef struct a_tenb { unsigned char c[t_size]; } tenbyte;
*/





