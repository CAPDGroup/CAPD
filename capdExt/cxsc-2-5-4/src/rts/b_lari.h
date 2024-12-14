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

/* CVS $Id: b_lari.h,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : b_lari.h           Processor : C                   *
 *                                                                       *
 * Header File for Multiple Precision Standard Functions                 *
 * =====================================================                 *
 *                                                                       *
 * Include files :  macros.h  - Definitions for floating point arithmetic*
 *                  o_defs.h  - Definitions                              *
 *                  Lnames.h  - Names of routines and global variables   *
 *                                                                       *
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif

#undef max
#define max(a,b)                (a<b ? b : a)
#undef min
#define min(a,b)                (a>b ? b : a)

/*---------------------------------------------------------------------*
 | Definitions                                                         |
 | ===========                                                         |
 *---------------------------------------------------------------------*/

#ifdef Main
#define Reset_Name  Lroutine = NULL
#define Reset_Prec  Maxl = Lcurrprec
#else
#define Reset_Name  ;
#define Reset_Prec  Maxl = Currprec
#endif


/*---------------------------------------------------------------------*
 | Macros for argument check                                           |
 | =========================                                           |
 *---------------------------------------------------------------------*/

#define CHKNAN(x)    if (x->m[0] == 0) ERREXIT(NANDE,NANDE,0)


/*---------------------------------------------------------------------*
 | Macro for return from subroutine                                    |
 | ================================                                    |
 |                                                                     |
 |  Parameters :  RC         Return code                               |
 *---------------------------------------------------------------------*/

#ifdef Debug
#define EXIT(RC) { Reset_Name;\
           if (Ldebug >= 0)\
              printf("\n Leaving  Routine %s, rc = %d\n",LRoutine,RC);\
           Ldebug += 1;\
           return(RC); }
#else
#define EXIT(RC) { Reset_Name;  return(RC); }
#endif

#define RETURN(RC) { Reset_Prec;  EXIT(RC); }


/*---------------------------------------------------------------------*
 | Macro for immediate exit in case of errors                          |
 | ==========================================                          |
 |                                                                     |
 |  Parameters :  ERR        Error code                                |
 |                RC         Return code                               |
 |                DROP       Number of Pool Variables currently in use |
 *---------------------------------------------------------------------*/

#define ERREXIT(ERR,RC,DROP) { Lerror(ERR);  Lvardrop(DROP);  RETURN(RC); }


/*---------------------------------------------------------------------*
 | Macros for floating point standard functions                        |
 | ============================================                        |
 *---------------------------------------------------------------------*/

#ifdef NODSTD
#define B_ASIN(x)       (b_kasn(x))
#define B_ATAN(x)       (b_katn(x))
#define B_LOG_(x)       (b_klog(x))
#define B_SQRT(x)       (b_ksqt(x))
#else
#ifdef LINT_ARGS
extern double asin(double),atan(double),log(double),sqrt(double);
#else
extern double asin(),atan(),log(),sqrt();
#endif
#define B_ASIN(x)       asin(x)
#define B_ATAN(x)       atan(x)
#define B_LOG_(x)       log(x)
#define B_SQRT(x)       sqrt(x)
#endif


/*---------------------------------------*
 | Constants describing internal numbers |
 *---------------------------------------*/

#define LdBase    32           /* Number of dual digits of base          */
#define Minl      1            /* Minimal length of INTERN numbers       */
#define Minlerr   2            /* Minimal length for error estimates     */
#define Minlfl    3            /* Minimal length for double numbers      */
#define Minlflsqr 5            /* Minimal length for square of double's  */
#if SHORTABTYP
#define Smask     0x80000000  /* Mask for sign bit of type double       */
#define Emask     0x7FF00000  /* Mask for characteristic of type double */
#define EFUFfac   0x0000001F  /* Exponent if INTERN underflow factor    */
#define MFUFfac   0x7DF00000  /* Leading part of double underflow factor*/
#define MFIUFfac  0x01F00000  /* Leading part of inverse underflow factor*/
#define Mindbl    0x00100000  /* Leading 4 bytes of MINREAL            */
#define Mindble   (-20)        /* Exponent of (INTERN)(MINREAL)         */
#define Minsqrd   0x20000000  /* Leading 4 bytes of sqrt(MINREAL)      */
#else
#define Smask     0x80000000L  /* Mask for sign bit of type double       */
#define Emask     0x7FF00000L  /* Mask for characteristic of type double */
#define EFUFfac   0x0000001FL  /* Exponent if INTERN underflow factor    */
#define MFUFfac   0x7DF00000L  /* Leading part of double underflow factor*/
#define MFIUFfac  0x01F00000L  /* Leading part of inverse underflow factor*/
#define Mindbl    0x00100000L  /* Leading 4 bytes of MINREAL            */
#define Mindble   (-20)        /* Exponent of (INTERN)(MINREAL)         */
#define Minsqrd   0x20000000L  /* Leading 4 bytes of sqrt(MINREAL)      */
#endif
#define Minsqre   (-10)        /* Exponent of (INTERN)(sqrt(MINREAL))   */
#define Minsqrm   2            /* Mantissa of (INTERN)(sqrt(MINREAL))   */

#define Lpi_init  100          /* Initialized number of digits of pi    */


/*----------------*
 | Boolean Values |
 *----------------*/

#define true  1
#define false 0


/*----------------*
 | Error Codes :  |
 *----------------*/

#define CONVD     1003     /* Error during conversion to double         */
#define DUFLW     1006     /* Floating point underflow while computing  */
#define PEVAL     1001     /* Error during polynomial evaluation        */
#define RESUL     1002     /* Error during result computation           */
#define ASSGN      999     /* Error during assignment of result         */
#define ERRBD      998     /* Error during computation of error bound   */
#define EPERR     1004     /* Error during computation of error bound   */
#define WIDTH     1007     /* Error due to insufficient error bound     */
#define UNITS     1005     /* Error on computation of number of ulp's   */


/*-------------*
 | Roundings : |
 *-------------*/

#define up        '^'      /* Rounding upwards                          */
#define down      'v'      /* Rounding downwards                        */
#define nearest   'n'      /* Rounding nearest                          */
#define lbound    'l'      /* Value is lower bound                      */
#define ubound    'u'      /* Value is upper bound                      */
#define ibound    'i'      /* Value is rounded towards 0                */
#define obound    'o'      /* Value is rounded towards infinity         */
#define rounded   'r'      /* Value is rounded                          */
#define signed    's'      /* Value is rounded                          */
#define exact     'e'      /* Value is exact                            */
#define absolute  'a'      /* Error is absolute with correct sign       */


/*-----------------------*
 | Definition of Types : |
 *-----------------------*/

typedef struct intern  dynamic;    /* Short hand for type INTERN       */
typedef char           rounding;   /* Rounding specifications          */


/*----------------------------*
 | Number of pool variables : |
 *----------------------------*/

#define numvar    30           /* Number of Pool Variables            */



/*--------------*
 | Arithmetic : |
 *--------------*/

#define ADD_(IOP1,IOP2,IRES)       rc += b_badd(IOP1,IOP2,IRES)
#define SUB_(IOP1,IOP2,IRES)       rc += b_bsub(IOP1,IOP2,IRES)
#define MUL_(IOP1,IOP2,IRES,IRMD)  rc += b_bmul(IOP1,IOP2,IRES,IRMD)
#define DIV_(IOP1,IOP2,IRES)       rc += b_bdiv(IOP1,IOP2,IRES)
#define MULINT_(IOP,INT,IRES)      rc += b_bmun(IOP,INT,IRES)
#define DIVINT_(IOP,INT,IRES)      rc += b_bdvn(IOP,INT,IRES)
#define NEXT_(IOP,IRES)            rc += b_bnxt(IOP,IRES)
#define COPY_(IOP,IRES)            rc += b_bcpy(IOP,IRES)
#define SHIFT_(N,IOP,IRES)         rc += b_bshf(N,IOP,IRES)
#define DBLTODYN(FOP,IOP,IRCV)     IRCV = b_bcdi(FOP,IOP);\
                                   if (IRCV != ROUND) rc += IRCV;
/* #define DBLTODYN(FOP,IOP,IRCV)     rc += b_bcdi(FOP,IOP) */
#define DYNTODBL(IOP,FOP,IRCV)     IRCV = b_bcid(IOP,FOP);\
                                   if (IRCV != ROUND) rc += IRCV
#define NEW_(IRES)                 IRES = Lvarget();  poolvars++


/*------------*
 | Compares : |
 *------------*/

#define GT_(IOP1,IOP2)   (b_bcmp(IOP1,IOP2) > 0)
#define GE_(IOP1,IOP2)   (b_bcmp(IOP1,IOP2) >= 0)
#define LT_(IOP1,IOP2)   (b_bcmp(IOP1,IOP2) < 0)
#define LE_(IOP1,IOP2)   (b_bcmp(IOP1,IOP2) <= 0)
#define EQ_(IOP1,IOP2)   (b_bcmp(IOP1,IOP2) == 0)
#define NEQ_(IOP1,IOP2)  (b_bcmp(IOP1,IOP2) != 0)

#define GT_ABS_(IOP1,IOP2)   (b_bacm(IOP1,IOP2) > 0)
#define GE_ABS_(IOP1,IOP2)   (b_bacm(IOP1,IOP2) >= 0)
#define LT_ABS_(IOP1,IOP2)   (b_bacm(IOP1,IOP2) < 0)
#define LE_ABS_(IOP1,IOP2)   (b_bacm(IOP1,IOP2) <= 0)
#define EQ_ABS_(IOP1,IOP2)   (b_bacm(IOP1,IOP2) == 0)
#define NEQ_ABS_(IOP1,IOP2)  (b_bacm(IOP1,IOP2) != 0)


/*----------------------*
 | Standard Functions : |
 *----------------------*/

#define SQRT_(IARG)    if ((rc = sqrtve(IARG)  ) != 0) \
                          ERREXIT(0,rc,poolvars);
#define EXP_(IARG)     if ((rc = expve(IARG)   ) != 0) \
                          ERREXIT(0,rc,poolvars);
#define LN_(IARG)      if ((rc = lnve(IARG)    ) != 0) \
                          ERREXIT(0,rc,poolvars);
#define LN_T_(IARG)    if ((rc = lnvea(IARG)   ) != 0) \
                          ERREXIT(0,rc,poolvars);
#define SICO_(IARG)    if ((rc = sicovea(IARG) ) != 0) \
                          ERREXIT(0,rc,poolvars);
#define ARCSIN_(IARG)  if ((rc = arcsinve(IARG,&LPiOv2)) != 0) \
                          ERREXIT(0,rc,poolvars);
#define ARCTAN_(IARG)  if ((rc = arctanve(IARG,&LPiOv2)) != 0) \
                          ERREXIT(0,rc,poolvars);
#define SINH_(IARG)    if ((rc = sinhvea(IARG) ) != 0) \
                          ERREXIT(0,rc,poolvars);
#define ASSIGN_(IRES)  if ((rc = Lassign(IRES) ) != 0) \
                          ERREXIT(0,rc,poolvars);

#define EXP_ABS_(IARG) { sgn = IARG->s; rc = expve(IARG); IARG->s = sgn; \
                         if (rc != 0) ERREXIT(0,rc,poolvars); }



/*-------------*
 | Precision : |
 *-------------*/

#define PRECISION_(L)  Maxl = L

/*-------------------*
 | Improved Accuracy |
 *-------------------*/

#define INT_HPREC       1

#if INT_HPREC
                                /* test result of standard function f(x) */
                                /* in ASSIGN_() to satisfy restrictions  */
#define LESS_ABS_ARG      1     /* |f(x)| <= |x| */
#define GREATER_ABS_ARG   2     /* |f(x)| >= |x| */
#define LESS_ABS_ONE      4     /* |f(x)| <= |1| */
#define GREATER_ABS_ONE   8     /* |f(x)| >= |1| */
extern int b_case;
extern dynamic *b_farg;
#endif





