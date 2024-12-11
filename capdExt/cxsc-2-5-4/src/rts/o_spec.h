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

/* CVS $Id: o_spec.h,v 1.22 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_spec.h                              */
/*                                                              */
/****************************************************************/

#if CONVEX1_UNIX_C+IBM_RT_C+IBM_RS6000_C+HP_9000_C+SUN4_OS4_C+VAX_UNIX_C+\
    DEC_ULTRIX_C
#ifdef LINT_ARGS
>>>>>> LINT_ARGS must not be defined in p88rts.h <<<<<<
#endif
#else
#ifndef LINT_ARGS
>>>>>> LINT_ARGS must be defined in p88rts.h <<<<<<
#endif
#endif

#ifdef LINT_ARGS

/* ANSI_C is used in mathematical standard functions using      */
/* tenbyte ops. (see t_mach.h).                                 */
/* It is assumed, that only C compilers with prototypes are     */
/* ANSI C compatible.                                           */
 
#define ANSI_C
#endif

/*--------------------------------------------------------------*/
/* If standard functions asin(), atan(), log(), and sqrt()      */
/* are available for the IEEE 64-bit double format, then        */
/* the NODSTD macro need not be defined.                        */
/*--------------------------------------------------------------*/

#if ATARI_TURBO_C+IBM_370_C+VAX_VMS_C+VAX_UNIX_C+CONVEX2_UNIX_C+\
    IBM_3090_C+CRAY_UNIX_C+VP_EX_C+DEC_ALPHA_C+SUN4_OS5_GNU_C+CHAM_32_C+\
    CHAM_64_C
#define NODSTD          /* No Double StanDardfunctions exist    */
#endif
#ifdef DEC_ARITH
#define NODSTD          /* No Double StanDardfunctions exist    */
#endif

/*--------------------------------------------------------------*/
/* If long is larger then 32 bit, a_btyp can be specified by a  */
/* 32 bit unsigned word. This is the case in 64 bit system.     */
/* Then this constant has to be set.                            */
/*--------------------------------------------------------------*/
#define SHORTABTYP DEC_ALPHA_C+CHAM_64_C+GNU_X86_64

/* -------- C-runtime/machine properties ---------------------- */
/* These switches must be used carefully according to the       */
/* specifications of your C-compiler or machine.                */

/*    Arithmetic shift right macro.                             */
#if CONVEX2_UNIX_C + IBM_3090_C+CRAY_UNIX_C+SILICON_UNIX_C+\
    VP_EX_C+DEC_ALPHA_C+SUN4_OS5_GNU_C
#define B_ASHR(x,y)             (((x)>=0) ? (x)>>(y) : (~(~(x)>>(y))))
#else
#define B_ASHR(x,y)             ((x)>>(y))
#endif

/*    Assign an a_real value to an a_real variable.             */
/*    It may be necessary in some installation to copy a_real   */
/*    values byte-wise to their destination address.            */
#define R_ASSIGN(d,s)           (d) = (s)

/* 1 -Addition of (unsigned) 'a_btyp' values never causes an    */
/*    overflow exception.                                       */
/*    Any possible carry from the most significant 'a_btyp'     */
/*    position is ignored.                                      */
#define C_P_1                   1

/* 1 -Subtraction of (unsigned) 'a_btyp' values never causes    */
/*    an exception.                                             */
/*    Subtraction is done by adding the two's complement.       */
/*    Any possible borrow from the most significant 'a_btyp'    */
/*    position is ignored.                                      */
#define C_P_2                   1

/*    Special ratio of data length.                             */
/* 0 -Loops are used in order to handle 'a_btyp' arrays which   */
/*    are equivalenced with 'a_real' values.                    */
/* 1 -The ratio 'sizeof(a_real)/sizeof(a_btyp)' equals 2.       */
/*    Explicit code is generated for the elements of the        */
/*    'a_btyp' arrays which are equivalenced with 'a_real'      */
/*    values.                                                   */
#define C_P_3                   1

/*    Data type void.                                           */
/*    Some old versions of C-compilers do not know the data     */
/*    type 'void' when used with function declarations.         */
/* 0 -If 'void' is not available, then 'int' is used.           */
/* 1 -'void' is available.                                      */
#define C_P_4                   1

/*    32-bit unsigned integer data type.                        */
/* 0 -a 32-bit unsigned integer data type is used.              */
#define C_P_5                   0

/*    Compiler supports IEEE double data type and operations    */
/*    according to the IEEE standard.                           */
/* 0 -software simulation of IEEE operations are used.          */
#define C_P_6                   0

/*    Compiler supports variable argument lists.                */
/* 0 -IBM parameter passing conventions                         */
/* 1 -Standard parameter passing conventions                    */
#define DUMMY_IBMVA CONVEX1_UNIX_C+IBM_PS2_C+IBM_RT_C+IBM_RS6000_C+HP_9000_C
#if DUMMY_IBMVA+SUN4_OS4_C+VAX_UNIX_C+DEC_ULTRIX_C
#ifdef IBMva
#define C_P_7                   0
#else
>>>>>> IBMva must be defined in p88rts.h <<<<<<
#endif
#else
#ifdef IBMva
>>>>>> IBMva must not be defined in p88rts.h <<<<<<
#else
#define C_P_7                   1
#endif
#endif

/*    Number of bits forming a value of type a_intg.            */
#define C_P_8                   32

/* -------- Compiler Dependent Macros ------------------------- */

/*    Overloading of data type "void" if not supplied by        */
/*    compiler.                                                 */
#if C_P_4
#else
#define void                    int
#endif

/*    The first argument of an actual variable length argument  */
/*    list must be of type "int" and represents the actual      */
/*    number of supplied arguments.                             */
/*    e_args = last formal argument if variable length argument */
/*             list is used.                                    */
/*    e_list    = declaration of last formal argument.          */
/*    e_argc    = last named formal argument. (control variable)*/
/*    e_argv    = reference variable to argument list.          */
/*    e_argr    = value of reference variable for restore.      */
/*    e_open    = open the block which allows argument access.  */
/*    e_sarg    = save reference variable to argument list.     */
/*    e_rarg    = restore refernence variable to argument list. */
/*    e_ref(x)  = value of argument access of type "x".         */
/*    e_shut    = shut the block which allows argument access.  */
#ifdef LINT_ARGS
#if C_P_7
#define e_args          ...                     /* argument list        */
#define e_list                                  /* list declaration     */
#define e_open(name)    va_start(e_argv,name)   /* open access          */
#else
#define e_args          va_alist                /* argument list        */
#define e_list          va_dcl                  /* list declaration     */
#define e_open(name)    va_start(e_argv)        /* open access          */
#endif
#else
#define e_args          va_alist                /* argument list        */
#define e_list          va_dcl                  /* list declaration     */
#define e_open(name)    va_start(e_argv)        /* open access          */
#endif
#define e_ref(x)        va_arg(e_argv,x)        /* get argument         */
#define e_shut          va_end(e_argv)          /* close access         */
#define e_sarg          va_list e_argr=e_argv;  /* save  access         */
#define e_rarg          e_argv = e_argr         /* restore access       */





