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

/* CVS $Id: o_slct.h,v 1.28 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_slct.h                              */
/*                                                              */
/*      Description     : selecting C compiler                  */
/*                        defining string constants             */
/*                                                              */
/****************************************************************/

#include "cxscconf.h"

#if SUN4_FORTE
#undef SUN4_GNU_C
#define SUN4_GNU_C 1  
#endif

#define ATARI_TURBO_C   0
#define CHAM_32_C  0
#define CHAM_64_C  0
#define CONVEX1_UNIX_C  0
#define CONVEX2_UNIX_C  0
#define CRAY_UNIX_C     0
#define DEC_ULTRIX_C    0
#define GNU_C           0
#define HP_9000_C       0
#define IBM_3090_C      0
#define IBM_370_C       0
#define IBM_AT_BOR_C    0
#define IBM_AT_MS_C     0
#define IBM_AT_TURBO_C  0
#define IBM_EMX_C       0
#define IBM_OS2_ICC_C   0
#define IBM_PS2_C       0
#define IBM_RS6000_C    0
#define IBM_RT_C        0
#define IBM_WATCOM_C   0
#define SILICON_UNIX_C  0
#define SUN4_CPP_C      0
#define SUN4_OS5_GNU_C  0
#define T800_HELIOS_C   0
#define VAX_UNIX_C      0
#define VAX_VMS_C       0
#define VP_EX_C         0
#define ZORTECH_C       0

/*--------------------------------------------------------------*/
/* Installation Macros :                                        */
/* At most one of the following macros must be non-zero.        */
/*--------------------------------------------------------------*/
#ifdef DUMMY

                                /* TURBO C compiler for ATARI   */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #define NODSTD             */
                                /*   A_real is structure        */
#define ATARI_TURBO_C   0
                                /* Chameleon 32 Bit C for sim.  */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #define NODSTD             */
#define CHAM_32_C  0
                                /* Chameleon 64 Bit C for sim.  */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #define NODSTD             */
#define CHAM_64_C  0
                                /* CONVEX compiler for CONVEX C1*/
                                /*   IEEE_HARDWARE not supported*/
#define CONVEX1_UNIX_C  0
                                /* CONVEX compiler for CONVEX C2*/
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #define NODSTD             */
#define CONVEX2_UNIX_C  0
                                /* CRAY compiler for CRAY UNIX  */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #define NODSTD             */
                                /*   A_real is structure        */
#define CRAY_UNIX_C     0
                                /* OSF/1 compiler for DEC Alpha */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #define Using_lint_args    */
                                /*   #define INTEL      1       */
                                /*   #define Using 64 bit       */
                                /*   #define NODSTD             */
#define DEC_ALPHA_C    0
                                /* ULTRIX compiler for DEC      */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #define IBMva              */
                                /*   #undef  LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define DEC_ULTRIX_C    0
                                /* GNU_C compiler for DOS system*/
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define GNU_C           0
                                /* HP C compiler for HP9000     */
                                /*   IEEE_HARDWARE supported    */
                                /*   #define IBMva              */
                                /*   #undef  LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef  NODSTD             */
#define HP_9000_C       0
                                /* IBM C compiler for IBM 3090  */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #define NODSTD             */
#define IBM_3090_C      0
                                /* IBM C compiler for System/370*/
                                /*   IEEE_HARDWARE not supported*/
#define IBM_370_C       0
                                /* BORLAND compiler for IBM AT  */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #define IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define IBM_AT_BOR_C    0
                                /* MS C compiler for IBM AT     */
                                /*   IEEE_HARDWARE not supported*/
#define IBM_AT_MS_C     0
                                /* TURBO C compiler for IBM AT  */
                                /*   IEEE_HARDWARE not supported*/
#define IBM_AT_TURBO_C  0
                                /* EMX compiler for DOS/OS2     */
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef IBMva               */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef  NODSTD             */
#define IBM_EMX_C       0
                                /* GnuC compiler for LINUX/PC   */
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef IBMva               */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define IBM_LINUX_C       0
                                /* ICC compiler for IBM AT OS2  */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef IBMva               */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define IBM_OS2_ICC_C   0
                                /* IBM C compiler for IBM PS2   */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #define IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define IBM_PS2_C       0
                                /* IBM C compiler for IBM RS6000*/
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef  NODSTD             */
#define IBM_RS6000_C    0
                                /* IBM C compiler for IBM RT    */
                                /*   IEEE_HARDWARE supported    */
                                /*   #define IBMva              */
                                /*   #undef  LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef  NODSTD             */
#define IBM_RT_C        0
                                /* Watcom C++ V10.0 compiler    */
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef IBMva               */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define IBM_WATCOM_C   0
                                /* IRIX compiler for SG R4000   */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef NODSTD              */
#define SILICON_UNIX_C  0
                                /* OS4 C++ compiler for SUN4    */
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef NODSTD              */
#define SUN4_CPP_C      0
                                /* GNU C compiler for SUN4      */
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef IBMva               */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #undef NODSTD              */
#define SUN4_GNU_C      0
                                /* OS4 C compiler for SUN4      */
                                /*   IEEE_HARDWARE supported    */
                                /*   #define IBMva              */
                                /*   #undef LINT_ARGS           */
                                /*   #define INTEL      0       */
                                /*   #undef NODSTD              */
#define SUN4_OS4_C      0
                                /* Solaris 3.3 GnuC Compiler    */
						  /* Sun Super Sparc              */
                                /*   IEEE_HARDWARE supported    */
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #define NODSTD             */
#define SUN4_OS5_GNU_C  0
                                /* HELIOS C compiler for T800   */
                                /*   IEEE_HARDWARE supported    */
#define T800_HELIOS_C   0
                                /* UNIX C compiler for VAX      */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #define IBMva              */
                                /*   #undef LINT_ARGS           */
                                /*   #define INTEL      1       */
                                /*   #define NODSTD             */
                                /*   A_real is structure        */
#define VAX_UNIX_C      0
                                /* VMS C compiler for VAX       */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #define IBMva              */
                                /*   #undef LINT_ARGS           */
                                /*   #define INTEL      1       */
                                /*   #define NODSTD             */
                                /*   A_real is structure        */
#define VAX_VMS_C       0
                                /* C-EX compiler for VP600/20 UX*/
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef  IBMva              */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      0       */
                                /*   #define NODSTD             */
#define VP_EX_C         0
                                /* ZORTECH C++ compiler for DOS */
                                /*   IEEE_HARDWARE not supported*/
                                /*   #undef IBMva               */
                                /*   #define LINT_ARGS          */
                                /*   #define INTEL      1       */
                                /*   #undef  NODSTD             */
#define ZORTECH_C       0    

#endif

/* - - - - - - - - - - - - - - - - - - -*/
/* DEMO_VERSION                         */
/* Demo_version controls the pxsc demo  */
/* version. Non standard  I/O is not    */
/* supported. Product text is set       */
/* Files necessary to compile:          */
/* - f_assg.c, b_popt.c, p_init.c       */
/*   dont forget to set IEEE_HARDWARE   */
/*   correctly                          */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define DEMO_VERSION 
*/

/* - - - - - - - - - - - - - - - - - - -*/
/* TEST_VERSION                         */
/* Test_version marks version string    */
/* that only a test version is running  */
/* format 'T ??.??.?? date.             */
/* Set this flag only from commandline  */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define TEST_VERSION
*/
/* - - - - - - - - - - - - - - - - - - -*/
/* DEC_ARITH                            */
/* enables decimal version of pxsc.     */
/* Set this flag only from commandline  */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define DEC_ARITH   
*/
/* - - - - - - - - - - - - - - - - - - -*/
/* IEEE_HARDWARE                        */
/* If the selected machine/compiler/OS  */
/* supports the IEEE double data format */
/* and rounded arithmetic operations in */
/* hardware, then the definition of this*/
/* macro generates code for the support */
/* of hardware operations +,-,*,/ and   */
/* rounding modes nearest,up,down.      */
/* Requirements:                        */
/* - a_gets(),a_sets() must exist for   */
/*   the selected installation.         */
/* - a_real = C double = IEEE double    */
/*   (defined in p88rts.h)              */
/* - signaling NaN must be identified   */
/*   (defined in e_defs.h)              */
/* - INTEL flag must be set correctly   */
/*   (defined in o_syst.h)              */
/* - - - - - - - - - - - - - - - - - - -*/
#if IBM_RT_C
/*
#define IEEE_HARDWARE
*/
#endif

/*--------------------------------------------------------------*/
/* Flags that may be set according to the                       */
/* specification of the runtime environment.                    */
/*--------------------------------------------------------------*/

/* - - - - - - - - - - - - - - - - - - -*/
/* Runtime version text and copyright   */
/* statement.                           */
/* The following sequence of text is    */
/* displayed on runtime option -vn:     */
/*   (1) PRODUCT_TEXT                   */
/*   (2) VERSION_TEXT                   */
/*   (3) VERSION_SYS_TEXT               */
/*   (4) VERSION_EXT_TEXT               */
/*   (5) VERSION_ARI_TEXT               */
/*   (6) COPYRIGHT_TEXT                 */
/* - - - - - - - - - - - - - - - - - - -*/
#ifdef DEMO_VERSION
#define PRODUCT_TEXT   "PASCAL-XSC Demonstration, "
#else
#define PRODUCT_TEXT   "PASCAL-XSC Runtime "
#endif
#define VERSION_TEXT    "Version 1997-09-26"

#if ATARI_TURBO_C
#define VERSION_SYS_TEXT "(Atari)"
#else
#if CHAM_32_C
#define VERSION_SYS_TEXT "(Chameleon 32 C Compiler)"
#else
#if CHAM_64_C
#define VERSION_SYS_TEXT "(Chameleon 64 C Compiler)"
#else
#if CONVEX1_UNIX_C
#define VERSION_SYS_TEXT "(CONVEX C1)"
#else
#if CONVEX2_UNIX_C
#define VERSION_SYS_TEXT "(CONVEX C2)"
#else
#if CRAY_UNIX_C
#define VERSION_SYS_TEXT "(CRAY UNIX)"
#else
#if DEC_ALPHA_C
#define VERSION_SYS_TEXT "(DEC ALPHA)"
#else
#if DEC_ULTRIX_C
#define VERSION_SYS_TEXT "(DEC ULTRIX)"
#else
#if HP_9000_C
#define VERSION_SYS_TEXT "(HP 9000/xxx)"
#else
#if IBM_3090_C
#define VERSION_SYS_TEXT "(IBM 3090)"
#else
#if IBM_370_C
#define VERSION_SYS_TEXT "(IBM 370)"
#else
#if IBM_AT_BOR_C
#define VERSION_SYS_TEXT "(IBM PC BORLAND C)"
#else
#if IBM_AT_MS_C
#define VERSION_SYS_TEXT "(IBM PC MSC)"
#else
#if IBM_AT_TURBO_C
#define VERSION_SYS_TEXT "(IBM PC TURBO C)"
#else
#if IBM_EMX_C
#define VERSION_SYS_TEXT "(IBM EMX GNU C)"
#else
#if IBM_LINUX_C 
#define VERSION_SYS_TEXT "(PC LINUX GNU C)"
#else
#if IBM_OS2_ICC_C 
#define VERSION_SYS_TEXT "(IBM OS2 ICC C)"
#else
#if IBM_PS2_C
#define VERSION_SYS_TEXT "(IBM PS2)"
#else
#if IBM_RS6000_C
#define VERSION_SYS_TEXT "(IBM RS6000)"
#else
#if IBM_RT_C
#define VERSION_SYS_TEXT "(IBM 6150)"
#else
#if IBM_WATCOM_C
#define VERSION_SYS_TEXT "(IBM WATCOM C++ V10.0)"
#else
#if SILICON_UNIX_C
#define VERSION_SYS_TEXT "(Silicon Graphics IRIX C)"
#else
#if SUN4_CPP_C
#define VERSION_SYS_TEXT "(Sun4 C++)"
#else
#if SUN4_GNU_C
#define VERSION_SYS_TEXT "(Sun4 GNU C)"
#else
#if SUN4_OS4_C
#define VERSION_SYS_TEXT "(Sun4)"
#else
#if SUN4_OS5_GNU_C
#define VERSION_SYS_TEXT "(Sun4 Solaris 3.3 GNU C)"
#else
#if T800_HELIOS_C
#define VERSION_SYS_TEXT "(T800 HELIOS)"
#else
#if VAX_UNIX_C
#define VERSION_SYS_TEXT "(ULTRIX)"
#else
#if VAX_VMS_C
#define VERSION_SYS_TEXT "(VMS)"
#else
#if VP_EX_C
#define VERSION_SYS_TEXT "(VP 600/20 UX)"
#else
#if ZORTECH_C
#define VERSION_SYS_TEXT "(IBM PC ZORTECH C)"
#else
#define VERSION_SYS_TEXT " "
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#ifdef TEST_VERSION
#ifdef IEEE_HARDWARE
#define VERSION_EXT_TEXT	    "HWS T 1997-09-26"
#else
#define VERSION_EXT_TEXT	    "SWE T 1997-09-26"
#endif
#else
#ifdef IEEE_HARDWARE
#define VERSION_EXT_TEXT	    "HWS"
#else
#define VERSION_EXT_TEXT	    "SWE"
#endif
#endif

#ifdef DEC_ARITH
#define VERSION_ARI_TEXT "DEC"
#else
#define VERSION_ARI_TEXT "BIN"
#endif

#define COPYRIGHT_TEXT "(C) Copyright University of Karlsruhe 1991-1997\n"

/* - - - - - - - - - - - - - - - - - - -*/
/* Runtime messages are loaded from     */
/* MSG_FILE if this macro is defined.   */
/* Otherwise file "o_msg1.h" is included*/
/* by the runtime system.               */
/* - - - - - - - - - - - - - - - - - - -*/
#define MSG_FILE_ENABLED

/* - - - - - - - - - - - - - - - - - - -*/
/* No function names are pushed to the  */
/* trace back stack for all the runtime */
/* routines except for standard         */
/* functions if this macro name is      */
/* defined.                             */
/* - - - - - - - - - - - - - - - - - - -*/
#define RTSTRC_DISABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* No function names of standard        */
/* functions are pushed to the trace    */
/* back stack if this macro is defined. */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define STDTRC_DISABLE
*/

/* - - - - - - - - - - - - - - - - - - -*/
/* This macro enables the generation    */
/* of code for the command line option  */
/* PPOPT_COPYIB if it is defined.       */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define COPYIB_ENABLE
*/

/* - - - - - - - - - - - - - - - - - - -*/
/* This macro enables the generation    */
/* of code for the command line option  */
/* PPOPT_OPTION if it is defined.       */
/* - - - - - - - - - - - - - - - - - - -*/
#define OPTION_ENABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* Enables the generation of code for   */
/* the runtime option which toggles     */
/* normalized/denormalized generation   */
/* as successor/predecessor of 0.0 if   */
/* defined.                             */
/* - - - - - - - - - - - - - - - - - - -*/
#define NORMALIZE_ENABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* This macro enables the generation    */
/* of code for the command line options */
/* PPOPT_SYSDIR and PPOPT_USRDIR if     */
/* defined.                             */
/* - - - - - - - - - - - - - - - - - - -*/
#define DIRECT_ENABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* This macro enables the generation of */
/* code for the command line option     */
/* PPOPT_SIGNED if it is defined.       */
/* - - - - - - - - - - - - - - - - - - -*/
#define SIGNED_ENABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* This macro enables the generation of */
/* code for the command line option     */
/* PPOPT_TRACEC if it is defined.       */
/* - - - - - - - - - - - - - - - - - - -*/
#define TRACEC_ENABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* This macro enables the generation of */
/* code for the command line option     */
/* PPOPT_STATTB if it is defined.       */
/* - - - - - - - - - - - - - - - - - - -*/
#define STATTB_ENABLE

/* - - - - - - - - - - - - - - - - - - -*/
/* Output of decimal conversion of      */
/* a_real value in trace back parameter */
/* printing using r_writ() is done if   */
/* this macro is defined.               */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define TRAP_REAL
*/

/* - - - - - - - - - - - - - - - - - - -*/
/* Check for text file in binary read   */
/* and write routine. This status is    */
/* already checked by the PASCAL-XSC    */
/* compiler. Do not define this macro   */
/* in PASCAL-XSC context.               */
/* - - - - - - - - - - - - - - - - - - -*/
/*
#define CHECK_TEXT_FILE
*/

/*--------------------------------------------------------------*/
/* Boolean constants that may be set according to the           */
/* specification of the runtime environment.                    */
/*--------------------------------------------------------------*/

/* - - - - - - - - - - - - - - - - - - -*/
/* Continue after call to trap handler  */
/* (default value)                      */
/* FALSE = abort processing             */
/* TRUE  = continue processing          */
/* - - - - - - - - - - - - - - - - - - -*/
#define EXIT_CONT              FALSE


/*--------------------------------------------------------------*/
/* Integer constants that may be set according to the           */
/* specification of the runtime environment.                    */
/*--------------------------------------------------------------*/

/* - - - - - - - - - - - - - - - - - - -*/
/* Value returned to runtime environ-   */
/* ment when processing is aborted by   */
/* trap handler (default value).        */
/* Value must be greater than or equal  */
/* to 0 and less than or equal to 31.   */
/* - - - - - - - - - - - - - - - - - - -*/
#define EXIT_VALUE              1





