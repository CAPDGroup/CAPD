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

/* CVS $Id: e_defs.h,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_defs.h                              */
/*                                                              */
/****************************************************************/

/*    Macro for testing for a signaling NaN                     */
#ifdef IEEE_HARDWARE
#if IBM_RT_C+IBM_RS6000_C+SUN4_OS4_C+SUN4_GNU_C+SUN4_CPP_C+SUN4_OS5_GNU_C
#define SIGNAL_BIT              ((a_btyp)0x00080000L)
#define SIGNALING(a)            (!((a) & SIGNAL_BIT))
#define SET_QUIET(a)            ((a) | SIGNAL_BIT)
#define SET_SIGNAL(a)		((a) & ~SIGNAL_BIT)
#endif
#if HP_9000_C
#define SIGNAL_BIT              ((a_btyp)0x00080000L)
#define SIGNALING(a)            ((a) & SIGNAL_BIT)
#define SET_QUIET(a)            ((a) & ~SIGNAL_BIT)
#define SET_SIGNAL(a)		((a) | SIGNAL_BIT)
#endif
#if IBM_LINUX_C+IBM_EMX_C
#define SIGNAL_BIT              ((a_btyp)0x00080000L)
#define SIGNALING(a)            (!((a) & SIGNAL_BIT))
#define SET_QUIET(a)            ((a) | SIGNAL_BIT)
#define SET_SIGNAL(a)		((a) & ~SIGNAL_BIT)
#endif
#else	/* IEEE_HARDWARE */
#if SHORTABTYP
#define SIGNAL_BIT              ((a_btyp)0x00080000)
#else
#define SIGNAL_BIT              ((a_btyp)0x00080000L)
#endif
#define SIGNALING(a)            ((a) & SIGNAL_BIT)
#define SET_QUIET(a)            ((a) & ~SIGNAL_BIT)
#define SET_SIGNAL(a)		((a) | SIGNAL_BIT)
#endif

/*    Macro for testing for a mantissa used for infinity        */
#define MANT_INFINITY(a)        ((a)[0]==HIDDEN_BIT && (a)[1]==ZERO)

/*    Constants for trap handling routines.                     */
#if SHORTABTYP
#define E_CHNG                  ((a_btyp)0x00004000)
#define E_POPP                  ((a_btyp)0x00002000)
#define E_PUSH                  ((a_btyp)0x00001000)

#define E_ACTIVE                ((a_btyp)0x00000100)
#define E_BACK                  ((a_btyp)0x00000080)
#define E_CONT                  ((a_btyp)0x00000040)
#define E_DEFAULT               ((a_btyp)0x00000020)
#else
#define E_CHNG                  ((a_btyp)0x00004000L)
#define E_POPP                  ((a_btyp)0x00002000L)
#define E_PUSH                  ((a_btyp)0x00001000L)

#define E_ACTIVE                ((a_btyp)0x00000100L)
#define E_BACK                  ((a_btyp)0x00000080L)
#define E_CONT                  ((a_btyp)0x00000040L)
#define E_DEFAULT               ((a_btyp)0x00000020L)
#endif

/*    Value for non-existing exception handler.                 */
#ifdef LINT_ARGS
#define NO_TRAP                 ((void (*)(a_btyp,int,va_list))0)
#else
#define NO_TRAP                 ((void (*)())0)
#endif

/*------ Error characterizations -------------------------------*/
/*    Splitting of error condition code into error code and     */
/*    error characterization.                                   */
#if SHORTABTYP
#define E_MASK                  ((a_btyp)0xFFFFFF00)

#define E_ECNT                  ((a_btyp)0x00000080)
#define E_EMSG                  ((a_btyp)0x00000040)
#define E_ETBC                  ((a_btyp)0x00000020)
#define E_EXIT                  ((a_btyp)0x00000010)
#define E_EARG                  ((a_btyp)0x00000008)
#define E_EQIE                  ((a_btyp)0x00000004)
#define E_IEEE                  ((a_btyp)0x00000001)

/*------ Error groups ------------------------------------------*/
/*    Do not change the value order of error groups.            */
/*    Changing value order has side effect on file "e_data.c".  */
#define NO_ERROR                ((a_btyp)0x00000000)
#define INV_OP                  ((a_btyp)0x00000100)
#define DIV_BY_ZERO             ((a_btyp)0x00000a00)
#define OVERFLOW                ((a_btyp)0x00000b00)
#define UNDERFLOW               ((a_btyp)0x00000c00)
#define INEXACT                 ((a_btyp)0x00000d00)
#define ALLOCATION              ((a_btyp)0x00000e00)
#define I_O_ERROR               ((a_btyp)0x00001000)
#define I_O_BUFFER              ((a_btyp)0x00001100)
#define INV_ARG                 ((a_btyp)0x00001200)
#define INDEX_RANGE             ((a_btyp)0x00001300)
#else
#define E_MASK                  ((a_btyp)0xFFFFFF00L)

#define E_ECNT                  ((a_btyp)0x00000080L)
#define E_EMSG                  ((a_btyp)0x00000040L)
#define E_ETBC                  ((a_btyp)0x00000020L)
#define E_EXIT                  ((a_btyp)0x00000010L)
#define E_EARG                  ((a_btyp)0x00000008L)
#define E_EQIE                  ((a_btyp)0x00000004L)
#define E_IEEE                  ((a_btyp)0x00000001L)

/*------ Error groups ------------------------------------------*/
/*    Do not change the value order of error groups.            */
/*    Changing value order has side effect on file "e_data.c".  */
#define NO_ERROR                ((a_btyp)0x00000000L)
#define INV_OP                  ((a_btyp)0x00000100L)
#define DIV_BY_ZERO             ((a_btyp)0x00000a00L)
#define OVERFLOW                ((a_btyp)0x00000b00L)
#define UNDERFLOW               ((a_btyp)0x00000c00L)
#define INEXACT                 ((a_btyp)0x00000d00L)
#define ALLOCATION              ((a_btyp)0x00000e00L)
#define I_O_ERROR               ((a_btyp)0x00001000L)
#define I_O_BUFFER              ((a_btyp)0x00001100L)
#define INV_ARG                 ((a_btyp)0x00001200L)
#define INDEX_RANGE             ((a_btyp)0x00001300L)
#endif

/*------ Data type characterization ----------------------------*/
#define E_TCHR                  ((int)1)        /* char *            */
#define E_TDBL                  ((int)2)        /* double *          */
#define E_TDTP                  ((int)3)        /* dotprecision *    */
#define E_TINT                  ((int)4)        /* int *             */
#define E_TMLT                  ((int)5)        /* multiprecision *  */
#define E_TSTR                  ((int)6)        /* char *            */
#define E_TULT                  ((int)7)        /* ultraprecision *  */
#define E_TSTG                  ((int)8)        /* s_trng *          */
#define E_TLNG                  ((int)9)        /* a_long *          */
#define E_TRES                  ((int)0x100)    /* --- result        */
#define E_TMSG                  ((int)0x7e00)   /* --- message       */
                                /* identification of displayed type  */
#define E_TEXT(a)               ((a)*(2*E_TRES))

/*------ Error codes of float arithmetic -----------------------*/

#define DENOR 1         /* denormalized number converted                */
#define MINFI 2         /* -infinity detected                           */
#define NANDE 3         /* NAN detected                                 */
#define OFLOW 4         /* exponent overflow                            */
#define PINFI 5         /* +infinity detected                           */
#define ROUND 6         /* double value is rounded                      */
#define UFLOW 7         /* exponent underflow                           */
#define ZEROD 8         /* divide by zero                               */
#define RANGE 9         /* range error                                  */
#define ALLOC 10        /* allocation error                             */
#define NALLO 11        /* data not allocated                           */

/*------ Macros specifying push and pop in runtime -------------*/

#ifdef RTSTRC_DISABLE
#define E_TPUSH(a)
#define E_TPOPP(a)
#else
#define E_TPUSH(a) \
 {extern char *o_text[]; e_push((char *)a,(char *)o_text[6]);}
#define E_TPOPP(a)      e_popp();
#endif
#ifdef STDTRC_DISABLE
#define E_SPUSH(a)      a_intg flags; e_save(&flags);
#define E_SPOPP(a)      e_rest(flags);
#else
#define E_SPUSH(a) \
 a_intg flags; \
 {extern char *o_text[]; e_push((char *)a,(char *)o_text[6]);} \
 e_save(&flags);
#define E_SPOPP(a)      e_rest(flags); e_popp();
#endif





