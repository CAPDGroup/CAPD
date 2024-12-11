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

/* CVS $Id: b_defs.h,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_defs.h                              */
/*                                                              */
/****************************************************************/

/*    Number of bits generated in one cycle of the decimal-to-  */
/*    binary conversion.                                        */
/*    The product 'BIN_PACK_POWER*DEC_PACK_POWER' must be       */
/*    less than '2**B_LENGTH'.                                  */
#define BIN_PACK                8
#define BIN_PACK_POWER          256

/*    Number of bits forming a 'char' value.                    */
#define BITS_PER_CHAR           8
#define LOG_BITS_PER_CHAR       3
#define MOD_BITS_PER_CHAR       0x7

/*    Number of bits forming a 'a_btyp' value.                  */
/*    This number must be 32.                                   */
#define B_LENGTH                32

/*    Number of decimals generated in one cycle of the binary-  */
/*    to-decimal binary conversion.                             */
#define B2D_LOG10               4
#define B2D_POWER               10000

/*    Number of decimals packed for decimal-to-binary           */
/*    conversion.                                               */
/*    This number must be greater than sizeof(a_btyp).          */
/*    The product 'BIN_PACK_POWER*DEC_PACK_POWER' must be       */
/*    less than '2**B_LENGTH'.                                  */
#define DEC_PACK                7
#if SHORTABTYP
#define DEC_PACK_POWER          10000000
#else
#define DEC_PACK_POWER          10000000L
#endif

/*    Size of character array reserved for SYSTEM_DIRECTORY     */
/*    and USER_DIRECTORY including final '\0'.                  */
#define DIRNAME_LENGTH          64

/*    Number of decimal digits in exponent.                     */
#define ExpDigits               3

/*    Number of program parameters in PASCAL-XSC program        */
#define FPPLT                   16

/*    Number of bits forming an 'a_intg' value.                 */
/*    This number must be less than or equal to B_LENGTH.       */
#define I_LENGTH                C_P_8

/*    Log base 2 of B_LENGTH.                                   */
#define LOG_B_LENGTH            5

/*    MAXINT maximum signed integer value                       */
#if I_LENGTH-32
#if I_LENGTH-16
#define MAXINT                  ((a_intg)0)
#else
#define MAXINT                  ((a_intg)0x7FFF)
#endif
#else
#if SHORTABTYP
#define MAXINT                  ((a_intg)0x7FFFFFFF)
#else
#define MAXINT                  ((a_intg)0x7FFFFFFFL)
#endif
#endif

/*    SQRT_MAXINT square root of maximum signed integer value   */
#if I_LENGTH-32
#if I_LENGTH-16
#define SQRT_MAXINT             ((a_intg)0)
#else
#define SQRT_MAXINT             ((a_intg)0x00b5)
#endif
#else
#if SHORTABTYP
#define SQRT_MAXINT             ((a_intg)0x0000b504)
#else
#define SQRT_MAXINT             ((a_intg)0x0000b504L)
#endif
#endif

/*    MININT minimum signed integer value                       */
#if I_LENGTH-32
#if I_LENGTH-16
#define MININT                  ((a_intg)0)
#else
#define MININT                  ((a_intg)0x8000)
#endif
#else
#if SHORTABTYP
#define MININT                  ((a_intg)0x80000000)
#else
#define MININT                  ((a_intg)0x80000000L)
#endif
#endif

/*------ constants of type 'a_bool' ----------------------------*/

/*    TRUE must be one                                          */
#ifndef TRUE
#define TRUE                    ((a_bool)1)
#endif

/*    FALSE must be zero                                        */
#ifndef FALSE
#define FALSE                   ((a_bool)0)
#endif

/*------ constants of type 'a_btyp' ----------------------------*/

/*    Mask for upper half of a_btyp value.                      */
#if SHORTABTYP
#define HIGH_MASK               ((a_btyp)0xFFFF0000)

/*    Mask for lower half of a_btyp value.                      */
#define LOW_MASK                ((a_btyp)0x0000FFFF)

/*    Mask for least significant bit position of a_btyp value.  */
#define LSB                     ((a_btyp)0x00000001)

/*    Maximum possible a_btyp value.                            */
#define MAX_BASETYPE            ((a_btyp)0xFFFFFFFF)

/*    Mask for most significant bit position of a_btyp value.   */
#define MSB                     ((a_btyp)0x80000000)

/*    Inverse mask for most significant bit position of         */
/*    a_btyp value.                                             */
#define NOT_MSB                 ((a_btyp)0x7FFFFFFF)

/*    Minimum possible a_btyp value.                            */
#define ZERO                    ((a_btyp)0x00000000)
#else
#define HIGH_MASK               ((a_btyp)0xFFFF0000L)

/*    Mask for lower half of a_btyp value.                      */
#define LOW_MASK                ((a_btyp)0x0000FFFFL)

/*    Mask for least significant bit position of a_btyp value.  */
#define LSB                     ((a_btyp)0x00000001L)

/*    Maximum possible a_btyp value.                            */
#define MAX_BASETYPE            ((a_btyp)0xFFFFFFFFL)

/*    Mask for most significant bit position of a_btyp value.   */
#define MSB                     ((a_btyp)0x80000000L)

/*    Inverse mask for most significant bit position of         */
/*    a_btyp value.                                             */
#define NOT_MSB                 ((a_btyp)0x7FFFFFFFL)

/*    Minimum possible a_btyp value.                            */
#define ZERO                    ((a_btyp)0x00000000L)
#endif

/*------ constants of type 'char' ------------------------------*/

/*    end of line character                                     */
#define EOLN                    '\n'

/*    tabulator character                                       */
#define TAB                     '\t'

/* -------- Macro definitions --------------------------------- */

/*    Free allocated 'a_btyp' array.                            */
#define B_FREE(x)               free((char *)(x));

/*    Get upper half of 'a_btyp' value moved to lower half.     */
#define GETHIGH(x)              (((x) >> (B_LENGTH/2)) & LOW_MASK)

/*    Get lower half of 'a_btyp' value.                         */
#define GETLOW(x)               ((x) & LOW_MASK)

/*    Move lower half to upper half position.                   */
#define MOVEHIGH(x)             (((x) << (B_LENGTH/2)) & HIGH_MASK)

/*    Logical negation                                          */
#define NOT(b)                  (a_bool)(TRUE-(b))

/*    Constants for tenbyte arithmetic                          */
#define tCHARAC                 ((int)0x3fff)
#define tEXPO_ADJUST            ((int)0x3fff)
#define tEXPO_MAX               (((int)0x7ffe)-(tCHARAC))
#define tEXPO_MIN               (1-(tCHARAC))
#if SHORTABTYP
#define tFIRST_BIT              ((a_btyp)0x00800000)
#define tMASK                   ((a_btyp)0x000000FF)
#else
#define tFIRST_BIT              ((a_btyp)0x00800000L)
#define tMASK                   ((a_btyp)0x000000FFL)
#endif
#define tMANTL                  64
#define tMSB                    ((unsigned char)0x80)
#define tSHIFT                  BITS_PER_CHAR





