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

/* CVS $Id: r_defs.h,v 1.22 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_defs.h                              */
/*                                                              */
/*                   MINUS_SIGN / PLUS_SIGN                     */
/****************************************************************/

/*------ IEEE constants ----------------------------------------*/
/*    These constants are used in connection with the IEEE      */
/*    "double" data type.                                       */

/*    Buffer length used for "double" operations.               */
#define BSIZE                   (2*D_U_RATIO+1)

/*    Characteristic of IEEE value.                             */
#define CHARAC                  1023

/*    Character used as decimal point.                          */
#define XSC_DECIMAL_POINT           '.'

/*    Number of bits forming an IEEE value.                     */
/*    (sizeof(a_real)*BITS_PER_CHAR)                            */
#define DOUBLE_LENGTH           64

/*    Number of 'a_btyp' values forming an IEEE value.          */
#define D_U_RATIO               (DOUBLE_LENGTH/B_LENGTH)

/*    IEEE value classifications.                               */
#define E_CLS0                  0       /* signaling NaN        */
#define E_CLS1                  1       /* quiet NaN            */
#define E_CLS2                  2       /* -infinity            */
#define E_CLS3                  3       /* - normalized         */
#define E_CLS4                  4       /* -denormalized        */
#define E_CLS5                  5       /* -zero                */
#define E_CLS6                  6       /* +zero                */
#define E_CLS7                  7       /* +denormalized        */
#define E_CLS8                  8       /* +normalized          */
#define E_CLS9                  9       /* +infinity            */

/*    Exponent adjust for a "double" value in case of overflow  */
/*    or underflow.                                             */
#define EXPO_ADJUST             1536

/*    Mask for extracting exponent of IEEE value.               */
#if SHORTABTYP
#define EXPO_MASK               ((a_btyp)0x7FF00000)
#else
#define EXPO_MASK               ((a_btyp)0x7FF00000L)
#endif

/*    Maximum exponent of a normalized "double" value.          */
#define EXPO_MAX                1023

/*    Minimum exponent of a normalized "double" value.          */
#define EXPO_MIN                (-1022)

/*    Number of bit positions needed to shift the exponent      */
/*    of an IEEE value within a 'a_btyp' value most right.      */
#define EXPO_SHIFT              (B_LENGTH-ZERO_BITS)

/*    Character used as exponent delimiter.                     */
#define EXPONENT_E              'E'

/*    Mask for hidden mantissa bit if IEEE value is equivalenced*/
/*    with an array of type 'a_btyp'.                           */
#if SHORTABTYP
#define HIDDEN_BIT              ((a_btyp)0x00100000)
#else
#define HIDDEN_BIT              ((a_btyp)0x00100000L)
#endif

/*    Mask for mantissa bits in most significant 'a_btyp'       */
/*    value if IEEE value is equivalenced with an array of type */
/*    'a_btyp'.                                                 */
#if SHORTABTYP
#define MANT_HIGH               ((a_btyp)0x000FFFFF)
#else
#define MANT_HIGH               ((a_btyp)0x000FFFFFL)
#endif

/*    Number of bits forming the mantissa of a "double" value.  */
#define MANTL                   53

/*    Minimum of B_LENGTH/2 and ZERO_BITS-1                     */
#define MIN_FAC                 11

/*    Character used as minus sign.                             */
#define MINUS_SIGN              '-'

/*    Inverse mask for hidden mantissa bit if IEEE value is     */
/*    equivalenced with an array of type 'a_btyp'.              */
#if SHORTABTYP
#define NOT_HIDDEN_BIT          ((a_btyp)0xFFEFFFFF)
#else
#define NOT_HIDDEN_BIT          ((a_btyp)0xFFEFFFFFL)
#endif

/*    Character used as plus sign.                              */
#define PLUS_SIGN               '+'

/*    Mask for bit positions which nether hold a hidden or      */
/*    visible mantissa bit if IEEE value is equivalenced with   */
/*    an array of type 'a_btyp'.                                */
#if SHORTABTYP
#define SHFT_MASK               ((a_btyp)0xFFE00000)
#else
#define SHFT_MASK               ((a_btyp)0xFFE00000L)
#endif

/*    Number of bits used for exponent (=11) and sign (=1) in   */
/*    an IEEE value.                                            */
#define ZERO_BITS               12





