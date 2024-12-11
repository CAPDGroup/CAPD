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

/* CVS $Id: b_64bt.h,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

#ifndef _BIT64_
#define _BIT64_
/****************************************************************/
/*                                                              */
/*      Filename        : b_defs.h                              */
/*                                                              */
/****************************************************************/

/*    Number of bits forming a 'a_btyp' value.                  */
/*    This number must be 32.                                   */
#define B_LENGTH_64                64
#define LOG_B_LENGTH_64            6

/*    Number of bits forming an 'a_intg' value.                 */
/*    This number must be less than or equal to B_LENGTH.       */
#define I_LENGTH_64                64      

#define MAXINT_64                  ((a_intg)0x7FFFFFFFFFFFL)

#define MININT_64                  ((a_intg)0x800000000000L)

/*------ constants of type 'a_btyp' ----------------------------*/

/*    Mask for upper half of a_btyp value.                      */
#define HIGH_MASK_64               ((a_btyp_64)0xFFFFFFFF00000000ul)

/*    Mask for lower half of a_btyp value.                      */
#define LOW_MASK_64                ((a_btyp_64)0x00000000FFFFFFFFl)

/*    Mask for least significant bit position of a_btyp value.  */
#define LSB_64                     ((a_btyp_64)0x0000000000000001l)

/*    Maximum possible a_btyp value.                            */
#define MAX_BASETYPE_64            ((a_btyp_64)0xFFFFFFFFFFFFFFFFul)

/*    Mask for most significant bit position of a_btyp value.   */
#define MSB_64                     ((a_btyp_64)0x8000000000000000ul)

/*    Inverse mask for most significant bit position of         */
/*    a_btyp value.                                             */
#define NOT_MSB_64                 ((a_btyp_64)0x7FFFFFFFFFFFFFFFl)

#define GETHIGH_64(x)			(((x) >> (B_LENGTH_64/2)) & LOW_MASK_64)
#define GETLOW_64(x)			((x) & LOW_MASK_64)

/****************************************************************/
/*                                                              */
/*      Filename        : d_defs.h                              */
/*                                                              */
/****************************************************************/

/*    Number of 'a_btyp' values used for fraction part.         */
#define A_F_LENGTH_64              (A_F_BITS/B_LENGTH_64+1)

/*    Number of 'a_btyp' values used for integer part.          */
#define A_I_LENGTH_64              (A_I_BITS/B_LENGTH_64+1)

#define A_LENGTH_64                (A_I_LENGTH_64+A_F_LENGTH_64+A_START)

/*    Offset to least significant 'a_btyp'-digit of integer     */
/*    part.                                                     */
#define A_D_P_64                   (A_START+A_I_LENGTH_64-1)

/*    Offset to least significant 'a_btyp' digit of integer     */
/*    part if dotprecision is interpreted as character array.   */
#define B_D_P_64                   (A_D_P_64*sizeof(a_btyp_64))


/****************************************************************/
/*                                                              */
/*      Filename        : r_defs.h                              */
/*                                                              */
/*                   MINUS_SIGN / PLUS_SIGN                     */
/****************************************************************/

/*    Mask for extracting exponent of IEEE value.               */
#define EXPO_MASK_64               ((a_btyp_64)0x7FF0000000000000l)
/*    Number of bit positions needed to shift the exponent      */
/*    of an IEEE value within a 'a_btyp' value most right.      */
#define EXPO_SHIFT_64              (B_LENGTH_64-ZERO_BITS)

/*    Mask for hidden mantissa bit if IEEE value is equivalenced*/
/*    with an array of type 'a_btyp'.                           */
#define HIDDEN_BIT_64              ((a_btyp_64)0x0010000000000000l)

/*    Mask for mantissa bits in most significant 'a_btyp'       */
/*    value if IEEE value is equivalenced with an array of type */
/*    'a_btyp'.                                                 */
#define MANT_64               	((a_btyp_64)0x000FFFFFFFFFFFFFl)

/*    Inverse mask for hidden mantissa bit if IEEE value is     */
/*    equivalenced with an array of type 'a_btyp'.              */
#define NOT_HIDDEN_BIT_64          ((a_btyp_64)0xFFEFFFFFFFFFFFFFul)

/*    Character used as plus sign.                              */
#define PLUS_SIGN               '+'

/*    Mask for bit positions which nether hold a hidden or      */
/*    visible mantissa bit if IEEE value is equivalenced with   */
/*    an array of type 'a_btyp'.                                */
#define SHFT_MASK_64               ((a_btyp_64)0xFFE0000000000000ul)
#define D_U_RATIO_64			(sizeof(a_real)/sizeof(a_btyp_64))
#define BSIZE_64				(2*D_U_RATIO_64+1)
#define ASSERT(a) 


/*    Clear all bits of a 'a_btyp' array.                       */
/*    Actual argument s must not be an expression.              */
#if O_P_1
#define B_CLEAR_64(s,n)    { \
        (void)memset((char *)(s),(int)'\0',(size_t)((n)*sizeof(a_btyp_64)));}
#else
#define B_CLEAR_64(s,n) {size_t _; _=n; while(--_>=0) (s)[_]=ZERO;}
#endif

/*    Copy all bits of a 'a_btyp'/'char' array from left to right.     */
/*    Actual arguments d and s must not be expressions.                */
#if O_P_2
#define B_COPY_64(d,s,n)   { \
        (void)memcpy((char *)(d),(char *)(s),(size_t)((n)*sizeof(a_btyp_64)));}
#else
#define B_COPY_64(d,s,n) {size_t _; _=n; while(--_>=0) (d)[_]=(s)[_];}
#endif

/****************************************************************/
/*                                                              */
/*      Filename        : e_defs.h                              */
/*                                                              */
/****************************************************************/

/*    Macro for testing for a signaling NaN                     */
#define SIGNAL_BIT_64              ((a_btyp_64)0x0008000000000000ul)
#define SIGNALING_64(a)            ((a) & SIGNAL_BIT_64)
#define SET_QUIET_64(a)            ((a) & ~SIGNAL_BIT_64)
#define SET_SIGNAL_64(a)		((a) | SIGNAL_BIT_64)

/*    Macro for testing for a mantissa used for infinity        */
#define MANT_INFINITY_64(a)        ((a) == HIDDEN_BIT_64)

typedef long int a_intg_64;
#define ZERO_64 				0x0ul
#endif





