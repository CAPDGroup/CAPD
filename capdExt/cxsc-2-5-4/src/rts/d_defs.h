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

/* CVS $Id: d_defs.h,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_defs.h                              */
/*                                                              */
/****************************************************************/

/*    Offset to 'a_btyp' array holding position of first        */
/*    non-zero 'a_btyp'-digit.                                  */
#define A_BEGIN                 0

/*    Offset to 'a_btyp' array holding position of last         */
/*    non-zero 'a_btyp'-digit.                                  */
#define A_END                   1

/*    Offset to 'a_btyp' array holding sign.                    */
#define A_SIGN                  2

/*    Offset to 'a_btyp' array holding status.                  */
#define A_STATUS                3

/*    Offset to 'a_btyp' array holding guard digits.            */
#define A_LOWNAN                4

/*    Offset to 'a_btyp' array holding guard digits.            */
#define A_START                 5

/*    Status information : dotprecision is temporary result.    */
#define A_TEMPORARY             ((a_btyp)0x00000001)

/*    Status information : reserved for special temporaries.    */
/*    (Mask introduced for purposes in CSC-Project)             */
#define A_OWNTEMPORARY          ((a_btyp)0x00000002)

/*    Status information : value is +infinity.                  */
#define A_PINFINITY             ((a_btyp)0x00000004)

/*    Status information : value is -infinity.                  */
#define A_MINFINITY             ((a_btyp)0x00000008)

/*    Status information : value is a quiet NaN.                */
#define A_QUIETNAN              ((a_btyp)0x00000010)

/*    Status information : value is not a sum of -zero          */
#define A_MZERO                 ((a_btyp)0x00000020)

/*    Status information : value is not a sum of +zero          */
#define A_PZERO                 ((a_btyp)0x00000040)

/*    Number of guard digits (bits).                            */
#define A_GUARD_BITS            64

/*    Number of bits used for fractional part.                  */
#define A_F_BITS                (-2*EXPO_MIN+2*MANTL-2)

/*    Number of bits used for integer part.                     */
#define A_I_BITS                (A_GUARD_BITS+2*EXPO_MAX+1)

/*    Number of 'a_btyp' values used for fraction part.         */
#define A_F_LENGTH              (A_F_BITS/B_LENGTH+1)

/*    Number of 'a_btyp' values used for integer part.          */
#define A_I_LENGTH              (A_I_BITS/B_LENGTH+1)

/*    Number of 'a_btyp' values used for dotprecision value.    */
#ifndef DEC_ARITH
#define A_LENGTH                (A_I_LENGTH+A_F_LENGTH+A_START)
#else
#define A_LENGTH                (A_I_LENGTH+A_F_LENGTH+A_START+1)
#endif

/*    Offset to least significant 'a_btyp'-digit of integer     */
/*    part.                                                     */
#define A_D_P                   (A_START+A_I_LENGTH-1)

/*    Offset to least significant 'a_btyp' digit of integer     */
/*    part if dotprecision is interpreted as character array.   */
#define B_D_P                   (A_D_P*sizeof(a_btyp))

/*    Number of characters which can be stored in a             */
/*    dotprecision value.                                       */
#define BUFFERSIZE              (A_LENGTH*sizeof(a_btyp))

/*    Value of not allocated dotprecision variable.             */
#define NO_DOTPRECISION         ((dotprecision)NULL)





