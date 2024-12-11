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

/* CVS $Id: t_cnst.h,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/* ------------------------------------------------------------- */
/* Header Develop Konstanten fuer beide Versionen                */
/* ------------------------------------------------------------- */

#define EPS_CO                                            /* 2**-62 */ \
           EXTREAL(0x3F, 0xC1,                                         \
                   0x80, 0x00, 0x00, 0x00,                             \
                   0x00, 0x00, 0x00, 0x00) /* +2.168404344971009E-019 */
#define EPS_L                                             /* 2**-64 */ \
           EXTREAL(0x3F, 0xBF,                                         \
                   0x80, 0x00, 0x00, 0x00,                             \
                   0x00, 0x00, 0x00, 0x00) /* +5.421010862427522E-020 */

/* ------------------------------------------------------------- */
/* Header Konstanten fuer EmulationsVersion                      */
/* ------------------------------------------------------------- */

#define EPS_SIN                                                  \
           EXTREAL(0x3F, 0xC1,                                   \
                   0xEB, 0xB4, 0x82, 0x8D,                       \
                   0x7E, 0x43, 0x98, 0x00) /* +3.993000000000000E-019 */

#define EPS_SIN_HPrec                      \
           EXTREAL(0x3F, 0xC1,             \
                   0xA0, 0xDB, 0x09, 0x25, \
                   0xA4, 0x62, 0x90, 0x00) /* +2.725000000000000E-019 */

#define EPS_TAN                            \
           EXTREAL(0x3F, 0xC2,             \
                   0x9C, 0x57, 0x00, 0x1B, \
                   0x3C, 0xE4, 0xA0, 0x00) /* +5.297000000000000E-019 */

/* ------------------------------------------------------------- */
/* Wurzel                                                        */
/* ------------------------------------------------------------- */
#define EPS_SQRT_HPrec                     \
           EXTREAL(0x3F, 0xBE,             \
                   0x80, 0x00, 0x00, 0x00, \
                   0x00, 0x00, 0x00, 0x0B) /* +2.710505431213762E-020 */

#define EPS_SQRT                           \
           EXTREAL(0x3F, 0xBF,             \
                   0xC0, 0x00, 0x00, 0x00, \
                   0x00, 0x00, 0x00, 0x05) /* +8.131516293641284E-020 */

/* --------------------------------------------------------------- */
/* trigonometrische Umkehrfunktionen                               */
/* --------------------------------------------------------------- */
#define EPS_ATAN                           \
           EXTREAL(0x3F, 0xC2,             \
                   0x8E, 0x96, 0x01, 0x03, \
                   0xEA, 0x41, 0xF0, 0x00) /* +4.831000000000000E-019 */

#define EPS_ASIN                           \
           EXTREAL(0x3F, 0xC3,             \
                   0xB1, 0xAD, 0xD5, 0x7E, \
                   0x99, 0x64, 0xE0, 0x00) /* +1.204000000000000E-018 */

#define EPS_ACOS                           \
           EXTREAL(0x3F, 0xC3,             \
                   0xFF, 0xE4, 0xAA, 0xF8, \
                   0x3D, 0xDF, 0xF8, 0x00) /* +1.734000000000000E-018 */

#define EPS_ACOT                           \
           EXTREAL(0x3F, 0xC3,             \
                   0x99, 0x2E, 0x88, 0x01, \
                   0x8E, 0x73, 0x78, 0x00) /* +1.038000000000000E-018 */

/* ---------------------------------------------------------------- */
/* exp, ln                                                          */
/* ---------------------------------------------------------------- */
#define EPS_EXP                            \
           EXTREAL(0x3F, 0xC1,             \
                   0x8A, 0xC7, 0x4E, 0xB9, \
                   0xE3, 0x50, 0x40, 0x00) /* +2.351000000000000E-019 */

#define EPS_LN                             \
           EXTREAL(0x3F, 0xC1,             \
                   0x92, 0x67, 0xB9, 0x04, \
                   0x59, 0x36, 0x00, 0x00) /* +2.480200000000000E-019 */

#define EPS_LN_REL                         \
           EXTREAL(0x3F, 0xC2,             \
                   0xAC, 0x4C, 0x6A, 0x02, \
                   0xEF, 0x3B, 0xC8, 0x00) /* +5.837700000000000E-019 */

#define EPS_LN_ABS                         \
           EXTREAL(0x3F, 0xC2,             \
                   0xA0, 0x2D, 0x40, 0xAD, \
                   0x47, 0xDC, 0x30, 0x00) /* +5.427000000000000E-019 */

/* -------------------------------------------------------------- */
/* hyperbolische                                                  */
/* -------------------------------------------------------------- */
#define EPS_SINH                           \
           EXTREAL(0x3F, 0xC2,             \
                   0xD5, 0xD5, 0xAC, 0x19, \
                   0xD9, 0x5A, 0x28, 0x00) /* +7.245000000000000E-019 */

#define EPS_COSH                           \
           EXTREAL(0x3F, 0xC1,             \
                   0xD0, 0x23, 0x67, 0xCE, \
                   0xD0, 0xF2, 0x88, 0x00) /* +3.526000000000000E-019 */

#define EPS_TANH                           \
           EXTREAL(0x3F, 0xC2,             \
                   0xF6, 0x18, 0x25, 0x93, \
                   0x06, 0x4C, 0x38, 0x00) /* +8.338000000000000E-019 */

#define EPS_COTH                           \
           EXTREAL(0x3F, 0xC2,             \
                   0xE4, 0xA6, 0xAD, 0x51, \
                   0xBC, 0xCF, 0x30, 0x00) /* +7.747000000000000E-019 */

/* ------------------------------------------------------------- */
/* Area                                                          */
/* ------------------------------------------------------------- */
#define EPS_ASINH                          \
           EXTREAL(0x3F, 0xC5,             \
                   0x8C, 0x90, 0x6E, 0xBE, \
                   0xD6, 0xB1, 0xB8, 0x00) /* +3.810000000000000E-018 */

#define EPS_ACOSH                          \
           EXTREAL(0x3F, 0xC5,             \
                   0x89, 0xCC, 0x13, 0xDF, \
                   0x5D, 0x8D, 0xF8, 0x00) /* +3.735000000000000E-018 */

#define EPS_ATANH                          \
           EXTREAL(0x3F, 0xC4,             \
                   0xDC, 0x53, 0xF1, 0xFD, \
                   0x4E, 0x5F, 0x90, 0x00) /* +2.986000000000000E-018 */

#define EPS_ACOTH                          \
           EXTREAL(0x3F, 0xC4,             \
                   0xFE, 0x6A, 0xE0, 0xE7, \
                   0x74, 0xBB, 0xD8, 0x00) /* +3.448000000000000E-018 */

/* ------------------------------------------------------------ */
/* Pow                                                          */
/* ------------------------------------------------------------ */
#define EPS_POW                            \
           EXTREAL(0x3F, 0xC4,             \
                   0xE2, 0x00, 0x8B, 0x92, \
                   0x53, 0xC2, 0xD0, 0x00) /* +3.062900000000000E-018 */





