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

/* CVS $Id: t_mach.h,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

/**************************************************************
**                                                           **
**   mach.h                                                  **
**      Maschinen- und Betriebssystemabhaengige              **
**      Parameter fuer die ieee-Bibliothek                   **
**                                                           **
**   Datum: 08.02.1991 cb                                    **
**                                                           **
**************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_name.h"
#else
#include "t_name.h"
#endif


/* LSBFIRST: 1, falls least significant byte (LSB) auf niedriger */
/*           Speicheradresse liegt, sonst 0                      */
/*           1 fuer Intel-Prozessoren (80386)                    */
/*           0 fuer Motorola-Prozessoren (680x0)                 */
/*           0 fuer HP 9000                                      */
/*           0 fuer IBM RT/PC                                    */

#define LSBFIRST        INTEL

/* ANSI_STRINGIZE: 1, falls der Praeprozessor das ANSI-Stringize- */
/*                 Token (#) kennt, sonst 0                       */
/*                 1 fuer IBM MS C und Turbo C,                   */
/*                 1 fuer Atari Turbo C,                          */
/*                 0 fuer AIX PS/2 cc                             */
/*                 0 fuer HP 9000 HPUX cc                         */
/*                 0 fuer IBM RT/PC AIX cc                        */

#define ANSI_STRINGIZE 0

/* ANSI_C: definiert, falls der C-Compiler dem ANSI-Standard      */
/*         entspricht, sonst undefiniert                          */
/*         definiert fuer IBM MS C und Turbo C                    */
/*         definiert fuer Atari Turbo C                           */
/*         definiert fuer AIX PS/2 cc                             */
/*         nicht definiert fuer IBM PC/RT AIX cc                  */
/*         nicht definiert fuer HP 9000 HPUX cc                   */
/*         nicht definiert fuer SUN cc                            */

/* #define ANSI_C */

/* LINT_ARGS: definiert, falls Typpruefung in den Argumentlisten  */
/*            moeglich ist, sonst undefiniert                     */
/*            definiert fuer IBM MS C und Turbo C                 */
/*            definiert fuer Atari Turbo C                        */
/*            definiert fuer AIX PS/2 cc                          */
/*            nicht definiert fuer IBM PC/RT AIX cc               */
/*            nicht definiert fuer HP 9000 HPUX cc                */

/* #ifndef LINT_ARGS */
/* #define LINT_ARGS */     /* LINT_ARGS already set in p88rts.h  */
/* #endif            */





