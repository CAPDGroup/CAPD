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

/* CVS $Id: o_user.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : o_user.c                              */
/*                                                              */
/*      Description     : Definition of global objects          */
/*                        which may be changed by the user.     */
/*                                                              */
/*      Note            : All variables must be initialized.    */
/*                                                              */
/****************************************************************/

#include <stdio.h>
/* 
 * some system don't define NULL in stdio.h. In many cases, 
 * especially ANSI C Compiler, this is done in stddef.
 */
#ifndef NULL
#include <stddef.h>
#endif

/*--------------------------------------------------------------*/
/*   o_sdir : Path name of "fixed system directory".            */
/*   ======                                                     */
/*            The length of this variable is restricted         */
/*            to 63 readable characters.                        */
/*--------------------------------------------------------------*/
char o_sdir[64] = "/usr/pxsc/sys/";

/*--------------------------------------------------------------*/
/*   o_udir : Path name of "fixed user directory".              */
/*   ======                                                     */
/*            The length of this variable is restricted         */
/*            to 63 readable characters.                        */
/*--------------------------------------------------------------*/
char o_udir[64] = "/pxsc/";

/*--------------------------------------------------------------*/
/*   o_text : Array of pointers to character sequences used by  */
/*   ======   the runtime system.                               */
/*--------------------------------------------------------------*/
char *o_text[] = 
{

        /*------------------------------------------------------*/
        /*    Header string of runtime error messages.          */
        /*------------------------------------------------------*/
"--- ",                                                  /* [0] */

        /*------------------------------------------------------*/
        /*    Reserved.                                         */
        /*------------------------------------------------------*/
"",                                                      /* [1] */

        /*------------------------------------------------------*/
        /*    Reserved.                                         */
        /*------------------------------------------------------*/
"",                                                      /* [2] */

        /*------------------------------------------------------*/
        /*    Name of environment variable holding              */
        /*    "home directory" of the user.                     */
        /*------------------------------------------------------*/
"HOME",                                                  /* [3] */

        /*------------------------------------------------------*/
        /*    Name of environment variable holding              */
        /*    "user directory" of PASCAL-XSC user.              */
        /*------------------------------------------------------*/
"PXSC_USR",                                              /* [4] */

        /*------------------------------------------------------*/
        /*    Name of environment variable holding              */
        /*    "system directory" of PASCAL-XSC system.          */
        /*------------------------------------------------------*/
"PXSC_SYS",                                              /* [5] */

        /*------------------------------------------------------*/
        /*    Name of environment displayed when exceptions     */
        /*    are traced within the runtime system.             */
        /*------------------------------------------------------*/
"runtime system",                                        /* [6] */

        /*------------------------------------------------------*/
        /*    Name of file containing informations about the    */
        /*    runtime system.                                   */
        /*------------------------------------------------------*/
"info.txt",                                              /* [7] */

        /*------------------------------------------------------*/
        /*    Name of file containing message texts displayed   */
        /*    in case of exceptions.                            */
        /*------------------------------------------------------*/
"o_msg1.h",                                              /* [8] */

        /*------------------------------------------------------*/
        /*    Prompting message displayed for entering the      */
        /*    file name of an input file for a local file       */
        /*    variable.                                         */
        /*------------------------------------------------------*/
"\nEnter filename for input : ",                         /* [9] */

        /*------------------------------------------------------*/
        /*    Prompting message displayed for entering the      */
        /*    file name of an output file for a local file      */
        /*    variable.                                         */
        /*------------------------------------------------------*/
"\nEnter filename for output : ",                       /* [10] */

        /*------------------------------------------------------*/
        /*    Prompting message displayed for entering the      */
        /*    file name of a file associated with a file        */
        /*    variable listed in the program parameter list.    */
        /*------------------------------------------------------*/
"\nEnter program parameter : ",                         /* [11] */

        /*------------------------------------------------------*/
        /*    Runtime options.                                  */
        /*------------------------------------------------------*/
"|?|-?|-info|-help|-h|/h|/help|",                       /* [12] */
"-pr",                                                  /* [13] */
"-cc",                                                  /* [14] */
"-ib",                                                  /* [15] */
"-info",                                                /* [16] */
"-nn",                                                  /* [17] */
"-pp",                                                  /* [18] */
"-sd",                                                  /* [19] */
"-ud",                                                  /* [20] */
"-sz",                                                  /* [21] */
"-tr",                                                  /* [22] */
"-tb",                                                  /* [23] */
"-tf",                                                  /* [24] */
"-vn",                                                  /* [25] */
"-ieee",                                                /* [26] */
"d",                                                    /* [27] */
"i",                                                    /* [28] */
"o",                                                    /* [29] */
"u",                                                    /* [30] */
"x",                                                    /* [31] */

        /*------------------------------------------------------*/
        /*    Device name associated with console input.        */
        /*------------------------------------------------------*/
"",                                                     /* [32] */

        /*------------------------------------------------------*/
        /*    Device name associated with console output.       */
        /*------------------------------------------------------*/
"",                                                     /* [33] */

        /*------------------------------------------------------*/
        /*    Text displayed for boolean value FALSE.           */
        /*------------------------------------------------------*/
"FALSE",                                                /* [34] */

        /*------------------------------------------------------*/
        /*    Text displayed for boolean value TRUE.            */
        /*------------------------------------------------------*/
"TRUE ",                                                /* [35] */

        /*------------------------------------------------------*/
        /*    Text displayed for real value infinity.           */
        /*------------------------------------------------------*/
"infinity",                                             /* [36] */

        /*------------------------------------------------------*/
        /*    Text displayed for real value signaling NaN.      */
        /*------------------------------------------------------*/
"sNaN",                                                 /* [37] */

        /*------------------------------------------------------*/
        /*    Text displayed for real value quiet NaN.          */
        /*------------------------------------------------------*/
"qNaN",                                                 /* [38] */

        /*------------------------------------------------------*/
        /*    Name of directory for storing temporary files.    */
        /*                                                      */
        /*    There must be a trailing path delimiter character.*/
        /*------------------------------------------------------*/
"",                                                     /* [39] */

        /*------------------------------------------------------*/
        /*    Base name of temporary files.                     */
        /*                                                      */
        /*    Total number of characters must be exactly 10.    */
        /*    Leading character may be chosen arbitrarily.      */
        /*    There must be 5 digit positions.                  */
        /*    The extension must have exactly 4 characters.     */
        /*------------------------------------------------------*/
"t00000.tmp",                                           /* [40] */

        /*------------------------------------------------------*/
        /*    Reserved.                                         */
        /*------------------------------------------------------*/
""                                                      /* [41] */
};

/*--------------------------------------------------------------*/
/*   o_errr : File pointer variable specifying the device for   */
/*   ======   displaying of error messages.                     */
/*--------------------------------------------------------------*/
FILE *o_errr = NULL;

/*--------------------------------------------------------------*/
/*   o_pmti : File pointer variable specifying the device for   */
/*   ======   reading file names when prompted by the runtime   */
/*            system.                                           */
/*--------------------------------------------------------------*/
FILE *o_pmti = NULL;

/*--------------------------------------------------------------*/
/*   o_pmto : File pointer variable specifying the device for   */
/*   ======   displaying the prompting message by the runtime   */
/*            system.                                           */
/*--------------------------------------------------------------*/
FILE *o_pmto = NULL;

/*--------------------------------------------------------------*/
/*   o_user : This routine is called immediately after starting */
/*   ======   an executable PASCAL-XSC program.                 */
/*            The initialization of file pointers must be done  */
/*            in this routine.                                  */
/*--------------------------------------------------------------*/
void o_user()
        {
        o_errr = stderr;
        o_pmti = stdin;
        o_pmto = stdout;
        }






