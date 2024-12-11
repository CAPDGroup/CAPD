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

/* CVS $Id: f_rint.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rint.c                              */
/*                                                              */
/*      Entry           : a_intg f_rint(device,i)               */
/*                        FILE *device;                         */
/*                        a_intg *i;                            */
/*                                                              */
/*      Arguments       : device - text file to be read         */
/*                        i - integer value read                */
/*                            at input: i holds window character*/
/*                                                              */
/*      Return value    : last character read                   */
/*                                                              */
/*      Description     : perform PASCAL integer read.          */
/*                                                              */
/*      Note            : BS PASCAL Standard 6.9.1 c)           */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local a_intg f_rint(FILE *device,a_intg *i)
#else
local a_intg f_rint(device,i)

FILE *device;
a_intg *i;
#endif
        {
        int c;
        a_bool sign;
        char ch;

        E_TPUSH("f_rint")

        sign = FALSE;
        c = (int)*i;

        /* ignore leading blanks and linefeed characters        */
        while (c==' ' || c==EOLN) c = fgetc(device);

        /* check for sign                                       */
        if (c=='+') c = fgetc(device);
        else if (c=='-')
           {
           c = fgetc(device);
           sign = TRUE;
           }

        /* handle first decimal digit                           */
        if (!isdigit(c))
           {
           if (c==EOF)
              {
              e_trap(I_O_ERROR,2,E_TMSG,20);
              }
           else
              {
              ungetc(c,device);
              ch = c;
              e_trap(I_O_ERROR,4,E_TMSG,21,E_TCHR+E_TEXT(10),&ch);
              }
           E_TPOPP("f_rint")
           return(c);
           }

        *i = c-'0';
        c = fgetc(device);

        /* read decimal digits                                  */
        while (isdigit(c))
           {

           /* check for integer overflow                        */
           if (MAXINT/10-*i<=0)
              {
              if ((MAXINT/10-*i<0) || (c-'0'>MAXINT%10))
                 {
                 if (sign) *i = -*i;
                 e_trap(OVERFLOW,2,E_TMSG,15);
                 do c = fgetc(device);
                 while (isdigit(c));
                 if (c!=EOF) ungetc(c,device);
                 E_TPOPP("f_rint")
                 return(c);
                 }
              }

           /* add digit                                         */
           *i = (*i*10)+c-'0';
           c = fgetc(device);
           }

        if (sign) *i = -*i;

        if (c!=EOF) ungetc(c,device);

        E_TPOPP("f_rint")
        return(c);
        }





