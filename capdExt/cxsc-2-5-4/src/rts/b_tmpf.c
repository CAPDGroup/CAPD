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

/* CVS $Id: b_tmpf.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tmpf.c                              */
/*                                                              */
/*      Entry           : void b_tmpf(name,length)              */
/*                        char *name;                           */
/*                        a_intg length;                        */
/*                                                              */
/*      Arguments       : name = space for temporary filename   */
/*                        length = length of reserved space     */
/*                                                              */
/*      Description     : Generate a temporary filename         */
/*                        Characters 2 thru 6 of TEMP_BASE      */
/*                        are assumed to be digits.             */
/*                                                              */
/*                   2. mv str... in front of do-> start at t0..*/
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local void b_tmpf(char *name,a_intg length)
#else
local void b_tmpf(name,length)

char *name;
a_intg length;
#endif

        {
        int len,pos;
        FILE *fp;

        E_TPUSH("b_tmpf")

        if ((pos = strlen(o_text[39]))+strlen(o_text[40])>=length)
           {
           e_trap(I_O_BUFFER,2,E_TMSG,30);
           E_TPOPP("b_tmpf")
           return;
           }

        (void)strcpy(name,o_text[39]);
        (void)strcat(name,o_text[40]);
        do
           {
           /* test for existing file                            */
           if ((fp = fopen(name,"r"))!=NULL) fclose(fp);
           else break;

           /* generate next temporary filename                  */
           for (len=pos+5;len>pos;len--)
              {
              if (name[len]=='9') name[len] = '0';
              else
                 {
                 name[len]++;
                 break;
                 }
              }
           }
        while (fp!=NULL);

        E_TPOPP("b_tmpf")
        return;
        }





