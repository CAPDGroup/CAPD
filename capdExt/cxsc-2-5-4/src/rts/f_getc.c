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

/* CVS $Id: f_getc.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_getc.c                              */
/*                                                              */
/*      Entry           : void f_getc(desc)                     */
/*                        f_text *desc;                         */
/*                                                              */
/*      Arguments       : desc   - descriptor of device         */
/*                                                              */
/*      Description     : get a character from a checked        */
/*                        device.                               */
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
local void f_getc(f_text *desc)
#else
local void f_getc(desc)

f_text *desc;
#endif
        {
        int ch; /* !!! must be int !!! */
        size_t i;

        E_TPUSH("f_getc")

        if (desc->text)
           {
           if ((ch = fgetc(desc->fp))==EOF)
              {
              desc->win.ch[0] = ' ';
              desc->eof = TRUE;
              desc->eoln = FALSE;
              }
           else if (ch=='\n')
              {
              desc->win.ch[0] = ' ';
              desc->eoln = TRUE;
              }
           else
              {
              desc->win.ch[0] = ch;
              desc->eoln = FALSE;
              }
           }
        else
           {
           for (i=0;i<desc->ellen;i++)
              {
              if ((ch = fgetc(desc->fp))==EOF)
                 {
                 desc->eof = TRUE;
                 if (i)
                    e_trap(I_O_ERROR,4,E_TMSG,20,
                           E_TSTR+E_TEXT(8),desc->name);
                 break;
                 }
              desc->win.ch[i] = ch;
              }
           }

        E_TPOPP("f_getc")
        return;
        }





