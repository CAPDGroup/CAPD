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

/* CVS $Id: f_srse.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_srse.c                              */
/*                                                              */
/*      Entry           : void f_srse(desc,device)              */
/*                        f_text *desc;                         */
/*                        s_trng device;                        */
/*                                                              */
/*      Arguments       : desc   - descriptor of device         */
/*                        device - name of device               */
/*                          ""   = standard input (text only)   */
/*                          file = filename used for fopen      */
/*                                                              */
/*      Description     : reset PASCAL device.                  */
/*                                                              */
/*                   name argument removed from argument list   */
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
local void f_srse(f_text *desc,s_trng device)
#else
local void f_srse(desc,device)

f_text *desc;
s_trng device;
#endif
        {
        char *filename = NULL;

        E_TPUSH("f_srse")

        /* put final newline character to textfile if opened for output */
        if (desc->asgd==TRUE && desc->fp!=NULL && desc->text==TRUE &&
            desc->outf==TRUE && desc->err==FALSE && desc->eoln==FALSE)
           f_putc((a_char)'\n',desc);

        if (desc->stdo==FALSE && desc->stdi==FALSE && desc->fp!=NULL)
           fclose(desc->fp);
        desc->fp = NULL;

        desc->infl = TRUE;
        desc->outf = desc->err = FALSE;

        /* open standard input                          */
        if (device.clen<=0)
           {
           if (desc->temp)
              {
              remove(desc->name);
              desc->temp = FALSE;
              }

           if (desc->text==TRUE)
              {
              desc->stdi = TRUE;
              desc->name[0] = '\0';
              }
           else
              {
              e_trap(I_O_ERROR,2,E_TMSG,33);
              E_TPOPP("f_srse")
              return;
              }
           }

        /* use filename                                 */
        else
           {
           if (desc->temp)
              {
              remove(desc->name);
              desc->temp = FALSE;
              }

           if (device.clen>=f_fnsz-1)
              {
              e_trap(I_O_BUFFER,6,E_TMSG,29,E_TMSG,30,
                                  E_TSTG+E_TEXT(8),device);
              E_TPOPP("f_srse")
              return;
              }
           filename = &desc->name[0];
           memcpy(filename,device.ptr,device.clen);
           filename[device.clen] = '\0';

           desc->stdi = FALSE;
           }

        desc->eof = desc->eoln = desc->stdo = FALSE;
        desc->asgd = TRUE;

        if (desc->stdi)
           {
           desc->fp = stdin;
           desc->eoln = TRUE;
           desc->win.ch[0] = ' ';
           }
        else if ((desc->fp =
                  fopen(filename,((desc->text) ? "r" : "rb")))==NULL)
           {
           e_trap(I_O_ERROR,4,E_TMSG,31,E_TSTR+E_TEXT(8),filename);
           desc->err = TRUE;
           }
        else
           f_getc(desc);

        if (device.tmp) s_free(&device);

        E_TPOPP("f_srse")
        }





