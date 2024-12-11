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

/* CVS $Id: a_exit.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_exit.c                              */
/*                                                              */
/*      Entries         : void a_exit(p)                        */
/*                        a_intg p;                             */
/*                                                              */
/*      Arguments       : exit value passed to environment      */
/*                                                              */
/*      Description     : Exit to environment.                  */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_VOID f_ftop;
#endif

#ifdef LINT_ARGS
local void a_exit(a_intg p)
#else
local void a_exit(p)

a_intg p;
#endif
        {
        register f_text *desc;

        E_TPUSH("a_exit")

        /* close files */
        while (f_ftop!=NULL)
           {
           desc = (f_text *)f_ftop;

           /* put newline character to textfile if opened for output */
           if (desc->asgd==TRUE && desc->fp!=NULL && desc->text==TRUE &&
               desc->outf==TRUE && desc->err==FALSE && desc->eoln==FALSE)
              f_putc((a_char)'\n',desc);

           if (desc->stdi==FALSE && desc->stdo==FALSE && desc->fp!=NULL)
              {
              fclose(desc->fp);

              /* remove temporary file */
              if (desc->temp) remove(desc->name);
              }
           desc->temp = FALSE;
           desc->fp = NULL;

           f_ftop = desc->next;
           }

        exit((int)p);

        E_TPOPP("a_exit")
        return;
        }





