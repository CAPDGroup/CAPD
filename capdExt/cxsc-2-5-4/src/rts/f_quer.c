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

/* CVS $Id: f_quer.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_quer.c                              */
/*                                                              */
/*      Entry           : void f_quer(desc,status,name)         */
/*                        f_text *desc;                         */
/*                        a_intg *status;                       */
/*                        s_trng *name;                         */
/*                                                              */
/*      Arguments       : desc   - file descriptor              */
/*                        status - status information           */
/*                             1 - file assigned                */
/*                             2 - file open                    */
/*                             4 - ... with existing pp-name    */
/*                             8 - ... for reading              */
/*                            16 - ... from standard I/O device */
/*                            32 - ... with end-of-line status  */
/*                            64 - ... with end-of-file status  */
/*                           128 - ... with error status        */
/*                        name   - assigned file name           */
/*                                                              */
/*      Description     : Query file descriptor.                */
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
local void f_quer(f_text *desc,a_intg *status,s_trng *name)
#else
local void f_quer(desc,status,name)

f_text *desc;
a_intg *status;
s_trng *name;
#endif
        {

#define ASSIGNED_STATUS         1
#define OPEN_STATUS             2
#define PP_NAME_STATUS          4
#define READING_STATUS          8
#define STANDARD_STATUS        16
#define END_OF_LINE_STATUS     32
#define END_OF_FILE_STATUS     64
#define ERROR_STATUS          128

        *status = 0;

        if (desc->asgd)
           {
           (*status) += ASSIGNED_STATUS;
           if (desc->org!=NULL) (*status) += PP_NAME_STATUS;
           if (desc->outf || desc->infl)
              {
              (*status) += OPEN_STATUS;
              if (desc->infl) (*status) += READING_STATUS;
              if (desc->stdi || desc->stdo) (*status) += STANDARD_STATUS;
              if (desc->eoln) (*status) += END_OF_LINE_STATUS;
              if (desc->eof) (*status) += END_OF_FILE_STATUS;
              if (desc->err) (*status) += ERROR_STATUS;
              }

           s_asta(name,(a_char *)desc->name,(a_intg)strlen(desc->name));
           }
        else
           s_free(name);

        return;
        }





