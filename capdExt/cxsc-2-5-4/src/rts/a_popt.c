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

/* CVS $Id: a_popt.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_popt.c                              */
/*                                                              */
/*      Entry           : void a_popt(desc,option)              */
/*                        f_text *desc;                         */
/*                        s_trng option;                        */
/*                                                              */
/*      Arguments       : desc   - file descriptor              */
/*                        option - name of option               */
/*                                                              */
/*      Description     : Check for runtime option and process  */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern d_otpr b_cp__;
#endif

#ifdef LINT_ARGS
local void a_popt(f_text *desc,s_trng option)
#else
local void a_popt(desc,option)

f_text *desc;
s_trng option;
#endif

        {
        int i,k;
        char *buffer;

        E_TPUSH("a_popt")

        if (b_text(desc,FALSE))
           {
           for (i=0;i<option.clen;)
              {

              /* search for start of option */
              while (i<option.clen)
                 if (option.ptr[i]!=' ') break; else i++;

              /* search for end of option */
              k = i+1;
              while (k<option.clen)
                 if (option.ptr[k]==' ') break; else k++;

              /* process option if any */
              if (i<option.clen)
                 {
                 buffer = (char *)memcpy((char *)b_cp__,option.ptr+i,k-i);
                 buffer[k-i] = '\0';
                 (void)b_popt(desc->fp,buffer);
                 i = k+1;
                 }
              }
           }

        if (option.tmp) s_free(&option);

        E_TPOPP("a_popt")
        return;
        }





