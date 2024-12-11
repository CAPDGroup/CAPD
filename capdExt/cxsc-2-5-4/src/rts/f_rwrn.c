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

/* CVS $Id: f_rwrn.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rwrn.c                              */
/*                                                              */
/*      Entry           : void f_rwrn(desc,spec)                */
/*                        f_text *desc;                         */
/*                        a_intg spec;                          */
/*                                                              */
/*      Arguments       : desc   - descriptor of device         */
/*                        spec   - specification                */
/*                           0 = stdin                          */
/*                           1 = stdout                         */
/*                           2 = stderr                         */
/*                           3 = stdcon                         */
/*                           4 = stdprn                         */
/*                           5 = stdrdr                         */
/*                           6 = stdpun                         */
/*                           8 = stdtmp                         */
/*                           9 = stdorg                         */
/*                                                              */
/*      Description     : rewrite PASCAL device according to    */
/*                        specification                         */
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
local void f_rwrn(f_text *desc,a_intg spec)
#else
local void f_rwrn(desc,spec)

f_text *desc;
a_intg spec;
#endif
        {
        E_TPUSH("f_rwrn")

        /* put final newline character to textfile if opened for output */
        if (desc->asgd==TRUE && desc->fp!=NULL && desc->text==TRUE &&
            desc->outf==TRUE && desc->err==FALSE && desc->eoln==FALSE)
           f_putc((a_char)'\n',desc);

        if (desc->asgd &&
            desc->stdo==FALSE && desc->stdi==FALSE && desc->fp!=NULL)
           fclose(desc->fp);
        desc->fp = NULL;

        desc->eoln = desc->outf = TRUE;
        desc->infl = desc->err = FALSE;

        switch((int)spec)
           {
           case 1 :     /* assign stdout to TEXT file           */

           case 2 :     /* assign stderr to TEXT file           */

                        if (desc->text==FALSE)
                           {
                           e_trap(I_O_ERROR,2,E_TMSG,33);
                           E_TPOPP("f_rwrn")
                           return;
                           }

                        if (desc->asgd && desc->temp)
                           {
                           remove(desc->name);
                           desc->temp = FALSE;
                           }
                        desc->stdo = TRUE;
                        break;


#if IBM_RT_C+IBM_RS6000_C+IBM_PS2_C+HP_9000_C+SUN4_OS4_C+SUN4_CPP_C
           case 3 :     /* assign TTYOUT   to TEXT file         */

                        if (desc->asgd && desc->temp)
                           {
                           remove(desc->name);
                           desc->temp = FALSE;
                           }

                        /* get path name of controlling terminal*/
                        /* length of path name is defined by    */
                        /* L_ctermid in <stdio.h>.              */
                        (void)ctermid(desc->name);

                        desc->stdo = FALSE;
                        break;
#endif

           case 8 :     /* assign temporary attribute to file   */

                        if (desc->temp==FALSE)
                           {
                           b_tmpf(desc->name,f_fnsz);
                           desc->temp = TRUE;
                           }
                        desc->stdo = FALSE;
                        break;

           case 9 :     /* assign command line association      */

                        if (desc->asgd && desc->temp)
                           {
                           remove(desc->name);
                           desc->temp = FALSE;
                           }

                        if (desc->org==NULL || *desc->org=='\0')
                           {
                           if (desc->text==FALSE)
                              {
                              e_trap(I_O_ERROR,2,E_TMSG,33);
                              E_TPOPP("f_rwrn")
                              return;
                              }
                           desc->stdo = TRUE;
                           }
                        else
                           {
                           (void)strcpy(desc->name,desc->org);
                           desc->stdo = FALSE;
                           }

                        break;

#if ATARI_TURBO_C+IBM_AT_MS_C+IBM_AT_TURBO_C
           case 3 :     /* assign TTYOUT   to TEXT file         */
#endif
           case 0 :     /* assign stdin to TEXT tile            */
           case 4 :     /* assign printer to TEXT file          */
           case 5 :     /* assign reader to file                */
           case 6 :     /* assign puncher to file               */
           default :    e_trap(I_O_ERROR,6,E_TMSG,43,E_TINT,&spec,
                               E_TSTR+E_TEXT(8),desc->name);
                        E_TPOPP("f_rwrn")
                        return;
           }

        desc->stdi = FALSE;
        desc->eof = desc->asgd = TRUE;

        if (desc->stdo)
           {
           if (spec==2)
              desc->fp = stderr;
           else
              desc->fp = stdout;
           }
        else if ((desc->fp =
                 fopen(desc->name,((desc->text) ? "w" : "wb")))==NULL)
           {
           e_trap(I_O_ERROR,4,E_TMSG,32,E_TSTR,desc->name);
           desc->err = TRUE;
           }

        E_TPOPP("f_rwrn")
        return;
        }





