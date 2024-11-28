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

/* CVS $Id: f_rstn.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rstn.c                              */
/*                                                              */
/*      Entry           : void f_rstn(desc,spec)                */
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
/*      Description     : reset PASCAL device according to      */
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
local void f_rstn(f_text *desc,a_intg spec)
#else
local void f_rstn(desc,spec)

f_text *desc;
a_intg spec;
#endif
        {
        E_TPUSH("f_rstn")

        /* put final newline character to textfile if opened for output */
        if (desc->asgd==TRUE && desc->fp!=NULL && desc->text==TRUE &&
            desc->outf==TRUE && desc->err==FALSE && desc->eoln==FALSE)
           f_putc((a_char)'\n',desc);

        if (desc->stdo==FALSE && desc->stdi==FALSE && desc->fp!=NULL)
           fclose(desc->fp);
        desc->fp = NULL;

        desc->infl = TRUE;
        desc->outf = desc->err = FALSE;

        switch((int)spec)
           {
           case 0 :     /* assign stdin to TEXT file            */

                        if (desc->text==FALSE)
                           {
                           e_trap(I_O_ERROR,2,E_TMSG,33);
                           E_TPOPP("f_rstn")
                           return;
                           }

                        if (desc->asgd && desc->temp)
                           {
                           remove(desc->name);
                           desc->temp = FALSE;
                           }
                        desc->stdi = TRUE;
                        break;

#if IBM_RT_C+IBM_RS6000_C+IBM_PS2_C+SUN4_OS4_C+SUN4_CPP_C+HP_9000_C
           case 3 :     /* assign TTYIN    to TEXT file         */

                        if (desc->text==FALSE)
                           {
                           e_trap(I_O_ERROR,2,E_TMSG,33);
                           E_TPOPP("f_rstn")
                           return;
                           }

                        if (desc->asgd && desc->temp)
                           {
                           remove(desc->name);
                           desc->temp = FALSE;
                           }

                        /* get path name of controlling terminal*/
                        /* length of path name is defined by    */
                        /* L_ctermid in <stdio.h>.              */
                        if ((desc->fp = fopen(ctermid(desc->name),"r"))==
                            NULL)
                           {
                           e_trap(I_O_ERROR,4,E_TMSG,31,
                                              E_TSTR+E_TEXT(8),desc->name);
                           E_TPOPP("f_rstn")
                           return;
                           }

                        desc->stdi = FALSE;
                        break;
#endif

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
                              E_TPOPP("f_rstn")
                              return;
                              }
                           desc->stdi = TRUE;
                           }
                        else
                           {
                           (void)strcpy(desc->name,desc->org);
                           desc->stdi = FALSE;
                           if ((desc->fp =
                                fopen(desc->org,((desc->text) ? "r" : "rb"))
                               )==NULL)
                              {
                              e_trap(I_O_ERROR,4,E_TMSG,31,
                                     E_TSTR+E_TEXT(8),desc->org);
                              E_TPOPP("f_rstn")
                              return;
                              }
                           }

                        break;

#if ATARI_TURBO_C+IBM_AT_MS_C+IBM_AT_TURBO_C
           case 3 :     /* assign TTYIN    to TEXT file         */
#endif
           case 1 :     /* assign stdout to TEXT file           */
           case 2 :     /* assign stderr to TEXT file           */
           case 4 :     /* assign printer to TEXT file          */
           case 5 :     /* assign reader to file                */
           case 6 :     /* assign puncher to file               */
           case 8 :     /* assign temporary attribute to file   */
           default :    e_trap(I_O_ERROR,6,E_TMSG,43,
                               E_TINT,&spec,E_TSTR+E_TEXT(8),desc->name);
                        E_TPOPP("f_rstn")
                        return;
           }

        desc->eof = desc->eoln = desc->stdo = FALSE;
        desc->asgd = TRUE;

        if (desc->stdi)
           {
           desc->fp = stdin;
           desc->eoln = TRUE;
           desc->win.ch[0] = ' ';
           desc->name[0] = '\0';
           }
        else
           f_getc(desc);

        E_TPOPP("f_rstn")
        return;
        }





