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

/* CVS $Id: f_rset.c,v 1.21 2014/01/30 17:24:07 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_rset.c                              */
/*                                                              */
/*      Entry           : void f_rset(desc,name,device)         */
/*                        f_text *desc;                         */
/*                        a_char *name;                         */
/*                        a_char *device;                       */
/*                                                              */
/*      Arguments       : desc   - descriptor of device         */
/*                        name   - PASCAL file variable name    */
/*                        device - name of device               */
/*                          NULL = reuse file descriptor        */
/*                          ""   = standard input (text only)   */
/*                          file = filename used for fopen      */
/*                                                              */
/*      Description     : reset PASCAL device.                  */
/*                                                              */
/*                   prompting optional for local file variable */
/*                   f_pmti, f_pmto replace stdin, stderr       */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern f_text f_pmti;
extern f_text f_pmto;
extern a_bool f_pptf;
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local void f_rset(f_text *desc,a_char *name,a_char *device)
#else
local void f_rset(desc,name,device)

f_text *desc;
a_char *name;
a_char *device;
#endif
        {
        char *p;
        int i,ch;       /* !!! must be int !!! */

        E_TPUSH("f_rset")

        /* put final newline character to textfile if opened for output */
        if (desc->asgd==TRUE && desc->fp!=NULL && desc->text==TRUE &&
            desc->outf==TRUE && desc->err==FALSE && desc->eoln==FALSE)
           f_putc((a_char)'\n',desc);

        if (desc->stdo==FALSE && desc->stdi==FALSE && desc->fp!=NULL)
           fclose(desc->fp);
        desc->fp = NULL;

        desc->infl = TRUE;
        desc->outf = desc->err = FALSE;

        /* reuse descriptor                                     */
        if (device==NULL)
           {

           /* file assigned to descriptor                       */
           if (desc->asgd)
              {
              if (desc->stdo)
                 desc->stdi = TRUE;
              else if (desc->stdi==FALSE)
                 device = (a_char *)&desc->name[0];
              }

           /* no device assigned to descriptor                  */

           /* prompt for filename if program parameter or       */
           /* requested by runtime flag f_pptf                  */
           else if (desc->pp || f_pptf)
              {
new_file_name:
              fprintf(f_pmto.fp,"(PASCAL file variable %s) %s",
                      (char *)name,o_text[9]);
              device = (a_char *)&desc->name[0];

              /* length of input must be less than f_fnsz       */
              for (i=1;i<f_fnsz;i++)
                 if ((ch = fgetc(f_pmti.fp))=='\n') break;
                 else *device++ = ch;
              *device = '\0';

              desc->asgd = TRUE;
              device = (a_char *)&desc->name[0];

              if (i==1)
                 {
                 if (desc->text==TRUE)
                    {
                    desc->stdi = TRUE;
                    desc->name[0] = '\0';
                    }
                 else
                    {
                    e_trap(NO_ERROR+E_EMSG+E_ECNT,6,E_TMSG,68,E_TMSG,33,
                           E_TSTR+E_TEXT(9),name);
                    goto new_file_name;
                    }
                 }
              else if (i>=f_fnsz)
                 {
                 while ((ch = fgetc(f_pmti.fp))!='\n') { /* empty */ }
                 e_trap(I_O_BUFFER,6,E_TMSG,30,
                        E_TSTR+E_TEXT(9),name,E_TSTR+E_TEXT(8),device);
                 E_TPOPP("f_rset")
                 return;
                 }
              }

           /* using temporary file causes error                 */
           else
              {
              e_trap(I_O_ERROR,4,E_TMSG,38,E_TSTR+E_TEXT(9),name);
              E_TPOPP("f_rset")
              return;
              }
           }

        /* open standard input                          */
        else if (*device=='\0')
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
              e_trap(I_O_ERROR,4,E_TMSG,33,E_TSTR+E_TEXT(9),name);
              E_TPOPP("f_rset")
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

           p = &desc->name[0];
           for (i=1;(*p++ = device[i-1])!='\0' && i<f_fnsz;i++)
              { /* empty */ }
           *p = '\0';

           desc->stdi = FALSE;

           if (i>=f_fnsz)
              {
              e_trap(I_O_BUFFER,6,E_TMSG,30,
                     E_TSTR+E_TEXT(9),name,E_TSTR+E_TEXT(8),device);
              E_TPOPP("f_rset")
              return;
              }
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
                  fopen((char *)device,((desc->text) ? "r" : "rb")))==NULL)
           {
           e_trap(I_O_ERROR,6,E_TMSG,31,
                  E_TSTR+E_TEXT(9),name,E_TSTR+E_TEXT(8),device);
           desc->err = TRUE;
           }
        else
           f_getc(desc);

        E_TPOPP("f_rset")
        return;
        }





