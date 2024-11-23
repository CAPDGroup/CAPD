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

/* CVS $Id: f_assg.c,v 1.22 2014/02/27 15:29:42 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : f_assg.c                              */
/*                                                              */
/*      Entry           : void f_assg(desc,name,len)            */
/*                        f_text *desc;                         */
/*                        char *name;                           */
/*                        size_t len;                           */
/*                                                              */
/*      Arguments       : desc   - file descriptor              */
/*                        name   - name of PASCAL file variable */
/*                        len    - length of file components    */
/*                                                              */
/*      Description     : Assign command line argument to       */
/*                        file descriptor.                      */
/*                                                              */
/*                   f_pmto.fp used in fprintf()                */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern int f_argc;
extern char **f_argv;
extern char *f_pplt[];
extern int f_apos;
extern a_VOID f_ftop;
extern a_bool f_ppmt;
extern a_bool f_pppl;
extern a_bool f_pppd;
extern a_bool f_pptf;
extern f_text f_inpu;
extern f_text f_outp;
extern f_text f_pmti;
extern f_text f_pmto;
extern char *o_text[];
#endif

#ifdef DEMO_VERSION
          /* e_trap(I_O_ERROR,6,E_TMSG,31,E_TSTR+E_TEXT(9),name,\*/
#define DEMOMSG() {fprintf(stderr,\
		"\nPASCAL-XSC demo-version supports ONLY standard I/O !\n\n");\
           e_trap(I_O_ERROR+E_EQIE+E_EXIT,6,E_TMSG,31,E_TSTR+E_TEXT(9),name,\
                  E_TSTR+E_TEXT(8),&desc->name[0]);}

#endif

#ifdef LINT_ARGS
local void f_assg(f_text *desc,char *name,size_t len)
#else
local void f_assg(desc,name,len)

f_text *desc;
char *name;
size_t len;
#endif
        {
        int ch,i,k;  /* !!! must be int !!! */
        char *pos,*nme;

        E_TPUSH("f_assg")

        /* insert file descriptor in linked list */
        desc->next = f_ftop;
        f_ftop = (a_VOID)desc;

        /* display file variable names of program parameter list only */
        if (f_pppd==TRUE && f_pppl==FALSE)
           {
           fprintf(f_pmto.fp,"%s ",name);

           /* save file variable names of program parameter list */
           for (i=0;i<FPPLT;i++)
              {
              if (f_pplt[i]==NULL)
                 {
                 f_pplt[i] = name;
                 if (i<FPPLT) f_pplt[i+1] = NULL;
                 break;
                 }
              }

           E_TPOPP("f_assg")
           return;
           }

        /* initialize file descriptor                           */
        desc->eof  = FALSE;
        desc->eoln = FALSE; 
        desc->err  = FALSE; 
        desc->temp = FALSE;
        desc->infl = FALSE; 
        desc->outf = FALSE;
        desc->stdi = FALSE;
        desc->stdo = FALSE; 
        desc->err  = FALSE;
        desc->asgd = TRUE;
        desc->pp = f_pppl;

        /* reset standard input file variable */
        if (desc==&f_inpu)
           {
           desc->eoln = desc->text = desc->infl = desc->stdi = TRUE;
           desc->fp = stdin;
           desc->ellen = 1;
           desc->win.ch[0] = ' ';

           /* standard input device not listed in program parameter list */
           if (NOT(f_pppl))
              {
              E_TPOPP("f_assg")
              return;
              }
           }

        /* rewrite standard output file variable */
        else if (desc==&f_outp)
           {
           desc->text = desc->outf = desc->stdo = TRUE;
           desc->fp = stdout;
           desc->ellen = 1;

           /* standard output device not listed in program parameter list */
           if (NOT(f_pppl))
              {
              E_TPOPP("f_assg")
              return;
              }
           }

        /* file variable is not standard input/output variable */
        else
           {
#ifdef DEMO_VERSION
           /* no io allowed except of standard io */
		 DEMOMSG();
#else
           /* file variable is a TEXT file if len is 0 on input ... */
           if (len==0)
              {
              len = 1;
              desc->text = TRUE;
              }

           /* ... otherwise its a binary file   */
           else
              desc->text = FALSE;

           desc->fp = NULL;
           desc->ellen = len;

           /* no file name assigned to file variable    */
           /* since not in program parameter list       */
           if (NOT(f_pppl))
              {
              desc->asgd = FALSE;
              E_TPOPP("f_assg")
              return;
              }
#endif /* demo version */
           }

        /* check program parameters for keyword assignment      */
        for (i=1;i<f_argc;i++)
           {

           /* program parameter has keyword assignment          */
           if ((pos = strchr(f_argv[i],'='))!=NULL)
              {

              /* convert keyword to upper characters            */
              for (nme=f_argv[i];nme<pos;nme++)
                  *nme = toupper(*nme);

              /* keyword identified                             */
              if (strncmp(f_argv[i],name,pos-f_argv[i])==0)
                 {

                 /* copy assigned name                          */
                 pos++;
                 desc->org = pos;
                 nme = &desc->name[0];
                 for (k=1;k<f_fnsz;k++)
                    {
                    if (*pos)
                       *nme++ = *pos++;
                    else
                       break;
                    }
                 *nme = '\0';

                 /* file name is too long */
                 if (k==f_fnsz)
                    {
                    e_trap(I_O_BUFFER,6,E_TMSG,30,E_TSTR+E_TEXT(9),name,
                           E_TSTR+E_TEXT(8),&desc->name[0]);
                    }

                 /* standard input file variable is used */
                 else if (desc==&f_inpu)
                    {

                    /* assign standard input device if empty string */
                    if (desc->name[0]=='\0')
                       desc->fp = stdin;
                    else if ((desc->fp = fopen(desc->name,"r"))==NULL)
                       e_trap(I_O_ERROR,6,E_TMSG,31,E_TSTR+E_TEXT(9),name,
                              E_TSTR+E_TEXT(8),&desc->name[0]);

                    /* assign a file name */
                    else
                       {
                       desc->stdi = FALSE;
                       f_getc(desc);
                       }
                    }

                 /* standard output file variable is used */
                 else if (desc==&f_outp)
                    {

                    /* assign standard output device if empty string */
                    if (desc->name[0]=='\0')
                       desc->fp = stdout;
                    else if ((desc->fp = fopen(desc->name,"w"))==NULL)
                       {
                       e_trap(I_O_ERROR,6,E_TMSG,32,E_TSTR+E_TEXT(9),name,
                              E_TSTR+E_TEXT(8),&desc->name[0]);
                       }

                    /* assign a file name */
                    else
                       desc->stdo = FALSE;
                    }

                 /* standard I/O device assigned to file variable */
                 else if (desc->name[0]=='\0')
                    {
                    
                    /* No file name assigned to file variable implies       */
                    /* standard input or standard output device for text    */
                    /* file variables and temporary file for binary file    */
                    /* variable.                                            */
                    if (desc->text)
                       {
                       desc->stdi = desc->stdo = TRUE;
                       }
                    else
                       {
                       desc->asgd = FALSE;
                       desc->temp = TRUE;
                       }
                    }

                 E_TPOPP("f_assg")
                 return;
                 }
              }
           }

        /* update actual position for next positional parameter */
        while (f_apos<f_argc)
           {

           /* ignore command line argument with keyword assignment */
           if (strchr(f_argv[f_apos],'=')!=NULL) f_apos++;

           /* assign filename to file variable */
           else
              {
              desc->org = pos = f_argv[f_apos];
              nme = &desc->name[0];

              for (k=1;*pos!='\0' && k<f_fnsz;k++)
                 *nme++ = *pos++;
              *nme = '\0';

              /* file name too long */
              if (k==f_fnsz)
                 {
                 e_trap(I_O_BUFFER,6,E_TMSG,30,E_TSTR+E_TEXT(9),name,
                        E_TSTR+E_TEXT(8),&desc->name[0]);
                 }

              /* standard input file variable is used */
              else if (desc==&f_inpu)
                 {

                 /* open file for standard input variable */
                 if ((desc->fp = fopen(desc->name,"r"))==NULL)
                    e_trap(I_O_ERROR,6,E_TMSG,31,E_TSTR+E_TEXT(9),name,
                           E_TSTR+E_TEXT(8),&desc->name[0]);
                 else
                    {
                    desc->stdi = FALSE;
                    f_getc(desc);
                    }
                 }

              /* standard output file variable is used */
              else if (desc==&f_outp)
                 {

                 /* open file for standard output variable */
                 if ((desc->fp = fopen(desc->name,"w"))==NULL)
                    {
                    e_trap(I_O_ERROR,6,E_TMSG,32,E_TSTR+E_TEXT(9),name,
                           E_TSTR+E_TEXT(8),&desc->name[0]);
                    }
                 else
                    desc->stdo = FALSE;
                 }
              f_apos++;
              E_TPOPP("f_assg")
              return;
              }
           }

        desc->org = NULL;

        /* no prompting for standard input and output */
        if (desc==&f_inpu || desc==&f_outp)
           {
           E_TPOPP("f_assg")
           return;
           }

#ifdef DEMO_VERSION
           /* no io allowed except of standard io */
		 DEMOMSG();
#else
        /* prompt for program parameter if not input/output     */
        /* and prompting is switched on                         */
        if (f_ppmt)
           {

           /* display prompt */
           fprintf(f_pmto.fp,"(PASCAL file variable %s) %s",name,o_text[11]);

           /* get file name from prompting input device */
           nme = &desc->name[0];
           for (i=1;i<f_fnsz;i++)
              {
              if ((ch = fgetc(f_pmti.fp))=='\n')
                 break;
              else
                 *nme++ = ch;
              }
           *nme = '\0';

           /* No file name assigned to file variable implies       */
           /* standard input or standard output device for text    */
           /* file variables and temporary file for binary file    */
           /* variable.                                            */
           if (i==1)
              {
              if (desc->text)
                 {
                 desc->stdi = desc->stdo = TRUE;
                 }
              else
                 {
                 desc->asgd = FALSE;
                 desc->temp = TRUE;
                 }
              }

           /* file name too long */
           else if (i>=f_fnsz)
              {
              while ((ch = fgetc(f_pmti.fp))!='\n') { /* empty */ }
              e_trap(I_O_BUFFER,6,E_TMSG,30,E_TSTR,name,
                     E_TSTR,&desc->name[0]);
              }
           }

        /* no file name assigned to file variable in */
        /* program parameter list */
        else
           desc->asgd = FALSE;

#endif /* demo version */
        E_TPOPP("f_assg")
        return;
        }





