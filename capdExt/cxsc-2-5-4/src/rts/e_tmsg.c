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

/* CVS $Id: e_tmsg.c,v 1.23 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_tmsg.c                              */
/*                                                              */
/*      Entry           : void e_tmsg(msgid)                    */
/*                        int msgid;                            */
/*                                                              */
/*      Arguments       : msgid  - message identification       */
/*                                                              */
/*      Description     : Message handling.                     */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern char *e_head;
extern f_text f_errr;
extern char *o_text[];
#endif

#ifdef MSG_FILE_ENABLED
static char *e_mfld = NULL;
#else
static e_mtyp e_mfld[] =
{
/*----------------------------------------------*/
/* Include file name must be referencable via   */
/* o_text[8].                                   */
#ifdef AIX
#include "/u/p88c/runtime/o_msg1.h"
#else
#include "o_msg1.h"
#endif
/*----------------------------------------------*/

/* The first component of the last item         */
/* must be zero to indicate "END OF LIST".      */
{   0,  "" }
};
#endif

#ifdef LINT_ARGS
local void e_tmsg(int msgid)
#else
local void e_tmsg(msgid)
int msgid;
#endif
        {
#ifdef MSG_FILE_ENABLED
        f_text desc;
        s_trng name;
        char ch;
        char *q;
        size_t i=0,k=0;
        int item=0;
        char nr[8];

        if (e_mfld==NULL)
           {
           desc.text = TRUE;
           name.ptr = o_text[8];
           name.suba = TRUE;
           name.fix = name.tmp = FALSE;
           name.alen = name.clen = strlen(o_text[8]);
           if (b_op88(&desc,name,90)==0)
              {
              fprintf(f_errr.fp,"%se_tmsg : Unable to open ",e_head);
              fprintf(f_errr.fp,"message file \"%s\".\n",o_text[8]);
              return;
              }
           else
              {
              ch = fgetc(desc.fp);
              while(!feof(desc.fp))
                 {
                 if (ch=='{') item = 1;
                 if (item)
                    {
                    if (i>=k)
                       {
                       k += 256;
                       if ((q = (char*) malloc(k+1))==NULL)
                          {
                          fprintf(f_errr.fp,"%se_tmsg : Allocation error\n",
                                  e_head);
                          return;
                          }
                       if (e_mfld!=NULL)
                          {
                          memcpy(q,e_mfld,i);
                          B_FREE(e_mfld);
                          }
                       e_mfld = q;
                       }

                    switch(item)
                       {
                       case 1: if (isdigit((int)ch) || ch==',' || ch=='{')
                                  {
                                  e_mfld[i++] = ch;
                                  }
                               if (ch==',') item++;
                               break;
                       case 2: if (ch=='"')
                                  {
                                  e_mfld[i++] = ch;
                                  item++;
                                  }
                               break;
                       case 3: if (ch=='"') item++;
                               e_mfld[i++] = ch;
                               break;
                       case 4: if (ch=='}')
                                  {
                                  e_mfld[i++] = ch;
                                  item = 0;
                                  }
                               break;
                       }
                    }
                 ch = fgetc(desc.fp);
                 }
              e_mfld[i] = '\0';
              }
           fclose(desc.fp);
           }

        sprintf(nr,"%d",msgid);
        q = e_mfld;
        while ((q = strchr(q,'{'))!=NULL)
           {
           q++;
           if (memcmp(nr,q,strlen(nr))==0)
              {
              q += strlen(nr)+2;

              /* header displayed if "\r" is not the first character  */
              if (*q!='\\' || q[1]!='r') fprintf(f_errr.fp,"%s",e_head);
              else q += 2;

              /* quotation mark delimits the displayed string */
              while (*q!='"')
                 {

                 /* consider "\r","\n" */
                 if (*q=='\\')
                    {
                    q++;
                    if (*q=='r')
                       {
                       fprintf(f_errr.fp,"\n");
                       if (q[1]!='"')
                          fprintf(f_errr.fp,"%*s",(int)strlen(e_head)," ");
                       }
                    else if (*q=='n')
                       {
                       fprintf(f_errr.fp,"\n");
                       if (q[1]!='"') fprintf(f_errr.fp,"%s",e_head);
                       }
                    }

                 /* display string character */
                 else
                    fprintf(f_errr.fp,"%c",*q);
                 q++;
                 }
              return;
              }
           }
#else
        int i;  /* !!! must be int !!! */
        char *q;

        for (i=0;e_mfld[i].msgid;i++)
           {
           if (e_mfld[i].msgid==msgid)
              {
              q = e_mfld[i].text;
              if (*q)
                 {
                 if (*q!='\r') fprintf(f_errr.fp,"%s",e_head);
                 else q++;
                 while (*q)
                    {
                    if (*q=='\r')
                       {
                       fprintf(f_errr.fp,"\n");
                       if (q[1])
                          fprintf(f_errr.fp,"%*s",strlen(e_head)," ");
                       }
                    else if (*q=='\n')
                       {
                       fprintf(f_errr.fp,"\n");
                       if (q[1]) fprintf(f_errr.fp,"%s",e_head);
                       }
                    else
                       fprintf(f_errr.fp,"%c",*q);
                    q++;
                    }
                 }
              return;
              }
           }
#endif
        }





