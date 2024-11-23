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

/* CVS $Id: e_push.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_push.c                              */
/*                                                              */
/*      Entries         : void e_push(function,filename)        */
/*                        char *function;                       */
/*                        char *filename;                       */
/*                                                              */
/*      Arguments       : function - name of function           */
/*                        filename - name of file               */
/*                                                              */
/*      Description     : push function name and file name      */
/*                        to trace back stack                   */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern bentry *e_btop,*e_bstk;
extern int e_line;
extern char *e_head;
#ifdef TRACEC_ENABLE
extern a_bool f_pptr;
extern int e_tlvl;
#endif
#ifdef STATTB_ENABLE
extern a_bool f_pptb;
#endif
extern f_text f_errr;
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local void e_push(char *function,char *filename)
#else
local void e_push(function,filename)

char *function;
char *filename;
#endif
        {
#ifdef TRACEC_ENABLE
        int i;

        if (f_pptr)
#ifdef STATTB_ENABLE
#ifndef RTSTRC_DISABLE
        if (NOT(f_pptb) || (char *)filename!=o_text[6])
#endif
#endif
           {
           fprintf(f_errr.fp,"%s",e_head);
           for (i=0;i<e_tlvl;i++)
              fprintf(f_errr.fp,"%c",((i%5) ? '.' : '+'));
           fprintf(f_errr.fp,"%s in %s entered.\n",
                   (char *)function,(char *)filename);
           e_tlvl++;
           }
#endif

        /* stack not yet generated                              */
        if (e_bstk==(bentry *)NULL) {
            /* possible pointer allignment problems             */
            e_bstk = e_btop = (bentry *)malloc(sizeof(bentry));
#ifdef HEAP_CHECK
b_geth((a_char *)&e_bstk,(a_char *)e_bstk,(a_char *)"e_push");
#endif
            if (e_btop==(bentry *)NULL)
                {
                fprintf(f_errr.fp,"%se_push : ",e_head);
                fprintf(f_errr.fp,"Insufficient virtual storage\n");
                return;
                }
            else
               {
               e_btop->pred = e_btop->succ = (bentry *)NULL;
               e_btop->line = e_line;
               e_btop->function = (char *)function;
               e_btop->filename = (char *)filename;
               }
            return;
            }

        /* actually no stack entries                            */
        if (e_btop==(bentry *)NULL) e_btop = e_bstk;
        else if (e_btop->succ==(bentry *)NULL)
            {
            /* possible pointer allignment problems             */
            e_btop->succ = (bentry *)malloc(sizeof(bentry));
            if (e_btop->succ==(bentry *)NULL)
                {
                fprintf(f_errr.fp,"%se_push : ",e_head);
                fprintf(f_errr.fp,"Insufficient virtual storage\n");
                return;
                }
#ifdef HEAP_CHECK
b_geth((a_char *)&e_btop->succ,(a_char *)e_btop->succ,(a_char *)"e_push");
#endif
            e_btop->succ->pred = e_btop;
            e_btop->succ->succ = (bentry *)NULL;
            e_btop = e_btop->succ;
            }
        else
            e_btop = e_btop->succ;

        e_btop->line = e_line;
        e_btop->function = (char *)function;
        e_btop->filename = (char *)filename;

        return;
        }





