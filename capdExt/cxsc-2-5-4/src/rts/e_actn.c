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

/* CVS $Id: e_actn.c,v 1.23 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_actn.c                              */
/*                                                              */
/*      Entry           : void e_actn(action,code,function)     */
/*                        a_btyp action,code;                   */
/*                        void (*function)( );                  */
/*                                                              */
/*      Arguments       : code - error code                     */
/*                        action - action to be performed       */
/*                        function - trap function              */
/*                                                              */
/*      Description     : Action stack handler.                 */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern aentry *e_astk;
extern aentry e_anor;
extern char *e_head;
extern f_text f_errr;
#endif

#ifdef LINT_ARGS
local void e_actn(a_btyp action,a_btyp code,
                  void (*function)(a_btyp,int,va_list))
#else
local void e_actn(action,code,function)

a_btyp action;
a_btyp code;
void (*function)();
#endif
        {
        aentry *pnew,*p,*pold;
        a_btyp op;

        E_TPUSH("e_actn");

        op = code & E_MASK;

        /* Initialize exception handler stack                   */
        if (e_astk==(aentry *)NULL)
           {
#if VAX_VMS_C
           e_data();
#endif
           e_astk = &e_anor;
           }

/* This case has been withdrawn. >>>
        if (e_astk==(aentry *)NULL)
           {
           if (action & (E_CHNG | E_POPP))
              {
              fprintf(f_errr.fp,"%se_actn : No exception ",e_head);
              fprintf(f_errr.fp,"handler found ");
              fprintf(f_errr.fp,"for trap code %08.8lx\n",code);
              E_TPOPP("e_actn")
              return;
              }
           else if (action & E_PUSH)
              {
              p = e_astk = (aentry *)malloc(sizeof(aentry));
              if (p==(aentry *)NULL)
                 {
                 e_trap(ALLOCATION,2,E_TMSG,41);
                 E_TPOPP("e_actn")
                 return;
                 }
#ifdef HEAP_CHECK
b_geth((a_char *)&e_astk,(a_char *)e_astk,(a_char *)"e_actn");
#endif
              p->succ = NULL;
              p->stat = FALSE;
              }
           else
              {
              fprintf(f_errr.fp,"%se_actn : No action specified ",e_head);
              fprintf(f_errr.fp,"for trap code %08.8lx\n",code);
              E_TPOPP("e_actn")
              return;
              }
           }
        else
 <<< This case has been withdrawn. */
           {
           pold = (aentry *)NULL;
           p = e_astk;
           while ((p->code & E_MASK)<op) {
              if (p->succ==(aentry *)NULL) break;
              pold = p;
              p = p->succ;
              }

           if ((p->code & E_MASK)!=op)
              {
              if (action & (E_CHNG | E_POPP))
                 {
                 fprintf(f_errr.fp,"%se_actn : No exception ",e_head);
                 fprintf(f_errr.fp,"handler found ");
#if GNU_X86_64
                 fprintf(f_errr.fp,"for trap code %8.8x\n",code);
#else
                 fprintf(f_errr.fp,"for trap code %8.8lx\n",code);
#endif
                 E_TPOPP("e_actn")
                 return;
                 }
              }

           if (action & E_POPP)
              {
              if (p->stat)
                 {
                 fprintf(f_errr.fp,"%se_actn : Static exception ",e_head);
                 fprintf(f_errr.fp,"handler may not be ");
#if GNU_X86_64
                 fprintf(f_errr.fp,"removed for trap code %8.8x\n",code);
#else
                 fprintf(f_errr.fp,"removed for trap code %8.8lx\n",code);
#endif
                 E_TPOPP("e_actn")
                 return;
                 }
              if (pold==(aentry *)NULL)
                 {
                 e_astk = p->succ;
#ifdef HEAP_CHECK
b_freh((a_char *)&e_astk,(a_char *)p,(a_char *)"e_actn");
#endif
                 }
              else
                 {
                 pold->succ = p->succ;
#ifdef HEAP_CHECK
b_freh((a_char *)&pold->succ,(a_char *)p,(a_char *)"e_actn");
#endif
                 }
              free((char *)p);
              E_TPOPP("e_actn")
              return;
              }

           if (action & E_PUSH)
              {
              /* activate exception handler if inactive special */
              if ((p->code & E_MASK)==op && !(p->active)
                  && !(p->def))
                 {
                 p->active = TRUE;
                 E_TPOPP("e_actn")
                 return;
                 }

              /* possible pointer allignment problems           */
              pnew = (aentry *)malloc(sizeof(aentry));
              if (pnew==(aentry *)NULL)
                 {
                 e_trap(ALLOCATION,2,E_TMSG,41);
                 E_TPOPP("e_actn")
                 return;
                 }
#ifdef HEAP_CHECK
b_geth((a_char *)&pold->succ,(a_char *)pnew,(a_char *)"e_actn");
#endif
              pnew->succ = pold->succ;
              pold->succ = pnew;
              p = pnew;
              p->stat = FALSE;
              }
           else if (!(action & E_CHNG))
              {
              fprintf(f_errr.fp,"%se_actn : No exception handler ",e_head);
#if GNU_X86_64
              fprintf(f_errr.fp,"found for trap code %8.8x\n",code);
#else
              fprintf(f_errr.fp,"found for trap code %8.8lx\n",code);
#endif
              E_TPOPP("e_actn")
              return;
              }
           }

        p->code = code;
        if (function!=NO_TRAP) p->action = function;
        p->active = (action & E_ACTIVE) ? TRUE : FALSE;
        p->back = (action & E_BACK) ? TRUE : FALSE;
        p->cont = (action & E_CONT) ? TRUE : FALSE;
        p->def = (action & E_DEFAULT) ? TRUE : FALSE;

        E_TPOPP("e_actn")
        }





