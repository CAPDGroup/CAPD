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

/* CVS $Id: e_trap.c,v 1.23 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_trap.c                              */
/*                                                              */
/*      Entry           : void e_trap(code,e_argc,e_args)       */
/*                        a_btyp code;                          */
/*                        int e_argc;                           */
/*                        e_list                                */
/*                                                              */
/*      Arguments       : code - error code                     */
/*                        e_argc - number of arguments in list  */
/*                        e_args - variable argument list       */
/*                                                              */
/*      Description     : Trap handler.                         */
/*                                                              */
/*                   remove e_done                              */
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
#if C_P_7
local void e_trap(a_btyp code,int e_argc,...)
#else
local void e_trap(code,e_argc,e_args)
a_btyp code;
int e_argc;
e_list
#endif
#else
local void e_trap(code,e_argc,e_args)
a_btyp code;
int e_argc;
e_list
#endif

        {
        va_list e_argv;
        a_btyp op;
        aentry *ent;
        a_bool found = FALSE;

        e_open(e_argc);

        /* split code value                                     */
        op = code & E_MASK;

        /* Initialize exception handler stack                   */
        if ((ent = e_astk)==(aentry *)NULL)
           {
#if VAX_VMS_C
           e_data();
#endif
           ent = e_astk = &e_anor;
           }

        /* search for exception handler                         */
        while (ent!=(aentry *)NULL)
        if ((ent->code & E_MASK)==op)
           {

           /* active exception handler found                    */
           if (ent->active)
              {

              /* default setting for message/traceback          */
              if (!(code & E_EMSG))
                 {
                 code |= (ent->back) ? E_ETBC+E_EMSG+E_EARG : E_EMSG;
                 }

              /* default setting for exit/continue              */
              if (!(code & (E_ECNT+E_EXIT)))
                 {
                 code |= (ent->cont) ? E_ECNT : E_EXIT;
                 }

              /* invoke  exception handler                      */
              if (ent->action!=NO_TRAP)
                 {
                 (*ent->action)(code,e_argc,e_argv);
                 found = TRUE;
                 }
              else
                 {
                 fprintf(f_errr.fp,"%se_trap : No exception ",e_head);
                 fprintf(f_errr.fp,"handler defined ");
#if GNU_X86_64
                 fprintf(f_errr.fp,"for trap code %8.8x\n",code);
#else
                 fprintf(f_errr.fp,"for trap code %8.8lx\n",code);
#endif
                 }

              /* inactivate special exception handler           */
              if (!(ent->def)) ent->active = FALSE;

              break;
              }
           ent = ent->succ;
           }
        else if ((ent->code & E_MASK)<op)
           ent = ent->succ;
        else
           break;
        /* end of while                                         */

        if (NOT(found))
           {
           e_back(f_errr.fp);
           fprintf(f_errr.fp,"%se_trap : No active exception ",e_head);
           fprintf(f_errr.fp,"handler found ");
#if GNU_X86_64
           fprintf(f_errr.fp,"for error code %8.8x\n",code);
#else
           fprintf(f_errr.fp,"for error code %8.8lx\n",code);
#endif
           }

        e_shut;
        }





