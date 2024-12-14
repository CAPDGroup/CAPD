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

/* CVS $Id: e_popp.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_popp.c                              */
/*                                                              */
/*      Entries         : void e_popp                           */
/*                                                              */
/*      Description     : pop trace back entry from stack       */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern bentry *e_btop;
extern int e_line;
#ifdef TRACEC_ENABLE
extern a_bool f_pptr;
extern int e_tlvl;
#endif
extern f_text f_errr;
extern char *e_head;
#ifdef STATTB_ENABLE
extern a_bool f_pptb;
#endif
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local void e_popp(void)
#else
local void e_popp()
#endif
        {
#ifdef TRACEC_ENABLE
        int i;

        if (f_pptr)
#ifdef STATTB_ENABLE
#ifndef RTSTRC_DISABLE
        if (NOT(f_pptb) || e_btop->filename!=o_text[6])
#endif
#endif
           {
           e_tlvl--;
           fprintf(f_errr.fp,"%s",e_head);
           for (i=0;i<e_tlvl;i++)
              fprintf(f_errr.fp,"%c",((i%5) ? '.' : '+'));
           fprintf(f_errr.fp,"%s in %s terminated.\n",
                   e_btop->function,e_btop->filename);
           }
#endif

        if (e_btop==(bentry *)NULL) return;

        /* restore global variable e_line                       */
        e_line = e_btop->line;

        /* pop top element from traceback stack                 */
        e_btop = e_btop->pred;

        return;
        }





