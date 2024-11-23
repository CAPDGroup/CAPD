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

/* CVS $Id: e_bmsg.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_bmsg.c                              */
/*                                                              */
/*      Entries         : void e_bmsg(device)                   */
/*                        FILE *device;                         */
/*                                                              */
/*      Arguments       : device - output device for traceback  */
/*                                 message                      */
/*                                                              */
/*      Description     : print trace back message              */
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
extern char *e_head;
extern int e_line;
extern char *o_text[];
#endif

#ifdef LINT_ARGS
local void e_bmsg(FILE *device)
#else
local void e_bmsg(device)

FILE *device;
#endif
        {
        bentry *ent;

        if ((ent = e_btop)==(bentry *)NULL)
           {
           fprintf(device,"%se_bmsg : No items in trace back stack ",e_head);
           fprintf(device,"available.\n");
           return;
           }

        fprintf(device,"%sERROR",e_head);
        if (e_line>0) fprintf(device," at line %d",e_line);
#ifndef RTSTRC_DISABLE
        do
           {
           if (ent->filename==o_text[6])
              ent = ent->pred;
           else
              {
              if (ent->filename!=NULL)
                 fprintf(device," in '%s'",ent->filename);
              break;
              }
           }
        while (ent!=(bentry *)NULL);
#else
        if (ent->filename!=NULL)
           fprintf(device," in '%s'",ent->filename);
#endif
        fprintf(device,"\n");

        return;
        }





