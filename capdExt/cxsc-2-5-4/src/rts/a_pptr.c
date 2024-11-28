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

/* CVS $Id: a_pptr.c,v 1.21 2014/01/30 17:24:02 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : a_pptr.c                              */
/*                                                              */
/*      Entries         : void a_2pop()                         */
/*                      : void a_2psh(a_VOID new_d,a_VOID old)  */
/*                      : void a_cpsh(a_VOID new_d,a_VOID old)  */
/*                                                              */
/*      Description     : manage stack entries on stack         */
/*                        referenced by a_ptop                  */
/*                                                              */
/*                   free()                                     */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern struct a_pptr *a_ptop;
#endif

#ifdef LINT_ARGS
local void a_2pop(void)
#else
local void a_2pop()
#endif
        {
        a_VOID p;

        /* Neues Display schon installiert? */
        if ( a_ptop->new_d!=NULL )
           {  /* Nein ! */
           p = a_ptop->new_d;
           a_ptop->new_d = NULL;
           }
        else
           {
           /* Ja. Altes Display reinstallieren */
           a_pptr *temp = a_ptop;

           p = a_ptop->old;
           a_ptop = a_ptop->prev;

           /* Stack ruecksetzen                */
#ifdef HEAP_CHECK
b_freh((a_char *)&temp,(a_char *)temp,(a_char *)"a_pptr");
#endif
           free((char*) temp);
           }

        for (; p!=NULL;  p = *(a_VOID*)p)
           **((a_VOID**)p+1) = p;

        return;
        }

#ifdef LINT_ARGS
void a_2psh(a_VOID new_d, a_VOID old)
#else
void a_2psh(new_d, old)
a_VOID new_d, old;
#endif
        {
        a_pptr *temp;

        temp = (a_pptr *) malloc (sizeof(*temp));
#ifdef HEAP_CHECK
b_freh((a_char *)&temp,(a_char *)temp,(a_char *)"a_pptr");
#endif
        temp->new_d = new_d;
        temp->old = old;
        temp->prev = a_ptop;
        a_ptop = temp;

        return;
        }

#ifdef LINT_ARGS
void a_cpsh(a_VOID new_d, a_VOID old)
#else
void a_cpsh(new_d, old)
a_VOID new_d, old;
#endif
        {
        if (new_d!=NULL) a_2psh (new_d, old);

        return;
        }






