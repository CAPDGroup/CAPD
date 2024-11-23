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

/* CVS $Id: b_heap.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_heap.c                              */
/*                                                              */
/*      Entry           : void b_geth(ptr,addr,where)           */
/*      Entry           : void b_freh(ptr,addr,where)           */
/*      Entry           : void b_tmph(ptr)                      */
/*      Entry           : void b_varh(ptr,addr)                 */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#define local static
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#undef local
#define local
#endif

typedef struct listelement { a_char *ptr;
                             a_char *addr;
                             a_char *where;  } listelement;

#define HEAP_SIZE       1000
static listelement list[HEAP_SIZE];
static int used = 0;
static int count = 0;

extern f_text f_errr;

/****************************************************************/
/*                                                              */
/*      Filename        : b_heap.c                              */
/*                                                              */
/*      Entry           : void b_geth(ptr,addr,where)           */
/*                        a_char *ptr;                          */
/*                        a_char *addr;                         */
/*                        a_char *where;                        */
/*                                                              */
/*      Arguments       : ptr    - address of variable holding  */
/*                                 heap address                 */
/*                        addr   - address of allocated heap    */
/*                        size   - size of allocated heap       */
/*                        where  - routine where allocation     */
/*                                 is done                      */
/*                                                              */
/*      Description     : the address of allocated heap is      */
/*                        stored in a list.                     */
/*                                                              */
/*      Note            : This routine is used for debugging    */
/*                        purposes only.                        */
/*                        Up to HEAP_SIZE addresses may be      */
/*                        stored.                               */
/*                                                              */
/****************************************************************/

#ifdef LINT_ARGS
local void b_geth(a_char *ptr,a_char *addr,a_char *where)
#else
local void b_geth(ptr,addr,where)

a_char *ptr;
a_char *addr;
a_char *where;
#endif

        {
        int i;
        int unused;

        /* search for ptr */
        unused = used;
        for (i=0;i<used;i++)
           {

           /* the variable has already been assigned a heap address */
           if (list[i].ptr==ptr)
              {
      fprintf(f_errr.fp,"--------------------------------------\n");
      fprintf(f_errr.fp,"--- Reassignment of allocated heap to\n");
      fprintf(f_errr.fp,"--- variable at address: %p\n",ptr);
      fprintf(f_errr.fp,"--- Previous allocation\n");
      fprintf(f_errr.fp,"---    in routine '%s'\n",list[i].where);
      fprintf(f_errr.fp,"---    heap addr = %p\n",list[i].addr);
      fprintf(f_errr.fp,"--- Actual allocation\n");
      fprintf(f_errr.fp,"---    in routine '%s'\n",where);
      fprintf(f_errr.fp,"---    heap addr = %p\n",addr);
      fprintf(f_errr.fp,"--------------------------------------\n");
              exit(0);
              }

           /* unused position found in list */
           else if (unused==used && list[i].where==NULL)
              {
              unused = i;
              }
           }

        /* no space left in list */
        if (unused==HEAP_SIZE)
           {
           fprintf(f_errr.fp,"-----------------------------------\n");
           fprintf(f_errr.fp,"--- Insufficient HEAP_CHECK space\n");
           fprintf(f_errr.fp,"-----------------------------------\n");
           }

        /* new list element */
        else
           {
           list[unused].ptr = ptr;
           list[unused].addr = addr;
           list[unused].where = where;
           count++;

fprintf(stdout,"(%3d) : inserted element(%3d) = %p %p %s\n",
count,unused,ptr,addr,where);

           if (unused==used) used++;
           }

        return;
        }

/****************************************************************/
/*                                                              */
/*      Filename        : b_heap.c                              */
/*                                                              */
/*      Entry           : void b_freh(ptr,addr,where)           */
/*                        a_char *ptr;                          */
/*                        a_char *addr;                         */
/*                        a_char *where;                        */
/*                                                              */
/*      Arguments       : ptr    - address of variable holding  */
/*                                 heap address                 */
/*                        addr   - address of allocated heap    */
/*                        where  - routine where free is done   */
/*                                                              */
/*      Description     : the address of allocated heap is      */
/*                        removed from a list.                  */
/*                                                              */
/*      Note            : This routine is used for debugging    */
/*                        purposes only.                        */
/*                                                              */
/****************************************************************/

#ifdef LINT_ARGS
local void b_freh(a_char *ptr,a_char *addr,a_char *where)
#else
local void b_freh(ptr,addr,where)

a_char *ptr;
a_char *addr;
a_char *where;
#endif

        {
        int i;

        /* variable already freed */
        if (addr==NULL)
           {
           return;
           }

        /* search for addr */
        for (i=0;i<used;i++)
           {

           /* addr found */
           if (list[i].addr==addr)
              {

fprintf(stdout,"(%3d) :  removed element(%3d) = %p %p %s\n",
count,i,ptr,addr,where);

              /* remove list element */
              list[i].ptr = list[i].addr = list[i].where = NULL;
              count--;
              if (i==used-1)
                 {
                 do
                    used--;
                 while (used>0 && list[used-1].where==NULL);
                 }
              return;
              }
           }

        /* ptr/addr not found */
        fprintf(f_errr.fp,"-----------------------------------\n");
        fprintf(f_errr.fp,"--- Attempt to free unknown heap\n");
        fprintf(f_errr.fp,"--- in routine '%s'\n",where);
        fprintf(f_errr.fp,"---    Heap address = %p\n",addr);
        fprintf(f_errr.fp,"---    Variable address = %p\n",ptr);
        fprintf(f_errr.fp,"-----------------------------------\n");
        exit(0);

        return;
        }

/****************************************************************/
/*                                                              */
/*      Filename        : b_heap.c                              */
/*                                                              */
/*      Entry           : void b_tmph(ptr)                      */
/*                        a_char *ptr;                          */
/*                                                              */
/*      Arguments       : ptr    - address of variable holding  */
/*                                 heap address                 */
/*                                                              */
/*      Description     : the address of a variable is          */
/*                        removed from a list.                  */
/*                                                              */
/*      Note            : This routine is used for debugging    */
/*                        purposes only.                        */
/*                                                              */
/****************************************************************/

#ifdef LINT_ARGS
local void b_tmph(a_char *ptr)
#else
local void b_tmph(ptr)

a_char *ptr;
#endif

        {
        int i;

        /* search for ptr */
        for (i=0;i<used;i++)
           {

           /* addr found */
           if (list[i].ptr==ptr)
              {

fprintf(stdout,"(%3d) :     return value(%3d) = %p\n",count,i,ptr);

              /* mark list element to be a return value */
              list[i].ptr = NULL;
              return;
              }
           }

        return;
        }

/****************************************************************/
/*                                                              */
/*      Filename        : b_heap.c                              */
/*                                                              */
/*      Entry           : void b_varh(ptr,addr)                 */
/*                        a_char *ptr;                          */
/*                        a_char *addr;                         */
/*                                                              */
/*      Arguments       : ptr    - address of variable holding  */
/*                                 heap address                 */
/*                        addr   - address of allocated heap    */
/*                                                              */
/*      Description     : the address of a variable is          */
/*                        stored in a list.                     */
/*                                                              */
/*      Note            : This routine is used for debugging    */
/*                        purposes only.                        */
/*                                                              */
/****************************************************************/

#ifdef LINT_ARGS
local void b_varh(a_char *ptr,a_char *addr)
#else
local void b_varh(ptr,addr)

a_char *ptr;
a_char *addr;
#endif

        {
        int i;

        /* search for addr */
        for (i=0;i<used;i++)
           {

           /* addr found */
           if (list[i].addr==addr)
              {

fprintf(stdout,"(%3d) :         variable(%3d) = %p %p\n",
count,i,ptr,addr);

              /* mark list element to be a variable */
              list[i].ptr = ptr;
              return;
              }
           }

        return;
        }





