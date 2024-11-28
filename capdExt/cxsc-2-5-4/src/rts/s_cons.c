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

/* CVS $Id: s_cons.c,v 1.21 2014/01/30 17:24:13 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : s_cons.c                              */
/*                                                              */
/*      Entries         : a_VOID s_cons(s,e_argc,e_args)        */
/*                        s_etof s;                             */
/*                        a_char *e_argc;                       */
/*                        e_list                                */
/*                                                              */
/*      Arguments       : s = set to be constructed             */
/*                        e_argc = control string               */
/*                        e_args = variable length argument list*/
/*                                                              */
/*      Function value  : pointer to result set                 */
/*                                                              */
/*      Description     : Construct set.                        */
/*                        The control string "e_argc" consists  */
/*                        of a sequence of "1" or "2":          */
/*                        "1" : add one element to the set      */
/*                        "2" : add a contiguos range of        */
/*                              elements to the set             */
/*                        The element and the element range     */
/*                        are formed by the sequence of variable*/
/*                        length arguments "e_args".            */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
#if C_P_7
local a_VOID s_cons(s_etof s,a_char *e_argc,...)
#else
local a_VOID s_cons(s,e_argc,e_args)
s_etof s;
a_char *e_argc;
e_list
#endif
#else
local a_VOID s_cons(s,e_argc,e_args)
s_etof s;
a_char *e_argc;
e_list
#endif

        {
        va_list e_argv;
        register a_intg item1,item2;

        E_TPUSH("s_cons")

        e_open(e_argc);

        /* clear set */
        (void)memset(&s[0],0x00,s_SIZE);

        while (*e_argc)
           {
           switch(*e_argc++)
              {
              case '1'  :
                          item1 = e_ref(a_intg);
                          (void)s_ins1(s,item1);
                          break;
              case '2'  :
                          item1 = e_ref(a_intg);
                          item2 = e_ref(a_intg);
                          (void)s_ins2(s,item1,item2);
                          break;
              default   : ;
              }
           }

        e_shut;

        E_TPOPP("s_cons")
        return((a_VOID)&s[0]);
        }





