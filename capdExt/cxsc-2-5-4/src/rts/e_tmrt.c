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

/* CVS $Id: e_tmrt.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : e_tmrt.c                              */
/*                                                              */
/*      Entry           : void e_tmrt(e_argc,e_argv,print)      */
/*                        int e_argc;                           */
/*                        va_list e_argv;                       */
/*                        a_bool print;                         */
/*                                                              */
/*      Arguments       : e_argc - number of arguments          */
/*                        e_argv - reference to arguments       */
/*                        print - leading messages printed      */
/*                                                              */
/*      Description     : Print leading messages, but do not    */
/*                        print argument list.                  */
/*                                                              */
/*      Note            : e_argc must be even.                  */
/*                        e_argv is a list of pairs representing*/
/*                        integer values for the data type      */
/*                        specification and pointers to the     */
/*                        objects.                              */
/*                        Displaying stops if one result        */
/*                        argument is detected or the list ends.*/
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
extern a_VOID e_rptr;
extern int e_rtyp;
#endif

#ifdef LINT_ARGS
local void e_tmrt(int e_argc,va_list e_argv,a_bool print)
#else
local void e_tmrt(e_argc,e_argv,print)
int e_argc;
va_list e_argv;
a_bool print;
#endif
        {
        int i,type;    /* !!! must be int !!! */

        e_rtyp = 0;
        e_rptr = NULL;

        for (i=0;i<e_argc;)
           {
           /* possible pointer allignment problems              */
           type = e_ref(int);

           /* Process message */
           if (type==E_TMSG && i==0)
              {
              /* possible pointer allignment problems           */
              type = e_ref(int);
              if (print) e_tmsg(type);
              e_argc -= 2;
              continue;
              }

           switch (type & ~(E_TRES+E_TMSG))
              {
              case E_TCHR:
                 /* possible pointer allignment problems        */
                 e_rptr = e_ref(char *);
                 break;
              case E_TDBL:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(a_real *);
                 break;
              case E_TLNG:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(a_long *);
                 break;
              case E_TDTP:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(dotprecision *);
                 break;
              case E_TINT:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(a_intg *);
                 break;
              case E_TMLT:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(multiprecision *);
                 break;
              case E_TSTR:
                 /* possible pointer allignment problems        */
                 e_rptr = e_ref(char *);
                 break;
/*
              case E_TULT:
                 e_rptr = e_ref(ultraprecision *);
                 break;
*/
              case E_TSTG:
                 /* possible pointer allignment problems        */
                 e_rptr = (char *)e_ref(s_trng *);
                 break;
              }
           i += 2;

           /* Display stops if result */
           if ((type & E_TRES)!=0)
              {
              e_rtyp = type ^ E_TRES;
              return;
              }
           }

        return;
        }





