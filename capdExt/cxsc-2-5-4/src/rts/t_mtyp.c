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

/* CVS $Id: t_mtyp.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
#include <stdio.h>

/*--------------------------------------------------------------*
 | message exception typ and function name                      |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
void msg_exctyp(int typ, const char *name)
#else
void msg_exctyp(typ, name)
      int   typ;
const char  *name;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
void msg_exctyp(int typ, char *name)
#else
void msg_exctyp(typ, name)
int   typ;
char  *name;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   char  *s;

   fprintf(MsgOut,"ieee math error ");
   exc_to_a(typ, &s);
   fprintf(MsgOut, "%s ", s);
   fprintf(MsgOut,"in function %s\n",name);
} /* msg_exctyp() */

/*---------------------------------------------------------------------*
 | message exception arguments                                         |
 *---------------------------------------------------------------------*/
/* BEGIN OF msg_exc ---->
#ifdef ANSI_C
#ifdef LINT_ARGS
void msg_exc(const char *name, const ExtReal *arg)
#else
void msg_exc(name, arg)
const char     *name;
const ExtReal  *arg;
#endif
#else
#ifdef LINT_ARGS
void msg_exc(char *name, ExtReal *arg)
#else
void msg_exc(name, arg)
char     *name;
ExtReal  *arg;
#endif
#endif
{
   fprintf(MsgOut, "%s >", name);
   fprinte(MsgOut, arg);
   fprintf(MsgOut,"<\n");
}   <---- END-OF  msg_exc() */





