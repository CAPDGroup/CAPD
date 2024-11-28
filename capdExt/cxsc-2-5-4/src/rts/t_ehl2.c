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

/* CVS $Id: t_ehl2.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
#include <stdio.h>

/*--------------------------------------------------------------*
 | Exception Handle PunktFunktionen mit 2 Argumenten            |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int exc_handle_2(int fct, int exctype, const ExtReal *arg1,
                 const ExtReal *arg2, ExtReal *res)
#else
int exc_handle_2(fct, exctype, arg1, arg2, res)
      int      fct;        /* FunktionsNummer                   */
      int      exctype;    /* ExceptionNummer                   */
const ExtReal  *arg1;
const ExtReal  *arg2;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int exc_handle_2(int fct, int exctype, ExtReal *arg1,
                 ExtReal *arg2, ExtReal *res)
#else
int exc_handle_2(fct, exctype, arg1, arg2, res)
int      fct;              /* FunktionsNummer                   */
int      exctype;          /* ExceptionNummer                   */
ExtReal  *arg1;
ExtReal  *arg2;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   struct ieee_exception exc;
   IExtReal              ia1;
   IExtReal              ia2;
   IExtReal              ir;
   char                  *name;  /* NamensString der Funktion   */

   /* --- kopiere, falls arg==res --- */
   copyee(res, &ir.u);

   /* --- hole Standard Antwort und FunktionsName --- */
   default_handle_2(fct, exctype, arg1, arg2, &ir.u);
   excfct_to_a(fct, &name);

   /* --- matherr vorbereiten --- */
   copyee(arg1,   &ia1.u);
   copyee(&Zero,  &ia1.l);
   copyee(arg2,   &ia2.u);
   copyee(&Zero,  &ia2.l);
   copyee(&Zero,  &ir.l);

   exc.type = exctype;
   exc.name = name;
   exc.arg1 = &ia1;
   exc.arg2 = &ia2;
   exc.res  = &ir;

   /* --- Aufruf AnwenderFunktion --- */
   if (TakeOver==ieee_matherr(&exc)) {

      /* --- Uebernahme vom Anwender gesetzter Werte, Ende --- */
      copyee(&(exc.res->u), res);
      return exc.type;
   }

   /* --- sonst Ausgabe der unveraenderten FunktionsParameter --- */
   msg_exctyp(exc.type, name);
/*   msg_exc("arg1", arg1);
   msg_exc("arg2", arg2);
   msg_exc("res ", &ir.u); */

   /* --- kopiere default Antwort nach res (erst hier, falls arg==res) --- */
   copyee(&ir.u, res);

   return exctype;
} /* except_2() */





