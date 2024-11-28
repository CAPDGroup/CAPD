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

/* CVS $Id: t_merr.c,v 1.21 2014/01/30 17:24:17 cxsc Exp $ */


#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_ieee.h"
#else
#include "t_ieee.h"
#endif

#ifdef LINT_ARGS
#ifdef ANSI_C
static int err_ext_to_int(const ExtReal *arg, ExtReal *res);
#else  /* NOT ANSI_C */
static int err_ext_to_int(ExtReal *arg, ExtReal *res);
#endif /* ANSI_C */
#else  /* NOT LINT_ARGS */
static int err_ext_to_int();
#endif /* LINT_ARGS */

/*-------------------------------------------------------------------*
 | ieee err                                                          |
 | vom Anwender erweiterbar                                          |
 |                                                                   |
 | falls return == TakeOver wird der Wert fuer x->res als            |
 | Funktionswert und                                                 |
 | x->type als Rueckgabewert der Funktion uebernommen                |
 | es tritt dann kein Standard Fehlerbehandlung in Aktion            |
 | x->type = -1 ==> x->res wird direkt ohne Fehlermeldung als        |
 |                  Funktionswert verwendet                          |
 | x->type = 0  ==> nicht unbedingt sinnvoll, unterbrochene Funktion |
 |                  wird mit dem Originalargument ausgefuehrt        |
 | x->type > 0  ==> es erfolgt eine Fehlermeldung, in den makehead-  |
 |                  Rahmen erfolgt Programmabbruch                   |
 *-------------------------------------------------------------------*/
#ifdef LINT_ARGS
int ieee_matherr(struct ieee_exception *x)
#else
int ieee_matherr(x)
struct ieee_exception *x;
#endif /* LINT_ARGS */
{
/*   printf("von  ieee_matherr()\n"); */

   switch(x->type) {
   case DOMAIN:
      if(0==strcmp(x->name, "extreal_to_int"))
         return err_ext_to_int(&(x->arg1->u), &(x->res->u));
      break;
   case SING:
      break;
   case OVER_FLOW:
      break;
   case PLOSS:
      break;
   case TLOSS:
      break;
   case UNDER_FLOW:
      if(0==strcmp(x->name, "exp"))
      {
        copyee(&Zero, &(x->res->u));
        x->type = -1;
        return(TakeOver);
      }
      break;
   default:
      break;
   }

/*   printf("Ende ieee_matherr()\n\n"); */
   /* --- keine eigene Massnahme getroffen --- */
   return NoTakeOver;
} /* ieee_err() */

/*--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
static int err_ext_to_int(const ExtReal *arg, ExtReal *res)
#else
static int err_ext_to_int(arg, res)
const ExtReal *arg;
      ExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
static int err_ext_to_int(ExtReal *arg, ExtReal *res)
#else
static int err_ext_to_int(arg, res)
ExtReal *arg;
ExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   extern ExtReal IntMax;
   extern ExtReal IntMin;

   if(SGNE(arg)==NEG)
      copyee(&IntMin, res);
   else
      copyee(&IntMax, res);
   return TakeOver;
} /* err_ext_to_int() */





