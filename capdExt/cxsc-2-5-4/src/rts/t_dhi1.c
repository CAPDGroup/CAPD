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

/* CVS $Id: t_dhi1.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */

/*      ISin,ICos cases considered */

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | Default Handle IntervallFunktionen mit einem Argument        |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int default_handle_i1(int fct, int exc, const IExtReal *arg, IExtReal *res)
#else
int default_handle_i1(fct, exc, arg, res)
      int      fct;
      int      exc;
const IExtReal *arg;
      IExtReal *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int default_handle_i1(int fct, int exc, IExtReal *arg, IExtReal *res)
#else
int default_handle_i1(fct, exc, arg, res)
int      fct;
int      exc;
IExtReal *arg;
IExtReal *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   switch(fct) {
   case  ITan:
   case  ICot:
   case  ICoth:
      if(exc==ExcISing) {
         copyee(&PInfty, &(res->u));
         copyee(&MInfty, &(res->l));
      }
      else {
         copyee(&Zero, &(res->u));
         copyee(&Zero, &(res->l));
      }
      break;
   case  IAtanh:
      if(exc==DOMAIN||exc==ExcISing) {
         if(-1!=cmpee(&(arg->u), &One))
            copyee(&PInfty, &(res->u));
         if(-1!=cmpee(&(arg->l), &One))
            copyee(&PInfty, &(res->l));
         if(1!=cmpee(&(arg->u), &MinusOne))
            copyee(&MInfty, &(res->u));
         if(1!=cmpee(&(arg->l), &MinusOne))
            copyee(&MInfty, &(res->l));
      }
      else {
         copyee(&Zero, &(res->u));
         copyee(&Zero, &(res->l));
      }
      break;
   case  IAcoth:
      if(exc==DOMAIN||exc==ExcISing) {
         if(1!=cmpabsee(&(arg->u), &One))
            copyee(&PInfty, &(res->u));
         if(1!=cmpabsee(&(arg->l), &One))
            copyee(&MInfty, &(res->l));
      }
      else {
         copyee(&Zero, &(res->u));
         copyee(&Zero, &(res->l));
      }
      break;
   case  IExp:
      if(exc==OVER_FLOW) {
         copyee(&PInfty, &(res->u));
         copyee(&PInfty, &(res->l));
      }
      else {
         copyee(&Zero, &(res->u));
         copyee(&Zero, &(res->l));
      }
      break;
   case  ICos:
   case  ISin:
      copyee(&One, &(res->u));
      copyee(&MinusOne, &(res->l));
      break;
   default:
      copyee(&Zero, &(res->u));
      copyee(&Zero, &(res->l));
      break;
   } /* switch(fct) */

   /* --- kein Fehler moeglich --- */
   return NoErr;
} /* default_handle_i1() */





