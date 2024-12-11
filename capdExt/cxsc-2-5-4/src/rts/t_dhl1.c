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

/* CVS $Id: t_dhl1.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */


#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

/*--------------------------------------------------------------*
 | Default Handle PunktFunktionen mit einem Argument            |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int default_handle_1(int fct, int exc, const ExtReal *arg, ExtReal *res)
#else
int default_handle_1(fct, exc, arg, res)
      int      fct;
      int      exc;
const ExtReal  *arg;
      ExtReal  *res;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int default_handle_1(int fct, int exc, ExtReal *arg, ExtReal *res)
#else
int default_handle_1(fct, exc, arg, res)
int      fct;
int      exc;
ExtReal  *arg;
ExtReal  *res;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   switch(fct) {
   case  Exp:
      if(exc==OVER_FLOW)
         copyee(&PInfty, res);
      else
         copyee(&Zero, res);
      break;
   case  EToI:
      if(exc==OVER_FLOW)
         copyee(SGNE(arg)==POS ? &IntMax : &IntMin, res);
      else
         copyee(&Zero, res);
      break;
   case  EToL:
      switch(exc){
      case ExcPInf:
         copyee(&PInfty, res);
         break;
      case OVER_FLOW:
         copyee(&LongRealMax, res);
         break;
      case UNDER_FLOW:
         copyee(&LongRealDenormMin, res);
         break;
      case ExcMInf:
         copyee(&MInfty, res);
         break;
      default:
         copyee(&Zero, res);
         break;
      }
      break;
   default:
      copyee(&Zero, res);
      break;
   } /* switch(fct) */

   /* --- kein Fehler moeglich --- */
   return NoErr;
} /* default_handle_1() */





