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

/* CVS $Id: t_chk1.c,v 1.21 2014/01/30 17:24:15 cxsc Exp $ */


#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif

#define ReturnErr(A)    \
   if (NoErr!=(exc=A)) return exc

/*--------------------------------------------------------------*
 | Pruefe Argument                                              |
 *--------------------------------------------------------------*/
#ifdef ANSI_C
#ifdef LINT_ARGS
int chk_1(int fct, const ExtReal *arg)
#else
int chk_1(fct, arg)
int fct;                /* FunktionsNummer */
const ExtReal *arg;
#endif /* LINT_ARGS */
#else  /* NOT ANSI_C */
#ifdef LINT_ARGS
int chk_1(int fct, ExtReal *arg)
#else
int chk_1(fct, arg)
int fct;                /* FunktionsNummer */
ExtReal *arg;
#endif /* LINT_ARGS */
#endif /* ANSI_C */
{
   int         exc;

   /* --- ArgumentPruefung nach Funktionen --- */
   switch(fct) {
   case Sin  :
   case Cos  :
   case Tan  :
   case Cot  :
   case ISin :
   case ICos :
   case ITan :
   case ICot :
      ReturnErr (chk_extreal(arg, DBExtReal));
/* |arg|>2**127*/
      ReturnErr ((1==(cmpabsee(arg, &TwoPow127))?TLOSS:NoErr));
/* |arg|>2**63 */
      ReturnErr ((1==(cmpabsee(arg, &TwoPow63 ))?PLOSS:NoErr));
         return NoErr;
   case Asin :
   case Acos :
   case IAsin:
   case IAcos:
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((1==(cmpabsee(arg, &One))?DOMAIN:NoErr));/* |arg|>1*/
         return NoErr;
   case Atan :
   case Acot :
   case IAtan:
   case IAcot:
      ReturnErr (chk_extreal(arg, DBExtReal));
         return NoErr;

   case Sqrt :
   case ISqrt:
      ReturnErr (chk_extreal(arg, DBPos0ExtReal));
         return NoErr;
   case Sqrtm1:
      ReturnErr (chk_extreal(arg, DBPos0ExtReal));
      ReturnErr ((-1!=(cmpee(arg, &MaxArgSqrtm1))?DOMAIN:NoErr));
         return NoErr;
   case Sinh :
   case Cosh :
   case ISinh:
   case ICosh:
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((1==(cmpee(arg, &MaxArgHyp))?OVER_FLOW:NoErr));
      ReturnErr ((-1==(cmpee(arg, &MinArgHyp))?OVER_FLOW:NoErr));
         return NoErr;
   case Exp  :
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((1==(cmpee(arg, &MaxArgExp))?OVER_FLOW:NoErr));
      ReturnErr ((-1==(cmpee(arg, &MinArgExp))?UNDER_FLOW:NoErr));
         return NoErr;
   case IExp :
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((1==(cmpee(arg, &MaxArgExp))?OVER_FLOW:NoErr));
         return NoErr;
   case Coth :
   case ICoth:
      ReturnErr (chk_extreal(arg, DBPosExtReal|DBNegExtReal));
         return NoErr;
   case Tanh :
   case ITanh:
      ReturnErr (chk_extreal(arg, DBExtReal));
         return NoErr;
   case Expm1:
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((1==(cmpee(arg, &Ln2))?DOMAIN:NoErr));
         return NoErr;           /* Underflow unmoeglich        */

   case Ln   :
   case ILn  :
      ReturnErr (chk_extreal(arg, DBPosExtReal));
         return NoErr;
   case Lnp1 :
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((-1!=(cmpabsee(arg, &MaxArgLnp1))?DOMAIN:NoErr));
         return NoErr;
   case Asinh:
   case IAsinh:
      ReturnErr (chk_extreal(arg, DBExtReal));
         return NoErr;
   case Acosh:
   case IAcosh:
      ReturnErr (chk_extreal(arg, DBPosExtReal));
      ReturnErr ((-1==(cmpabsee(arg, &One))?DOMAIN:NoErr));/* |arg|<1*/
         return NoErr;
   case Atanh:
   case IAtanh:
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((-1!=(cmpabsee(arg, &One))?DOMAIN:NoErr));/* |arg|>=1*/
         return NoErr;
   case Acoth:
   case IAcoth:
      ReturnErr (chk_extreal(arg, DBExtReal));
      ReturnErr ((1!=(cmpabsee(arg, &One))?DOMAIN:NoErr));/* |arg|<=1*/
         return NoErr;

   case Pow:
   case IPow:
      ReturnErr (chk_extreal(arg, DBExtReal));
         return NoErr;

   case Div:                    /* 010291 cb */
   case IDiv:
      ReturnErr (chk_extreal(arg, DBExtReal));
         return NoErr;

   /* ------------------------------------------------------------- */
   case EToI :
      ReturnErr (chk_extreal(arg, DBExtReal));
/* |arg|>INT_MAX */
      ReturnErr ((1==(cmpabsee(arg, &IntMax))?OVER_FLOW:NoErr));
         return NoErr;

   case EToL :
      ReturnErr (chk_extreal(arg, DBExtReal));
/* |arg|>DBL_MAX */
      ReturnErr (1==cmpabsee(arg, &LongRealMax)?OVER_FLOW:NoErr);
         return NoErr;

   default: fprintf(MsgOut, "Fehler ArgumentCheck!\n");
         return ExcUnknown;
   } /* switch(fct) */

} /* chk_1() */





