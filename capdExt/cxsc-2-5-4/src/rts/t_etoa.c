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

/* CVS $Id: t_etoa.c,v 1.21 2014/01/30 17:24:16 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : t_etoa.c                              */
/*                                                              */
/****************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/tbyte/t_defs.h"
#else
#include "t_defs.h"
#endif
#include <string.h>

/*--------------------------------------------------------------*
 | exception to ascii string                                    |
 *--------------------------------------------------------------*/
#ifdef LINT_ARGS
void exc_to_a(int exc, char **str)
#else
void exc_to_a(exc, str)
int   exc;
char  **str;
#endif /* LINT_ARGS */
{
   static char s[20];
   /*
   char s1[20];
   */

   switch(exc) {
   case DOMAIN:     *str = "argument domain";       break;
   case SING:       *str = "singularity";           break;
   case OVER_FLOW:  *str = "overflow";              break;
   case UNDER_FLOW: *str = "underflow";             break;
   case PLOSS:      *str = "partial loss of precision"; break;
   case TLOSS:      *str = "total loss of precision";   break;
   case ExcNoI:     *str = "no interval";           break;
   case ExcUnknown: *str = "unknown";               break;
   case ExcPInf:    *str = "+Infinity";             break;
   case ExcMInf:    *str = "-Infinity";             break;
   case ExcPZero:   *str = "+Zero";                 break;
   case ExcMZero:   *str = "-Zero";                 break;
   case ExcPNorm:   *str = "+Normal";               break;
   case ExcMNorm:   *str = "-Normal";               break;
   case ExcPDenorm: *str = "+Denorm";               break;
   case ExcMDenorm: *str = "-Denorm";               break;
   case ExcPNAN:    *str = "not a number (+)";      break;
   case ExcMNAN:    *str = "not a number (-)";      break;
   case ExcInvalid: *str = "invalid";               break;
   case ExcISing:   *str = "interval singularity";  break;
   case ExcDBZ:     *str = "division by zero";      break;
   case ExcDIZ:     *str = "division by an interval containing zero"; break;
   default:
   /* itoa(exc, s1, 10);
      strcpy(s, "exc no. ");
      strcat(s, s1); */
      sprintf(s, "exc no. %d", exc);
      *str = s;
      break;
   } /* switch(exc) */

   return;
} /* exc_to_a() */

/*--------------------------------------------------------------*
 | FunktionsNummer nach ascii string                            |
 *--------------------------------------------------------------*/
#ifdef LINT_ARGS
void excfct_to_a(int fct, char **str)
#else
void excfct_to_a(fct, str)
int   fct;
char  **str;
#endif /* LINT_ARGS */
{
   static char s[20];
   /*
   char s1[20];
   */

   switch(fct) {
   case Sqrt:       *str = "sqrt";   break;
   case Sqrtm1:     *str = "sqrtm1"; break;
   case ISqrt:      *str = "isqrt";  break;

   case Exp :       *str = "exp";    break;
   case IExp:       *str = "iexp";   break;
   case Expm1:      *str = "expm1";  break;
   case Ln  :       *str = "ln";     break;
   case ILn :       *str = "iln";    break;
   case Lnp1:       *str = "lnp1";   break;

   case Sin  :      *str = "sin";    break;
   case Cos  :      *str = "cos";    break;
   case Tan  :      *str = "tan";    break;
   case Cot  :      *str = "cot";    break;
   case ISin :      *str = "isin";   break;
   case ICos :      *str = "icos";   break;
   case ITan :      *str = "itan";   break;
   case ICot :      *str = "icot";   break;

   case Asin:       *str = "asin";   break;
   case Acos:       *str = "acos";   break;
   case Atan:       *str = "atan";   break;
   case Acot:       *str = "acot";   break;
   case IAsin:      *str = "iasin";  break;
   case IAcos:      *str = "iacos";  break;
   case IAtan:      *str = "iatan";  break;
   case IAcot:      *str = "iacot";  break;

   case Sinh :      *str = "sinh";   break;
   case Cosh :      *str = "cosh";   break;
   case Tanh :      *str = "tanh";   break;
   case Coth :      *str = "coth";   break;
   case ISinh:      *str = "isinh";  break;
   case ICosh:      *str = "icosh";  break;
   case ITanh:      *str = "itanh";  break;
   case ICoth:      *str = "icoth";  break;

   case Asinh :      *str = "asinh";   break;
   case Acosh :      *str = "acosh";   break;
   case Atanh :      *str = "atanh";   break;
   case Acoth :      *str = "acoth";   break;
   case IAsinh:      *str = "iasinh";  break;
   case IAcosh:      *str = "iacosh";  break;
   case IAtanh:      *str = "iatanh";  break;
   case IAcoth:      *str = "iacoth";  break;

   case Pow:         *str = "pow";     break;
   case IPow:        *str = "ipow";    break;

   case EToI:   *str = "extreal_to_int";       break;
   case EToL:   *str = "extreal_to_longreal";  break;
   case LToE:   *str = "longreal_to_extreal";  break;

   case Div:    *str = "r_div";     break;
   case IDiv:   *str = "i_div"; break;

   default:
   /* itoa(fct, s1, 10);
      strcpy(s, "fct no. ");
      strcat(s, s1); */
      sprintf(s, "fct no. %d", fct);
      *str = s;
      break;
   } /* switch(fct) */

   return;
} /* excfct_to_a() */





