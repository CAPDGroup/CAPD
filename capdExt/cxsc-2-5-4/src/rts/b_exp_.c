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

/* CVS $Id: b_exp_.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lexp               Processor : C                   *
 *                                                                       *
 * Exponential Function for Multiple Precision Arithmetic                *
 * ======================================================                *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Argument Range :    All nonnegative numbers                           *
 * Result Range :      Nonnegative numbers                               *
 *                                                                       *
 * Used Functions :                                                      *
 *    expve           Computation of function value and error bound      *
 *    Lassign         Assignment of function value and number of ulp's   *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhF, LhE        Used for debugging only                            *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#define  Name        "Lexp"
#define  Main
        
#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

static char  *function = Name;
        
/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_exp_(dynamic *xi,dynamic *ri)
#else
int b_exp_(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Lexp"
{
   extern dynamic Lone;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec;
   extern char    *Lroutine;
        
   int  rc;
#define poolvars  0

#ifdef Debug
   extern dynamic   LhF, LhE;
   extern int       Ldebug;
#endif
        
#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level         */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;                    /* Save precision setting */
   rc = 0;
        
        
  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/
        
   if (xi->z == 1) {
      if ((COPY_(&Lone,ri)) != 0)       /* ri = exp(0) = 1      */
         ERREXIT(rc,rc,0);
      ri->r = 0;                        /* Result exact         */
      EXIT(0);
      }

        
  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/
        
   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling       */
      }
        
   rc  = expve(xi);             /* Function value and error bound   */
   if (rc != 0)
      {
      if (rc != UFLOW)
         ERREXIT(0,RESUL,0)     /* Error message and handling       */
      else {
         ri->z = 1;             /* Set result 0                 */
         ri->s = 0;             /* Set sign to positive         */
         ri->r = 1;             /* Set rounding bit to 1 ulp    */
         ERREXIT(0,0,0);        /* Error message and handling   */
         }
      }
        
#if INT_HPREC
   b_case = (xi->s) ? LESS_ABS_ONE : GREATER_ABS_ONE;
#endif

   ASSIGN_(ri);                 /* Assign result and ulp's          */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lexp 139) LhF    = ");  Lprinti(&LhF);
      printf("\n  (Lexp 140) LhE    = ");  Lprinti(&LhE);
      printf("\n  (Lexp 141) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif
        
   RETURN(0);
}
        
        
        





