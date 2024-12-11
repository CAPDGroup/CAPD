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

/* CVS $Id: b_tan_.c,v 1.22 2014/01/30 17:24:05 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Ltan               Processor : C                   *
 *                                                                       *
 * Tangent Function for Multiple Precision Arithmetic                    *
 * ==================================================                    *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Argument Range :    All numbers                                       *
 * Result Range :      All numbers                                       *
 *                                                                       *
 * Used Functions :                                                      *
 *    sicovea         Computation of function value and error bound      *
 *    Lassign         Assignment of function value and number of ulp's   *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhF, LhD        Used for debugging only                            *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Ldebug          Flag for printing additional information           *
 *    Lversion        Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#define  Name        "Ltan"
#define  Main

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

#define Lguard  2

static char  *function = Name;

/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/

#ifdef LINT_ARGS
int b_tan_(dynamic *xi,dynamic *ri)
#else
int b_tan_(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Ltan"
{
#if INT_HPREC
   extern dynamic Lone;
#endif
   extern dynamic   LhF, LhD, LhE, Leps;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec;
   extern char      *Lroutine;
   extern a_btyp LhI;  /* entier(x/(pi/4)) mod 16 when calling sin(x) */

   int    rc;
#define poolvars  0

#ifdef Debug
   extern int       Ldebug;
#endif

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level          */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif



  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;
   LhI = 0;

  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z == 1) {
      ri->z = 1;                /* ri = tan(0) = 0              */
      ri->r = 0;
      EXIT(0);
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling   */
      }

   SICO_(xi);           /* Compute function value and error bound   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Ltan 132) LhF = ");  Lprinti(&LhF);
      printf("\n  (Ltan 133) LhD = ");  Lprinti(&LhD);
      }
#endif

   Maxl = Lcurrprec + Lguard;
   DIV_(&LhF,&LhD,&LhF);

   Leps.e = 1 - Maxl;
   Maxl = Minlerr;
   SHIFT_(1,&LhE,&LhE);
   NEXT_(&LhE,&LhE);
   ADD_(&LhE,&Leps,&LhE);
   NEXT_(&LhE,&LhE);

   if (rc != 0) {
      ERREXIT(RESUL,RESUL,0);   /* Error message and handling   */
      }

#if INT_HPREC
   if (LT_ABS_(xi,&Lone))
      {
      b_farg = xi;
      b_case = GREATER_ABS_ARG;
      }
#endif

   ASSIGN_(ri);         /* Assign tangent and set number of ulp's   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Ltan 155) tan = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





