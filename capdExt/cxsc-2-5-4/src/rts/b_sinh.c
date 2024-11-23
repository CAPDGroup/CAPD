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

/* CVS $Id: b_sinh.c,v 1.22 2014/01/30 17:24:05 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lsinh              Processor : C                   *
 *                                                                       *
 * Hyperbolic Sine Function for Multiple Precision Arithmetic            *
 * ==========================================================            *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Argument Range :    Restricted only due to overflow                   *
 * Result Range :      All nonnegative internal numbers                  *
 *                                                                       *
 * Used Functions :                                                      *
 *    expve           Intern Exponential Function                        *
 *    sinhvea         Intern Hyperbolic Sine Function for |x|<0.5        *
 *    Lginit          Initialization of global constants                 *
 *                    Intern arithmetic                                  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhD, LhE, LhF                                        (Variables)   *
 *    Lone              ( = 1 )                            (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 * Dependencies on other INTERN Functions :                              *
 *    expve        Computed function value and error bound are returned  *
 *                 by the global variables p (LhF) and errexp (LhE).     *
 *    sinhvea      Computed function value and error bound are returned  *
 *                 by the global variables p (LhF) and errexp (LhE).     *
 *                                                                       *
 *************************************************************************/

#define  Name        "Lsinh"
#define  Main
        
#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

static char  *function = Name;
        
/****************************************************************
 * Redefinition of Names of Intermediate Variables Used         *
 ****************************************************************/
        
#define errexp    LhE   /* Rounding error    ( == expve !!! )     */
#define val       LhF   /* Function value    ( == expve !!! )     */
#define errsinh   LhE   /* Rounding error    ( == sinhvea !!! )   */
/*#define val     LhF *//* Function value    ( == sinhvea !!! )   */
#define err       LhE   /* Rounding error                         */
#define oneovexp  LhD   /* Value of 1/exp(xi)                     */
        
/*****************************************************************
 * Constants for Computation of sinh(x)                          *
 *****************************************************************/
        
#define Lguard    2
#define poolvars  0
        
static a_btyp  mbds[3]   = { 0x80000000L,0x00000000L,0x00000000L };
static a_btyp  m5o2[3]   = { 0x00000002L,0x80000001L,0x00000000L };
        
static dynamic  bdapp   = { 0, 0, 0, 0, -1, Minl, &mbds[0] };
static dynamic  l5o2    = { 0, 0, 0, 0,  0, Minlerr, &m5o2[0] };

/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_sinh(dynamic *xi,dynamic *ri)
#else
int b_sinh(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Lsinh"
{
   extern dynamic   Lone, val, err, oneovexp;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern rounding  Lrnd;
   extern char      *Lroutine;
        
   int       sgn, rc;

#ifdef Debug
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
   Lcurrprec = Maxl;                    /* Save precision setting       */
   rc = 0;
        
        
  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/

        
  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/
        
   if (xi->z == 1) {
      if ((rc = b_bini(ri)) == 0) {         /* ri = sinh(0) = 0  */
         ri->r = 0;                        /* Result exact      */
         EXIT(0);
         }
       else
         ERREXIT(rc,rc,0);      /* Error message and handling    */
      }
        
  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/
        
   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling    */
      }

        
        
  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/
        
   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */
        
        
        
  /**********************************************************************
   * Computation of the Function Value                                  *
   **********************************************************************/


  /*----------------*
   | Initialization |
   *----------------*/

   Maxl = Lcurrprec+Lguard;                /* Increase accuracy         */

   if GE_ABS_(xi,&bdapp) {                 /* Use exponential function  */


  /**********************************************************************
   * Computation using exp(xi)                                          *
   **********************************************************************/

  /* Computation of exp(|xi|) using internal return.
     Function value is returned in variable &val ( == LhF ),
     rounding error is returned in variable &errexp ( == LhE ).         */

      EXP_ABS_(xi);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lsinh 206) exp  = ");
         Lprinti(&val);
         }
#endif


  /*---------------------------------*
   | Computation of the Final Result |
   *---------------------------------*/

  /* The final result is computed by squaring the value for the reduced
     argument NSQR times, where NSQR is the number of halvings during
     computation of the reduced argument.                               */

      DIV_(&Lone,&val,&oneovexp);  /* sinh(xi) =                    */
      SUB_(&val,&oneovexp,&val);   /*     [exp(xi)-1/exp(xi)]/2     */
      DIVINT_(&val,2,&val);
      val.s = xi->s;                       /* Set sign of result        */

      l5o2.e = 1 - Maxl;
      Maxl = Minlerr;
      ADD_(&errexp,&l5o2,&err);            /* add division error        */
      if (err.r != 0) NEXT_(&err,&err);
      Lrnd = rounded;

      if (rc != 0) {
         ERREXIT(RESUL,RESUL,0);           /* Error message and handling*/
         }
      }

   else {
        
  /**********************************************************************
   * Computation using Taylor approximation                             *
   **********************************************************************/
        
  /* Computation of sinh(xi) using internal return.
     Function value is returned in variable &val ( == LhF ),
     rounding error is returned in variable &errsinh ( == LhE ).        */
        
      SINH_(xi);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lsinh 250) app  = ");
         Lprinti(&val);
         }
#endif
        
      if (errsinh.r != 0) NEXT_(&errsinh,&err);
      if (rc != 0) {
         ERREXIT(RESUL,RESUL,0);           /* Error message and handling*/
         }
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lsinh 263) |err| = ");  Lprinti(&err);
      printf("\n  (Lsinh 264) val   = ");  Lprinti(&val);
      }
#endif
        
        
        
  /*----------------------------------------------------------*
   | Assign result and get number of units for an outer bound |
   *----------------------------------------------------------*/
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lsinh 280) Maxl = %d\n",Maxl);
      }
#endif

   ASSIGN_(ri);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lsinh 288) LhF    = ");  Lprinti(&LhF);
      printf("\n  (Lsinh 289) LhE    = ");  Lprinti(&LhE);
      printf("\n  (Lsinh 290) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





