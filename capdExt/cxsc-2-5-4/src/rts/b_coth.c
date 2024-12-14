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

/* CVS $Id: b_coth.c,v 1.22 2014/01/30 17:24:03 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lcoth              Processor : C                   *
 *                                                                       *
 * Hyperbolic Cotangent Function for Multiple Precision Arithmetic       *
 * ===============================================================       *
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
 *    Lassign         Assignment of funtion value and number of ulp's    *
 *    Lvarget         Fetching an unused INTERN variable                 *
 *    Lvardrop        Returning a no longer used INTERN variable         *
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
 *                 by the global variables p (LhF) and errsinh (LhE).    *
 *                                                                       *
 *************************************************************************/

#define  Name        "Lcoth"
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
        
#define errexp    LhE           /* Rounding error ( == expve !!! )  */
#define val       LhF           /* Function value ( == expve !!! )  */
#define errsinh   LhE           /* Rounding error ( == sinhvea !!! )*/
/*#define val       LhF  */     /* Function value ( == sinhvea !!! )*/
#define err       LhE           /* Rounding error                   */
#define oneovexp  LhD           /* Value of 1/exp(xi)               */
#define dummy     LhD           /* Dummy variable (not further used)*/

#define sum       cosh          /* exp(|xi|) + exp(-|xi|)           */
#define diff      errcosh       /* exp(|xi|) - exp(-|xi|)           */

/********************************************************************
 * Constants for Computation of coth(x)                             *
 ********************************************************************/
        
#define Lguard   2
        
static a_btyp  mbds[3]   = { 0x80000000L,0x00000000L,0x00000000L };
static a_btyp  m7o2[3]   = { 0x00000003L,0x80000001L,0x00000000L };
static a_btyp  m385[3]   = { 0x00000003L,0xD9D5C53DL,0x00000000L };
static a_btyp  mln2[3]   = { 0xB17217F7L,0xD1CF79ACL,0x00000000L };
        
static dynamic  bdapp   = { 0, 0, 0, 0, -1, Minl, &mbds[0] };
static dynamic  l7o2    = { 0, 0, 0, 0,  0, Minlerr, &m7o2[0] };
static dynamic  l385    = { 0, 0, 0, 0,  0, Minlerr, &m385[0] };
static dynamic  ln2     = { 0, 0, 0, 0, -1, Minl, &mln2[0] };
        
/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_coth(dynamic *xi,dynamic *ri)
#else
int b_coth(xi,ri)
dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Lcoth"
{
   extern dynamic   Lone, Leps, val, err, oneovexp;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern rounding  Lrnd;
   extern char      *Lroutine;
        
   dynamic  *cosh, *errcosh;

   int      sgn, rc, poolvars = 0;

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
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;
        
        
  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/
        
        
  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/
        
   if (xi->z == 1) {
      ERREXIT(OFLOW,OFLOW,0);   /* Error message and handling   */
      }
        
  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/
        
   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling   */
      }

        
        
  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/
        
   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS    */
        
        
        
  /**************************************************************
   * Computation of the Function Value                          *
   **************************************************************/
        
        
  /*----------------*
   | Initialization |
   *----------------*/
        
   Maxl = Minlerr;
        
        
  /**************************************************************
   * Check and handle case coth(x)=sgn(x)                       *
   **************************************************************/
        
   MULINT_(&ln2,(B_LENGTH/2)*Lcurrprec,&dummy);
   ADD_(&dummy,&ln2,&dummy);
   if GE_ABS_(xi,&dummy) {
      sgn  = xi->s;
      Maxl = Lcurrprec;
      Leps.e = -Maxl;
      COPY_(&Lone,ri);
      if (rc != 0)
         ERREXIT(rc,rc,0);      /* Error message and handling   */
      ri->s = sgn;
      ri->r = 1;

      RETURN(0);
      }
        
        
  /**************************************************************
   * Computation of exp(xi)                                     *
   **************************************************************/
        
  /* Computation of exp(xi) using internal return.
     Function value is returned in variable &val ( == LhF ),
     rounding error is returned in variable &errexp ( == LhE ). */
        
   Maxl = Lcurrprec+Lguard;                /* Increase accuracy */
        
   EXP_ABS_(xi);

   cosh = Lvarget();            /* Get new pool variable        */
   errcosh = Lvarget();         /* Get new pool variable        */
   poolvars = 2;                /* Number of used pool variables*/

   if GE_ABS_(xi,&bdapp) {      /* Use exponential function     */

        
  /**************************************************************
   * Computation of coth(xi) using exp(xi)                      *
   **************************************************************/
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lcoth 246) exp  = ");
         Lprinti(&val);
         }
#endif
        
        
  /*---------------------------------*
   | Computation of the Final Result |
   *---------------------------------*/
        
  /* The final result is computed by squaring the value for the reduced
     argument NSQR times, where NSQR is the number of halvings during
     computation of the reduced argument.                       */
        
      DIV_(&Lone,&val,&oneovexp);    /* coth(xi) =           */
      SUB_(&val,&oneovexp,diff);     /* [exp(xi)+1/exp(xi)]  */
      ADD_(&val,&oneovexp,sum);      /* /[exp(xi)-1/exp(xi)] */
      DIV_(sum,diff,&val);
      val.s = xi->s;                       /* Set sign of result */
        
      l385.e = 1 - Maxl;
      Maxl = Minlerr;
      SHIFT_(1,&errexp,&errexp);           /* double error of exp(xi) */
      ADD_(&errexp,&l385,&err);            /* add division error      */
      if (err.r != 0) NEXT_(&err,&err);
      }
        
   else {
        
  /**************************************************************
   * Computation using sinh(xi) and cosh(xi)                    *
   **************************************************************/
        
  /* Computation of cosh(xi) using exp(|xi|).                   */
        
      DIV_(&Lone,&val,cosh);    /* cosh(xi) =                */
      ADD_(cosh,&val,cosh);     /*     [exp(xi)+1/exp(xi)]/2 */
      SHIFT_(-1,cosh,cosh);
        
      COPY_(&errexp,errcosh);

      if (rc != 0) {
         ERREXIT(288,288,poolvars);
         }
        
        
  /* Computation of sinh(xi) using internal return.
     Function value is returned in variable &val ( == LhF ),
     rounding error is returned in variable &errsinh ( == LhE ). */
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lcoth 298) arg  = ");
         Lprinti(xi);
         }
#endif

      SINH_(xi);
        
      DIV_(cosh,&val,&val);     /* coth(xi) = cosh(xi)/sinh(xi)  */
        
      l7o2.e = 1 - Maxl;
      Maxl = Minlerr;
      ADD_(&errsinh,errcosh,&err);
      ADD_(&err,&l7o2,&err);    /* add division error   */
      if (err.r != 0) NEXT_(&err,&err);
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lcoth 315) app  = ");
         Lprinti(&val);
         }
#endif
        
      }
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lcoth 324) |err| = ");  Lprinti(&err);
      printf("\n  (Lcoth 325) val   = ");  Lprinti(&val);
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,RESUL,0);   /* Error message and handling   */
      }
        
        
        
  /*----------------------------------------------------------*
   | Assign result and get number of units for an outer bound |
   *----------------------------------------------------------*/
        
   Maxl = Lcurrprec;
        
   Lrnd = rounded;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n (Lcoth 350) Maxl = %d\n",Maxl);
      }
#endif
        
   ASSIGN_(ri);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n (Lcoth 358) LhF    = ");  Lprinti(&LhF);
      printf("\n (Lcoth 359) LhE    = ");  Lprinti(&LhE);
      printf("\n (Lcoth 360) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   Lvardrop(poolvars);  /* Return used variables to pool */
   RETURN(0);
}





