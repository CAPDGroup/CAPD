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

/* CVS $Id: b_cosh.c,v 1.22 2014/01/30 17:24:03 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lcosh              Processor : C                   *
 *                                                                       *
 * Hyperbolic Cosine Function for Multiple Precision Arithmetic          *
 * ============================================================          *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Argument Range :    Restricted only due to overflow                   *
 * Result Range :      All internal numbers not less 1                   *
 *                                                                       *
 * Used Functions :                                                      *
 *    expve     Computation of exponential function with error bound     *
 *              Intern arithmetic                                        *
 *    Lginit    Initialization of global constants                       *
 *    Lassign   Assignment of function value and error bound             *
 *    Lvarget   Fetching an unused INTERN variable                       *
 *    Lvardrop  Returning a no longer used INTERN variable               *
 *              Intern arithmetic                                        *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    Lone              ( = 1 )                            (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 * Dependencies on other INTERN Functions :                              *
 *    expve        Computed function value and error bound are returned  *
 *                 by the global variables p (LhF) and err (LhE).        *
 *                                                                       *
 *************************************************************************/

#define  Name        "Lcosh"
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
        
#define err       LhE   /* Rounding error    ( == expve !!! )   */
#define p         LhF   /* Polynomial value  ( == expve !!! )   */
#define dummy     LhD   /* Dummy variable (not further used)    */
        
        
/****************************************************************
 * Constants for Computation of cosh(x)                         *
 ****************************************************************/

#define Lguard   2

static a_btyp  m5o2[3]   = { 0x00000002L,0x80000001L,0x00000000L };

static dynamic   l5o2    = { 0, 0, 0, 0,  0, Minlerr, &m5o2[0] };
        
/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_cosh(dynamic *xi,dynamic *ri)
#else
int b_cosh(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Lcosh"
{
   extern dynamic   Lone, p, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern rounding  Lrnd;
   extern char      *Lroutine;
        
   dynamic   *expo2;
   int       sgn, rc, poolvars;
   unsigned int ulp;    /* !!! must be int !!! */

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
      if (COPY_(&Lone,ri) == 0) {          /* ri = cosh(0) = 1  */
         ri->r = 0;                        /* Result exact      */
         EXIT(0);
         }
       else
         ERREXIT(rc,rc,0);      /* Error message and handling   */
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
        
   Maxl = Lcurrprec + Lguard;   /* Set appropriate accuracy  */
   poolvars = 0;


  /*------------------------*
   | Computation of exp(xi) |
   *------------------------*/
        
  /* Computation of exp(xi) using internal return.
     Function value is returned in variable &p ( == LhF ),
     rounding error is returned in variable &err ( == LhE ).    */

   EXP_ABS_(xi);
        
   expo2 = Lvarget();
   poolvars = 1;

   SHIFT_(-1,&p,expo2);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lcosh 205) exp   = ");  Lprinti(&p);
      printf("\n  (Lcosh 206) exp/2 = ");  Lprinti(expo2);
      }
#endif
        
        
  /*---------------------------------*
   | Computation of the Final Result |
   *---------------------------------*/
        
   DIV_(&Lone,&p,&dummy);         /* cosh(xi) =            */
   ADD_(&p,&dummy,&p);            /* [exp(xi)+1/exp(xi)]/2 */
   SHIFT_(-1,&p,&p);
        
   l5o2.e = 1 - Maxl;
   Maxl = Minlerr;
   ADD_(&err,&l5o2,&err);       /* add division error   */
   if (err.r != 0) NEXT_(&err,&err);

   if (rc != 0) {
      ERREXIT(RESUL,RESUL,0);   /* Error message and handling   */
      }
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lcosh 230) cosh = ");  Lprinti(&p);
      printf("\n  (Lcosh 231) err  = ");  Lprinti(&err);
      }
#endif


  /*---------------------------*
   | Restore current precision |
   *---------------------------*/

   Lrnd = rounded;              /* Value is upper bound          */
   ASSIGN_(ri);                 /* Assign result and set ulp's   */

   if LT_(ri,&Lone) {           /* Enforce cosh(xi) >= 1 :       */
      ulp = ri->r;              /*    Save number of ulp's       */
      COPY_(&Lone,ri);          /*    ri = 1                     */
      ri->r = ulp - 1;          /*    Set corrected ulp's        */
      }
   if LT_(ri,expo2) {           /* Enforce cosh(xi) >= exp(xi)/2 */
      ulp = ri->r;              /*    Save number of ulp's       */
      COPY_(expo2,ri);          /*    ri = exp(|xi|)/2           */
      ri->r = max(ulp-1,1);     /*    Set corrected ulp's        */
      }

   if (rc != 0) {
      ERREXIT(RESUL,RESUL,poolvars);    /* Error message and handling */
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lcosh 260) LhF    = ");  Lprinti(&LhF);
      printf("\n  (Lcosh 261) LhE    = ");  Lprinti(&LhE);
      printf("\n  (Lcosh 262) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   Lvardrop(poolvars);          /* Return used variables to pool*/
   RETURN(0);
}





