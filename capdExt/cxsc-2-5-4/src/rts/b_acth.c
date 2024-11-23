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

/* CVS $Id: b_acth.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Larcoth            Processor : C                   *
 *                                                                       *
 * Inverse Hyperbolic Cotangent for Multiple Precision Arithmetic        *
 * ==============================================================        *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Argument Range :    |xi| < 1                                          *
 * Result Range :      All nonzero numbers                               *
 *                                                                       *
 * Used Functions :                                                      *
 *    On computation of function values the computed value and the error *
 *    bound are returned by the global variables LhF and LhE.   However, *
 *    in many cases these names are redefined to a symbolic name related *
 *    to the called function.                                            *
 *                                                                       *
 *    name      macro   description                                      *
 *    --------  ------   ------------------------------------------------*
 *    sqrtve    SQRT_    Intern Square Root Function                     *
 *    lnve      LN_      Intern Logarithm Function                       *
 *    lnvea     LN_T_    Intern Logarithm Function (arguments near 1)    *
 *    expve     EXP_     Intern Exponential Function                     *
 *    Lassign   ASSIGN_  Assignment of function value and number of ulp's*
 *    Lginit    (none)   Initialization of global constants              *
 *                    Intern arithmetic (macro names see file b_lari.h)  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhF, LhE, LhD   Used for debugging only                            *
 *    Lone              ( = 1 )                            (Constant)    *
 *    Leps              ( = current rounding error )       (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#define  Name        "Larcoth"
#define  Main

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

static char  *function = Name;

/************************************************************************
 * Redefinition of Names of Intermediate Variables Used                 *
 ************************************************************************/

#define Lres      LhF           /* Function value                       */
#define err       LhE           /* Rounding error                       */
#define expy      LhF           /* Function value exp(y)                */
#define errexp    LhE           /* Rounding error of exp(y)             */
#define sqty      LhF           /* Function value sqrt(y)               */
#define dummy     LhD           /* Dummy variable (not further used)    */


/************************************************************************
 * Constants for Computation of arcoth(x)                               *
 ************************************************************************/

#define Lguard    2
#define M0        4
#define ebdb      1
#define invexp    dvsr
#define Lhlp      dvnd

static a_btyp  macc[3]   = { 0x00100000L,0x00000000L,0x00000000L };
static a_btyp  mbdb[3]   = { 0x00000001L,0x55555555L,0x00000000L };
static a_btyp  maerr[3]  = { 0x00000001L,0x877C5693L,0x00000000L };
static a_btyp  mbda[3]   = { 0x00000014L,0x55555556L,0x00000000L };
static a_btyp  m5o2[3]   = { 0x00000002L,0x80000001L,0x00000000L };

static dynamic  acc   = { 0, 0, 0, 0, 0, Minl, &macc[0] };
static dynamic  bdb   = { 0, 0, 0, 0, 0, Minlerr, &mbdb[0] };
static dynamic  bda   = { 0, 0, 0, 0, 0, Minlerr, &mbda[0] };
static dynamic  aerr  = { 0, 0, 0, 0, 0, Minlerr, &maerr[0] };
static dynamic  l5o2  = { 0, 0, 0, 0, 0, Minlerr, &m5o2[0] };

#ifdef Debug
static a_btyp  mbds[3]   = { 0x00004000,0x00000000,0x00000000 };
static dynamic  bds   = { 0, 0, 0, 0, 0, Minl, &mbds[0] };
#endif

/************************************************************************
 * Algorithm                                                            *
 ************************************************************************/

#ifdef LINT_ARGS
int b_acth(dynamic *xi,dynamic *ri)
#else
int b_acth(xi,ri)
dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larcoth"
{
   extern dynamic   Lone, Lres, err, dummy;
   extern a_btyp    Maxl;
   extern a_btyp    Lcurrprec, Lgiflag;
   extern char      *Lroutine;
   extern rounding  Lrnd;
   extern a_real    *r_1o2_, *r_one_;

   dynamic   *y, *Res, *Arg, *sqt, *sinhy, *coshy, *dvsr, *dvnd;
   int       rc, poolvars;
   a_btyp    prec;
   a_real    yd;

#ifdef Debug
   extern int       Ldebug;
#endif

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;


  /**********************************************************************
   * Check and Handle Range Error and Special Cases                     *
   **********************************************************************/


  /*-------------*
   | Check Range |
   *-------------*/

   if (xi->z == 1) {
      ERREXIT(RANGE,RANGE,0);   /* Error message and handling   */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling   */
      }


  /*-------------*
   | Check Range |
   *-------------*/

   if LE_ABS_(xi,&Lone) {
      ERREXIT(RANGE,RANGE,0);   /* Error message and handling   */
      }



  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcoth 207) Maxl = %d\n",Maxl);
      }
#endif


   Arg = Lvarget();
   y   = Lvarget();
   poolvars = 2;
   Maxl = xi->l;
   COPY_(xi,Arg);                        /* Arg = |xi|  */
   Arg->s = 0;
   Maxl = Lcurrprec;

   if GE_(Arg,&bda) {

     /*-----------------------------------------------------------------*
      | Computation of arcoth(xi) using logarithm                       |
      *-----------------------------------------------------------------*/

     /* Compute for |xi| >= 61/3
            arcoth(xi) = artanh(1/xi)         (using lnvea)     */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 231) |xi|<=61/3");
         }
#endif

      Maxl += 2;

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 243) Maxl = %d",Maxl);
         }
#endif

      DIV_(&Lone,Arg,y);
      LN_T_(y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 252) Maxl = %d",Maxl);
         }
#endif
      SHIFT_(-1,&Lres,&Lres);

      aerr.e = 1 - Maxl;
      ADD_(&err,&aerr,&err);

      }

    else if LT_(Arg,&bdb) {

     /* Compute for |xi| > 3/4
            arcoth(xi) = ln([xi+1]/[xi-1])/2
            (using lnve)      */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 269) |xi|>3/4");
         }
#endif

      dvnd = Lvarget();
      dvsr = Lvarget();
      poolvars += 2;

      Maxl += 1;

      ADD_(Arg,&Lone,dvnd);
      SUB_(Arg,&Lone,dvsr);
      Maxl += 1;
      DIV_(dvnd,dvsr,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 290) Maxl = %d",Maxl);
         }
#endif

      if (rc != 0) ERREXIT(RESUL,294,poolvars);
      LN_(y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 299) Maxl = %d",Maxl);
         }
#endif

      SHIFT_(-1,&Lres,&Lres);

      aerr.e = 1 - Maxl;
      ADD_(&err,&aerr,&err);

      }

    else {

     /*-----------------------------------------------------------------*
      | Computation of an Initial Approximation                         |
      *-----------------------------------------------------------------*/

     /* An initial approximation is computed by taking the floating point
        logarithm of the argument xi.           */

      Maxl = Minlfl;
#ifdef Debug
      if GT_(Arg,&bds) {
         printf("\n  (Larcoth 326) This case should not have occurred !");
         rc = 0;
         COPY_(Arg,y);                     /* y = |xi|          */
         }
       else {
#endif
              /* y = ln([xi+1]/[xi-1])/2 */
      rc = b_bcid(Arg,&yd,(a_intg)0);
      if (rc==ROUND) rc = 0;
/* Cordes      rc += b_bcdi(log((yd+1)/(yd-1))/2,y);            */
      rc += b_bcdi(r_muln(*r_1o2_,
                   B_LOG_(r_divn(r_addn(yd,*r_one_),r_subn(yd,*r_one_)))
                   ),
                &y,(a_intg)0);
#ifdef Debug
         }
#endif

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 341) y = arcoth(|xi|) = ");
         Lprinti(y);
         }
#endif


     /*-----------------------------------------------------------------*
      | Newton Iteration Steps for gaining the required accuracy        |
      *-----------------------------------------------------------------*/

     /* The accuracy of the initial approximation is increased by a number
        of Newton steps until the relative error bound is less than 1% of
        the initial accuracy setting.                              */

      acc.e = -Lcurrprec;

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 363) acc = ");
         Lprinti(&acc);
         }
#endif

      Maxl = Minlerr;
      sqt  = Lvarget();
      poolvars += 1;

      MUL_(Arg,Arg,sqt,&dummy);
      SUB_(sqt,&Lone,sqt);
      if (rc != 0) ERREXIT(RESUL,374,poolvars);
      SQRT_(sqt);
      COPY_(&sqty,sqt);

      COPY_(&acc,&err);
      NEXT_(&err,&err);
      err.s = 0;

      sinhy  = Lvarget();
      coshy  = Lvarget();
      invexp = Lvarget();
      Lhlp   = Lvarget();
      poolvars += 4;

      for (prec=Maxl=M0; GT_(&err,&acc);
           prec=Maxl=min(2*Maxl,Lcurrprec+2)) {

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (Larcoth 392) y   = ");  Lprinti(y);
            printf("\n  (Larcoth 393) Maxl = %d",Maxl);
            }
#endif

         if (rc != 0) ERREXIT(RESUL,397,poolvars);
         EXP_(y);               /* Lhlp := exp(y) - |xi|        */
         DIV_(&Lone,&expy,invexp);
         SUB_(&expy,invexp,sinhy);
         SHIFT_(-1,sinhy,sinhy);
         ADD_(&expy,invexp,coshy);
         SHIFT_(-1,coshy,coshy);
         MUL_(Arg,sinhy,Lhlp,&dummy);
         SUB_(coshy,Lhlp,Lhlp);

         Maxl = Minlerr;
         NEXT_(Lhlp,&dummy);    /* determine error bound &err   */
         DIV_(Lhlp,sqt,&err);
         DIV_(&err,y,&err);
         err.s = 0;

         Maxl = prec;
         MUL_(Lhlp,sinhy,Lhlp,&dummy);  /* y := y + [cosh(y) */
         ADD_(y,Lhlp,y);                /*   -Arg*sinh(y)]*sinh(y) */

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (Larcoth 419) expy  = ");  Lprinti(&expy);
            printf("\n  (Larcoth 420) sinhy = ");  Lprinti(sinhy);
            printf("\n  (Larcoth 421) coshy = ");  Lprinti(coshy);
            printf("\n  (Larcoth 422) Lhlp  = ");  Lprinti(Lhlp);
            printf("\n  (Larcoth 423) dummy = ");  Lprinti(&dummy);
            printf("\n  (Larcoth 424) y   = ");  Lprinti(y);
            printf("\n  (Larcoth 425) err = ");  Lprinti(&err);
            }
#endif
         }

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 432) Iterations performed !");
         printf("\n");
         }
#endif

      Maxl = Lcurrprec+2;
      Res = Lvarget();
      poolvars += 1;
      rc = 0;
      COPY_(y,Res);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcoth 445) Res  = ");  Lprinti(Res);
         printf("\n");
         }
#endif


     /*-----------------------------------------------------------------*
      | Determine number of units to be added for an upper bound        |
      *-----------------------------------------------------------------*/

      if (rc != 0) ERREXIT(RESUL,459,poolvars);
      EXP_(Res);                        /* compute exp(Res)     */
      DIV_(&Lone,&expy,invexp);
      SUB_(&expy,invexp,sinhy);
      SHIFT_(-1,sinhy,sinhy);
      ADD_(&expy,invexp,coshy);
      SHIFT_(-1,coshy,coshy);
      MUL_(Arg,sinhy,Lhlp,&dummy);

      l5o2.e = 1-Maxl;
      l5o2.s = err.s;
      Maxl = Minlerr;
      ADD_(&err,&l5o2,&err);    /* err = err(cosh y) = err(sinh y)  */
      NEXT_(&err,&err);         /*     = err(exp y) + 2.5*b^(1-l)   */
      ADD_(Lhlp,coshy,&dummy);
      NEXT_(&dummy,&dummy);
      MUL_(&err,&dummy,&err,&dummy);    /* err = total absolute error   */
      NEXT_(&err,&err);

      SUB_(Lhlp,coshy,Lhlp);
      NEXT_(Lhlp,Lhlp);
      err.s = Lhlp->s;          /* add absolute error of difference */
      ADD_(Lhlp,&err,Lhlp);
      NEXT_(Lhlp,&err);

      DIV_(&err,sqt,&err);
      NEXT_(&err,&err);
      DIV_(&err,y,&err);
      NEXT_(&err,&err);

      Lrnd = rounded;

      Maxl = Lcurrprec+2;
      COPY_(Res,&Lres);

      }

   Lres.s = xi->s;              /* Set sign of result   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcoth 500) Maxl = %d\n",Maxl);
      printf("\n  (Larcoth 501) Lres = ");  Lprinti(&Lres);
      printf("\n");
      }

   if (Ldebug >= 0) printf("\n");
#endif

   Lvardrop(poolvars);

   if (rc != 0) {
      ERREXIT(RESUL,511,0);   /* Error handling         */
      }

   ASSIGN_(ri);         /* Assign result and set number of ulp's    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcoth 518) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larcoth 519) LhE    = ");  Lprinti(&err);
      printf("\n  (Larcoth 520) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





