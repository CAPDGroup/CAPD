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

/* CVS $Id: b_acsh.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/**************************************************************************
 *                                                                        *
 * Descriptive Name : Larcosh            Processor : C                    *
 *                                                                        *
 * Inverse Hyperbolic Cosine for Multiple Precision Arithmetic            *
 * ===========================================================            *
 *                                                                        *
 * Include files :  b_lari.h  - definitions for INTERN standard functions *
 *                                                                        *
 * Function value : int       - 0                                         *
 *                                                                        *
 * Used Number Base :  2**32                                              *
 *                                                                        *
 * Argument Range :    Argument >= 1                                      *
 * Result Range :      Result >= 0                                        *
 *                                                                        *
 * Used Functions :                                                       *
 *    On computation of function values the computed value and the error  *
 *    bound are returned by the global variables LhF and LhE.   However,  *
 *    in many cases these names are redefined to a symbolic name related  *
 *    to the called function.                                             *
 *                                                                        *
 *    name      macro   description                                       *
 *    --------  ------   -------------------------------------------------*
 *    sqrtve    SQRT_    Intern Square Root Function                      *
 *    lnve      LN_      Intern Logarithm Function                        *
 *    lnvea     LN_T_    Intern Logarithm Function (arguments near 1)     *
 *    expve     EXP_     Intern Exponential Function                      *
 *    Lassign   ASSIGN_  Assignment of function value and number of ulp's *
 *    Lginit    (none)   Initialization of global constants               *
 *                       Intern arithmetic (macro names see file b_lari.h)*
 *                                                                        *
 *                                                                        *
 * Used Global INTERN Variables and Constants :                           *
 *    LhF, LhE, LhD                                                       *
 *    Lone              ( = 1 )                            (Constant)     *
 *    Leps              ( = current rounding error )       (Constant)     *
 *                                                                        *
 *                                                                        *
 * Used UNSIGNED Global Variables :                                       *
 *    Maxl            Length of INTERN Variables                          *
 *    Lcurrprec       Initial Length of INTERN Variables                  *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables  *
 *    Ldebug          Flag for printing additional information            *
 *                                                                        *
 **************************************************************************/

#define  Name        "Larcosh"
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

#define Lres      LhF       /* Function value                           */
#define err       LhE       /* Rounding error                           */
#define expy      LhF       /* Function value exp(y)                    */
#define errexp    LhE       /* Rounding error of exp(y)                 */
#define sqty      LhF       /* Function value sqrt(y)                   */
#define errsqt    LhE       /* Rounding error of sqrt(y)                */
#define dummy     LhD       /* Dummy variable (not further used)        */

/************************************************************************
 * Constants for Computation of arcosh(x)                               *
 ************************************************************************/

#define Lguard    2
#define M0        4
#define ebdb      1
#define invexp    dvsr
#define summ2     dvsr
#define summ1     dvnd
#define serr      sqt

static a_btyp  macc[3]   = { 0x00100000L,0x00000000L,0x00000000L };
static a_btyp  mbdb[3]   = { 0x00000001L,0x40000000L,0x00000000L };
static a_btyp  mbdt[3]   = { 0x00000001L,0x0051EB85L,0x00000000L };
static a_btyp  maerr[3]  = { 0x9CED9168L,0x72B020C5L,0x00000000L };
static a_btyp  m5o2[3]   = { 0x00000002L,0x80000001L,0x00000000L };

static dynamic  acc   = { 0, 0, 0, 0, 0, Minl, &macc[0] };
static dynamic  bdb   = { 0, 0, 0, 0, 0, Minlerr, &mbdb[0] };
static dynamic  bdt   = { 0, 0, 0, 0, 0, Minlerr, &mbdt[0] };
static dynamic  aerr  = { 0, 0, 0, 0, 0, Minlerr, &maerr[0] };
static dynamic  l5o2  = { 0, 0, 0, 0, 0, Minlerr, &m5o2[0] };


/***********************************************************************
 * Algorithm                                                           *
 ***********************************************************************/

#ifdef LINT_ARGS
int b_acsh(dynamic *xi,dynamic *ri)
#else
int b_acsh(xi,ri)
dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larcosh"
{
   extern dynamic   Lone, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;
   extern rounding  Lrnd;
   extern a_real      Fln2, Flnb, *r_one_;

   dynamic   *y, *Res, *sqt, *sinhy, *coshy, *dvsr, *dvnd, *Lhlp;
   int       rc, poolvars;
   a_intg    argexp;
   a_btyp    prec;
   a_real    yd;

#ifdef Debug
   extern int       Ldebug;
#endif

#ifdef Debug
   Ldebug -= 1;             /* Diminish Debug Level                     */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;        /* Save precision setting                   */
   rc = 0;


  /**********************************************************************
   * Check and Handle Range Error and Special Cases                     *
   **********************************************************************/


  /*-------------*
   | Check Range |
   *-------------*/

   if (xi->z == 1) {
      ERREXIT(RANGE,RANGE,0);   /* Error message and handling           */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling           */
      }


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   if EQ_(xi,&Lone) {
      if ((rc = b_bini(ri)) == 0) {    /* ri = arcosh(1) = 0            */
         ri->r = 0;                    /* Result exact                  */
         EXIT(0);
         }
       else
         ERREXIT(rc,rc,0);             /* Error message and handling    */
      }


  /*-------------*
   | Check Range |
   *-------------*/

   if LT_(xi,&Lone) {
      ERREXIT(RANGE,RANGE,0);   /* Error message and handling           */
      }



  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();           /* Initialization of GLOBALS     */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcosh 222) Maxl = %d\n",Maxl);
      }
#endif


   dvnd = Lvarget();
   dvsr = Lvarget();
   y    = Lvarget();
   sqt  = Lvarget();
   poolvars = 4;


   if LE_(xi,&bdb) {                   /* save exponent of xi           */

 /*---------------------------------------------------------------------*
  | Computation of arcosh(xi) using logarithm                           |
  *---------------------------------------------------------------------*/

 /* Compute for |xi| <= 5/4
            arcosh(xi) = ln{xi + sqrt([xi+1]/[xi-1])}
            (using lnve)
                       = artanh{sqrt([xi+1]/[xi-1])/xi}
            (using lnvea)      */

      rc = 0;
      ADD_(xi,&Lone,dvnd);
      Maxl = xi->l + 1;
      SUB_(xi,&Lone,dvsr);
      Maxl = Lcurrprec + 2;
      MUL_(dvnd,dvsr,y,&dummy);
      if (rc != 0) ERREXIT(RESUL,254,poolvars);
      SQRT_(y);
      DIVINT_(&err,5,serr);

      if LE_(xi,&bdt) {
     /* arcosh(xi) = artanh{sqrt([xi+1]/[xi-1])/xi}
            (using lnvea)      */

         DIV_(&sqty,xi,y);
         if (rc != 0) ERREXIT(RESUL,262,poolvars);
         LN_T_(y);
         Maxl += 1;
         SHIFT_(-1,&Lres,&Lres);
         Maxl -= 1;
         }

       else {
     /* arcosh(xi) = ln{xi + sqrt([xi+1]/[xi-1])}
            (using lnve)       */

         ADD_(&sqty,xi,y);
         if (rc != 0) ERREXIT(RESUL,273,poolvars);
         LN_(y);
         }

      aerr.e = - Maxl;
      Maxl = Minlerr;
      ADD_(&err,&aerr,&err);
      ADD_(&err,serr,&err);

      }

    else {

 /*------------------------------------------------------------------*
  | Computation of an Initial Approximation                          |
  *------------------------------------------------------------------*/

 /* An initial approximation is computed by taking the floating point
        logarithm of the argument xi.                                   */

      Maxl = Minlfl;

      if ((argexp = xi->e) >= ebdb) {  /* save exponent of xi           */
         xi->e = 0;                    /* set exponent of xi to 0       */
         rc = b_bcid(xi,&yd,(a_intg)0); /* y = ln(mant(xi)) + ln2 +   */
                                    /*   + expo(xi)*lnb        */
         if (rc==ROUND) rc = 0;
/*Cordes rc += b_bcdi(log(yd)+Fln2+argexp*Flnb,y);      */
         rc += b_bcdi(r_addn(r_addn(
                             B_LOG_(yd),Fln2
                             ),
                      r_muln(r_flot(argexp),Flnb)
                      ),
               &y,(a_intg)0);
         xi->e = argexp;               /* restore exponent of xi        */
         }
       else {
         rc = b_bcid(xi,&yd,(a_intg)0); /* y = ln(xi) + ln(1 + sqrt(1 + */
                                    /*   + 1/(xi*xi)))             */
         if (rc==ROUND) rc = 0;
/*Cordes rc += b_bcdi(log(yd)+log(1+sqrt(1-1/(yd*yd))),y);      */
         rc += b_bcdi(r_addn(B_LOG_(yd),
                      B_LOG_(r_addn(*r_one_,B_SQRT(
                         r_subn(*r_one_,r_divn(*r_one_,r_muln(yd,yd)))
                         )))
                      ),
               &y,(a_intg)0);
         }

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcosh 314) y = arcosh(xi) = ");
         Lprinti(y);
         }
#endif


 /*---------------------------------------------------------------------*
  | Newton Iteration Steps for gaining the required accuracy            |
  *---------------------------------------------------------------------*/

 /* The accuracy of the initial approximation is increased by a number
        of Newton steps until the relative error bound is less than 1% of
        the initial accuracy setting.                                   */

      acc.e = -Lcurrprec;

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcosh 336) acc = ");
         Lprinti(&acc);
         }
#endif

      Maxl = Minlerr;

      MUL_(xi,xi,sqt,&dummy);
      SUB_(sqt,&Lone,sqt);
      if (rc != 0) ERREXIT(RESUL,345,poolvars);
      SQRT_(sqt);
      COPY_(&sqty,sqt);

      COPY_(&acc,&err);
      NEXT_(&err,&err);
      err.s = 0;

      sinhy = Lvarget();
      coshy = Lvarget();
      Lhlp  = Lvarget();
      poolvars += 3;

      for (prec=Maxl=M0; GT_(&err,&acc);
           prec=Maxl=min(2*Maxl,Lcurrprec+2)) {

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (Larcosh 362) y   = ");  Lprinti(y);
            printf("\n  (Larcosh 363) Maxl = %d",Maxl);
            }
#endif

         if (rc != 0) ERREXIT(RESUL,367,poolvars);
         EXP_(y);                   /* Lhlp := exp(y) - xi              */
         DIV_(&Lone,&expy,invexp);
         SUB_(&expy,invexp,sinhy);
         SHIFT_(-1,sinhy,sinhy);
         ADD_(&expy,invexp,coshy);
         SHIFT_(-1,coshy,coshy);
         SUB_(coshy,xi,Lhlp);

         Maxl = Minlerr;
         NEXT_(Lhlp,&dummy);        /* determine error bound &err       */
         MUL_(sqt,coshy,summ2,&dummy);
         MUL_(xi,sinhy,summ1,&dummy);
         ADD_(summ1,summ2,dvsr);
         ADD_(xi,coshy,dvnd);
         DIV_(dvnd,dvsr,&err);
         MUL_(&err,Lhlp,&err,&dummy);
         DIV_(&err,y,&err);
         err.s = 0;

         Maxl = prec;
         DIV_(Lhlp,sinhy,&dummy);  /* y := y - [cosh(y)-xi]/sinh(y)  */
         SUB_(y,&dummy,y);

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (Larcosh 393) expy  = ");  Lprinti(&expy);
            printf("\n  (Larcosh 394) sinhy = ");  Lprinti(sinhy);
            printf("\n  (Larcosh 395) coshy = ");  Lprinti(coshy);
            printf("\n  (Larcosh 396) Lhlp  = ");  Lprinti(Lhlp);
            printf("\n  (Larcosh 397) dummy = ");  Lprinti(&dummy);
            printf("\n  (Larcosh 398) y   = ");  Lprinti(y);
            printf("\n  (Larcosh 399) err = ");  Lprinti(&err);
            }
#endif
         }

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcosh 406) Iterations performed !");
         printf("\n");
         }
#endif

      Maxl = Lcurrprec+2;
      Res = Lvarget();
      poolvars += 1;
      rc  = copyii(y,Res);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcosh 418) Res  = ");  Lprinti(Res);
         printf("\n");
         }
#endif


 /*---------------------------------------------------------------------*
  | Determine number of units to be added for an upper bound            |
  *---------------------------------------------------------------------*/

      if (rc != 0) ERREXIT(RESUL,432,poolvars);
      EXP_(Res);                       /* compute exp(Res)              */
      DIV_(&Lone,&expy,invexp);
      ADD_(&expy,invexp,coshy);
      SHIFT_(-1,coshy,coshy);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcosh 440) coshy  = ");  Lprinti(coshy);
         printf("\n  (Larcosh 441) errexp = ");  Lprinti(&errexp);
         printf("\n");
         }
#endif

      l5o2.e = 1-Maxl;
      l5o2.s = err.s;
      Maxl = Minlerr;
      SUB_(&expy,invexp,sinhy);
      SHIFT_(-1,sinhy,sinhy);
      SUB_(coshy,xi,Lhlp);
      NEXT_(Lhlp,Lhlp);

      ADD_(&err,&l5o2,&err);        /* err = err(cosh y)                */
      NEXT_(&err,&err);             /*     = err(exp y) + 2.5*b^(1-l)   */

      MUL_(&err,coshy,&err,&dummy); /* absolute error = err(cosh)*cosh  */
      NEXT_(&err,&err);
      err.s = Lhlp->s;              /* add absolute error of difference */
      ADD_(Lhlp,&err,Lhlp);
      NEXT_(Lhlp,Lhlp);

      MUL_(sqt,coshy,summ2,&dummy);
      MUL_(xi,sinhy,summ1,&dummy);
      ADD_(summ1,summ2,dvsr);
      ADD_(xi,coshy,dvnd);
      NEXT_(dvnd,dvnd);
      DIV_(dvnd,dvsr,&err);
      NEXT_(&err,&err);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larcosh 473) err   = ");  Lprinti(&err);
         printf("\n  (Larcosh 474) Lhlp  = ");  Lprinti(Lhlp);
         printf("\n");
         }
#endif

      MUL_(&err,Lhlp,&err,&dummy);
      NEXT_(&err,&err);
      DIV_(&err,y,&err);
      NEXT_(&err,&err);

      Lrnd = rounded;

      Maxl = Lcurrprec+2;
      COPY_(Res,&Lres);

      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcosh 493) Maxl = %d\n",Maxl);
      printf("\n  (Larcosh 494) Lres = ");  Lprinti(&Lres);
      printf("\n");
      }
#endif

   Lvardrop(poolvars);

   if (rc != 0) {
      ERREXIT(RESUL,502,0); /* Error handling                           */
      }

   ASSIGN_(ri);             /* Assign result and set number of ulp's    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcosh 509) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larcosh 510) LhE    = ");  Lprinti(&err);
      printf("\n  (Larcosh 511) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





