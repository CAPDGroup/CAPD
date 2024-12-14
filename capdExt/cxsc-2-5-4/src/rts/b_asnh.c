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

/* CVS $Id: b_asnh.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Larsinh            Processor : C                   *
 *                                                                       *
 * Inverse Hyperbolic Sine for Multiple Precision Arithmetic             *
 * =========================================================             *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Argument Range :    All numbers                                       *
 * Result Range :      All numbers                                       *
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
 *    expve     EXP_     Intern Exponential Function                     *
 *    sinhvea   SINH_    Intern Hyperbolic Sine Function                 *
 *    Lassign   ASSIGN_  Assignment of function value and number of ulp's*
 *    Lginit    (none)   Initialization of global constants              *
 *                    Intern arithmetic (macro names see file b_lari.h)  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhF, LhE, LhD                                        (Variables)   *
 *    Lone              ( = 1 )                            (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#define  Name        "Larsinh"
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
#define sqty      LhF           /* Function value sqrt(y)               */
#define expy      LhF           /* Function value exp(y)                */
#define errexp    LhE           /* Rounding error of exp(y)             */
#define errsinh   LhE           /* Rounding error of sinh(y)            */
#define dummy     LhD           /* Dummy variable (not further used)    */

/************************************************************************
 * Constants for Computation of arsinh(x)                               *
 ************************************************************************/

#define Lguard    2
#define M0        4
#define ebdb      1
#define invexp    dvsr
#define summ2     dvsr
#define summ1     dvnd

static a_btyp  macc[3]   = { 0x00000001L,0x00000000L,0x00000000L };
static a_btyp  mbes[3]   = { 0x80000000L,0x00000000L,0x00000000L };
static a_btyp  mbds[3]   = { 0x00010000L,0x00000000L,0x00000000L };
static a_btyp  m5o2[3]   = { 0x00000002L,0x80000001L,0x00000000L };

static dynamic  acc   = { 0, 0, 0, 0, 0, Minl, &macc[0] };
static dynamic  bds   = { 0, 0, 0, 0,-1, Minl, &mbds[0] };
static dynamic  bes   = { 0, 0, 0, 0,-1, Minl, &mbes[0] };
static dynamic  besl  = { 0, 0, 0, 0,-2, Minl, &macc[0] };
static dynamic  l5o2  = { 0, 0, 0, 0, 0, Minlerr, &m5o2[0] };

/************************************************************************
 * Algorithm                                                            *
 ************************************************************************/

#ifdef LINT_ARGS
int b_asnh(dynamic *xi,dynamic *ri)
#else
int b_asnh(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larsinh"
{
   extern dynamic   Lone, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;
   extern rounding  Lrnd;
   extern a_real      Fln2, Flnb, *r_one_;

   dynamic   *y, *Arg, *Res, *sqt, *sinhy, *coshy, *dvsr, *dvnd, *Lhlp;
   int       rc, poolvars;
   a_btyp    prec;
   a_real    yd;

#ifdef Debug
   extern int       Ldebug;
#endif

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level                 */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting               */
   rc = 0;

  /**********************************************************************
   * Check and Handle Range Error and Special Cases                     *
   **********************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z == 1) {
      if ((rc = b_bini(ri)) == 0) {        /* ri = arsinh(0) = 0        */
         ri->r = 0;                        /* Result exact              */
         EXIT(0);
         }
       else
         ERREXIT(rc,rc,0);                 /* Error message and handling*/
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling           */
      }



  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 202) Maxl = %d\n",Maxl);
      printf("\n  (Larsinh 203) xi   = ");  Lprinti(xi);
      }
#endif


  /**********************************************************************
   * Computation of an Initial Approximation                            *
   **********************************************************************/

  /* An initial approximation is computed by taking the floating point
     logarithm of the argument xi.                                      */

   Arg = Lvarget();
   y   = Lvarget();
   poolvars = 2;
   Maxl = xi->l;
   COPY_(xi,Arg);
   Arg->s = 0;
   if (rc != 0) ERREXIT(CONVD,225,poolvars);

   Maxl = Minlfl;

   if LT_(Arg,&bds)
      COPY_(Arg,y);                        /* y = |xi|                  */
    else if (Arg->e >= ebdb) {
      Arg->e = 0;                          /* set exponent of Arg to 0  */
      rc = b_bcid(Arg,&yd,(a_intg)0); /* y = ln(mant(Arg)) + ln2 + */
                                          /*   + expo(xi)*lnb          */
      if (rc==ROUND) rc = 0;
/* Cordes      rc += b_bcdi(log(yd)+Fln2+xi->e*Flnb,y); */
      rc += b_bcdi(r_addn(r_addn(B_LOG_(yd),Fln2),
                   r_muln(r_flot((a_intg)xi->e),Flnb)
                   ),
                 &y,(a_intg)0);

      Arg->e = xi->e;                      /* restore exponent of Arg   */
      }
    else {
      rc = b_bcid(Arg,&yd,(a_intg)0);  /* y = ln(Arg) +     */
                                       /*   + ln(1 + sqrt(1/(Arg*Arg))) */
      if (rc==ROUND) rc = 0;
/* Cordes      rc += b_bcdi(log(yd)+log(1+sqrt(1+1/(yd*yd))),y);        */
      rc += b_bcdi(r_addn(B_LOG_(yd),
                   B_LOG_(r_addn(*r_one_,B_SQRT(r_addn(*r_one_,
                                           r_divn(*r_one_,r_muln(yd,yd))
                                           ))
                       ))
                   ),
                &y,(a_intg)0);
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 246) y = arsinh(|xi|) = ");  Lprinti(y);
      }
#endif


  /**********************************************************************
   * Newton Iteration Steps for gaining the required accuracy           *
   **********************************************************************/

  /* The accuracy of the initial approximation is increased by a number
     of Newton steps until the relative error bound is less than 1% of
     the initial accuracy setting.                                      */

   acc.e = -Lcurrprec;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 267) acc = ");
      Lprinti(&acc);
      }
#endif

   Maxl = Minlerr;
   Lhlp = Lvarget();
   sqt  = Lvarget();
   poolvars += 2;

   MUL_(Arg,Arg,sqt,&dummy);
   ADD_(sqt,&Lone,sqt);
   if (rc != 0) ERREXIT(RESUL,279,poolvars);
   SQRT_(sqt);
   COPY_(&sqty,sqt);

   COPY_(&acc,&err);
   NEXT_(&err,&err);
   err.s = 0;

   sinhy = Lvarget();
   coshy = Lvarget();
   dvsr  = Lvarget();
   dvnd  = Lvarget();
   poolvars += 4;

   for (prec=Maxl=M0; GT_(&err,&acc);
        prec=Maxl=min(2*Maxl,Lcurrprec+2)) {
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Larsinh 296) y   = ");  Lprinti(y);
         printf("\n  (Larsinh 297) Maxl = %d",Maxl);
         }
#endif

      if (rc != 0) ERREXIT(RESUL,301,poolvars);
      EXP_(y);                          /* Lhlp := exp(y) - Arg         */
      DIV_(&Lone,&expy,invexp);
      ADD_(&expy,invexp,coshy);
      SHIFT_(-1,coshy,coshy);
      besl.e = 2 - prec;
      if GE_(y,&besl) {
         SUB_(&expy,invexp,sinhy);
         SHIFT_(-1,sinhy,sinhy);
         }
       else {
         if (rc != 0) ERREXIT(RESUL,312,poolvars);
         SINH_(y);
         COPY_(&Lres,sinhy);
         }
      if (EQ_(sinhy,Arg) & (Maxl >= Lcurrprec+2)) break;
      SUB_(sinhy,Arg,Lhlp);

      Maxl = Minlerr;                     /* determine error bound &err */
      MUL_(sqt,sinhy,summ2,&dummy);
      MUL_(Arg,coshy,summ1,&dummy);
      ADD_(summ1,summ2,dvsr);
      ADD_(Arg,sinhy,dvnd);
      DIV_(dvnd,dvsr,&err);
      if (Lhlp->z == 1) {
         COPY_(&Lone,&dummy);
         dummy.e = 1 - prec;
         }
       else
         {
         NEXT_(Lhlp,&dummy);
         }
      MUL_(&err,&dummy,&err,&dummy);
      DIV_(&err,y,&err);
      err.s = 0;

      Maxl = prec;
      DIV_(Lhlp,coshy,&dummy); /* y := y - [sinh(y)-Arg]/cosh(y) */
      SUB_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
/*       printf("\n  (Larsinh 341) expy  = ");  Lprinti(&expy);  */
         printf("\n  (Larsinh 342) sinhy = ");  Lprinti(sinhy);
         printf("\n  (Larsinh 343) coshy = ");  Lprinti(coshy);
         printf("\n  (Larsinh 344) Lhlp  = ");  Lprinti(Lhlp);
         printf("\n  (Larsinh 345) dummy = ");  Lprinti(&dummy);
         printf("\n  (Larsinh 346) y   = ");  Lprinti(y);
         printf("\n  (Larsinh 347) err = ");  Lprinti(&err);
         }
#endif
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 354) Iterations performed !");
      printf("\n");
      }
#endif

   Maxl = Lcurrprec+2;
   Res = Lvarget();
   poolvars += 1;
   COPY_(y,Res);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 366) Res  = ");  Lprinti(Res);
      printf("\n");
      }
#endif


  /**********************************************************************
   * Determine number of units to be added for an upper bound           *
   **********************************************************************/

   if (rc != 0) ERREXIT(RESUL,380,poolvars);

   if GE_(Res,&bes) {
      EXP_(Res);                        /* compute exp(Res)             */
      DIV_(&Lone,&expy,invexp);
      SUB_(&expy,invexp,sinhy);
      SHIFT_(-1,sinhy,sinhy);

      l5o2.e = 1-Maxl;                  /* err = err(sinh y)            */
      l5o2.s = err.s;                   /*   = err(exp y) + 2.5*b^(1-l) */
      Maxl = Minlerr;
      NEXT_(&err,&err);
      ADD_(&err,&l5o2,&err);
      }
    else {
      SINH_(Res);
      COPY_(&Lres,sinhy);
      Maxl = Minlerr;
      COPY_(&errsinh,Lhlp);
      if (rc != 0) ERREXIT(RESUL,399,poolvars);
      EXP_(Res);                        /* compute exp(Res)             */
      DIV_(&Lone,&expy,invexp);
      COPY_(Lhlp,&err);                 /* err = err(sinh y)            */
      }

   NEXT_(&err,&err);                    /* err = err(sinh y) * sinh y   */
   MUL_(&err,sinhy,&err,&dummy);
   NEXT_(&err,&err);

   Maxl = Lcurrprec+2;
   COPY_(Res,&Lres);
   Lres.s = xi->s;                      /* Set sign of result           */

   Maxl = Minlerr;
   ADD_(&expy,invexp,coshy);
   SHIFT_(-1,coshy,coshy);
   SUB_(sinhy,Arg,Lhlp);
   err.s = Lhlp->s;             /* add absolute error of difference */
   ADD_(Lhlp,&err,Lhlp);
   NEXT_(Lhlp,Lhlp);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 423) expy  = ");  Lprinti(&expy);
      printf("\n  (Larsinh 424) sinhy = ");  Lprinti(sinhy);
      printf("\n  (Larsinh 425) coshy = ");  Lprinti(coshy);
      printf("\n  (Larsinh 426) Lhlp  = ");  Lprinti(Lhlp);
      }
#endif

   MUL_(sqt,sinhy,summ2,&dummy);
   MUL_(Arg,coshy,summ1,&dummy);
   ADD_(summ1,summ2,dvsr);
   ADD_(Arg,sinhy,dvnd);
   NEXT_(dvnd,dvnd);
   DIV_(dvnd,dvsr,&err);
   NEXT_(&err,&err);
   MUL_(&err,Lhlp,&err,&dummy);
   NEXT_(&err,&err);
   DIV_(&err,y,&err);
   NEXT_(&err,&err);

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 446) dummy = ");  Lprinti(&dummy);
      printf("\n  (Larsinh 447) err = ");  Lprinti(&err);
      }
#endif

   Lrnd = rounded;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 455) Maxl = %d\n",Maxl);
      printf("\n  (Larsinh 456) xi   = ");  Lprinti(xi);
      printf("\n  (Larsinh 457) Lres = ");  Lprinti(&Lres);
      printf("\n  (Larsinh 458) err  = ");  Lprinti(&err);
      printf("\n");
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,464,0);     /* Error message and handling    */
      }

#if INT_HPREC
   b_farg = xi;
   b_case = LESS_ABS_ARG;
#endif

   ASSIGN_(ri);         /* Assign result and set number of ulp's    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larsinh 471) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larsinh 472) LhE    = ");  Lprinti(&err);
      printf("\n  (Larsinh 473) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





