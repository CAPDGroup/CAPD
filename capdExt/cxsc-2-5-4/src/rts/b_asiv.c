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

/* CVS $Id: b_asiv.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/************************************************************************
 *                                                                      *
 * Descriptive Name : arcsinve           Processor : C                  *
 *                                                                      *
 * Inverse Internal Sine for Multiple Precision Arithmetic              *
 * =======================================================              *
 *                                                                      *
 *                                                                      *
 * Function value : int       - 0                                       *
 *                                                                      *
 *                                                                      *
 * Used Number Base :  2**32                                            *
 *                                                                      *
 *                                                                      *
 * Argument Range :    |Arg| <= 1                                       *
 * Result Range :      |Res| <= pi/2 (+eps)                             *
 *                                                                      *
 *                                                                      *
 * Used Functions :                                                     *
 *    On computation of function values the computed value and the error*
 *    bound are returned by the global variables LhF and LhE.   However,*
 *    in many cases these names are redefined to a symbolic name related*
 *    to the called function.                                           *
 *                                                                      *
 *    name      macro   description                                     *
 *    --------  ------   -----------------------------------------------*
 *    sqrtve    SQRT_    Intern Square Root Function                    *
 *    sicovea   SICO_    Intern Trigonometric Sine and Cosine Function  *
 *    Lginit    (none)   Initialization of global constants             *
 *                   Intern arithmetic (macro names see file b_lari.h)  *
 *                                                                      *
 *                                                                      *
 * Used Global INTERN Variables and Constants :                         *
 *    LhF, LhE, LhD                                        (Variables)  *
 *    Lone              ( = 1 )                            (Constant)   *
 *                                                                      *
 *                                                                      *
 * Used UNSIGNED Global Variables :                                     *
 *    Maxl            Length of INTERN Variables                        *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables*
 *    Ldebug          Flag for printing additional information          *
 *                                                                      *
 ************************************************************************/


#undef  Main
#undef  Name
#undef  LRoutine
#undef  poolvars

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

/************************************************************************
 * Redefinition of Names of Intermediate Variables Used                 *
 ************************************************************************/

#define Lres      LhF           /* Function value                       */
#define err       LhE           /* Rounding error                       */
#define sqty      LhF           /* Function value sqrt(y)               */
#define dummy     LhD           /* Dummy variable (not further used)    */

/************************************************************************
 * Constants for Computation of arcsin(x)                               *
 ************************************************************************/

#define Lguard    2
#define PiMore    20
#define M0        4
#define summ2     dvsr
#define summ1     dvnd
#define erri      dvnd

static a_btyp  macc[3]   = { 0x00000001L,0x00000000L,0x00000000L };
static a_btyp  msqth[3]  = { 0xB504F333L,0xF9DE6484L,0x597D89B4L };

static dynamic  acc   = { 0, 0, 0, 0, 0, Minl, &macc[0] };
static dynamic  Lsqth = { 0, 0, 0, 0,-1, Minlfl, &msqth[0] };

/************************************************************************
 * Algorithm                                                            *
 ************************************************************************/

#ifdef LINT_ARGS
int b_asiv(dynamic *xi,dynamic *LPiOv2)
#else
int b_asiv(xi,LPiOv2)

dynamic *xi;
dynamic *LPiOv2;
#endif
#define  LRoutine    "arcsinve"
{
   extern dynamic   Lone, LPiov4, Lres, err, dummy;
   extern a_btyp    Maxl;
   extern a_btyp    Lgiflag;
   extern rounding  Lrnd;

   dynamic   *y, *Arg, *Res, *sqt, *siny, *cosy, *dvsr, *dvnd, *Lhlp;
   int       rc, poolvars;
   a_btyp    Currprec, prec;
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

   Currprec = Maxl;             /* Save precision setting       */
   rc = 0;


  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 635) Maxl = %d\n",Maxl);
      printf("\n  (arcsinve 636) xi   = ");  Lprinti(xi);
      }
#endif


  /**********************************************************************
   * Computation of an Initial Approximation                            *
   **********************************************************************/

  /* An initial approximation is computed by taking the floating point
     logarithm of the argument xi.                              */

   Arg = Lvarget();
   y   = Lvarget();
   poolvars = 2;
   Maxl = xi->l;
   COPY_(xi,Arg);
   Arg->s = 0;
   if (rc != 0) ERREXIT(CONVD,658,poolvars);

   if GT_(Arg,&Lsqth) {                 /* Argument reduction    */
      Maxl = Currprec + Lguard;
      SUB_(&Lone,Arg,y);
      ADD_(&Lone,Arg,&dummy);
      MUL_(y,&dummy,y,&dummy);
      if (rc != 0) ERREXIT(CONVD,665,poolvars);
      SQRT_(y);
      COPY_(&sqty,Arg);
      if (LPiOv2->m == NULL || LPiOv2->e < Currprec+Lguard+1) {
         Maxl++;
         SHIFT_(1,&LPiov4,LPiOv2);
         }
      }
    else {
      if (LPiOv2->m == NULL) {
         Maxl = Minlerr;
         SHIFT_(1,&LPiov4,LPiOv2);
         }
      }

   Maxl = Minlfl;

   if (Arg->e >= EFUFfac) {
      rc = b_bcid(Arg,&yd,(a_intg)0);   /* y = arcsin(Arg)   */
      if (rc==ROUND) rc = 0;
      rc += b_bcdi(B_ASIN(yd),&y,(a_intg)0);
      }
    else
      COPY_(Arg,y);                        /* y = |xi|          */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 692) y = arcsin(|xi|) = ");  Lprinti(y);
      }
#endif


  /**********************************************************************
   * Newton Iteration Steps to gain the required accuracy               *
   **********************************************************************/

  /* The accuracy of the initial approximation is increased by a number
     of Newton steps until the relative error bound is less than 1% of
     the initial accuracy setting.                               */

   acc.e = -Currprec;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 713) acc = ");  Lprinti(&acc);
      }
#endif

   Maxl = Minlerr;
   Lhlp = Lvarget();
   sqt  = Lvarget();
   poolvars += 2;

   if LE_ABS_(xi,&Lsqth) {              /* sqt = v(sqrt(1 - Arg^2)) */
      MUL_(Arg,Arg,sqt,&dummy);
      SUB_(&Lone,sqt,sqt);
      if (rc != 0) ERREXIT(725,725,poolvars);
      SQRT_(sqt);
      if (err.s == 1) {
         MUL_(&sqty,&err,Lhlp,&dummy);
         NEXT_(Lhlp,Lhlp);
         SUB_(&sqty,Lhlp,sqt);
         }
       else
         COPY_(&sqty,sqt);
      }
    else {
      COPY_(xi,sqt);            /* = |xi| if Arg=sqrt(1-xi^2)   */
      sqt->s = 0;
      }

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arcsinve 742) sqt = ");  Lprinti(sqt);
         }
#endif

   COPY_(&acc,&err);
   NEXT_(&err,&err);
   err.s = 0;

   siny = Lvarget();
   cosy = Lvarget();
   dvsr = Lvarget();
   dvnd = Lvarget();
   poolvars += 4;
   COPY_(&Lone,erri);

   for (prec=Maxl=M0 ;;
        prec=Maxl=min(2*Maxl,min(max(M0,-2*err.e),Currprec)+Lguard)) {

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arcsinve 762) y   = ");  Lprinti(y);
         printf("\n  (arcsinve 763) err = ");  Lprinti(&err);
         printf("\n  (arcsinve 764) Maxl = %d",Maxl);
         }
#endif

      if (rc != 0) ERREXIT(RESUL,768,poolvars);
      SICO_(y);                 /* LhF = sin(y), LhD = cos(y)  */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arcsinve 773) err_sc = ");  Lprinti(&err);
         }
#endif

      if (LT_(erri,&acc) || (EQ_(&LhF,Arg) && (Maxl>=Currprec+Lguard)))
         break;
      COPY_(&LhF,siny);
      COPY_(&LhD,cosy);
      SUB_(siny,Arg,Lhlp);

      Maxl = Minlerr;                    /* determine error bound &err  */
      MUL_(sqt,siny,summ2,&dummy);
      MUL_(Arg,cosy,summ1,&dummy);
      ADD_(summ1,summ2,dvsr);
      ADD_(Arg,siny,dvnd);
      if (Lhlp->z == 1) {
         MUL_(&err,siny,&dummy,&LhF);
         }
       else
         NEXT_(Lhlp,&dummy);
      DIV_(dvnd,dvsr,&err);
      MUL_(&err,&dummy,&err,&dummy);
      MUL_(&err,LPiOv2,&err,&dummy);
      DIV_(&err,y,&err);
      err.s = 0;
      COPY_(&err,erri);

      Maxl = prec;
      DIV_(Lhlp,cosy,&dummy);  /* y := y - [sin(y)-Arg]/cos(y)*/
      SUB_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arcsinve 806) siny  = ");  Lprinti(siny);
         printf("\n  (arcsinve 807) cosy  = ");  Lprinti(cosy);
         printf("\n  (arcsinve 808) Lhlp  = ");  Lprinti(Lhlp);
         printf("\n  (arcsinve 809) dummy = ");  Lprinti(&dummy);
         }
#endif
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 816) Iterations performed !");
      printf("\n");
      }
#endif

   Maxl = Currprec + Lguard;
   Res = Lvarget();
   poolvars += 1;
   COPY_(y,Res);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 828) Res  = ");  Lprinti(Res);
      printf("\n");
      }
#endif


  /**********************************************************************
   * Determine number of units to be added for an upper bound           *
   **********************************************************************/

   if (rc != 0) ERREXIT(RESUL,842,poolvars);

   Maxl = Minlerr;
   err.s = 0;                   /* err = absolute error of sin(y)   */

   MUL_(&LhD,&err,Lhlp,cosy);   /* cosy = v|cos(y)*(1-err)|         */
   SUB_(&LhD,Lhlp,cosy);

   MUL_(&LhF,&err,Lhlp,&dummy); /* Lhlp = sin(y)*err                */

   Maxl = Currprec + Lguard;
   SUB_(&LhF,Lhlp,siny);        /* siny = v|sin(y)*(1-err)|         */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 857) siny = ");  Lprinti(siny);
      printf("\n  (arcsinve 858) cosy = ");  Lprinti(cosy);
      }
#endif

   Maxl = Minlerr;
   NEXT_(Lhlp,Lhlp);            /* dvnd = ^|Arg + ^sin(y)*(1+err)|  */
   ADD_(&LhF,Lhlp,dvnd);
   NEXT_(dvnd,dvnd);
   ADD_(Arg,siny,dvnd);
   NEXT_(dvnd,dvnd);

   MUL_(sqt,siny,summ2,&dummy); /* dvsr = v|Arg*cos(y)+sqt*sin(y)|  */
   MUL_(Arg,cosy,summ1,&dummy);
   ADD_(summ1,summ2,dvsr);

   DIV_(dvnd,dvsr,&dummy);      /* err = |dvnd / dvsr|              */
   NEXT_(&dummy,&dummy);
   dummy.s = 0;

   SUB_(Arg,siny,siny);         /* Lhlp = ^|Arg - sin(y)|           */
   if (siny->z != 1)
      COPY_(siny,Lhlp);
   NEXT_(Lhlp,Lhlp);
   siny->s = 0;

   MUL_(&dummy,Lhlp,&err,&dummy);
   NEXT_(&err,&err);
   MUL_(&err,LPiOv2,&err,&dummy);
   NEXT_(&err,&err);
   DIV_(&err,y,&err);
   NEXT_(&err,&err);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 892) Lhlp = ");  Lprinti(Lhlp);
      }
#endif

   Maxl = Currprec + Lguard;
   COPY_(Res,&Lres);
   Lres.s = xi->s;              /* Set sign of result               */

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 904) dummy = ");  Lprinti(&dummy);
      printf("\n  (arcsinve 905) err = ");  Lprinti(&err);
      }
#endif

   Lrnd = rounded;

   if (rc != 0) {
      ERREXIT(RESUL,912,0);     /* Error message and handling    */
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arcsinve 917) LhF    = ");  Lprinti(&Lres);
      printf("\n  (arcsinve 918) LhE    = ");  Lprinti(&err);
      }
#endif

   RETURN(0);
}





