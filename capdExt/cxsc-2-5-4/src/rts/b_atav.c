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

/* CVS $Id: b_atav.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : arctanve           Processor : C                   *
 *                                                                       *
 * Inverse Internal Sine for Multiple Precision Arithmetic               *
 * =======================================================               *
 *                                                                       *
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 *                                                                       *
 * Argument Range :    |Arg| <= 1                                        *
 * Result Range :      |Res| <= pi/2 (+eps)                              *
 *                                                                       *
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
 *    sicovea   SICO_    Intern Trigonometric Sine and Cosine Function   *
 *    Lginit    (none)   Initialization of global constants              *
 *                    Intern arithmetic (macro names see file b_lari.h)  *
 *                                                                       *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhF, LhE, LhD                                        (Variables)   *
 *    Lone              ( = 1 )                            (Constant)    *
 *                                                                       *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

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
 * Constants for Computation of arctan(x)                               *
 ************************************************************************/

#define Lguard    2
#define PiMore    20
#define M0        4
#define erri      dvnd

static a_btyp  macc[3]   = { 0x00000001,0x00000000,0x00000000 };

static dynamic  acc   = { 0, 0, 0, 0, 0, Minl, &macc[0] };

#ifdef Version
static int   arctanve_sccsid = true;
#endif

/************************************************************************
 * Algorithm                                                            *
 ************************************************************************/

#ifdef LINT_ARGS
int b_atav(dynamic *xi,dynamic *LPiOv2)
#else
int b_atav(xi,LPiOv2)

dynamic *xi;
dynamic *LPiOv2;
#endif
#define  LRoutine    "arctanve"
{
   extern dynamic   Lone, LPiov4, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lgiflag;
   extern rounding  Lrnd;

   dynamic   *y, *Arg, *Res, *sqt, *siny, *cosy, *dvnd, *Lhlp;
   int       rc, poolvars;
   a_btyp    Currprec, prec;
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

   Currprec = Maxl;             /* Save precision setting               */
   rc = 0;


  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 867) Maxl = %d\n",Maxl);
      printf("\n  (arctanve 868) xi   = ");  Lprinti(xi);
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
   if (rc != 0) ERREXIT(CONVD,890,poolvars);

   if GT_(Arg,&Lone) {                  /* Argument reduction           */
      Maxl = Currprec + Lguard;
      DIV_(&Lone,Arg,Arg);
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

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 916) Arg = ");  Lprinti(Arg);
      }
#endif

   Maxl = Minlfl;

   if (Arg->e >= EFUFfac) {
      rc = b_bcid(Arg,&yd,(a_intg)0);       /* y = arctan(Arg)   */
      if (rc==ROUND) rc = 0;
      rc += b_bcdi(B_ATAN(yd),&y,(a_intg)0);
      }
    else
      {
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 918) COPY_ ");
      }
#endif
      COPY_(Arg,y);                        /* y = |xi|          */
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 919) y = arctan(|xi|) = ");  Lprinti(y);
      }
#endif


  /**********************************************************************
   * Newton Iteration Steps to gain the required accuracy               *
   **********************************************************************/

  /* The accuracy of the initial approximation is increased by a number
     of Newton steps until the relative error bound is less than 1% of
     the initial accuracy setting.                                      */

   acc.e = -Currprec;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 940) acc = ");  Lprinti(&acc);
      }
#endif

   Maxl = Minlerr;
   Lhlp = Lvarget();
   sqt  = Lvarget();
   poolvars += 2;

   MUL_(Arg,Arg,sqt,&dummy);            /* sqt = v(sqrt(1 + Arg^2)) */
   ADD_(&Lone,sqt,sqt);
   if (rc != 0) ERREXIT(951,951,poolvars);
   SQRT_(sqt);
   if (err.s == 1) {
      MUL_(&sqty,&err,Lhlp,&dummy);
      NEXT_(Lhlp,Lhlp);
      SUB_(&sqty,Lhlp,sqt);
      }
    else
      COPY_(&sqty,sqt);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 963) sqt = ");  Lprinti(sqt);
      }
#endif

   COPY_(&acc,&err);
   NEXT_(&err,&err);
   err.s = 0;

   siny = Lvarget();
   cosy = Lvarget();
   dvnd = Lvarget();
   poolvars += 4;
   COPY_(&Lone,erri);

   for (prec=Maxl=M0 ;;
        prec=Maxl=min(2*Maxl,min(max(M0,-2*err.e),Currprec)+Lguard)) {

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arctanve 982) y   = ");  Lprinti(y);
         printf("\n  (arctanve 983) err = ");  Lprinti(&err);
         printf("\n  (arctanve 984) Maxl = %d",Maxl);
         }
#endif

      if (rc != 0) ERREXIT(RESUL,988,poolvars);
      SICO_(y);                          /* LhF = sin(y), LhD = cos(y)  */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arctanve 993) err_sc = ");  Lprinti(&err);
         }
#endif

      COPY_(&LhF,siny);
      COPY_(&LhD,cosy);
      MUL_(Arg,cosy,Lhlp,&dummy);        /* Lhlp = sin(y) - Arg*cos(y)  */
      if (LT_(erri,&acc) || (EQ_(siny,Lhlp) && (Maxl>=Currprec+Lguard)))
         break;
      SUB_(siny,Lhlp,Lhlp);

      Maxl = Minlerr;                    /* determine error bound &err  */
      if (Lhlp->z == 1) {
         MUL_(&err,siny,&dummy,erri);
         SHIFT_(1,&dummy,&dummy);
         }
       else
         NEXT_(Lhlp,&dummy);
      DIV_(&dummy,sqt,&err);
      MUL_(&err,LPiOv2,&err,&dummy);
      DIV_(&err,y,&err);
      err.s = 0;
      COPY_(&err,erri);

      Maxl = prec;
           /* y:=y-cos(y)*[sin(y)-Arg*cos(y)] */
      MUL_(Lhlp,cosy,&LhF,&dummy);
      SUB_(y,&LhF,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (arctanve 1023) siny = ");  Lprinti(siny);
         printf("\n  (arctanve 1024) cosy = ");  Lprinti(cosy);
         printf("\n  (arctanve 1025) Lhlp = ");  Lprinti(Lhlp);
         printf("\n  (arctanve 1026) LhF  = ");  Lprinti(&LhF);
         }
#endif
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 1033) Iterations performed !");
      printf("\n");
      }
#endif

   Maxl = Currprec + Lguard;
   Res = Lvarget();
   poolvars += 1;
   COPY_(y,Res);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 1045) Res  = ");  Lprinti(Res);
      printf("\n");
      }
#endif


  /**********************************************************************
   * Determine number of units to be added for an upper bound           *
   **********************************************************************/

   if (rc != 0) ERREXIT(RESUL,1059,poolvars);

   Maxl = Minlerr;
   NEXT_(&dummy,&dummy);        /* dummy = x*cos(y) (trailing part) */
   dummy.s = 0;

   err.s = 0;                   /* err = absolute error of sin(y)   */

   SUB_(Lhlp,siny,dvnd);    /* dvnd = ^|x*cos(y) - sin(y)|      */
   NEXT_(dvnd,dvnd);
   dvnd->s = 0;
   ADD_(dvnd,&dummy,dvnd);      /* add trailing part of |x*cos| */
   NEXT_(dvnd,dvnd);
   Lhlp->s = 0;             /*    use ^|x*cos(y)| from now on   */
   NEXT_(Lhlp,Lhlp);
   MUL_(Lhlp,&err,Lhlp,&dummy); /* multiply |x*cos| by rel. err. */
   NEXT_(Lhlp,Lhlp);
   ADD_(dvnd,Lhlp,dvnd);        /* add error of |x*cos(y)|       */
   NEXT_(dvnd,dvnd);
   MUL_(siny,&err,Lhlp,&dummy); /* multiply sin(y) by rel. err.  */
   NEXT_(Lhlp,Lhlp);
   Lhlp->s = 0;
   ADD_(dvnd,Lhlp,dvnd);        /* add error of |sin(y)|         */
   NEXT_(dvnd,dvnd);

   DIV_(dvnd,sqt,&err);         /* err = [|x*cos(y) - sin(y)| /     */
   NEXT_(&err,&err);            /*  sqrt(1+Arg^2)]*(pi/2)/y   */
   MUL_(&err,LPiOv2,&err,&dummy);
   NEXT_(&err,&err);
   DIV_(&err,y,&err);
   NEXT_(&err,&err);

   Maxl = Currprec + Lguard;
   COPY_(Res,&Lres);
   Lres.s = xi->s;              /* Set sign of result           */

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 1099) dummy = ");  Lprinti(&dummy);
      printf("\n  (arctanve 1100) err = ");  Lprinti(&err);
      }
#endif

   Lrnd = rounded;

   if (rc != 0) {
      ERREXIT(RESUL,1107,0);    /* Error message and handling   */
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (arctanve 1112) LhF    = ");  Lprinti(&Lres);
      printf("\n  (arctanve 1113) LhE    = ");  Lprinti(&err);
      }
#endif

   RETURN(0);
}





