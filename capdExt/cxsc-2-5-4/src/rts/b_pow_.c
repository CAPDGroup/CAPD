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

/* CVS $Id: b_pow_.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lpower             Processor : C                   *
 *                                                                       *
 * Square Root Function for Multiple Precision Arithmetic                *
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

#define  Name        "Lpower"
#define  Main

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
#define dummy     LhD           /* Dummy variable (not further used)    */

/************************************************************************
 * Constants for Computation of Power Function                          *
 ************************************************************************/

#define  Lguard    3
#undef   BIT0
#define  BIT0      1

#define  exponent  n
#define  lnx       y
#define  power     x
#define  MAXDIG    0x7FFFFFFFL

static a_btyp  mbapp[3]  = { 0x18000000L,0x00000000L,0x00000000L };

static dynamic   bdapp   = { 0, 0, 0, 0, -1, Minl, &mbapp[0] };

static char  *function = Name;

/************************************************************************
 * Algorithm                                                            *
 ************************************************************************/

#ifdef LINT_ARGS
int b_pow_(dynamic *xi,dynamic *yi,dynamic *ri)
#else
int b_pow_(xi,yi,ri)

dynamic *xi;
dynamic *yi;
dynamic *ri;
#endif
#define  LRoutine    "Lpower"
{
   extern dynamic   Lone, LhF, LhE, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;
   extern rounding  Lrnd;

   dynamic    *y, *x, *errln;
   a_btyp   Guarddig, ulp, ulppwr, n;
   int        rc, poolvars, ysgn, xsgn, rsgn;

#ifdef Debug
   extern int       Ldebug;
#endif

  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level                 */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;                    /* Save precision setting       */
   rsgn = 0;
   xsgn = xi->s;
   rc = 0;


  /*----------------------------*
   | Check and Handle Case yi=0 |
   *----------------------------*/

   if (yi->z != 0) {
      if (xi->z != 0)
         ERREXIT(RANGE,RANGE,0);        /* 0^0 undefined        */
      if ((COPY_(&Lone,ri)) != 0)       /* ri = xi^0 = 1        */
         ERREXIT(rc,rc,0);
      ri->r = 0;                        /* Result exact         */
      EXIT(0);
      }


  /*----------------------------------------------------*
   | Check and Handle Undefined Mantissa of Exponent yi |
   *----------------------------------------------------*/

   if (yi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);           /* Error message and handling   */
      }


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z != 0) {
      if (yi->s != 0)
         ERREXIT(RANGE,RANGE,0);        /* 0^(-|yi|) undefined  */
      ri->z = 1;                        /* 0^y = 0              */
      ri->r = 0;                        /* Result exact         */
      EXIT(0);
      }


  /*------------------------------------------------*
   | Check and Handle Undefined Mantissa of Base xi |
   *------------------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);           /* Error message and handling   */
      }

   power = Lvarget();
   poolvars = 1;

   Maxl = xi->l;
   COPY_(xi,x);                         /* power == x           */
   x->s = 0;
   Maxl = Lcurrprec;


  /*-------------------------------------*
   | Check and Handle Entier Exponent yi |
   *-------------------------------------*/

   if (yi->e < 0)  goto EXPLN;

   if (yi->l > yi->e+1)                 /* Fraction digits present      */
      for (n=yi->e+1; n<yi->l; n++)
         if (yi->m[n] != 0)  goto EXPLN;

   rsgn = (int)(yi->m[yi->e] & BIT0)*xi->s;  /* Sign of Result       */

   xsgn = 0;                            /* Accept negative argument     */
   if (yi->e > 0)  goto EXPLN;

   exponent = yi->m[0];
   if (exponent >= MAXDIG)
      Maxl = Lcurrprec + 3;
    else
      Maxl = Lcurrprec + 2;

   if ((exponent & BIT0) != 0) {
      COPY_(x,&Lres);
      }
    else
      COPY_(&Lone,&Lres);

   ulppwr = 1;
   ulp    = 0;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpower 227) exponent = %8.8x",exponent);
      }
#endif

   for (exponent>>=1; exponent!=0; exponent>>=1) {
      MUL_(power,power,power,&dummy);
      ulppwr *= 2;                      /* (2^i)-1 ulp from power       */
      if ((exponent & BIT0) != 0) {
                                /* Multiply by current power of xi  */
         MUL_(&Lres,power,&Lres,&dummy);
         ulp += ulppwr;         /* (2^i)-1 ulp from power + 1 ulp   */
         }                      /*   due to multiplication          */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lpower 241) exponent = %8.8x",exponent);
         printf("\n  (Lpower 242) ulppwr = %8.8x,  ulp = %8.8x",ulppwr,ulp);
         }
#endif

      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpower 250) Product = ");  Lprinti(&Lres);
      printf("\n  (Lpower 251)  ulp = %8.8x",ulp);
      }
#endif

   COPY_(&Lone,&err);
   err.m[0] = ulp+1;    /* 1 additional ulp     */
   err.e = 1-Maxl;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpower 261) err     = ");  Lprinti(&err);
      }
#endif

   if (yi->s != 0) {
      DIV_(&Lone,&Lres,&Lres);  /* Compute outer bound              */
      NEXT_(&Lres,&Lres);
      err.m[0]++;       /* 1 additional ulp due to division */
      Lrnd = obound;            /* Computed Value outer bound       */
      }
    else
      Lrnd = ibound;            /* Computed Value is truncated      */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n Lres = ");  Lprinti(&Lres);
      printf("\n err  = ");  Lprinti(&err);
      }
#endif

   goto RESULT;


  /*-----------------------------------*
   | Check and Handle Negative Base xi |
   *-----------------------------------*/

EXPLN:
   if (xsgn != 0) {
      ERREXIT(RANGE,RANGE,poolvars); /* Error message and handling  */
      }


  /*---------------------------------------------------*
   | Check and Handle Different Ranges for Exponent YI |
   *---------------------------------------------------*/

   y = Lvarget();
   poolvars += 1;

   SUB_(x,&Lone,y);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpower 305) base-1 = ");  Lprinti(y);
      }
#endif

   Guarddig = max(yi->e,0) + Lguard;
   ysgn = y->s;
   y->s = 0;

   if LE_(y,&bdapp) {           /* Approxmation if |y| > 0.09375 */

     /*------------------------------------------*
      | Initializations for Taylor Approximation |
      *------------------------------------------*/

      y->s = ysgn;              /* y = (1 - base) / (1 + base)   */
      Maxl = min(Lcurrprec,x->l) + 1;
      ADD_(x,&Lone,&dummy);     /* EXACT sum                     */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lpower 325) Maxl = %d",Maxl);
         printf("\n  (Lpower 326) base+1 = ");  Lprinti(&dummy);
         }
#endif

      Maxl = Lcurrprec + Guarddig + 2;
      DIV_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lpower 335) Maxl = %d",Maxl);
         printf("\n  (Lpower 336) (base-1)/(base+1) = ");  Lprinti(y);
         }
#endif

      if (rc != 0) ERREXIT(PEVAL,340,poolvars);

      LN_T_(y);         /* Compute function value and error bound   */
      }
    else {
      Maxl = Lcurrprec + Guarddig;
      LN_(x);           /* Compute function value and error bound   */
      }

   Maxl = Lres.l;
   COPY_(&Lres,lnx);

   Maxl = Lcurrprec + Guarddig;
   MUL_(yi,&Lres,y,&dummy);             /* y = yi*ln(|xi|)      */

   Maxl = Minlerr;
   errln = Lvarget();
   poolvars += 1;
   err.s = 0;                   /* Absolute error of Argument of    */
   NEXT_(&err,&err);            /* Exponential Function :           */
   COPY_(&Lone,&dummy);   /*  yi*ln(|xi|)*[err(ln)+eps(Maxl)] */
   dummy.e = 1 - Lcurrprec - Guarddig;
   ADD_(&err,&dummy,&err);
   NEXT_(&err,&err);
   NEXT_(y,&dummy);
   MUL_(&err,&dummy,errln,&dummy);
   errln->s = 0;

   Maxl = Lcurrprec + 1;
   rc  = expve(y);              /* Function value and error bound   */

   if (rc != 0)
      {
      if (rc != UFLOW) {
         ERREXIT(0,RESUL,poolvars); /* Error message and handling */
         }
      else {
         ri->z = 1;             /* Set result 0                     */
         ri->r = 1;             /* Set rounding bit to 1 ulp        */
         ERREXIT(0,UFLOW,poolvars); /* Error message and handling   */
         }
      }

   Maxl = Minlerr;
   ADD_(&err,errln,&err);
   NEXT_(&err,&err);
   Lrnd = rounded;

RESULT:
   ASSIGN_(ri);                         /* Assign result and ulp's  */
   ri->s = rsgn;
   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpower 393) LhF    = ");  Lprinti(&LhF);
      printf("\n  (Lpower 394) LhE    = ");  Lprinti(&LhE);
      printf("\n  (Lpower 395) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





