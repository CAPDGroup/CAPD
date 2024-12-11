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

/* CVS $Id: b_lnva.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : lnvea              Processor : C                   *
 *                                                                       *
 * Internal Function ln((1+x)/(1-x)) for Multiple Precision Arithmetic   *
 * ===================================================================   *
 *                                                                       *
 * Compute the function value of log((1-x)/(1+x)) using its Taylor series*
 * The argument supplied is the fraction (1-x)/(1+x).                    *
 * This routine should be called only in case of |x-1| < 3/32.           *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Arguments :                                                           *
 *      dynamic  *Larg          Argument                Parameter        *
 *      a_btyp   Maxl           Accuracy                global           *
 *                                                                       *
 * Results   :                                                           *
 *      dynamic  Lres           Result                  global ( == LhF )*
 *      dynamic  err            Rounding Error          global ( == LhE )*
 *                                                                       *
 * Argument Range :    All positive numbers                              *
 * Result Range :      Positive numbers                                  *
 *                                                                       *
 * Used Functions :                                                      *
 *    Lginit          Initialization of global constants                 *
 *    r_succ          Increasing a double real number by 1 ulp           *
 *    r_pred          Decreasing a double real number by 1 ulp           *
 *                    Intern arithmetic                                  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    Lone              ( = 1 )                            (Constant)    *
 *    Leps              ( = current rounding error )       (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 * Dependencies :                                                        *
 *    Number Base       assumed 2^32                                     *
 *                                                                       *
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

/*********************************************************************
 * Redefinition of Names of Intermediate Variables Used              *
 *********************************************************************/

#define Lres      LhF           /* Function value                    */
#define err       LhE           /* Rounding error                    */
#define dummy     LhD           /* Dummy variable (not further used) */

/*********************************************************************
 * Constants for Computation of ln(x)                                *
 *********************************************************************/

#define t         Larg
#define Lpd       p

static a_btyp  mepsd[3]  = { 0x000000CEL,0x38E38E38L,0xE38E38E3L };
static a_btyp  mt2fac[3] = { 0x2AC10000L,0x00000000L,0x00000000L };
static a_btyp  m5o2[3]   = { 0x00000002L,0x80000000L,0x00000000L };

static dynamic  Lepsd   = { 0, 0, 0, 0,  0, Minlerr, &mepsd[0] };
static dynamic  t2fac   = { 0, 0, 0, 0, -1, Minlerr, &mt2fac[0] };
static dynamic  l5o2    = { 0, 0, 0, 0,  0, Minlerr, &m5o2[0] };

/*********************************************************************
 * Algorithm                                                         *
 *********************************************************************/

#ifdef LINT_ARGS
int b_lnva(dynamic *t)
#else
int b_lnva(t)

dynamic *t;
#endif
#define  LRoutine    "lnvea"
{
   extern dynamic   Lres, LMinreal, Lmindbl, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lgiflag;
   extern a_btyp    EUFfac;
   extern a_real      FUFfac, FIUFfac, *r_mone, *r_one_;
   extern rounding  Lrnd;

#ifdef Debug
   extern int       Ldebug;
#endif

   dynamic   *t2, *p, *Lrhs;
   int       rc, poolvars, n, nscaled;
   a_btyp    Currprec;
   a_real    t2d, pd, epsd;

  /*******************************************************************
   * Initialization of Global Constants and Storage Locations        *
   *******************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level          */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

   Currprec = Maxl;             /* Save Maxl                     */


  /**************************************************************
   * Computation of the Degree of the Approximation Polynomial  *
   **************************************************************/

   /* The degree N of the polynomial is computed using standard (double)
      real arithmetic (IEEE) with adding an rounding error of 1 digit of
      the last place for each operation according to the formula
           t^(2N) <= 0.5*B^(1-L)*(2N+3)*3712/9  ,
      where Currprec=Maxl is the used number of mantissa digits with Maxl
      being the number of mantissa digits when the routine is entered.
      If the inequality above holds, the approximation error EAPP is bound
      by 0.5*B^(1-L).
      The test is started with N=0, i.e. the inequality
           t^2 <= 0.5*B^(1-L)*3*3712/9   .
      The loop counter n is containing the value N+1.           */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnvea 168) Currprec = %d",Currprec);
      printf("\n  (lnvea 169) Arg = ");   Lprinti(t);
      }
#endif

   t2 = Lvarget();              /* Get new pool variable         */
   poolvars = 1;                /* Number of used pool variables */

   nscaled = 0;
   Lepsd.e = 1 - Currprec;      /* eps = (1856/9)*b^(1-Currprec) */
   rc = 0;
   MUL_(t,t,t2,&dummy);         /* t2 = t^2 (low precision)      */
   if (rc != 0) ERREXIT(CONVD,CONVD,poolvars);
   n = 0;
   Maxl = Minl;

   if GE_(&Lepsd,&LMinreal) {
      while GT_(&Lmindbl,&Lepsd) {
         Lepsd.e += (int)EUFfac;
         nscaled++;
         }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnvea 196) Lepsd = ");   Lprinti(&Lepsd);
      printf("\n  (lnvea 197) Number of scaling operations : %d",nscaled);
      printf("\n  (lnvea 198) t2 = ");   Lprinti(t2);
      }
#endif

      n = b_bcid(&Lepsd,&epsd,(a_intg)0);  /* Conversion to double */
      if (n==ROUND) n = 0;
      if (n!=0) {
/*Cordes pd = -1.0;     */
            R_ASSIGN(pd,*r_mone);
            goto CONTINUE;
            }

      NEXT_(t2,&dummy);
      n = b_bcid(&dummy,&t2d,(a_intg)0); /* Conversion to double  */
      if (n==ROUND || n==UFLOW) {
                                       /* upper bound by adding 1 ulp */
/*Cordes t2d = Faddeps(t2d);    */     /*   if conversion is not exact*/
         R_ASSIGN(t2d,r_succ(t2d));    /*   if conversion is not exact*/
         n = 0;
         }
      rc += n;

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (lnvea 219) dummy = ");   Lprinti(&dummy);
         printf("\n  (lnvea 220) t2d = %g",t2d);
         }
#endif

      if (rc != 0) {
/*Cordes pd = -1.0;    */
         R_ASSIGN(pd,*r_mone);
         rc = 0;
         goto CONTINUE;
         }

/*Cordes pd = t2d;         */
      R_ASSIGN(pd,t2d);

#ifdef Debug
      if (Ldebug >= 0 && pd == 0.0) {
         printf("\n  (lnvea 234) Underflow of pd");
         }
#endif

/* Cordes      if (pd == 0.0) {
   Cordes         pd = -1.0;
   Cordes         }
   Cordes      }
   Cordes    else
   Cordes      pd = -1.0;  */

      if (r_sign(pd)==0) {
         R_ASSIGN(pd,*r_mone);
         }
      }
    else
      R_ASSIGN(pd,*r_mone);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnvea 247) Returncode conversion Lepsd->epsd : %d",rc);
      printf("\n  (lnvea 248) t2d = %g,    pd = %g",t2d,pd);
      printf("\n  (lnvea 249) epsd = 2^(%d) = ",UL*Lepsd.e);
      Lprinti(&Lepsd); printf(" = %g",epsd);
      }
#endif

CONTINUE:

   p = Lvarget();       /* Get new variable              */
   poolvars++;
/* Cordes if (pd>0) {         */
   if (r_sign(pd)>0) {
/* Cordes if ((nscaled > 0) & (epsd > 1.0)) {      */
      if ((nscaled > 0) & (r_gt(epsd,*r_one_))) {
/* Cordes epsd *= FIUFfac;      */
         R_ASSIGN(epsd,r_mulu(epsd,FIUFfac));
         nscaled--;
         }
/* Cordes for (n=1; (nscaled > 0) | (pd >= (2*(n++)+1)*epsd);) {   */
      for (n=1; (nscaled > 0) |
                (r_ge(pd,r_mulu(r_flot((a_intg)(2*(n++)+1)),epsd)));) {
/* Cordes pd   = Faddeps(pd*t2d);   */  /* pd   = ^{pd*t2d}  */
         R_ASSIGN(pd,r_mulu(pd,t2d));   /* pd   = ^{pd*t2d} */
/* Cordes if ((nscaled > 0) & (pd < 1.0)) {   */
         if ((nscaled > 0) & (r_lt(pd,*r_one_))) {
/* Cordes   pd *= FUFfac;     */
            R_ASSIGN(pd,r_mulu(pd,FUFfac));
            nscaled--;
            }

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (lnvea 272) epsd=%g,  pd=%g,  t2d=%g",epsd,pd,t2d);
            }
#endif

         }
      }
    else {
      Maxl = Minlerr;           /* Set appropriate accuracy     */
      Lrhs = Lvarget();         /* Get new variable Lrhs        */
      poolvars++;
      COPY_(t2,Lpd);
      COPY_(&Lepsd,Lrhs);
      for (n=1; (MULINT_(Lrhs,2*(n++)+1,&dummy),GE_(Lpd,&dummy));) {
         MUL_(Lpd,t2,Lpd,&dummy); /* pd   = ^{pd*t2d}           */
         if (dummy.z==0) NEXT_(Lpd,Lpd);
         }
      Lvardrop(1);                /* Drop variable Lhrs         */
      poolvars--;
      }
   n--;

#ifdef Debug
   if (Ldebug >= 0) printf("\n  (lnvea 294) Polynomial degree =%d",n);
#endif


  /*-----------------------------------*
   | Computation of the Rounding Error |
   *-----------------------------------*/

   /* The rounding error ERR is computed using (short, i.e. 'Minlerr' digits)
      INTERN numbers.                                            */

   Maxl = Minlerr;              /* Set appropriate accuracy  */

   MUL_(t2,&t2fac,&err,&dummy); /* err=(2.5+t2fac*t2)*b^(1-Maxl) */
   if (err.r != 0) NEXT_(&err,&err);
   ADD_(&err,&l5o2,&err);
   if (err.r != 0) NEXT_(&err,&err);
   NEXT_(&err,&err);            /* quadratic terms compensation  */
   err.e += (1 - Currprec);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnvea 321) et = ");  Lprinti(&err);
      }
#endif

   if (rc != 0) ERREXIT(EPERR,EPERR,poolvars);


  /*----------------------------------------*
   | Evaluation of Approximation Polynomial |
   *----------------------------------------*/

  /* Evaluation of the Taylor series of the exponential function for the
     reduced argument T using Horner's scheme.
     The coefficients N! are generated during evaluation by using N*T
     instead of T for the general Horner's scheme.              */

   Maxl = Currprec;             /* Accuracy polynomial evaluation */

   DIVINT_(t2,2*(n--)+1,p);     /*  p = t*t/(2n+1)               */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnvea 348) p = ");  Lprinti(p);
      printf("\n  (lnvea 349) t = ");  Lprinti(t);
      }
#endif

   for (; n>0; n--) {
       MUL_(p,t2,p,&dummy);     /* p = p*t^2 + t*t/(2n+1)    */
       DIVINT_(t2,2*n+1,&dummy);
       ADD_(p,&dummy,p);

#ifdef Debug
       if (Ldebug >= 0) {
          printf("\n  (lnvea 360) p = ");  Lprinti(p);
          }
#endif

       }

   Maxl++;
   SHIFT_(1,t,t2);              /*     t2 = 2*t                  */
   Maxl--;

   MUL_(p,t2,p,&dummy);         /*     p = 2*t + 2*t*p           */
   ADD_(p,t2,&Lres);
   Lrnd = rounded;

   if (rc != 0) ERREXIT(PEVAL,PEVAL,poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnvea 383) Lres = ");  Lprinti(&Lres);
      }
#endif

   Lvardrop(poolvars);
   RETURN(rc);
}





