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

/* CVS $Id: b_snhv.c,v 1.22 2014/01/30 17:24:05 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : sinhvea            Processor : C                   *
 *                                                                       *
 * Hyperbolic Sine Function for Multiple Precision Arithmetic for |x|<0.5*
 * ======================================================================*
 *                                                                       *
 * Compute the function value of the hyperbolic sine function using its  *
 * Taylor series.                                                        *
 * This routine should be called only in case of an argument not greater *
 * than 0.5 .                                                            *
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
 *    sqrtve          Intern Square Root Function                        *
 *    Lginit          Initialization of global constants                 *
 *    r_succ          Increasing a double real number by 1 ulp           *
 *    r_pred          Decreasing a double real number by 1 ulp           *
 *                    Intern arithmetic                                  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    Lone              ( = 1 )                            (Constant)    *
 *    Lmindbl           ( = 0x0010000000000000 )           (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 * Dependencies on other INTERN Functions :                              *
 *    sqrtve       Computed function value and error bound are returned  *
 *                 by the global variables ? (LhF) and sqrterr (LhE).    *
 *                                                                       *
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

/****************************************************************
 * Redefinition of Names of Intermediate Variables Used         *
 ****************************************************************/
        
#define res       LhF   /* Result variable                          */
#define err       LhE   /* Error of result                          */
#define sqrtval   LhF   /* Function value    ( == sqrtve !!! ) LhF  */
#define sqrterr   LhE   /* Rounding error    ( == sqrtve !!! ) LhE  */
#define dummy     LhD   /* Dummy variable (not further used)        */
        
        
/****************************************************************
 * Constants for Computation of sinh(x)                         *
 ****************************************************************/
        
#define Lguard    2

#define c2        t     /* Value of cosh(xi)                        */
#define Lpd       errsh /* Intern Variable replacing pd (Underflow) */
#define Lrhs      sh    /* Right hand side for degree computation   */
        
static a_btyp  mbds[3]   = { 0x18000000L,0x00000000L,0x00000000L };
static a_btyp  mefaca[3] = { 0x00000002L,0xBC28F5C3L,0x00000000L };
static a_btyp  mefacf[3] = { 0x00000001L,0x770A3D71L,0x00000000L };
static a_btyp  mepsd[3]  = { 0x00000002L,0xFC1A6314L,0x00000000L };
static a_btyp  m3[3]     = { 0x00000003L,0x00418938L,0x00000000L };
static a_btyp  mxfac[3]  = { 0xAC84B5DCL,0xC63F1413L,0x00000000L };
        
static dynamic  bdapp   = { 0, 0, 0, 0, -1, Minl, &mbds[0] };
static dynamic  errfacf = { 0, 0, 0, 0,  0, Minlerr, &mefacf[0] };
static dynamic  erra    = { 0, 0, 0, 0,  0, Minlerr, &mefaca[0] };
static dynamic  Lepsd   = { 0, 0, 0, 0,  0, Minlerr, &mepsd[0] };
static dynamic  l3      = { 0, 0, 0, 0,  0, Minl, &m3[0] };
static dynamic  xfac    = { 0, 0, 0, 0, -1, Minl, &mxfac[0] };

/* Cordes : unused
static dynamic  bdtsqr  = { 0, 0, 0, 0, 15, Minl, &mbds[1] };
*/
        
/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_snhv(dynamic *xi)
#else
int b_snhv(xi)

dynamic *xi;
#endif
#define  LRoutine    "sinhvea"
{
   extern dynamic   Lone, Lmindbl, LMinreal, LFminsq, sqrtval,
                    sqrterr, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  EUFfac, Lgiflag;
   extern a_real      FUFfac, FIUFfac, *r_one_, *r_mone;
   extern rounding  Lrnd;

#ifdef Debug
   extern int       Ldebug;
#endif
        
   dynamic   *errsh, *sh, *t2, *t;
   int       i, j, n, nscaled, rc, poolvars;
   a_btyp    Currprec, nsqr, u, L;
   a_real    xd, pd, epsd;
        
  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/
        
   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level          */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

   Currprec = Maxl;             /* Save precision setting        */
        
        
        
  /**************************************************************
   * Computation of the Function Value                          *
   **************************************************************/


  /*----------------*
   | Initialization |
   *----------------*/
        
   t = Lvarget();               /* Get new pool variable         */
   poolvars = 1;                /* Number of used pool variables */
   Maxl = xi->l + 1;            /* Set appropriate accuracy      */
   rc = 0;
   if ((COPY_(xi,t)) != 0) {    /* t = |xi|                      */
      ERREXIT(rc,rc,poolvars);  /* Error message and handling    */
      }
   t->s = 0;
        
   Maxl = Currprec;             /* Reset Maxl for expve          */
        
        
  /*----------------------------------------*
   | Computation using Taylor approximation |
   *----------------------------------------*/
        
  /* Argument reduction is performed by halving the argument until the
     absolute value  of the result is not greater 0.09375.
     The formula used for recursive argument reduction and computation
     of the final result is
          sinh x = 2*sinh(x/2)*cosh(x/2)
                 = 2*sinh(x/2)*sqrt(1 + sinh(x/2)*sinh(x/2)) .
     The function value for the original argument is computed by using
     this formula according to the number of halvings.

     Increasing the current mantissa length by 1, halving is performed
     without any rounding error except possible discarding of mantissa
     digits and/or bits exceeding the specified accuracy.       */
        
   nsqr=0;
   if GT_(t,&bdapp) {           /* Reduction if |xi| > 0.09375   */
      for (u=t->m[0]; (u & MSB) == 0; u <<= 1) nsqr++;
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 207) nsqr = %d",nsqr);
         printf("\n  (sinhvea 208) u    = %d",u);
         }
#endif
      if (nsqr < B_LENGTH-1) {
         u <<= 1;
         if ((u & MSB) != 0) nsqr--;
         }
       else
         if ((t->m[1] & MSB) != 0) nsqr--;
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 219) nsqr = %d",nsqr);
         printf("\n  (sinhvea 220) u    = %d",u);
         }
#endif
      nsqr = max(0,3-nsqr);     /* t->e = 0 since bdapp<=t<=0.5  */

/* nsqr += Currprec / 2; */      nsqr += Currprec / 8;
      if ((SHIFT_(-nsqr,t,t)) != 0) {
         ERREXIT(rc,rc,poolvars); /* Error message and handling    */
         }
      }
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 233) Number of halvings = %d",nsqr);
         printf("\n  (sinhvea 234) Reduced Argument = "); Lprinti(t);
         }
#endif


  /*-----------------------------------------------------------*
   | Computation of the Degree of the Approximation Polynomial |
   *-----------------------------------------------------------*/
        
   /* The degree N of the polynomial is computed using standard (double)
      real arithmetic (IEEE) with adding an rounding error of 1 digit of
      the last place for each operation according to the formula
           T^(N+1) <= 0.5*B^(1-L)*(N+1)!  ,
      where L=Maxl+Lguard is the used number of mantissa digits with Maxl
      being the number of mantissa digits when the routine is entered.
      If the inequality above holds, the approximation error EAPP is bound
      by 0.5*B^(1-L).
      The test is started with N=0, i.e. the inequality
           T^2<=3*B^(1-L)/1.0051   .
      The loop counter n is containing the value N+1.   */
        
   L = Currprec + Lguard;       /* Number of digits used         */
   u = max(nsqr+19,0) / 20;     /* Number of additional digits   */
   nscaled = 0;
   Lepsd.e = 1 - L - u;         /* epsd = 2.984777*b^(1-L-u) */
   Maxl = L+u;
   t2 = Lvarget();              /* Value of t*t,  (== c2)        */
   poolvars += 1;               /* Number of used pool variables */
   rc = 0;
   MUL_(t,t,t2,&dummy);         /* t2 = t^2 (low precision)  */
   if (rc != 0) {
      ERREXIT(rc,rc,poolvars);  /* Error message and handling    */
      }
   Maxl = Minl;
        
   if (GE_(&Lepsd,&LMinreal) & GE_(t,&LFminsq)) {
      while GT_(&Lmindbl,&Lepsd) {
         Lepsd.e += EUFfac;
         nscaled++;
         }

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 282) Number of scaling operations : %d",
                 nscaled);
         }
#endif

      rc = b_bcid(&Lepsd,&epsd,(a_intg)0); /* Conversion to double */
      if (rc==ROUND) rc = 0;
      NEXT_(t2,&dummy);
      i = b_bcid(&dummy,&xd,(a_intg)0); /* Conversion to double */
      if (i==ROUND) {       /*   upper bound by adding 1 ulp */
/* Cordes xd = Faddeps(xd);     */   /*   if conversion is not exact  */
            R_ASSIGN(xd,r_succ(xd)); /*   if conversion is not exact  */
            i = 0;
            }
      rc += i;

      if (rc != 0)
/*Cordes pd = -1.0;         */  /* Errors occurred during conv   */
         R_ASSIGN(pd,*r_mone);  /* Errors occurred during conv   */
      else
/*Cordes pd = xd;           */  /* pd=xd  */
         R_ASSIGN(pd,xd);       /* pd=xd  */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 303) Returncode conversion Lepsd->epsd : %d",
                 rc);
         }
#endif

      }
    else
/*Cordes pd = -1.0;   */
      R_ASSIGN(pd,*r_mone);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sinhvea 313) epsd = 2^(%d) = ",B_LENGTH*Lepsd.e);
      Lprinti(&Lepsd);
      printf(" = %g",epsd);
      printf("\n  (sinhvea 316) xd = %g,    pd = %g",xd,pd);
      }
#endif

   Lpd  = Lvarget();                       /*    (== errsh)     */
   Lrhs = Lvarget();                       /*    (== sh)        */
   poolvars += 2;

/*Cordes if (pd>0) {  */
   if (r_sign(pd)>0) {
/*Cordes for (n=1; (nscaled > 0) | (pd >= epsd); n++) {    */
      for (n=1; (nscaled > 0) | (r_ge(pd,epsd)); n++) {
/*Cordes pd = Faddeps(pd*xd); */        /* pd   = ^{pd*xd*xd} */
         R_ASSIGN(pd,r_mulu(pd,xd));    /* pd   = ^{pd*xd*xd} */
         i = 2*n*(2*n+1);
/*Cordes epsd = Fsubeps(epsd*i);  */    /* epsd = v{epsd*2n*(2n+1)}      */
         R_ASSIGN(epsd,r_muld(epsd,r_flot((a_intg)i)));
                                        /* epsd = v{epsd*2n*(2n+1)} */
/*Cordes if ((nscaled > 0) & (epsd > 1.0)) {    */
         if ((nscaled > 0) & (r_gt(epsd,*r_one_))) {
/* Cordes   epsd *= FIUFfac;   */
            R_ASSIGN(epsd,r_mulu(epsd,FIUFfac));
            nscaled--;
            }
/* Cordes if ((nscaled > 0) & (pd < 1.0)) {   */
         if ((nscaled > 0) & (r_lt(pd,*r_one_))) {
/* Cordes   pd *= FUFfac;      */
            R_ASSIGN(pd,r_mulu(pd,FUFfac));
            nscaled--;
            }
#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sinhvea 339) %g %g",epsd,pd);
            }
#endif
         }
      n--;
      }
    else {
      Maxl = Minlerr;           /* Set appropriate accuracy      */
      COPY_(t2,Lpd);
      COPY_(&Lepsd,Lrhs);
      for (n=1; GE_(Lpd,Lrhs); n++) {
         MUL_(Lpd,t2,Lpd,&dummy); /* pd   = ^{pd*xd*xd} */
         if (dummy.z==0) NEXT_(Lpd,Lpd);
         i = 2*n*(2*n+1);
         MULINT_(Lrhs,i,Lrhs);  /* epsd = v{epsd*2n*(2n+1)}      */
         }
      }
        
        
  /*-----------------------------------*
   | Computation of the Rounding Error |
   *-----------------------------------*/
        
   /* The rounding error ERR is computed using (short, i.e. 2 digits)
      INTERN numbers.                                            */
        
   Maxl = Minlerr;              /* Set appropriate accuracy      */
        
   MUL_(t2,&xfac,errsh,&dummy); /* errsh = {3.001 +              */
   if (errsh->r != 0) NEXT_(errsh,errsh);  /*  0.6739*x^2}*b^(1-L)  */
   ADD_(&l3,errsh,errsh);
   if (errsh->r != 0) NEXT_(errsh,errsh);
   NEXT_(errsh,errsh);          /* quadratic terms compensation  */
   errsh->e += (1 - L - u);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sinhvea 381) et = ");
      Lprinti(errsh);
      }
#endif

   if (rc != 0) {
      ERREXIT(EPERR,387  ,poolvars);    /* Error message and handling    */
      }

#ifdef Debug
   if (Ldebug & pd == 0.0) {
      printf("\n  (expve 371) Underflow of pd");
      ERREXIT(DUFLW,393  ,poolvars);    /* Error message and handling    */
      }

   if ((Ldebug >= 0) & (n == 99)) {
      printf("\n  (sinhvea 397) Maximal degree reached");
      printf("\n  (sinhvea 398) Approximation polynomial degree %d",n);
      n = min(n,100);
      }
#endif

        
  /*----------------------------------------*
   | Evaluation of Approximation Polynomial |
   *----------------------------------------*/
        
  /* Evaluation of the Taylor series of the exponential function for the
     reduced argument T using Horner's scheme.
     The coefficients N! are generated during evaluation by using N*T
     instead of T for the general Horner's scheme.            */
        
   Maxl = L + u;                /* Accuracy polynomial evaluation*/

   rc = 0;
   MUL_(t,t,t2,&dummy);         /*  t2 = t*t = t^2               */
   j = 2*n*(2*n+1);
   DIVINT_(t2,j,sh);            /*  sh = t^2/(2n*(2n+1))         */
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sinhvea 427) sh = ");  Lprinti(sh);
      printf("\n  (sinhvea 428) t  = ");  Lprinti(t);
      }
#endif
        
   for (i=n-1; i>0; i--) {
       ADD_(&Lone,sh,sh);       /*     sh=(sh+1)*t^2/(2i*(2i+1)) */
       j = 2*i*(2*i+1);
       MUL_(sh,t2,sh,&dummy);
       DIVINT_(sh,j,sh);

#ifdef Debug
       if (Ldebug >= 0) {
          printf("\n  (sinhvea 440) sh = ");  Lprinti(sh);
          }
#endif
       }
   ADD_(&Lone,sh,sh);           /*     sh = (sh+1)*t             */
   MUL_(sh,t,sh,&dummy);
   if (rc != 0) {
      ERREXIT(PEVAL,447  ,poolvars); /* Error message and handling    */
      }
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sinhvea 452) sh = ");  Lprinti(sh);
      printf("\n  (sinhvea 453) t  = ");  Lprinti(t);
      }
#endif

        
  /*---------------------------------*
   | Computation of the Final Result |
   *---------------------------------*/
        
  /* The final result is computed by squaring the value for the reduced
     argument NSQR times, where NSQR is the number of halvings during
     computation of the reduced argument.                        */
        
   rc = 0;
   erra.e = 1 - L;              /* Set exponent of 2.735*b^(1-l) */
   for (i=1; i<=nsqr; i++) {
      Maxl = L;                 /* Accuracy for computation      */
      MUL_(sh,sh,c2,&dummy);    /* sinh(2x) = 2*sinh(x)*sqrt(1+  */
      ADD_(c2,&Lone,c2);        /*              sinh(x)*sinh(x)) */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 479) c2    = ");  Lprinti(c2);
         }
#endif

      SQRT_(c2);                /* val: sqrtval, error: sqrterr  */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 487) Maxl  = %d",Maxl);
         printf("\n  (sinhvea 488) sqt   = ");  Lprinti(&sqrtval);
         }
#endif

      MUL_(sh,&sqrtval,sh,&dummy);
      SHIFT_(1,sh,sh);

      Maxl = Minlerr;           /* Accuracy for error estimates  */
      MUL_(errsh,&errfacf,errsh,&dummy);   /* errsh = 2.735*errsh  */
      ADD_(errsh,&erra,errsh);  /*       + 1.465*b^(1-l) */
      sqrterr.s = 0;
      ADD_(errsh,&sqrterr,errsh); /*       + |sqrterr|             */
      NEXT_(errsh,errsh);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sinhvea 504) |err| = ");  Lprinti(errsh);
         printf("\n  (sinhvea 505) sh    = ");  Lprinti(sh);
         }
#endif
      }
   sh->s = xi->s;               /* Set sign of result            */
   erra.e = -1;                 /* Reset exponent of 5/2         */
   Maxl = L;                    /* Accuracy for computation      */
   COPY_(sh,&res);              /* Resturn result                */

   Maxl = Minlerr;
   MUL_(errsh,sh,&err,&dummy);
   if (err.r != 0) NEXT_(&err,&err);
   Lrnd = ibound;

   if (rc != 0) ERREXIT(EPERR,519  ,poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sinhvea 523) |err| = ");  Lprinti(&err);
      printf("\n  (sinhvea 524) res   = ");  Lprinti(&res);
      }
#endif

   Lvardrop(poolvars);          /* Return used variables to pool */
   RETURN(rc);
}





