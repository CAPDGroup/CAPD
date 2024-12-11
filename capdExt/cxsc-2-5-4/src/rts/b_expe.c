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

/* CVS $Id: b_expe.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : expve              Processor : C                   *
 *                                                                       *
 * Exponential Function for Multiple Precision Arithmetic                *
 * ======================================================                *
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
 *    Lginit          Initialization of global constants                 *
 *    r_succ          Increasing a double real number by 1 ulp           *
 *    r_pred          Decreasing a double real number by 1 ulp           *
 *                    Intern arithmetic                                  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    Lone              ( = 1 )                            (Constant)    *
 *    Leps              ( = current rounding error )       (Constant)    *
 *    Lmindbl           ( = 0x0010000000000000 )           (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    LFunits         Distance to Upper Bound in ulps                    *
 *    Lintern         Flag for call by other standard functions          *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
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

#define p         LhF   /* Function value                           */
#define err       LhE   /* Rounding error                           */
#define dummy     LhD   /* Dummy variable (not further used)        */

/********************************************************************
 * Constants for Computation of exp(x)                              *
 ********************************************************************/

/*
#undef  UFLOW
#define UFLOW     0
*/

#define Lguard    2

#define errsqr    Lpd

static a_btyp  mbds[3]   = { 0x18000000L,0x00000000L,0x00000000L};
static a_btyp  m3o2[3]   = { 0x00000001L,0x80000000L,0x00000000L};
static a_btyp  m17o16[3] = { 0x00000001L,0x10000000L,0x00000000L};

static dynamic  bdapp   = { 0, 0, 0, 0, -1, Minl, &mbds[0] };
static dynamic  l3o2    = { 0, 0, 0, 0,  0, Minlerr, &m3o2[0] };
static dynamic  l17o16  = { 0, 0, 0, 0,  0, Minlerr, &m17o16[0] };

/* Cordes : unused
static dynamic  bdtsqr  = { 0, 0, 0, 0, 15, Minl, &mbds[1] };
*/

/********************************************************************
 * Algorithm                                                        *
 ********************************************************************/

#ifdef LINT_ARGS
int b_expe(dynamic *Larg)
#else
int b_expe(Larg)
dynamic *Larg;
#endif
#define  LRoutine    "expve"
{
   extern dynamic   Lone, Leps, Lmindbl, LMinreal, LFminsq, p, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  EUFfac, Lgiflag;
   extern a_real      FUFfac, FIUFfac, *r_one_, *r_mone;
   extern rounding  Lrnd;

#ifdef Debug
   extern int       Ldebug;
#endif

   dynamic   *t, *Lpd, *Lepsd;
   int       i, n, nscaled, rc, poolvars;
   a_intg    nsqr;
   a_btyp    Currprec, u, L;
   a_real    xd, pd, epsd;

  /******************************************************************
   * Initialization of Global Constants and Storage Locations       *
   ******************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level          */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

   Currprec = Maxl;             /* Save precision setting        */

  /***************************************************************
   * Reduction of Argument                                       *
   ***************************************************************/

  /* Argument reduction is performed by halving the argument until the
     absolute value  of the result is not greater 0.09375.
     The formula used for recursive argument reduction and computation
     of the final result is
          e^x = (e^(x/2))^2 .
     The function value for the original argument is computed by squaring
     the value for the reduced argument according to the number of
     halvings.
     Increasing the current mantissa length by 1, halving is performed
     without any rounding error except possible discarding of mantissa
     digits and/or bits exceeding the specified accuracy.       */

   t = Lvarget();
   poolvars = 1;
   Maxl = Larg->l + 1;          /* Set appropriate accuracy      */
   rc = 0;
   if ((COPY_(Larg,t)) != 0)    /* t = |Larg|                    */
      ERREXIT(rc,rc,poolvars);  /* Error message and handling    */
   t->s = 0;
   nsqr=0;
   if GT_(t,&bdapp) {           /* Reduction if |Larg| > 0.09375 */
      for (u=t->m[0]; (u & MSB) == 0; u <<= 1) nsqr++;
      if (nsqr < B_LENGTH-1) {
         u <<= 1;
         if ((u & MSB) != 0) nsqr--;
         }
       else
         if ((t->m[1] & MSB) != 0) nsqr--;
      nsqr = (t->e+1)*B_LENGTH - nsqr + 3;
      if (nsqr>=39) {
         if (Larg->s == 0)                 /* Result overflow            */
            ERREXIT(OFLOW,OFLOW,poolvars)  /* Error message and handling */
          else {                           /* Result underflow :         */
            p.z = 1;                       /* ri = exp(-oo) = 0          */
            p.r = 1;                       /* Set rounding bit           */
            p.s = 0;                       /* Set positive sign          */
            ERREXIT(0,UFLOW,poolvars);     /* Error message and handling */
            }
         }
/* nsqr += Currprec / 2; */      nsqr += Currprec / 4;
      if ((SHIFT_(-nsqr,t,t)) != 0)
         ERREXIT(rc,rc,poolvars);          /* Error message and handling */
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 194) Number of halvings = %d",nsqr);
      /* printf("\n  (expve 195) Reduced Argument = "); Lprinti(t); */
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
      The test is started with N=1, i.e. the inequality T^2<=B^(1-L).
      The loop counter n is containing the value N+1.           */

   Maxl = Minl;
   L = Currprec + Lguard;       /* Number of digits used        */
   u = max(nsqr+9,0) / 32;      /* Additional guard digits      */
   nscaled = 0;
   Leps.e = 1 - L - u;          /* eps = b^(1-L-u)   */

   if (GE_(&Leps,&LMinreal) & GE_(t,&LFminsq)) {
      while GT_(&Lmindbl,&Leps) {
         Leps.e += EUFfac;
         nscaled++;
         }

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (expve 234) Number of scaling operations : %d",
                nscaled);
         }
#endif

      rc = b_bcid(&Leps,&epsd,(a_intg)0); /* Conversion to double */
      if (rc==ROUND) rc = 0;
      i = b_bcid(t,&xd,(a_intg)0);  /* Conversion to double         */
      if (i==ROUND) {               /* upper bound by adding 1 ulp  */
/*Cordes xd = Faddeps(xd);      */  /*   if conversion is not exact */
         R_ASSIGN(xd,r_succ(xd));   /*   if conversion is not exact */
         i = 0;
         }
      else rc += i;

      if (rc != 0)
/*Cordes pd = -1.0;             */  /* Errors occurred during conv */
         R_ASSIGN(pd,*r_mone);      /* Errors occurred during conv */
      else if (nscaled > 0) {
                                    /* pd=xd*xd with rounding error */
/* Cordes pd = Faddeps((xd*FUFfac)*xd); */
         R_ASSIGN(pd,r_mulu(r_mulu(xd,FUFfac),xd));
                                    /* pd=xd*xd with rounding error */
         Leps.e = 1 - L - u;        /* eps = b^(1-L-u) (restore) */
         nscaled--;
         }
       else
/*Cordes pd = Faddeps(xd*xd);   */  /* pd=xd*xd with rounding error */
          R_ASSIGN(pd,r_mulu(xd,xd));
                                    /* pd=xd*xd with rounding error */
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (expve 259) Returncode conversion Leps->epsd : %d",
                rc);
         }
#endif

      }
    else
/* Cordes      pd = -1.0;       */
      R_ASSIGN(pd,*r_mone);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 269) epsd = 2^(%d) = ",B_LENGTH*Leps.e);
      Lprinti(&Leps);
      printf(" = %g",epsd);
      printf("\n  (expve 271) xd = %g,    pd = %g",xd,pd);
      }
#endif

   Lpd = Lvarget();
   poolvars += 1;

/* Cordes   if (pd>0) { */
   if (r_sign(pd)>0) {
/*Cordes      for (n=3; (nscaled > 0) | (pd >= epsd); n++) {    */
      for (n=3; (nscaled > 0) | (r_ge(pd,epsd)); n++) {
/*Cordes pd   = Faddeps(pd*xd);         */ /* pd   = ^(pd*xd)   */
         R_ASSIGN(pd,r_mulu(pd,xd));         /* pd   = ^(pd*xd)   */
/*Cordes epsd = Fsubeps(epsd*n);        */ /* epsd = v(epsd*n)  */
         R_ASSIGN(epsd,r_muld(epsd,r_flot((a_intg)n)));
/*Cordes         if ((nscaled > 0) & (epsd > 1.0)) {    */
         if ((nscaled > 0) & (r_gt(epsd,*r_one_))) {
/* Cordes   epsd *= FIUFfac;        */
            R_ASSIGN(epsd,r_mulu(epsd,FIUFfac));
            nscaled--;
            }
/* Cordes         if ((nscaled > 0) & (pd < 1.0)) {     */
         if ((nscaled > 0) & (r_lt(pd,*r_one_))) {
/* Cordes   pd *= FUFfac;       */
            R_ASSIGN(pd,r_mulu(pd,FUFfac));
            nscaled--;
            }

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (expve 293) %g %g",epsd,pd);
            }
#endif

         }
      n--;
      }
    else {
      Lepsd = Lvarget();
      COPY_(t,Lpd);
      COPY_(&Leps,Lepsd);
      Maxl = Minlerr;              /* Set appropriate accuracy  */
      for (n=2; GE_(Lpd,Lepsd); n++) {
         MUL_(Lpd,t,Lpd,&dummy);   /* pd   = ^(pd*xd)   */
         if (dummy.z==0) NEXT_(Lpd,Lpd);
         MULINT_(Lepsd,n,Lepsd);   /* epsd = v(epsd*n)  */
         }
      Lvardrop(1);
      }


  /*-----------------------------------*
   | Computation of the Rounding Error |
   *-----------------------------------*/

   /* The rounding error ERR is computed using
      (short, i.e. 'Minlerr' digits)
      INTERN numbers.                         */

   Maxl = Minlerr;      /* Set appropriate accuracy      */

   MULINT_(t,3,&err);   /* err = 3*t  (polynomial eval)  */
   if (err.r != 0) NEXT_(&err,&err);

   ADD_(&l3o2,&err,&err);       /* et  = (1.5 + 3*t)*B^(1-L) */
   if (err.r != 0) NEXT_(&err,&err);
   NEXT_(&err,&err);            /* quadratic terms compensation  */
   err.e += (1 - L - u);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 337) et = ");  Lprinti(&err);
      }
#endif

   SHIFT_(nsqr,&err,&err);      /* err = 1.0625*2^nsqr*err   */
   NEXT_(&err,&err);
   MUL_(&err,&l17o16,&err,&dummy);      /*     = 2^nsqr*err*17/16    */
   if (err.r != 0) NEXT_(&err,&err);

   if (nsqr > 0) {
      SHIFT_(nsqr,&Lone,errsqr);        /* errsqr = (2^nsqr - 1)*eps */
      NEXT_(errsqr,errsqr);
      SUB_(errsqr,&Lone,errsqr);
      NEXT_(errsqr,errsqr);
      errsqr->e = 1 - L - u;            /* Number of digits result adapt */
      ADD_(&err,errsqr,&err);   /* err = 1.0625*2^nsqr*(eps+err) */
      if (err.r != 0) NEXT_(&err,&err); /*   + 2^nsqr*eps  */
      }

   NEXT_(&err,&err);            /* quadratic terms compensation  */
   Lvardrop(1);
   poolvars -= 1;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 362) err = ");  Lprinti(&err);
      }
#endif

   if (rc != 0)
      ERREXIT(EPERR,EPERR,poolvars);    /* Error message and handling   */

#ifdef Debug
   if (Ldebug >= 0 && pd == 0.0) {
      printf("\n  (expve 371) Underflow of pd");
      ERREXIT(DUFLW,372,poolvars);      /* Error message and handling   */
      }

   if (Ldebug >= 0 && n == 99) {
      printf("\n  (expve 376) Maximal degree reached");
      n = min(n,100);
      }
   if (Ldebug >= 0) {
      printf("\n  (expve 380) Approximation polynomial degree %d",n);
      }
#endif


  /*----------------------------------------*
   | Evaluation of Approximation Polynomial |
   *----------------------------------------*/

  /* Evaluation of the Taylor series of the exponential function for the
     reduced argument T using Horner's scheme.
     The coefficients N! are generated during evaluation by using N*T
     instead of T for the general Horner's scheme.              */

   Maxl = L + u;                /* Accuracy polynomial evaluation*/

   rc = 0;
   COPY_(t,&p);                 /*  p = t/n                      */
   DIVINT_(&p,n,&p);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 407) p = ");  Lprinti(&p);
      printf("\n  (expve 408) t = ");  Lprinti(t);
      }
#endif

   for (i=n-1; i>0; i--) {
       ADD_(&Lone,&p,&p);       /*     p = (p+1)*|Larg|/i        */
       MUL_(&p,t,&p,&dummy);
       DIVINT_(&p,i,&p);

#ifdef Debug
       if (Ldebug >= 0) {
          printf("\n  (expve 419) p = ");  Lprinti(&p);
          }
#endif

       }
   ADD_(&Lone,&p,&p);           /*     p = p+1                   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 428) p = ");  Lprinti(&p);
      printf("\n  (expve 429) t = ");  Lprinti(t);
      }
#endif

   if (rc != 0)
      ERREXIT(PEVAL,PEVAL,poolvars); /* Error message and handling */


  /*---------------------------------*
   | Computation of the Final Result |
   *---------------------------------*/

  /* The final result is computed by squaring the value for the reduced
     argument NSQR times, where NSQR is the number of halvings during
     computation of the reduced argument.               */

   Maxl = L + u;                /* Accuracy result adaptation    */

   rc = 0;
   for (i=1; i<=nsqr; i++)
      if ((MUL_(&p,&p,&p,&dummy)) != 0) break;

   if (rc != 0)
      {
      if (rc == OFLOW) {

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (expve 460) Overflow during result adaptation\n");
            }
#endif

         if (Larg->s == 0)                 /* Result overflow   */
            i  = OFLOW;                    /* Error message     */
         else {                            /* Result underflow  */
            p.z = 0;                       /* Set result 0      */
            p.r = 1;                       /* Set rounding bit  */
            i  = 0;                        /* Error message     */
            rc = UFLOW;                    /* Return code       */
            }
         ERREXIT(i,rc,poolvars); /* Error message and handling    */
         }
      else
         ERREXIT(RESUL,RESUL,poolvars); /* Error message and handling */
      }

   if (Larg->s == 1) {
      DIV_(&Lone,&p,&p);
      NEXT_(&p,&p);             /* Get upper bound for division */
      Maxl = Minlerr;
      ADD_(&err,&Leps,&err);    /* add division error           */
      if (err.r != 0) NEXT_(&err,&err);
      NEXT_(&err,&err);         /* quadratic terms compensation */
      Lrnd = ubound;            /* Value is upper bound         */
      }
    else
      Lrnd = lbound;            /* Value is lower bound         */

   if (rc != 0)  ERREXIT(EPERR,489,poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (expve 493) res = ");  Lprinti(&p);
      printf("\n  (expve 494) err = ");  Lprinti(&err);
      }
#endif

   Lvardrop(poolvars);
   RETURN(rc);
}





