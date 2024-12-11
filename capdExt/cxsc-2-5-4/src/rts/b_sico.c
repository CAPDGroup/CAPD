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

/* CVS $Id: b_sico.c,v 1.22 2014/01/30 17:24:05 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : sicovea              Processor : C                 *
 *                                                                       *
 * Sine and Cosine Function for Multiple Precision Arithmetic            *
 * ==========================================================            *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Argument Range :    Not restricted                                    *
 * Result Range :      [-1,1]                                            *
 *                                                                       *
 * Used Functions :                                                      *
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
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

/*******************************************************************
 * Subroutine for Decomposition of Dynamic to Entier and Fraction  *
 *******************************************************************/

#ifdef LINT_ARGS
static int GetFrac(multiprecision ix,multiprecision ir,a_btyp  *ep)
#else
static int GetFrac(ix,ir,ep)

dynamic *ix;
dynamic *ir;
a_btyp  *ep;
#endif
{
   dynamic   frac;
   a_btyp    bfp;
   int       rc;

   rc = 0;

   if (ix->e < 0) {     /* |xi| < 1 ==> fraction part = xi  */
      COPY_(ix,ir);
      *ep = 0;
      return(rc);
      }

   if (ix->l < ix->e+1) {       /* |xi| >= B ==> zero fraction part */
      rc = b_bini(ir);
      *ep = 0;
      return(rc);
      }

   bfp = ix->e+1;
   *ep = ix->m[ix->e] % 16;
   frac.l = ix->l - bfp;
   frac.m = &(ix->m[bfp]);

   for (bfp=0; (bfp<frac.l) & (frac.m[bfp]==0); bfp++) { /* empty */ }

   if (bfp < frac.l && frac.m[bfp] != 0) {
      frac.z = 0;
      frac.s = 0;
      frac.r = 0;
      frac.e = -1-bfp;
      frac.l -= bfp;
      frac.m = &(frac.m[bfp]);
      COPY_(&frac,ir);
      }
    else
      rc = b_bini(ir);

   return(rc);
}


/****************************************************************
 * Redefinition of Names of Intermediate Variables Used         *
 ****************************************************************/

#define p         LhF   /* Function value                           */
#define err       LhE   /* Rounding error                           */
#define dummy     LhD   /* Dummy variable (not further used)        */
#define sinvea    LhF   /* Function value of sin(arg)    ( Result ) */
#define cosvea    LhD   /* Function value of cos(arg)    ( Result ) */
                        /*    must be saved immediately by caller ! */


/****************************************************************
 * Constants for Computation of exp(x)                          *
 ****************************************************************/

#define M6ULP  0xFFFFFFFAL   /* Maximal last digit for sure entier part*/
#define MOD4   0x00000004L   /* Mask                                   */
#define MOD2   0x00000002L   /* Mask                                   */

#define Lguard    2     /* Number of guard digits for arithmetic    */
#define Rguard    4     /* Number of guard digits argument reduction*/
#define Cguard    20    /* Number of guard digits computation of Pi */

#define errsqr    Lpd
#define fac_el    fac

extern dynamic   LPiov4, L4ovPi;

extern a_btyp  LhI;


static a_btyp  mepsd[3]    = { 0xB4481CD8L,0x5689039CL,0x00000000L };
static a_btyp  m39o2[3]    = { 0x00000013L,0x80000000L,0x00000000L };
static a_btyp  mp675[3]    = { 0xACCCCCCDL,0x00000000L,0x00000000L };
static a_btyp  m2sqt2[3]   = { 0x00000002L,0xD413CCCFL,0xE7799212L };
static a_btyp  m2sqt2m1[3] = { 0x00000001L,0xD413CCCFL,0xE7799211L };

static dynamic  Lepsd    = { 0, 0, 0, 0,  0, Minlerr, &mepsd[0] };
static dynamic  l39o2    = { 0, 0, 0, 0,  0, Minlerr, &m39o2[0] };
static dynamic  lp675    = { 0, 0, 0, 0, -1, Minlerr, &mp675[0] };
static dynamic  L2sqt2   = { 0, 0, 0, 0,  0, Minlerr, &m2sqt2[0] };
static dynamic  L2sqt2m1 = { 0, 0, 0, 0,  0, Minlerr, &m2sqt2m1[0] };

static dynamic  PiOvFour = { 0, 0, 1, 0, 0, 0, NULL };
static dynamic  FourOvPi = { 0, 0, 1, 0, 0, 0, NULL };

/********************************************************************
 * Algorithm                                                        *
 ********************************************************************/

#ifdef LINT_ARGS
int b_sico(dynamic *Larg)
#else
int b_sico(Larg)

dynamic *Larg;
#endif
#define  LRoutine    "sicovea"
{
   extern dynamic   Lone, Lmindbl, LMinreal, LFminsq, p, err, dummy;

   extern a_btyp  Maxl;
   extern a_btyp  EUFfac, Lgiflag;
   extern a_real      FUFfac, FIUFfac, *r_one_, *r_mone;
   extern rounding  Lrnd;

#ifdef Debug
   extern int       Ldebug;
#endif

   dynamic   *sin, *cos, *t, *Lpd, *MSP, *LSP, *t2, *Lrhs, *dvnd, *fac;
   int       nhalve, nscaled, rc = 0, poolvars;
   a_intg    i, n, numdig;
   a_btyp    Currprec, j, k, u, L, ep, FOPl, *r;
   a_real    xd, pd, epsd;

  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level          */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
   if (Ldebug >= 0) {
      printf("\n  (sicovea 206) Maxl = %d",Maxl);
      }
#endif

   Currprec = Maxl;             /* Save precision setting        */

   u = 0;                       /* In case of using variable u      */


  /**************************************************************
   * Reduction of Argument                                      *
   **************************************************************/

  /* Argument reduction is performed by multiplying the argument by PI/4
     and decomposition of the result into entier part EP and the fraction
     part T.
     Since splitting off the entier part means loss of accuracy, the
     fraction part is checked to have the accuracy required.
     In case of unsufficient accuracy the argument reduction is continued
     with more digits of PI/4 (iterative process).

     Additionally, a second argument reduction may be performed by halving
     the current reduced argument N times. This reduction must not increase
     the rounding error, i.e. halving must be done with 1 additional digit.
     Since the result adaptation increases the rounding errors, sine and
     cosine of the twice reduced argument must be computed with higher
     accuracy.                                                        */

   FOPl = max(Larg->e,0) + Currprec + Rguard;
   if (FOPl > L4ovPi.l) {
      Maxl = FOPl + Cguard;
      rc += LpiGen();
      }

   PiOvFour.m = LPiov4.m;
   PiOvFour.e = LPiov4.e;
   PiOvFour.l = Currprec + Lguard;
   FourOvPi.m = L4ovPi.m;
   FourOvPi.e = L4ovPi.e;
   FourOvPi.l = FOPl;

#ifdef Debug
/*   if (Ldebug >= 0) {
      printf("\n  (sicovea 253) FourOvPi = "); Lprinti(&FourOvPi);
      printf("\n  (sicovea 254) L4ovPi   = "); Lprinti(&L4ovPi);
      } */
#endif

   t   = Lvarget();
   Lpd = Lvarget();
   MSP = Lvarget();
   LSP = Lvarget();
   poolvars = 4;
   rc = 0;
   ep = 0;                      /* Initialize entier part with 0    */

   if (GT_ABS_(Larg,&LPiov4)) {
      t->z = 1;
      t->e = 0;
      FOPl = FourOvPi.l;
      Maxl = max(max(Larg->e,0)+Currprec+Rguard,(FOPl+Larg->l)/2);

      for ( ;; ) {

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 276) Maxl = %d",Maxl);
            }
#endif

         MUL_(Larg,&FourOvPi,MSP,LSP);  /* MSP+LSP = |Arg*4/Pi|   */

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 284) FOPl = %d",FOPl);
            printf("\n  (sicovea 285) MSP = "); Lprinti(MSP);
            printf("\n  (sicovea 286) LSP = "); Lprinti(LSP);
            }
#endif

         MSP->s = 0;
         Maxl = MSP->l + max(t->e - MSP->e + 1,0);
         ADD_(t,MSP,t);                 /* t = trunc(|Arg*4/Pi|,l)      */
         if (t->e >= 0) {
            rc += GetFrac(t,t,&LhI);    /* Get entier and fraction part */
            ep = (ep + LhI)%16;
            }

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 300) t   = "); Lprinti(t);
            }
#endif


         Maxl = Currprec + Lguard;
         if (ep%2 == 1)
            SUB_(&Lone,t,Lpd);   /* 1-t in case odd entier part     */
          else
            COPY_(t,Lpd);        /*  t  otherwise                   */

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 313) Lpd = "); Lprinti(Lpd);
            }
#endif


         numdig = 1 + Currprec + Lguard + Larg->e - Lpd->e - FOPl;

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 322) numdig  = %d",numdig);
            }
#endif

         if (numdig<=0) {       /* Accuracy test                   */
            if (Lpd->l >= Maxl) {
               if (Lpd->e < -1 || Lpd->m[Maxl-1] < M6ULP)
                  break;        /* Required accuracy reached       */
               for ( n=(a_intg)Maxl-1, r=Lpd->m+n; --n>=0; )
                  if (*(--r)<MAX_BASETYPE) break;
               if (n>=0) break; /* Required accuracy reached       */
               }
            numdig = Lguard;
            }

         if (LSP->z==0) {
            LSP->s = 0;
            Maxl = LSP->l + max(t->e - LSP->e + 1,0);
            ADD_(t,LSP,t);
            }

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 343) Maxl = %d",Maxl);
            printf("\n  (sicovea 344) t   = "); Lprinti(t);
            }
#endif


         if (FOPl+numdig > L4ovPi.l) {
            Maxl = FOPl + numdig + Cguard;
            rc += LpiGen();
            }

         FourOvPi.m = &(L4ovPi.m[FOPl]);
         FourOvPi.e -= FourOvPi.l;
         FourOvPi.l = numdig;
         FOPl += numdig;

         for ( ; L4ovPi.m[FOPl] == 0; ) {
            FourOvPi.m = &(L4ovPi.m[++FOPl]);
            FourOvPi.l--;
            FourOvPi.e--;
            }
         if (FourOvPi.l <= 0) FourOvPi.z = 0;

         Maxl = max(numdig,(FOPl+Larg->l)/2);

#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 370) FOPl = %d,",FOPl);
            printf("    Larg->l = %d    Maxl = %d",Larg->l,Maxl);
            }
#endif

         }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 379) Reduced Argument = "); Lprinti(Lpd);
      }
#endif

      Maxl = Currprec + Lguard;
      MUL_(Lpd,&PiOvFour,t,&dummy);     /* t = <red arg>*pi/4   */
      }
    else {
      Maxl = Larg->l;
      COPY_(Larg,t);
      t->s = 0;
      }

   nhalve=0;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 396) Entier part = %d",ep);
      printf("\n  (sicovea 397) Series Argument  = "); Lprinti(t);
      printf("\n  (sicovea 398) Maxl = %d",Maxl);
      }
#endif


  /*-----------------------------------------------------------*
   | Computation of the Degree of the Approximation Polynomial |
   *-----------------------------------------------------------*/

   /* The degree N of the polynomial is computed using standard (double)
      real arithmetic (IEEE) with adding an rounding error of 1 digit of
      the last place for each operation according to the formula
           (T*T)^(N+1) <= 0.5*B^(1-L)*(2N+2)!/1.42  ,
      where L=Maxl+Lguard is the used number of mantissa digits with Maxl
      being the number of mantissa digits when the routine is entered.
      If the inequality above holds, the approximation error EAPP is bound
      by 0.5*B^(1-L).
      The test is started with N=1, i.e. the inequality T^2<=2K*B^(1-L).
      The loop counter n is containing the value N+1.           */

   L = Currprec + Lguard;          /* Number of digits used         */
/* u = max(nsqr+19,0) / 20;      *//* Number of additional digits   */
   nscaled = 0;
/* Lepsd.e = 1 - L - u;          *//* epsd = 0.704225*b^(1-L-u) */
/* Maxl = L+u;                   */
   Lepsd.e = 1 - L;                /* epsd = 0.704225*b^(1-L)   */
   Maxl = L;

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
         printf("\n  (sicovea 448) Number of scaling operations : %d",
                 nscaled);
         }
#endif

      rc = b_bcid(&Lepsd,&epsd,(a_intg)0); /* Conversion to double */
      if (rc==ROUND) rc = 0;
      NEXT_(t2,&dummy);
      i = b_bcid(&dummy,&xd,(a_intg)0); /* Conversion to double  */
      if (i==ROUND) {            /* upper bound by adding 1 ulp */
/* Cordes xd = Faddeps(xd);  */         /* if conversion is not exact  */
            R_ASSIGN(xd,r_succ(xd));    /* if conversion is not exact  */
            i = 0;
            }
      else rc += (int)i;

      if (rc != 0)
/*Cordes pd = -1.0;     */              /* Errors occurred during conv  */
         R_ASSIGN(pd,*r_mone);          /* Errors occurred during conv  */
       else
/*Cordes pd = xd;       */                 /* pd=xd             */
         R_ASSIGN(pd,xd);                  /* pd=xd             */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sicovea 469) Returncode conversion Lepsd->epsd : %d",
                rc);
         }
#endif

      }
    else
/* Cordes pd = -1.0;       */
      R_ASSIGN(pd,*r_mone);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 479) epsd = 2^(%d) = ",UL*Lepsd.e);
         Lprinti(&Lepsd);   printf(" = %g",epsd);
      printf("\n  (sicovea 481) xd = %g,    pd = %g",xd,pd);
      }
#endif

   Lpd  = Lvarget();                       /*     (== errsh)     */
   Lrhs = Lvarget();                       /*     (== sh)        */
   poolvars += 2;

/* Cordes if (pd>0) {         */
   if (r_sign(pd)>0) {
/* Cordes for (n=1; (nscaled > 0) | (pd >= epsd); n++) {   */
      for (n=1; (nscaled > 0) | (r_ge(pd,epsd)); n++) {
/*Cordes pd = Faddeps(pd*xd);      */   /* pd   = ^{pd*xd*xd}   */
         R_ASSIGN(pd,r_mulu(pd,xd));    /* pd   = ^{pd*xd*xd}   */
         i = 2*n*(2*n-1);
/* Cordes epsd = Fsubeps(epsd*i);  */   /* epsd = v{epsd*2n*(2n-1)} */
         R_ASSIGN(epsd,r_mulu(epsd,r_flot((a_intg)i)));
                                        /* epsd = v{epsd*2n*(2n-1)} */
/* Cordes if ((nscaled > 0) & (epsd > 1.0)) {   */
         if ((nscaled > 0) & (r_gt(epsd,*r_one_))) {
/* Cordes   epsd *= FIUFfac;   */
            R_ASSIGN(epsd,r_mulu(epsd,FIUFfac));
            nscaled--;
            }
/* Cordes if ((nscaled > 0) & (pd < 1.0)) {     */
         if ((nscaled > 0) & (r_lt(pd,*r_one_))) {
/* Cordes pd *= FUFfac;      */
            R_ASSIGN(pd,r_mulu(pd,FUFfac));
            nscaled--;
            }
#ifdef Debug
         if (Ldebug >= 0) {
            printf("\n  (sicovea 504) %g %g",epsd,pd);
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
         MUL_(Lpd,t2,Lpd,&dummy);  /* pd   = ^{pd*xd*xd}*/
         if (dummy.z==0) NEXT_(Lpd,Lpd);
         i = 2*n*(2*n-1);
         MULINT_(Lrhs,i,Lrhs); /* epsd = v{epsd*2n*(2n-1)}      */
         }
      }


  /*-----------------------------------*
   | Computation of the Rounding Error |
   *-----------------------------------*/

   /* The rounding error ERR is computed using (short, i.e. 'Minlerr' digits)
      INTERN numbers.                                       */

   Maxl = Minlerr;                      /* Set appropriate accuracy      */

   MUL_(t2,&lp675,&err,&dummy);         /* err = (.675*t^2 + */
   if (err.r != 0) NEXT_(&err,&err);
   ADD_(&l39o2,&err,&err);              /*        + 19.5)*B^(1-L)    */
   if (err.r != 0) NEXT_(&err,&err);
   NEXT_(&err,&err);                    /* quadratic terms compensation  */
   err.e += (1 - L - u);

   if (rc != 0)
      ERREXIT(EPERR,EPERR,poolvars);    /* Error message and handling    */

#ifdef Debug
   if (Ldebug >= 0 && pd == 0.0) {
      printf("\n  (sicovea 548) Underflow of pd");
      ERREXIT(DUFLW,DUFLW,poolvars);    /* Error message and handling    */
      }

   if (Ldebug >= 0) {
      printf("\n  (sicovea 553) err = ");  Lprinti(&err);
      printf("\n  (sicovea 554) Approximation polynomial degree %d",n);
      }
#endif


  /*----------------------------------------*
   | Evaluation of Approximation Polynomial |
   *----------------------------------------*/

  /* Evaluation of the Taylor series of the sine and cosine function for the
     reduced argument T using Horner's scheme.
     The coefficients (2N)! and (2N+1)! are generated during evaluation by
     using I*T instead of T for the general Horner's scheme, I=2N*(2N+1)
     or I=(2N-1)*2N for sine and cosine, respectively.               */

   Maxl = L + u;                   /* Accuracy polynomial evaluation*/

   sin = Lvarget();
   cos = Lvarget();
   poolvars += 2;
   (void)b_bini(sin);
   (void)b_bini(cos);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 583) t  = ");  Lprinti(t);
      printf("\n  (sicovea 584) t2 = ");  Lprinti(t2);
      }
#endif

   rc = 0;
   k  = 2*n;
   j  = k*(k+1);

#ifdef Debug
   if (Ldebug >= 0) printf("\n   j = %d",j);
#endif

   DIVINT_(t2,j,sin);              /*  sin = t^2/(2n*(2n+1))    */

   j  = k*(k-1);

#ifdef Debug
   if (Ldebug >= 0) printf("\n   j = %d",j);
#endif

   DIVINT_(t2,j,cos);           /*  cos = t^2/((2n-1)*2n)    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 608) sinN = ");  Lprinti(sin);
      printf("\n  (sicovea 609) cosN = ");  Lprinti(cos);
      }
#endif

   for (i=n-1; i>0; i--) {
       k = 2*i;
       SUB_(&Lone,sin,sin);     /* sin = (1+sin)*t^2/(2i*(2i+1))*/
       MUL_(sin,t2,sin,&dummy);
       j  = k*(k+1);
       DIVINT_(sin,j,sin);
       SUB_(&Lone,cos,cos);     /* cos = (1+cos)*t^2/((2i-1)*2i)*/
       MUL_(cos,t2,cos,&dummy);
       j  = k*(k-1);
       DIVINT_(cos,j,cos);

#ifdef Debug
       if (Ldebug >= 0) {
          printf("\n  (sicovea 626) sinI = ");  Lprinti(sin);
          printf("\n  (sicovea 627) cosI = ");  Lprinti(cos);
          }
#endif

       }
   SUB_(&Lone,sin,sin);         /* sin = t*(1-sinI)                 */
   MUL_(sin,t,sin,&dummy);
   SUB_(&Lone,cos,cos);         /* cos = 1-cosI                     */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 638) sin = ");  Lprinti(sin);
      printf("\n  (sicovea 639) cos = ");  Lprinti(cos);
      printf("\n  (sicovea 640) nhalve = %d",nhalve);
      }
#endif


  /*---------------------------------*
   | Computation of the Final Result |
   *---------------------------------*/

  /* The final result is computed by squaring the value for the reduced
     argument NSQR times, where NSQR is the number of halvings during
     computation of the reduced argument.               */

   fac = Lvarget();
   poolvars++;
   COPY_(&Lone,fac);

   if (rc != 0)
      ERREXIT(RESUL,662,poolvars); /* Error message and handling       */

   for (i=1; i<=nhalve; i++) {
      rc = 0;
      Maxl = L + u;                /* Accuracy result adaptation       */
      COPY_(sin,&LhF);
      MUL_(sin,cos,sin,&dummy); /* sin(2t) = 2*sin(t)*cos(t)    */
      SHIFT_(1,sin,sin);
      MUL_(&LhF,&LhF,&LhF,&dummy); /* cos(2t) = cos^2(t)-sin^2(t)*/
      MUL_(cos,cos,cos,&dummy);
      SUB_(cos,&LhF,cos);

      Maxl = Minlerr;
      MUL_(&L2sqt2,fac,fac,&dummy); /* fac = fac*2*sqrt(2)      */
      NEXT_(fac,fac);

      if (rc!=0) ERREXIT(RESUL,678,poolvars);
      }

   if (nhalve > 0) {
      dvnd = Lvarget();
      poolvars++;
      MUL_(&err,fac,&err,&dummy);  /* err = fac^n*err + */
      NEXT_(&err,&err);
      SUB_(fac,&Lone,dvnd);        /* +2*(fac^n-1)/(fac-1)*eps */
      NEXT_(dvnd,dvnd);
      DIV_(dvnd,&L2sqt2m1,fac_el);
      NEXT_(fac_el,fac_el);
      SHIFT_(1,fac_el,fac_el);
      NEXT_(fac_el,fac_el);
      fac_el->e += (1-(L+u));
      ADD_(&err,fac_el,&err);

      if (err.r != 0) NEXT_(&err,&err);
      NEXT_(&err,&err);           /* quadratic terms compensation  */
      }

   if (rc != 0)  ERREXIT(EPERR,EPERR,poolvars);

   Lrnd = rounded;                /* Value is rounded              */

   if (((ep+1) & MOD2) == 0) {
      COPY_(sin,&sinvea);
      COPY_(cos,&cosvea);
      }
    else {
      COPY_(cos,&sinvea);
      COPY_(sin,&cosvea);
      }

   if ((ep & MOD4) != 0) sinvea.s = 1;
   if (((ep+2) & MOD4) != 0) cosvea.s = 1;

   if (Larg->s != 0) {
      LhI = 15-ep;
      sinvea.s = (sinvea.s + 1) % 2;
      }
    else
      LhI = ep;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sicovea 722) sin = ");  Lprinti(&LhF);
      printf("\n  (sicovea 723) cos = ");  Lprinti(&LhD);
      printf("\n  (sicovea 724) err = ");  Lprinti(&err);
      printf("\n  (sicovea 725) entier part : %d",LhI);
      }
#endif

   Lvardrop(poolvars);
   RETURN(rc);
}





