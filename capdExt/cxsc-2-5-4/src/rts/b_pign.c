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

/* CVS $Id: b_pign.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : LpiGen             Processor : C                   *
 *                                                                       *
 * Function computing PI with current acurracy                           *
 * ===========================================                           *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Used Functions :                                                      *
 *    Lsqrt           Internal Square Root subroutine                    *
 *    Lginit          Initialization of global constants                 *
 *    Lintern         Flag for internal subroutines                      *
 *    r_succ          Increasing a double real number by 1 ulp           *
 *                    Intern arithmetic                                  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LPiov4                                               (Variable)    *
 *    L4ovPi                                               (Variable)    *
 *    Lone              ( = 1 )                            (Constant)    *
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

/****************************************************************
 * Redefinition of Names of Intermediate Variables Used         *
 ****************************************************************/

#define Lres      LhF   /* Computed approximation for pi            */
#define err       LhE   /* Total error pi                           */
#define dummy     LhD   /* Dummy variable for temporarily use only  */

/****************************************************************
 * Constants for Computation of pi                              *
 ****************************************************************/

#define Lguard    2
#define ip        an
#define epsI      epsa

static a_btyp mfour[3] = { 4, 0, 0 };
static a_btyp mtwo[3]  = { 2, 0, 0 };

static dynamic  Ltwo    = { 0, 0, 0, 0, 0, Minl, &mtwo[0] };
static dynamic  Lfour   = { 0, 0, 0, 0, 0, Minl, &mfour[0] };

/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_pign(void)
#else
int b_pign()
#endif
#define  LRoutine    "LpiGen"
{
   extern dynamic LPiov4, L4ovPi, Lone, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lgiflag, Lcurrprec;
   extern a_real  *r_one_;

   static a_btyp o0415[] = B_TYPINI( 0x3fda8f5cL, 0x28f5c290L );
   static a_btyp o0503[] = B_TYPINI( 0x3fe01893L, 0x74bc6a80L );
   static a_btyp o1001[] = B_TYPINI( 0x3ff00418L, 0x9374bc6bL );
   static a_btyp o1500[] = B_TYPINI( 0x3ff80000L, 0x00000000L );
   static a_btyp o1915[] = B_TYPINI( 0x3ffea3d7L, 0x0a3d70a5L );
   static a_btyp o4000[] = B_TYPINI( 0x40100000L, 0x00000000L );
   static a_btyp o5000[] = B_TYPINI( 0x40140000L, 0x00000000L );

   dynamic   *an, *bn, *pn, *Lha, *Lhb, *apperr, *rnderr, *errI, *upbound,
             *upboundI;
#ifdef Debug
   int       n;
#endif
   int       rc;
   a_btyp  ulpI, ulp, Currprec;
   int       poolvars;  /* !!! must be int !!! */
   a_real      eps, epsa, epsb, epsp;

#ifdef Debug
   extern int       Ldebug;
#endif

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level         */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*----------------*
   | Initialization |
   *----------------*/

   if (LPiov4.l == Lpi_init) {
      (void)b_bini(&LPiov4);
      (void)b_bini(&L4ovPi);
      }

        
   Currprec = Maxl;             /* Save precision setting       */
   rc = 0;
        

  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/
        
   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (LpiGen 144) Maxl = %d\n",Maxl);
      }
#endif


   an       = Lvarget();
   bn       = Lvarget();
   pn       = Lvarget();
   Lha      = Lvarget();
   Lhb      = Lvarget();
   apperr   = Lvarget();
   rnderr   = Lvarget();        /* Rounding error               */
   errI     = Lvarget();        /* Total error 4/pi             */
   upbound  = Lvarget();        /* Upper bound for 4*pi         */
   upboundI = Lvarget();        /* Upper bound for 4/pi         */
   poolvars = 10;
        
        
  /**************************************************************
   * Computation of the PI                                      *
   **************************************************************/

  /*----------------*
   | Initialization |
   *----------------*/
        
  /* The initial values for the iterative steps are computed and the first
     iteration step is performed in the interval [.5,2] .
     ......................................................     */
        
   Maxl += Lguard;                  /* Set appropriate accuracy      */
/* Lintern += 1;    */              /* Set flag for called routines  */
   ulp = Lcurrprec;
   if ((rc=Lsqrt(&Ltwo,an)) != 0) { /* a0 = sqrt(2)                  */
      Lcurrprec = ulp;
      ERREXIT(rc,183,poolvars);     /* Error message and handling    */
      }
   Lcurrprec = ulp;
/* Cordes   eps = an->r;        */      /* Error of root        */
   R_ASSIGN(eps,r_flot((a_intg)an->r)); /* Error of root        */
   ADD_(&Ltwo,an,pn);                   /* p0 = 2 + sqrt(2)     */
/* Cordes   epsp = Faddeps(1.0 + 0.415*eps);    */
/* possible allignment problem if a_real and a_btyp
   are not equally alligned */
   R_ASSIGN(epsp,r_addu(*r_one_,
                 r_mulu(*((a_real *)&o0415[0]),eps)));
   Maxl = Minlerr;
   DIVINT_(&Lone,10,apperr);
   Maxl = Currprec + Lguard;

   if (rc != 0) ERREXIT(RESUL,193,poolvars);
        
        
  /*-----------------*
   | Iteration Steps |
   *-----------------*/
        
  /* The accuracy of the initial approximation is increased by a number
     of Iteration steps until the relative error bound is less than 1% of
     the initial accuracy setting.                                     */
        
   ulp = Lcurrprec;
   rc  = Lsqrt(an,Lha);
   Lcurrprec = ulp;
        
   DIV_(&Lone,Lha,bn);                     /* b1 = 1/sqrt(a0)   */
/* Cordes   epsa = Faddeps(1.915 + 1.001*Lha->r + 0.503*eps);   */
/* possible allignment problem if a_real and a_btyp
   are not equally alligned */
   R_ASSIGN(epsa,r_addu(r_addu(*((a_real *)&o1915[0]),
        r_mulu(*((a_real *)&o1001[0]),r_flot((a_intg)Lha->r))),
        r_mulu(*((a_real *)&o0503[0]),eps)));

   ADD_(bn,Lha,an);             /* a1 = (sqrt(a0)+1/sqrt(a0))/2  */
   SHIFT_(-1,an,an);
/*Cordes epsb = Faddeps(1.0 + Lha->r + 0.503*eps); */
/* possible allignment problem if a_real and a_btyp
   are not equally alligned */
   R_ASSIGN(epsb,r_addu(r_addu(*r_one_,r_flot((a_intg)Lha->r)),
                      r_mulu(*((a_real *)&o0503[0]),eps)));
        
#ifdef Debug
   if (Ldebug >= 0) {
      n = 1;
      }
#endif

   for (;;) {
      ADD_(&Lone,bn,Lhb);       /* Lhb = 1 + bn                  */
      ADD_(&Lone,an,Lha);       /* Lha = 1 + an                  */
        
      MUL_(pn,bn,pn,&dummy);    /* pn = pn*bn*(1+an)/(1+bn)      */
      MUL_(pn,Lha,pn,&dummy);
      DIV_(pn,Lhb,pn);
/*Cordes epsp=Faddeps(Faddeps(Faddeps(epsp+epsa)+Faddeps(1.5*epsb))+5.0); */
/* possible allignment problem if a_real and a_btyp
   are not equally alligned */
      R_ASSIGN(epsp,r_addu(r_addu(r_addu(epsp,epsa),
                                r_mulu(*((a_real *)&o1500[0]),epsb)
                               ),
                         *((a_real *)&o5000[0])
                        )
           );

      Maxl = Minlerr;
      MUL_(apperr,apperr,apperr,&dummy);
      NEXT_(apperr,apperr);
      Maxl = Currprec + Lguard;
      if (apperr->e < -Currprec) break; /* Required accuracy reached    */
        
#ifdef Debug
      if (Ldebug >= 0) {
         n++;
         printf("\n  (LpiGen 243) n = %d",n);
         printf("\n  (LpiGen 244) pn = ");  Lprinti(pn);
         }
#endif
        
      ulp = Lcurrprec;
      rc += Lsqrt(an,Lha);                 /* Lha = sqrt(an)            */
      Lcurrprec = ulp;
/* Cordes      eps = Faddeps(0.503*epsa+Lha->r);        */
/* possible allignment problem if a_real and a_btyp
   are not equally alligned */
      R_ASSIGN(eps,r_addu(r_mulu(*((a_real *)&o0503[0]),epsa),
                        r_flot((a_intg)Lha->r)));
        
      ADD_(an,bn,bn);           /* bn = sqrt(an)*(1+bn)/(an+bn)  */
      DIV_(Lhb,bn,bn);
      MUL_(Lha,bn,bn,&dummy);
/* Cordes     epsb = Faddeps(Faddeps(Faddeps(epsb+epsa)+eps)+4.0); */
/* possible allignment problem if a_real and a_btyp
   are not equally alligned */
      R_ASSIGN(epsb,r_addu(r_addu(r_addu(epsb,epsa),eps),
                         *((a_real *)&o4000[0])
                        )
           );

        
      DIV_(&Lone,Lha,an);       /* an = (sqrt(an)+1/sqrt(an))/2  */
      ADD_(an,Lha,an);
      SHIFT_(-1,an,an);
/* Cordes      epsa = Faddeps(eps+1);   */
      epsa = r_addu(eps,*r_one_);
        
      }

/* Lintern -= 1;  */            /* Reset flag                    */

   DIV_(&Lfour,pn,ip);
        
   rc = b_bcdi(epsp,&rnderr,(a_intg)0);
   rnderr->e += (1-Maxl);
/* Cordes   retcode = b_bcdi(Faddeps(epsp+1),errI);     */
   rc += b_bcdi(r_addu(epsp,*r_one_),&errI,(a_intg)0);
   errI->e += (1-Maxl);
   Maxl = Minlerr;
   DIVINT_(apperr,3,apperr);
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (LpiGen 280)  pi%-2d     = ",n);  Lprinti(pn);
      printf("\n  (LpiGen 281)  rnderr  = ");  Lprinti(rnderr);
      printf("\n  (LpiGen 282)  apperr  = ");  Lprinti(apperr);
      printf("\n  (LpiGen 283)  1/pi    = ",n);  Lprinti(ip);
      printf("\n  (LpiGen 284)  rnderrI = ");  Lprinti(errI);
      }
#endif
        
        
  /*----------------------------------------------------------*
   | Determine number of units to be added for an upper bound |
   *----------------------------------------------------------*/
        
   ADD_(rnderr,apperr,&err);
   if (err.r == 1) NEXT_(&err,&err);
   COPY_(pn,&dummy);
   NEXT_(&dummy,&dummy);
   MUL_(&dummy,&err,&err,&dummy);
   if (err.r == 1) NEXT_(&err,&err);

   ADD_(errI,apperr,errI);
   if (errI->r == 1) NEXT_(errI,errI);
   COPY_(ip,&dummy);
   NEXT_(&dummy,&dummy);
   MUL_(&dummy,errI,errI,&dummy);
   if (errI->r == 1) NEXT_(errI,errI);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (LpiGen 314) err  = ");  Lprinti(&err);
      printf("\n  (LpiGen 315) errI = ");  Lprinti(errI);
      }
#endif


   Maxl = Currprec+1;
   SHIFT_(-2,pn,pn);
   Maxl = Currprec;
   SUB_(pn,&err,&LPiov4);
   SUB_(ip,errI,&L4ovPi);

   ADD_(pn,&err,upbound);
   if (upbound->r == 1) NEXT_(upbound,upbound);

   ADD_(ip,errI,upboundI);
   if (upboundI->r == 1) NEXT_(upboundI,upboundI);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (LpiGen 334) (pi/4)_l = ");  Lprinti(&LPiov4);
      printf("\n  (LpiGen 335) (pi/4)_u = ");  Lprinti(upbound);
      printf("\n  (LpiGen 336) (4/pi)_l = ");  Lprinti(&L4ovPi);
      printf("\n  (LpiGen 337) (4/pi)_u = ");  Lprinti(upboundI);
      }
#endif

        
   Maxl = Minlerr;

   SUB_(upbound,&LPiov4,&err);
   if (err.r == 1) NEXT_(&err,&err);
   err.e += Currprec;           /* ( = Currprec -1 - LPiov4.e )  */

   SUB_(upboundI,&L4ovPi,errI);
   if (errI->r == 1) NEXT_(errI,errI);
   errI->e += (Currprec-1);     /* ( = Currprec -1 - L4ovPi.e )  */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (LpiGen 354) err  = ");  Lprinti(&err);
      printf("\n  (LpiGen 355) errI = ");  Lprinti(errI);
      }
#endif

   for(ulp=0; compii(&err,&Lone)>=0 && ulp<4 ; ulp++)
      SUB_(&err,&Lone,&err);
   if (err.z == 0) ulp++;

   for(ulpI=0; compii(errI,&Lone)>=0 && ulpI<4 ; ulpI++)
      SUB_(errI,&Lone,errI);
   if (errI->z == 0) ulpI++;

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (LpiGen 371) ulp  = %d",ulp);
      printf("\n  (LpiGen 372) ulpI = %d",ulpI);
      }
#endif

   Maxl = Currprec;              /* Restore precision setting     */

   if ((ulp > 2) | (ulpI > 2)) {
      ERREXIT(WIDTH,379,0);      /* Error message and handling    */
      }

   LPiov4.r = (unsigned int)ulp; /* Return number of ulp's        */
   L4ovPi.r = (unsigned int)ulp;

   if (rc != 0) {
      ERREXIT(RESUL,386,0);      /* Error message and handling    */
      }
        
   RETURN(0);
}





