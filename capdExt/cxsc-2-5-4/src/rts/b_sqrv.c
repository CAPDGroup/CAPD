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

/* CVS $Id: b_sqrv.c,v 1.22 2014/01/30 17:24:05 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : sqrtve             Processor : C                   *
 *                                                                       *
 * Internal Square Root Function for Multiple Precision Arithmetic       *
 * ===============================================================       *
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
        
/*************************************************************************
 * Redefinition of Names of Intermediate Variables Used                  *
 *************************************************************************/
        
#define Lres      LhF           /* Function value                        */
#define err       LhE           /* Rounding error                        */
#define dummy     LhD           /* Dummy variable (not further used)     */
        
        
/****************************************************************
 * Constants for Computation of sqrt(x)                         *
 ****************************************************************/
        
#define Lguard    2
#define M0        4
        
static a_btyp  macc[3]   = { 0x00100000L,0x00000000L,0x00000000L };
static a_btyp  mtwo[3]   = { 0x00000002L,0x00000000L,0x00000000L };
        
static dynamic  acc   = { 0, 0, 0, 0, 0, Minl, &macc[0] };
static dynamic  Ltwo  = { 0, 0, 0, 0, 0, Minl, &mtwo[0] };
        
/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/
        
#ifdef LINT_ARGS
int b_sqrv(dynamic *Larg)
#else
int b_sqrv(Larg)

dynamic *Larg;
#endif
#define  LRoutine    "sqrtve"
{
   extern dynamic   Lone, Leps, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lgiflag;
   extern rounding  Lrnd;

#ifdef Debug
   extern int       Ldebug;
#endif
        
   dynamic   *y, *t, *Lhlp;
   int       rc, poolvars;
   a_intg    nexp, nmul2;
   a_btyp    Currprec, u;
   a_real    yd;
        
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
   * Argument Reduction                                         *
   **************************************************************/
        
  /* Argument reduction is performed by halving the argument until the
     reduced argument is part of the interval [.5,2] .
     digits and/or bits exceeding the specified accuracy.       */
        
   t = Lvarget();
   poolvars = 1;
        
   Maxl = Larg->l + 1;          /* Set appropriate accuracy      */
   rc = 0;
   if ((COPY_(Larg,t)) != 0) {  /* t = |arg|                     */
      ERREXIT(rc,157,poolvars); /* Error Exit                    */
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 162) t = ");
      Lprinti(t);
      printf(";   rc = %d",rc);
      }
#endif

   nmul2 = 0;
   for (u=t->m[0]; (u & MSB) == 0; u <<= 1) nmul2++;
   nmul2 = (B_LENGTH - nmul2) / 2;
   t->e = 0;
   if ((SHIFT_(-2*nmul2,t,t)) != 0) {
      ERREXIT(rc,173,poolvars);            /* Error Exit        */
      }
   nmul2 += (B_LENGTH/2)*(Larg->e % 2);
   nexp   = Larg->e / 2;
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 180) Number of halvings = %d",nmul2);
      printf("\n  (sqrtve 181) Reduced Argument = "); Lprinti(t);
      }
#endif
        
        
  /**************************************************************
   * Computation of an Initial Approximation                    *
   **************************************************************/
        
  /* An initial approximation is computed by taking the floating point
     square root of the reduced argument T.
     The error EPSo of this initial approximation Yo is computed according
     to the formula
         EPSo <= (Yo*Yo - T)/(T + min{T,Yo*Yo}  .               */
        
   Maxl = Minlfl;
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 204) t           = ");
      Lprinti(t);
      }
#endif

   y = Lvarget();
   poolvars += 1;
        
   if ((rc = b_bcid(t,&yd,(a_intg)0))==ROUND) rc = 0;

#ifdef Debug
   if (Ldebug >= 0) {
      printf(";   rc = %d",rc);
      }
#endif

   rc += b_bcdi(B_SQRT(yd),&y,(a_intg)0);
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 224) y = sqrt(t) = ");
      Lprinti(y);
      printf(";   rc = %d",rc);
      }
#endif
        
   Maxl = Minlflsqr;
   MUL_(y,y,&err,&dummy);
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 235) y*y = ");
      Lprinti(&err);
      printf(";   rc = %d",rc);
      }
#endif
        
   if (dummy.z == 0) NEXT_(&err,&err);
        
   Maxl = Minlerr;
   if GT_(&err,t)
      {
#ifdef Debug
   if (Ldebug>=0) printf("\n  (sqrtve 240) dummy = t << 1");
#endif
      SHIFT_(1,t,&dummy);
      }
   else
      {
#ifdef Debug
   if (Ldebug>=0) printf("\n  (sqrtve 241) dummy = t + y*y");
#endif
      ADD_(t,&err,&dummy);
      }
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 251) t + min{t,y*y} = ");
      Lprinti(&dummy);
      printf(";   rc = %d",rc);
      }
#endif
        
   SUB_(&err,t,&err);
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 261) y*y - t = ");
      Lprinti(&err);
      printf(";   rc = %d",rc);
      }
#endif
        
   DIV_(&err,&dummy,&err);
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 271) err = ");  Lprinti(&err);
      printf("\n  (sqrtve 272) nmul2 = %d,  rc = %d",nmul2,rc);
      }
#endif

   if (rc != 0) ERREXIT(PEVAL,276,poolvars);
   rc = 0;
        
        
  /**************************************************************
   * Newton Iteration Steps for gaining the required accuracy   *
   **************************************************************/
        
  /* The accuracy of the initial approximation is increased by a number
     of Newton steps until the relative error bound is less than 1% of
     the initial accuracy setting.                              */
        
   acc.e = -Currprec-1;
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 297) acc = ");
      Lprinti(&acc);
      }
#endif
        
   Leps.e = -Minlfl;
   Maxl = Minlerr;
   ADD_(&Lone,&Leps,&Lres);
   NEXT_(&Lres,&Lres);
   MUL_(&Lres,&Lres,&Lres,&dummy);
   NEXT_(&Lres,&Lres);
   if (err.s != 0) ADD_(&Lone,&err,&dummy);
              else SUB_(&Lone,&err,&dummy);
   DIV_(&Lres,&dummy,&Lres);
   NEXT_(&Lres,&Lres);
   SHIFT_(-1,&Lres,&Lres);
   NEXT_(&Lres,&Lres);
   err.s = 0;

   Lhlp = Lvarget();
        
   for (Maxl=M0; GT_(&err,&acc);
        Maxl=min(-2*(Leps.e-1),Currprec+2)) {
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sqrtve 322) Maxl = %d",Maxl);
         }
#endif
        
      DIV_(t,y,&Lres);
      ADD_(y,&Lres,y);
      SHIFT_(-1,y,y);
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sqrtve 332) y   = ");
         Lprinti(y);
         }
#endif

      Leps.e = 1 - Maxl;
      Maxl = Minlerr;
      MUL_(&err,&err,&err,&dummy);
      MUL_(&err,&Lres,&err,&dummy);
      ADD_(&Lres,&Leps,Lhlp);
      ADD_(Lhlp,&Ltwo,Lhlp);
      MUL_(Lhlp,&Leps,Lhlp,&dummy);
      ADD_(&err,Lhlp,&err);
        
#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (sqrtve 348) y   = ");  Lprinti(y);
         printf("\n  (sqrtve 349) err = ");  Lprinti(&err);
         }
#endif
        
      }
      Lvardrop(1);
        
        
  /**************************************************************
   * Computation of the Final Result                            *
   **************************************************************/
        
  /* The final result is computed by squaring the value for the reduced
     argument NSQR times, where NSQR is the number of halvings during
     computation of the reduced argument.                       */
        
   Maxl = Currprec + 2;
   SHIFT_(nmul2,y,y);
   y->e += nexp;
        
   if (rc != 0) ERREXIT(PEVAL,374,poolvars);
   rc = 0;

   Maxl = Currprec;
   COPY_(y,&Lres);
        

  /*----------------------------------------------------------*
   | Determine number of units to be added for an upper bound |
   *----------------------------------------------------------*/

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 392) err   = ");  Lprinti(&err);
      printf("\n  (sqrtve 393) dummy = ");  Lprinti(&dummy);
      }
#endif


   Maxl = Currprec+2;
/* Maxl = max(2*Currprec,Larg->l); */

   MUL_(&Lres,&Lres,&err,&dummy);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 405) err   = ");  Lprinti(&err);
      printf("\n  (sqrtve 406) dummy = ");  Lprinti(&dummy);
      }
#endif

   if ((dummy.z == 0) & GE_(&err,Larg)) NEXT_(&err,&err);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 414) err   = ");  Lprinti(&err);
      }
#endif

   if GT_(&err,t)
      {
#ifdef Debug
   if (Ldebug>=0) printf("\n  (sqrtve 423) dummy = Larg << 1");
#endif
      SHIFT_(1,Larg,&dummy);
      }
    else
       {
#ifdef Debug
   if (Ldebug>=0) printf("\n  (sqrtve 424) dummy = Larg + err");
#endif
      ADD_(Larg,&err,&dummy);
      }

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 425) err   = ");  Lprinti(&err);
      printf("\n  (sqrtve 426) dummy = ");  Lprinti(&dummy);
      }
#endif

   SUB_(Larg,&err,&err);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 434) err   = ");  Lprinti(&err);
      printf("\n  (sqrtve 435) dummy = ");  Lprinti(&dummy);
      }
#endif

   if (err.r == 1) NEXT_(&err,&err);

   Maxl = Minlerr;
   DIV_(&err,&dummy,&err);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 446) err   = ");  Lprinti(&err);
      }
#endif

   if (err.r == 1) NEXT_(&err,&err);

   if (rc != 0) ERREXIT(RESUL,452,poolvars);

   Lrnd = signed;

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (sqrtve 460) y   = ");  Lprinti(y);
      printf("\n  (sqrtve 461) err = ");  Lprinti(&err);
      printf("\n  (sqrtve 462) Lrnd = %d",Lrnd);
      }
#endif

   RETURN(rc);
}





