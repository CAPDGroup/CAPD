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

/* CVS $Id: b_lnve.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : lnve               Processor : C                   *
 *                                                                       *
 * Internal Logarithm Function for Multiple Precision Arithmetic         *
 * =============================================================         *
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

/********************************************************************
 * Redefinition of Names of Intermediate Variables Used             *
 ********************************************************************/
        
#define Lres      LhF   /* Function value                           */
#define err       LhE   /* Rounding error                           */
#define expy      LhF   /* Function value exp(y)                    */
#define errexp    LhE   /* Rounding error of exp(y)                 */
#define dummy     LhD   /* Dummy variable (not further used)        */
        
/********************************************************************
 * Constants for Computation of ln(x)                               *
 ********************************************************************/
        
#define Lguard    2
#define M0        4
#define Res       erry

static a_btyp  macc[3]   = { 0x00100000L,0x00000000L,0x00000000L };
        
static dynamic  acc     = { 0, 0, 0, 0,  0, Minl, &macc[0] };
        
/********************************************************************
 * Algorithm                                                        *
 ********************************************************************/
        
#ifdef LINT_ARGS
int b_lnve(dynamic *Larg)
#else
int b_lnve(Larg)

dynamic *Larg;
#endif
#define  LRoutine    "lnve"
{
   extern dynamic   Llnbase, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lgiflag;
   extern rounding  Lrnd;

#ifdef Debug
   extern int       Ldebug;
#endif
        
   dynamic   *y, *erry, *Lhlp;
   int       rc, poolvars;
   a_intg    argexp;
   a_btyp    Currprec, prec;
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

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnve 142) Larg = ");  Lprinti(Larg);
      printf("\n  (lnve 143) Maxl = %d",Maxl);
      }
#endif

  /***************************************************************
   * Computation of an Initial Approximation                     *
   ***************************************************************/
        
  /* An initial approximation is computed by taking the floating point
     logarithm of the argument Larg.                             */
        
   y = Lvarget();
   poolvars = 1;
   Maxl = Minlfl;

   if ((Larg->e > 0) | (Larg->e < -1)) {
      argexp = Larg->e;         /* save exponent of Larg         */
      Larg->e = 0;              /* set exponent of Larg to 0     */
      rc = b_bcid(Larg,&yd,(a_intg)0);
      if (rc==ROUND) rc = 0;
      rc += b_bcdi(B_LOG_(yd),&y,(a_intg)0);
                                /* y0 = ln(b*mant(Larg))    */
      Larg->e = argexp;
/* Cordes      MULINT_(&Llnbase,abs(argexp),&LhD);      */
      if (argexp >= 0)          /* y = ln(b*mant(Larg)) +        */
         {
         MULINT_(&Llnbase,argexp,&LhD);
         ADD_(y,&LhD,y);        /*       Larg.e*ln(b)            */
         }
       else
         {
         MULINT_(&Llnbase,-argexp,&LhD);
         SUB_(y,&LhD,y);
         }
      }
    else {
      rc = b_bcid(Larg,&yd,(a_intg)0);
      if (rc==ROUND) rc = 0;
      rc += b_bcdi(B_LOG_(yd),&y,(a_intg)0);
                                /* y = ln(b*mant(Larg))      */
      }
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnve 187) y = ln(Larg) = ");
      Lprinti(y);
      }
#endif

        
  /**************************************************************
   * Newton Iteration Steps for gaining the required accuracy   *
   **************************************************************/
        
  /* The accuracy of the initial approximation is increased by a number
     of Newton steps until the relative error bound is less than 1% of
     the initial accuracy setting.                              */
        
   acc.e = -Currprec;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnve 209) acc = ");
      Lprinti(&acc);
      }
#endif

   Maxl = Minlerr;
   erry = Lvarget();
   poolvars += 1;
   COPY_(&acc,erry);
   NEXT_(erry,erry);
   erry->s = 0;

   Lhlp = Lvarget();
   poolvars += 1;

   for (prec=Maxl=M0; GT_(erry,&acc);
        prec=Maxl=min(2*Maxl,Currprec+2)) {

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (lnve 228) Maxl = %d",Maxl);
         }
#endif

      if (rc != 0) ERREXIT(PEVAL,232,poolvars);
        
      EXP_(y);                  /* Lhlp := exp(y) - Larg            */
      SUB_(&expy,Larg,Lhlp);

      Maxl = Minlerr;
      NEXT_(Lhlp,&dummy);       /* determine error bound erry       */
      if GE_(Larg,&expy)
         DIV_(&dummy,&expy,erry);
       else
         DIV_(&dummy,Larg,erry);
      DIV_(erry,y,erry);
      erry->s = 0;

      Maxl = prec;
      DIV_(Lhlp,&expy,&dummy); /* y := y - [exp(y)-Larg]/exp(y) */
      SUB_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (lnve 252) y   = ");  Lprinti(y);
         printf("\n  (lnve 253) err = ");  Lprinti(erry);
         }
#endif

      }
        
#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnve 261) Iterations performed !");
      printf("\n");
      }
#endif

   if (rc != 0) ERREXIT(RESUL,266,poolvars);

   if ((rc = COPY_(y,Res)) != 0)  ERREXIT(rc,268,poolvars);

   Maxl = Currprec+2;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnve 274) Res  = ");  Lprinti(Res);
      printf("\n");
      }
#endif
        
        
  /**************************************************************
   * Determine number of units to be added for an upper bound   *
   **************************************************************/
        
   EXP_(Res);                   /* compute exp(Res)             */

   Maxl = Minlerr;
   MUL_(&expy,&errexp,&err,&dummy);
   NEXT_(&err,&err);
   SUB_(&expy,Larg,Lhlp);
   NEXT_(Lhlp,Lhlp);
   err.s = Lhlp->s;
   ADD_(Lhlp,&err,&err);
   NEXT_(&err,&err);
   if GE_(Larg,&expy)
      DIV_(&err,&expy,&err);
    else
      DIV_(&err,Larg,&err);
   NEXT_(&err,&err);
   DIV_(&err,y,&err);
   NEXT_(&err,&err);

/* Lrnd = signed; */
   Lrnd = rounded;

   Maxl = Currprec+2;
   COPY_(Res,&Lres);

   if (rc != 0)  ERREXIT(EPERR,317,poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (lnve 321) Lres = ");  Lprinti(&Lres);
      }
#endif

   Lvardrop(poolvars);
   RETURN(rc);
}





