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

/* CVS $Id: b_loga.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lloga              Processor : C                   *
 *                                                                       *
 * Logarithm to given base for Multiple Precision Arithmetic             *
 * =========================================================             *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Argument Range :    All positive numbers                              *
 * Base Range :        All positive numbers                              *
 * Result Range :      All numbers                                       *
 *                                                                       *
 * Used Functions :                                                      *
 *    lnve            Intern Logarithm Function                          *
 *    Lginit          Initialization of global constants                 *
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

#define  Name        "Lloga"
#define  Main

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
 * Constants for Computation of arsinh(x)                            *
 *********************************************************************/

#define Lguard    3

#define  rerr        "Argument 0 or negative"

static a_btyp  mbapp[3]  = { 0x18000000L,0x00000000L,0x00000000L };

static dynamic   bdapp   = { 0, 0, 0, 0, -1, Minl, &mbapp[0] };

static char  *function = Name;
static char  *Rerr     = rerr;

/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/

#ifdef LINT_ARGS
int b_loga(dynamic *xi,dynamic *bi,dynamic *ri)
#else
int b_loga(xi,bi,ri)

dynamic *xi;
dynamic *bi;
dynamic *ri;
#endif
#define  LRoutine    "Lloga"
{
   extern dynamic   dummy, Lone;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine, *Lerrmsg;
   extern dynamic   Lres, err;

#ifdef Debug
   extern int       Ldebug;
#endif

   dynamic    *y, *lny, *erry;
   int        rc, poolvars, ysgn;

  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level          */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting       */


  /*-----------------------------*
   | Check and Handle Case bi<=0 |
   *-----------------------------*/

   if (bi->s == 1 || bi->z == 1) {
      Lerrmsg = Rerr;       /* Error message for non positive argument  */
      ERREXIT(RANGE,143,0); /* Error message and handling               */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (bi->m[0] == 0) {
      ERREXIT(NANDE,152,0); /* Error message and handling               */
      }


  /*----------------------------*
   | Check and Handle Case bi=1 |
   *----------------------------*/

   rc = 0;
   if EQ_(bi,&Lone) {
      ERREXIT(RANGE,RANGE,0);
      }


  /*-----------------------------*
   | Check and Handle Case xi<=0 |
   *-----------------------------*/

   if (xi->s == 1 || xi->z == 1) {
      Lerrmsg = Rerr;       /* Error message for non positive argument  */
      ERREXIT(RANGE,172,0); /* Error message and handling               */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,181,0); /* Error message and handling               */
      }


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   if EQ_(xi,&Lone) {
      ri->z = 1;            /* ri = 0           */
      EXIT(rc);
      }


  /*----------------------*
   | Special case xi = bi |
   *----------------------*/

   if EQ_(xi,bi) {
      COPY_(&Lone,ri);      /* ri = 1           */
      EXIT(rc);         /* Error message and handling               */
      }


  /*---------------------------------------------------*
   | Check and Handle Different Ranges for Argument XI |
   *---------------------------------------------------*/

   y = Lvarget();
   poolvars = 1;

   rc = 0;
   SUB_(xi,&Lone,y);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lloga 217) Larg-1 = ");  Lprinti(y);
      }
#endif

   ysgn = y->s;
   y->s = 0;

   if LE_(y,&bdapp) {      /* Approxmation if |y| > 0.09375 */

     /*------------------------------------------*
      | Initializations for Taylor Approximation |
      *------------------------------------------*/

      y->s = ysgn;         /* y = (1 - Larg) / (1 + Larg)   */
      Maxl = min(Lcurrprec,xi->l) + 1;
      ADD_(xi,&Lone,&dummy);    /* EXACT sum             */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lloga 236) Maxl = %d",Maxl);
         printf("\n  (Lloga 237) Larg+1 = ");  Lprinti(&dummy);
         }
#endif

      Maxl = Lcurrprec + Lguard;
      DIV_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lloga 246) Maxl = %d",Maxl);
         printf("\n  (Lloga 247) (Larg-1)/(Larg+1) = ");  Lprinti(y);
         }
#endif

      if (rc != 0) ERREXIT(PEVAL,251,poolvars);

      LN_T_(y);         /* Compute function value and error bound   */
      }
    else
      LN_(xi);          /* Compute function value and error bound   */

   lny  = Lvarget();
   erry = Lvarget();
   poolvars += 2;
   Maxl = Lres.l;
   COPY_(&Lres,lny);
   Maxl = err.l;
   COPY_(&err,erry);
   erry->s = 0;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lloga 269) ln(xi) = ");  Lprinti(&Lres);
      printf("\n  (Lloga 270)        = ");  Lprinti(lny);
      printf("\n  (Lloga 271) err    = ");  Lprinti(&err);
      printf("\n  (Lloga 272)        = ");  Lprinti(erry);
      }
#endif


  /*-----------------------------------------------*
   | Check and Handle Different Ranges for Base BI |
   *-----------------------------------------------*/

   Maxl = Lcurrprec;
   SUB_(bi,&Lone,y);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lloga 286) Lbas-1 = ");  Lprinti(y);
      }
#endif

   ysgn = y->s;
   y->s = 0;

   if LE_(y,&bdapp) {              /* Approxmation if |y| > 0.09375 */

     /*------------------------------------------*
      | Initializations for Taylor Approximation |
      *------------------------------------------*/

      y->s = ysgn;                 /* y = (1 - Larg) / (1 + Larg)   */
      Maxl = min(Lcurrprec,bi->l) + 1;
      ADD_(bi,&Lone,&dummy);       /* EXACT sum                     */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lloga 305) Maxl = %d",Maxl);
         printf("\n  (Lloga 306) Larg+1 = ");  Lprinti(&dummy);
         }
#endif

      Maxl = Lcurrprec + Lguard;
      DIV_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lloga 315) Maxl = %d",Maxl);
         printf("\n  (Lloga 316) (Larg-1)/(Larg+1) = ");  Lprinti(y);
         }
#endif

      if (rc != 0) ERREXIT(PEVAL,320,poolvars);

      LN_T_(y);         /* Compute function value and error bound   */
      }
    else
      LN_(bi);          /* Compute function value and error bound   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lloga 329) ln(bi) = ");  Lprinti(&Lres);
      printf("\n  (Lloga 330) err    = ");  Lprinti(&err);
      }
#endif


   Maxl = Lcurrprec + Lguard;
   DIV_(lny,&Lres,&Lres); /* log_bi(xi) = ln(xi)/ln(bi)            */

   err.s = 0;             /* err = err(ln y) + err(ln bi) + eps(l) */
   ADD_(&err,erry,&err);
   NEXT_(&err,&err);
   COPY_(&Lone,y);
   y->e = 1 - Lcurrprec - Lguard;
   ADD_(&err,y,&err);
   NEXT_(&err,&err);


   ASSIGN_(ri); /* Assign result and set number of ulp's    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lloga 351) Lres   = ");  Lprinti(&Lres);
      printf("\n  (Lloga 352) err    = ");  Lprinti(&err);
      printf("\n  (Lloga 353) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   Lvardrop(poolvars);
   RETURN(0);
}





