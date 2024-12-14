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

/* CVS $Id: b_asin.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Larcsin            Processor : C                   *
 *                                                                       *
 * Inverse Sine and Cosine for Multiple Precision Arithmetic             *
 * =========================================================             *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Entries:                                                              *
 *    Function Call              Function               Formula          *
 *    -------------------------  ---------------------  ---------------- *
 *    Larcsin(Arg,Res)           Inverse Sine           Res=arcsin(Arg)  *
 *    Larccos(Arg,Res)           Inverse Cosine         Res=arccos(Arg)  *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Argument Range :    All numbers                                       *
 * Result Range :      All numbers                                       *
 *                                                                       *
 * Used Functions :                                                      *
 *    On computation of function values the computed value and the error *
 *    bound are returned by the global variables LhF and LhE.   However, *
 *    in many cases these names are redefined to a symbolic name related *
 *    to the called function.                                            *
 *                                                                       *
 *    name      macro    description                                     *
 *    --------  -------  ------------------------------------------------*
 *    arcsinve  ARCSIN_  Intern Inverse Sine                             *
 *    Lassign   ASSIGN_  Assignment of function value and number of ulp's*
 *    Lginit    (none)   Initialization of global constants              *
 *                     Intern arithmetic (macro names see file b_lari.h) *
 *                                                                       *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhF, LhE, LhD                                        (Variables)   *
 *    Lone              ( = 1 )                            (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#define  Name        "Larcsin"
#define  Main

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

static char  *function = Name;

/************************************************************************
 * Redefinition of Names of Intermediate Variables Used                 *
 ************************************************************************/

#define Lres      LhF           /* Function value                       */
#define err       LhE           /* Rounding error                       */
#define sqty      LhF           /* Function value sqrt(y)               */
#define dummy     LhD           /* Dummy variable (not further used)    */

/************************************************************************
 * Constants for Computation of arcsin(x)                               *
 ************************************************************************/

#define Lguard    2
#define PiMore    20
#define M0        4
#define summ2     dvsr
#define summ1     dvnd
#define erri      dvnd

static a_btyp  mtwo[3]   = { 0x00000002L,0x00000000L,0x00000000L };
static a_btyp  msqth[3]  = { 0xB504F333L,0xF9DE6484L,0x597D89B4L };

static dynamic  Ltwo  = { 0, 0, 0, 0, 1, Minl, &mtwo[0] };
static dynamic  Lsqth = { 0, 0, 0, 0,-1, Minlfl, &msqth[0] };

static dynamic  LPiOv2 = { 0, 0, 1, 0, 1, 0, NULL };

/************************************************************************
 * Algorithm arcsin                                                     *
 ************************************************************************/

#ifdef LINT_ARGS
int b_asin(dynamic *xi,dynamic *ri)
#else
int b_asin(xi,ri)
dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larcsin"
{
   extern dynamic Lone, LPiov4, Lres, err;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char    *Lroutine;

   int       rc;

#define  poolvars   0

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

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting               */
   rc = 0;


  /**********************************************************************
   * Check and Handle Range Error and Special Cases                     *
   **********************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z == 1) {
      if ((rc = b_bini(ri)) == 0) {        /* ri = arcsin(0) = 0        */
         ri->r = 0;                        /* Result exact              */
         EXIT(0);
         }
       else
         ERREXIT(rc,rc,0);      /* Error message and handling    */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling           */
      }


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   if EQ_ABS_(xi,&Lone) {
      if (LPiov4.l <= Maxl) {              /* ri = sgn(xi)*pi/2         */
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0);     /* Error message and handling*/
         }
/* Cordes : sign of xi is directly assigned to ri */
/*    prec = xi->s;          */            /* save sign of xi           */
      SHIFT_(1,&LPiov4,ri);
      if (rc != 0) ERREXIT(rc,rc,0);       /* Error message and handling*/
/* Cordes : sign of xi is directly assigned to ri */
/*    ri->s = prec;     */
      ri->s = xi->s;
      ri->r = 1;
      EXIT(0);
      }


  /*-------------*
   | Check Range |
   *-------------*/

   if GT_ABS_(xi,&Lone) {
      ERREXIT(RANGE,RANGE,0);   /* Error message and handling           */
      }



  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcsin 230) Maxl = %d\n",Maxl);
      printf("\n  (Larcsin 231) xi   = ");  Lprinti(xi);
      }
#endif

   ARCSIN_(xi);

   if GE_ABS_(xi,&Lsqth) {              /* Result adaptation            */
      Maxl = Lcurrprec + Lguard;
      Lres.s = 0;
      SUB_(&LPiOv2,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }

   Lres.s = xi->s;                      /* Set sign of result           */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcsin 251) Maxl = %d\n",Maxl);
      printf("\n  (Larcsin 252) xi   = ");  Lprinti(xi);
      printf("\n  (Larcsin 253) Lres = ");  Lprinti(&Lres);
      printf("\n  (Larcsin 254) err  = ");  Lprinti(&err);
      printf("\n");
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,260,0);                /* Error message and handling*/
      }

   ASSIGN_(ri);                 /* Assign result and set number of ulp's*/

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larcsin 267) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larcsin 268) LhE    = ");  Lprinti(&err);
      printf("\n  (Larcsin 269) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}
#undef  Name
#undef  LRoutine
#undef  poolvars

#define  Name        "Larccos"

static char  *functionc = (Name);

/************************************************************************
 * Algorithm arccos(x)                                                  *
 ************************************************************************/

#ifdef LINT_ARGS
int b_acos(dynamic *xi,dynamic *ri)
#else
int b_acos(xi,ri)
dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larccos"
{
   extern dynamic   Lone, LPiov4, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;

   int       rc;

#define  poolvars   0

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

   Lroutine = functionc;
   Lcurrprec = Maxl;            /* Save precision setting               */
   rc = 0;


  /**********************************************************************
   * Check and Handle Range Error and Special Cases                     *
   **********************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z == 1) {
      if (LPiov4.l <= Maxl) {              /* ri = arccos(0) = pi/2     */
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0);     /* Error message and handling*/
         }
      SHIFT_(1,&LPiov4,ri);
      if (rc != 0) ERREXIT(rc,rc,0);       /* Error message and handling*/
      ri->s = 0;
      ri->r = 1;
      EXIT(0);
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling           */
      }


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   if EQ_ABS_(xi,&Lone) {
      if (xi->s == 0) {
         if ((rc = b_bini(ri)) == 0) {      /* ri = arccos(1) = 0       */
            ri->r = 0;                     /* Result exact              */
            EXIT(0);
            }
          else
            ERREXIT(rc,rc,0);              /* Error message and handling*/
         }
       else {
         if (LPiov4.l <= Maxl) {           /* ri = arccos(-1) = pi      */
            Lcurrprec = Maxl;
            Maxl += PiMore;
            rc = LpiGen();
            Maxl = Lcurrprec;
            if (rc != 0) ERREXIT(0,rc,0);  /* Error message and handling*/
            }
         SHIFT_(2,&LPiov4,ri);
         if (rc != 0) ERREXIT(rc,rc,0);    /* Error message and handling*/
         ri->s = 0;
         ri->r = 1;
         EXIT(0);
         }
      }


  /*-------------*
   | Check Range |
   *-------------*/

   if GT_ABS_(xi,&Lone) {
      ERREXIT(RANGE,RANGE,0);   /* Error message and handling           */
      }



  /**********************************************************************
   * Initialization of Global Constants and Storage Locations           *
   **********************************************************************/

   if (Lgiflag==0) Lginit();               /* Initialization of GLOBALS */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larccos 435) Maxl = %d\n",Maxl);
      printf("\n  (Larccos 436) xi   = ");  Lprinti(xi);
      }
#endif

   ARCSIN_(xi);

   if LE_ABS_(xi,&Lsqth) {                 /* Result adaptation         */
      if (LPiOv2.l < Lcurrprec+Lguard) {   /* ri = pi/2 - arcsin(xi)    */
         if (LPiov4.l <= Lcurrprec+Lguard) {
            Maxl = Lcurrprec + Lguard + PiMore;
            rc = LpiGen();
            if (rc != 0) ERREXIT(0,rc,0);  /* Error message and handling*/
            }
         Maxl = Lcurrprec + Lguard;
         SHIFT_(1,&LPiov4,&LPiOv2);
         if (rc != 0) ERREXIT(rc,rc,0);    /* Error message and handling*/
         }
       else
         Maxl = Lcurrprec + Lguard;

      SUB_(&LPiOv2,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }
    else if (xi->s == 1) {            /* ri = pi-arcsin(sqrt(1-xi^2))   */
      if (LPiov4.l <= Lcurrprec+Lguard) {
         Maxl = Lcurrprec + Lguard + PiMore;
         rc = LpiGen();
         if (rc != 0) ERREXIT(0,rc,0);     /* Error message and handling*/
         }
      Maxl = Lcurrprec + Lguard;
      SHIFT_(2,&LPiov4,&dummy);
      if (rc != 0) ERREXIT(rc,rc,0);       /* Error message and handling*/

      Lres.s = 0;
      SUB_(&dummy,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }

   Lres.s = 0;                             /* Set sign of result        */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larccos 484) Maxl = %d\n",Maxl);
      printf("\n  (Larccos 485) xi   = ");  Lprinti(xi);
      printf("\n  (Larccos 486) Lres = ");  Lprinti(&Lres);
      printf("\n  (Larccos 487) err  = ");  Lprinti(&err);
      printf("\n");
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,493,0);                /* Error message and handling*/
      }

   ASSIGN_(ri);                 /* Assign result and set number of ulp's*/

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larccos 500) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larccos 501) LhE    = ");  Lprinti(&err);
      printf("\n  (Larccos 502) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





