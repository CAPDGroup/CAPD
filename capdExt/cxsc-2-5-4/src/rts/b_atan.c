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

/* CVS $Id: b_atan.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Larctan            Processor : C                   *
 *                                                                       *
 * Inverse Tangent and Cotangent for Multiple Precision Arithmetic       *
 * ===============================================================       *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                  math.h    - standard definitions for double reals    *
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Entries:                                                              *
 *    Function Call              Function               Formula          *
 *    -------------------------  ---------------------  ---------------- *
 *    Larctan(Arg,Res)           Inverse Tangent        Res=arctan(Arg)  *
 *    Larccot(Arg,Res)           Inverse Cotangent      Res=arccot(Arg)  *
 *    Latan2(X,Y,Res)            Inverse Tangent X/Y    Res=arctan(X/Y)  *
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
 *    arctanve  ARCTAN_  Intern Inverse Tangent                          *
 *    Lassign   ASSIGN_  Assignment of function value and number of ulp's*
 *    Lginit    (none)   Initialization of global constants              *
 *                    Intern arithmetic (macro names see file b_lari.h)  *
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

#define  Name        "Larctan"
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
 * Constants for Computation of arctan(x)                               *
 ************************************************************************/

#define Lguard    2
#define PiMore    20
#define M0        4
#define erri      dvnd

static a_btyp  mtwo[3]   = { 0x00000002,0x00000000,0x00000000 };

static dynamic  Ltwo  = { 0, 0, 0, 0, 1, Minl, &mtwo[0] };

static dynamic  LPiOv2 = { 0, 0, 1, 0, 1, 0, NULL };

/************************************************************************
 * Algorithm arctan                                                     *
 ************************************************************************/

#ifdef LINT_ARGS
int b_atan(dynamic *xi,dynamic *ri)
#else
int b_atan(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larctan"
{
   extern dynamic   Lone, LPiov4, Lres, err;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;

   int       rc;
/* Cordes : variable unused */
/* a_btyp  prec; */

#define  poolvars   0

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

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;


  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z == 1) {
      if ((rc = b_bini(ri)) == 0) {     /* ri = arctan(0) = 0   */
         ri->r = 0;                     /* Result exact         */
         EXIT(0);
         }
       else
         ERREXIT(rc,rc,0);              /* Error message and handling */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling   */
      }


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   if EQ_ABS_(xi,&Lone) {
      if (LPiov4.l < Maxl) {               /* ri = sgn(xi)*pi/4 */
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0); /* Error message and handling */
         }
/* Cordes : sign of xi assigned directly to ri */
/*    prec = xi->s;     */              /* save sign of xi      */
      COPY_(&LPiov4,ri);
      if (rc != 0) ERREXIT(rc,rc,0); /* Error message and handling */
/* Cordes : sign of xi assigned directly to ri */
/*    ri->s = prec;     */
      ri->s = xi->s;
      ri->r = 1;
      EXIT(0);
      }



  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();       /* Initialization of GLOBALS */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larctan 218) Maxl = %d\n",Maxl);
      printf("\n  (Larctan 219) xi   = ");  Lprinti(xi);
      }
#endif

   ARCTAN_(xi);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larctan 225) Lres = ");  Lprinti(&Lres);
      }
#endif

   if GT_ABS_(xi,&Lone) {       /* Result adaptation            */
      Maxl = Lcurrprec + Lguard;
      Lres.s = 0;
      SUB_(&LPiOv2,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }

   Lres.s = xi->s;                      /* Set sign of result   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larctan 239) Maxl = %d\n",Maxl);
      printf("\n  (Larctan 240) xi   = ");  Lprinti(xi);
      printf("\n  (Larctan 241) Lres = ");  Lprinti(&Lres);
      printf("\n  (Larctan 242) err  = ");  Lprinti(&err);
      printf("\n");
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,248,0);     /* Error message and handling   */
      }

#if INT_HPREC
   b_farg = xi;
   b_case = LESS_ABS_ARG;
#endif

   ASSIGN_(ri);         /* Assign result and set number of ulp's */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larctan 255) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larctan 256) LhE    = ");  Lprinti(&err);
      printf("\n  (Larctan 257) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}
#undef  Name
#undef  LRoutine
#undef  poolvars


#define  Name        "Larccot"

static char  *functionc = (Name);

/****************************************************************
 * Algorithm arccot(x)                                          *
 ****************************************************************/

#ifdef LINT_ARGS
int b_acot(dynamic *xi,dynamic *ri)
#else
int b_acot(xi,ri)

dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Larccot"
{
   extern dynamic   Lone, LPiov4, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;

   int       rc;
/*Cordes   a_btyp  prec;      */

#define  poolvars   0

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

   Lroutine = functionc;
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;


  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z == 1) {
      if (LPiov4.l <= Maxl) {           /* ri = arccot(0) = pi/2 */
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0); /* Error message and handling */
         }
      SHIFT_(1,&LPiov4,ri);
      if (rc != 0) ERREXIT(rc,rc,0);   /* Error message and handling */
      ri->s = 0;
      ri->r = 1;
      EXIT(0);
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling   */
      }


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   if EQ_ABS_(xi,&Lone) {
      if (LPiov4.l <= Maxl) {
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0);  /* Error message and handling   */
         }
/*Cordes    prec = xi->s;       */      /* save sign of xi(not used)    */
       if (xi->s == 0)
         COPY_(&LPiov4,ri);             /* ri = arccot(1) = pi/4        */
       else {
         SHIFT_(1,&LPiov4,&dummy);      /* ri = arccot(-1) = 3*pi/4     */
         ADD_(&LPiov4,&dummy,ri);
         }
      if (rc != 0) ERREXIT(rc,rc,0);    /* Error message and handling   */
      ri->s = 0;
      ri->r = 1;
      EXIT(0);
      }


  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larccot 410) Maxl = %d\n",Maxl);
      printf("\n  (Larccot 411) xi   = ");  Lprinti(xi);
      }
#endif

   ARCTAN_(xi);

   if LE_ABS_(xi,&Lone) {                  /* Result adaptation */
      if (LPiOv2.l < Lcurrprec+Lguard) {   /* ri = pi/2 - arctan(xi) */
         if (LPiov4.l <= Lcurrprec+Lguard) {
            Maxl = Lcurrprec + Lguard + PiMore;
            rc = LpiGen();
            if (rc != 0) ERREXIT(0,rc,0); /* Error message and handling */
            }
         Maxl = Lcurrprec + Lguard;
         SHIFT_(1,&LPiov4,&LPiOv2);
         if (rc != 0) ERREXIT(rc,rc,0); /* Error message and handling */
         }
       else
         Maxl = Lcurrprec + Lguard;

      SUB_(&LPiOv2,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }
    else if (xi->s == 1) {      /* ri = pi-arctan(sqrt(1-xi^2))  */
      if (LPiov4.l <= Lcurrprec+Lguard) {
         Maxl = Lcurrprec + Lguard + PiMore;
         rc = LpiGen();
         if (rc != 0) ERREXIT(0,rc,0); /* Error message and handling */
         }
      Maxl = Lcurrprec + Lguard;
      SHIFT_(2,&LPiov4,&dummy);
      if (rc != 0) ERREXIT(rc,rc,0); /* Error message and handling */

      Lres.s = 0;
      SUB_(&dummy,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }

   Lres.s = 0;                  /* Set sign of result   */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larccot 459) Maxl = %d\n",Maxl);
      printf("\n  (Larccot 460) xi   = ");  Lprinti(xi);
      printf("\n  (Larccot 461) Lres = ");  Lprinti(&Lres);
      printf("\n  (Larccot 462) err  = ");  Lprinti(&err);
      printf("\n");
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,468,0);     /* Error message and handling   */
      }

   ASSIGN_(ri);         /* Assign result and set number of ulp's*/

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Larccot 475) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Larccot 476) LhE    = ");  Lprinti(&err);
      printf("\n  (Larccot 477) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}
#undef  Name
#undef  LRoutine
#undef  poolvars


#define  Name        "Latan2"

static char  *function2 = (Name);

/****************************************************************
 * Algorithm atan2(x,y)                                         *
 ****************************************************************/

#ifdef LINT_ARGS
int b_atn2(dynamic *xi,dynamic *yi,dynamic *ri)
#else
int b_atn2(xi,yi,ri)

dynamic *xi;
dynamic *yi;
dynamic *ri;
#endif
#define  LRoutine    "Latan2"
{
   extern dynamic LPiov4, Lres, err, dummy;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine;

   dynamic   *y;
   int       rc;
/* Cordes : variable prec unused */
/*   a_btyp  prec;    */
   int  poolvars;       /* !!! must be int !!! */

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

   Lroutine = function2;
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;


  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=0 |
   *----------------------------*/

   if (xi->z != 0) {
      if (yi->z != 0)
         ERREXIT(RANGE,RANGE,0);        /* 0/0 undefined        */
      if (yi->m[0] == 0)
         ERREXIT(NANDE,NANDE,0);        /* Error message and handling */

      if (yi->s == 0)
         {
         if ((rc = b_bini(ri)) == 0) {  /* ri = arctan(0/|y|) = 0 */
            ri->r = 0;                  /* Result exact         */
            EXIT(0);
            }
          else
            ERREXIT(rc,rc,0);           /* Error message and handling */
         }

      if (LPiov4.l < Maxl) {
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0);  /* Error message and handling */
         }
      SHIFT_(2,&LPiov4,ri);             /* ri = arctan(0/(-|y|) = Pi  */
      if (rc != 0) ERREXIT(rc,rc,0);    /* Error message and handling */
      ri->r = 1;
      EXIT(0);
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling */
      }


  /*----------------------------*
   | Check and Handle Case yi=0 |
   *----------------------------*/

   if (yi->z != 0) {
      if (LPiov4.l < Maxl) {
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0);  /* Error message and handling */
         }
/* Cordes : sign of xi assigned directly to ri */
/*    prec = xi->s;     */          /* save sign of xi      */
      SHIFT_(1,&LPiov4,ri);         /* ri = sgn(y)*arcctg(0/|x|) = Pi/2 */
      if (rc != 0) ERREXIT(rc,rc,0); /* Error message and handling  */
      ri->r = 1;
/* Cordes : sign of xi assigned directly to ri */
/*    ri->s = prec;     */
      ri->s = xi->s;
      EXIT(0);
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (yi->m[0] == 0) {
      ERREXIT(NANDE,NANDE,0);   /* Error message and handling   */
      }


  /*-----------------------------*
   | Check and Handle Case xi=yi |
   *-----------------------------*/

   if EQ_ABS_(xi,yi) {
      if (LPiov4.l <= Maxl) {
         Lcurrprec = Maxl;
         Maxl += PiMore;
         rc = LpiGen();
         Maxl = Lcurrprec;
         if (rc != 0) ERREXIT(0,rc,0);  /* Error message and handling */
         }
/* Cordes : sign of xi assigned directly to ri */
/*    prec = xi->s;     */              /* save sign of xi      */
      if (yi->s == 0)
         COPY_(&LPiov4,ri);             /* ri=sgn(x)*arctan(1) =+-pi/4  */
       else {
         SHIFT_(1,&LPiov4,&dummy);      /* ri=sgn(x)*arccot(-1)=+-3*pi/4*/
         ADD_(&LPiov4,&dummy,ri);
         }
      if (rc != 0) ERREXIT(rc,rc,0);    /* Error message and handling   */
/* Cordes : sign of xi assigned directly to ri */
/*    ri->s = prec;     */
      ri->s = xi->s;
      ri->r = 1;
      EXIT(0);
      }



  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Latan2 672) Maxl = %d\n",Maxl);
      printf("\n  (Latan2 673) xi   = ");  Lprinti(xi);
      }
#endif

   y = Lvarget();
   poolvars = 1;

   Maxl = Lcurrprec + Lguard;

   if GE_ABS_(yi,xi)
      DIV_(xi,yi,y);
    else
      DIV_(yi,xi,y);

   Maxl = Lcurrprec;

   ARCTAN_(y);

   if GT_ABS_(xi,yi) {                  /* Result adaptation    */
      Maxl = Lcurrprec + Lguard;
      LPiOv2.s = xi->s;
      SUB_(&LPiOv2,&Lres,&Lres);
      LPiOv2.s = 0;

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }
    else if (yi->s != 0) {
      Maxl = Lcurrprec + Lguard;
      SHIFT_(2,&LPiov4,&dummy);
      dummy.s = xi->s;
      ADD_(&dummy,&Lres,&Lres);

      Ltwo.e = 1 - Maxl;
      ADD_(&err,&Ltwo,&err);
      NEXT_(&err,&err);
      }

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Latan2 716) Maxl = %d\n",Maxl);
      printf("\n  (Latan2 717) xi   = ");  Lprinti(xi);
      printf("\n  (Latan2 718) Lres = ");  Lprinti(&Lres);
      printf("\n  (Latan2 719) err  = ");  Lprinti(&err);
      printf("\n");
      }
#endif

   if (rc != 0) {
      ERREXIT(RESUL,725,0);     /* Error message and handling   */
      }

   ASSIGN_(ri);         /* Assign result and set number of ulp's*/

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Latan2 732) LhF    = ");  Lprinti(&Lres);
      printf("\n  (Latan2 733) LhE    = ");  Lprinti(&err);
      printf("\n  (Latan2 734) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   RETURN(0);
}





