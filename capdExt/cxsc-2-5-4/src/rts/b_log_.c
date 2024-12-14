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

/* CVS $Id: b_log_.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lln                Processor : C                   *
 *                                                                       *
 * Natural Logarithm for Multiple Precision Arithmetic                   *
 * ===================================================                   *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Argument Range :    All positive numbers                              *
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

#define  Name        "Lln"
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

#define Lguard    2

#define  rerr        "Argument 0 or negative"

static a_btyp  mbapp[3]  = { 0x18000000L,0x00000000L,0x00000000L };

static dynamic   bdapp   = { 0, 0, 0, 0, -1, Minl, &mbapp[0] };

static char  *function = Name;
static char  *Rerr     = rerr;

/****************************************************************
 * Algorithm                                                    *
 ****************************************************************/

#ifdef LINT_ARGS
int b_log_(dynamic *xi,dynamic *ri)
#else
int b_log_(xi,ri)
dynamic *xi;
dynamic *ri;
#endif
#define  LRoutine    "Lln"
{
   extern dynamic   dummy, Lone;
   extern a_btyp  Maxl;
   extern a_btyp  Lcurrprec, Lgiflag;
   extern char      *Lroutine, *Lerrmsg;

#ifdef Debug
   extern int       Ldebug;
   extern dynamic   Lres, err;
#endif

   dynamic    *y;
   int        rc, poolvars, ysgn;

  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level         */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /**************************************************************
   * Check and Handle Range Error and Special Cases             *
   **************************************************************/


  /*----------------------------*
   | Check and Handle Case xi=1 |
   *----------------------------*/

   rc = 0;
   if (compii(xi,&Lone) == 0) {
      ri->z = 1;                /* ri = 0                       */
      EXIT(rc);
      }

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting       */


  /*-----------------------------*
   | Check and Handle Case xi<=0 |
   *-----------------------------*/

   if (xi->s == 1 || xi->z == 1) {
      Lerrmsg = Rerr;       /* Error message for non positive argument  */
      ERREXIT(RANGE,153,0); /* Error message and handling               */
      }


  /*-------------------------------------*
   | Check and Handle Undefined Mantissa |
   *-------------------------------------*/

   if (xi->m[0] == 0) {
      ERREXIT(NANDE,162,0);     /* Error message and handling   */
      }


  /*-----------------------------------*
   | Check and Handle Different Ranges |
   *-----------------------------------*/

   y = Lvarget();
   poolvars = 1;

   rc = 0;
   SUB_(xi,&Lone,y);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lln 178) Larg-1 = ");  Lprinti(y);
      }
#endif

   ysgn = y->s;
   y->s = 0;

   if LE_(y,&bdapp) {      /* Approxmation if |y| > 0.09375 */

     /*------------------------------------------*
      | Initializations for Taylor Approximation |
      *------------------------------------------*/

      y->s = ysgn;              /* y = (1 - Larg) / (1 + Larg)  */
      Maxl = min(Lcurrprec,xi->l) + 1;
      ADD_(xi,&Lone,&dummy);    /* EXACT sum                    */

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lln 197) Maxl = %d",Maxl);
         printf("\n  (Lln 198) Larg+1 = ");  Lprinti(&dummy);
         }
#endif

      Maxl = Lcurrprec + Lguard;
      DIV_(y,&dummy,y);

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lln 207) Maxl = %d",Maxl);
         printf("\n  (Lln 208) (Larg-1)/(Larg+1) = ");  Lprinti(y);
         }
#endif

      if (rc != 0) ERREXIT(PEVAL,212,poolvars);

      LN_T_(y);         /* Compute function value and error bound   */
      }
    else
      {
      LN_(xi);          /* Compute function value and error bound   */
      }

   ASSIGN_(ri);         /* Assign result and set number of ulp's    */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lln 223) Lres   = ");  Lprinti(&Lres);
      printf("\n  (Lln 224) err    = ");  Lprinti(&err);
      printf("\n  (Lln 225) Result = ");  Lprinti(ri);
      printf("   (+ %d ulp)",ri->r);
      }
#endif

   Lvardrop(poolvars);
   RETURN(0);
}





