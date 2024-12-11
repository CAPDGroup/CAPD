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

/* CVS $Id: b_pi__.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lpi                Processor : C                   *
 *                                                                       *
 * Function returning PI with current acurracy                           *
 * ===========================================                           *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Used Functions :                                                      *
 *    LpiGen          Internal Routine to compute more digits of PI      *
 *    Lginit          Initialization of global constants                 *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#define  Name        "Lpi"
#define  Main

#define  PiMore      20

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

static char  *function = Name;

/************************************************************
 * Algorithm                                                *
 ************************************************************/

#ifdef LINT_ARGS
int b_pi__(dynamic *pi)
#else
int b_pi__(pi)

dynamic *pi;
#endif
#define  LRoutine    "Lpi"
{
   extern a_btyp  Maxl;
   extern a_btyp  Lgiflag;
   extern char      *Lroutine;
   extern dynamic   LPiov4;

#ifdef Debug
   extern int       Ldebug;
#endif

   int       rc;
   a_btyp  Lcurrprec;  /* Necessary while Lsqrt is used ! */
#define poolvars  0

#ifdef Debug
   Ldebug -= 1;         /* Diminish Debug Level         */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*----------------*
   | Initialization |
   *----------------*/

   Lroutine = function;
   Lcurrprec = Maxl;            /* Save precision setting       */
   rc = 0;


  /**************************************************************
   * Initialization of Global Constants and Storage Locations   *
   **************************************************************/

   if (Lgiflag==0) Lginit();    /* Initialization of GLOBALS     */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpi 123) Maxl = %d\n",Maxl);
      }
#endif

   if (Maxl > LPiov4.l) {
/*    Maxl = 2*Maxl; */
      Maxl += PiMore;
      rc = LpiGen();
      Maxl = Lcurrprec;

#ifdef Debug
      if (Ldebug >= 0) {
         printf("\n  (Lpi 135) Maxl = %d\n",Maxl);
         }
#endif

      if (rc != 0) {
         ERREXIT(0,rc,0);       /* Error message and handling    */
         }
      }


   SHIFT_(2,&LPiov4,pi);
   if (rc != 0) {
      ERREXIT(rc,rc,0);         /* Error message and handling    */
      }
   pi->r = 1;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lpi 153) Pi = ");  Lprinti(pi);
      }
#endif

   RETURN(0);
}





