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

/* CVS $Id: b_errr.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lerror             Processor : C                   *
 *                                                                       *
 * Routine for Printing Error Message and Handling Errors                *
 * ======================================================                *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Functions :                                                      *
 *    none                                                               *
 *                                                                       *
 * Used Global UNSIGNED Variables :                                      *
 *    Lgiflag         Flag for Initialization of Global INTERN Variables *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 * Used Global STRING Variables :                                        *
 *    Lroutine        Name of the Routine called by user                 *
 *                                                                       *
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

extern f_text f_errr;

#define  printe(x)     fprintf(f_errr.fp,x)
#define  printe2(x,y)  fprintf(f_errr.fp,x,y)

char  *Lroutine, *Lerrmsg;

static char  Unknown[] = "(unknown)";

#ifdef LINT_ARGS
int b_errr(a_btyp err)
#else
int b_errr(err)

a_btyp err;
#endif
#define  LRoutine    "Lerror"
{
   int   rc;

#ifdef Debug
   extern int       Ldebug;
#endif


   if (err == 0) return(0);

  /*--------------------------------------*
   | Print Version and modify Debug Level |
   *--------------------------------------*/

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level      */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif


  /*------------------------------------*
   | Check and Handle Undefined Routine |
   *------------------------------------*/

   if (Lroutine == NULL) {      /* Name of user called routine   */
      Lerrmsg = Unknown;        /*    is not known               */
      rc = NANDE;
      }
    else {
      rc = 0;
      }


  /*---------------------------*
   | Print appropriate Message |
   *---------------------------*/

   printe2("\n ***ERROR*** in Routine \"%s\" : ",Lroutine);
   switch ((int)err) {
      case DENOR   :  printe("Denormalized number converted");
                      break;
      case MINFI   :  printe("Minus infinity detected");
                      break;
      case NANDE   :  printe("NAN detected");
                      break;
      case OFLOW   :  printe("Exponent overflow");
                      break;
      case PINFI   :  printe("Plus infinity detected");
                      break;
      case ROUND   :  printe("Double value is rounded");
                      break;
      case UFLOW   :  printe("Exponent underflow");
                      break;
      case ZEROD   :  printe("Division by zero");
                      break;
      case RANGE   :  printe("Range error");
                      break;
      case ALLOC   :  printe("Allocation error");
                      break;
      case NALLO   :  printe("Data not allocated");
                      break;

      case ASSGN   :  printe("Assignment of result failed");
                      break;
      case ERRBD   :  printe("Determined error bound is invalid");
                      break;

      case PEVAL   :  printe("Error during polynomial evaluation");
                      break;
      case RESUL   :  printe("Error during result adaptation");
                      break;
      case CONVD   :  printe("Error during conversion to double");
                      break;
      case EPERR   :  printe("Error during computation of error bound");
                      break;
      case UNITS   :  printe("Error during computation of number of ulp's");
                      break;
      case DUFLW   :  printe("Floating point underflow during computation");
                      break;
      default      :  printe2("Return Code : %d",rc);
      }

   printe("\n");
   Lerrmsg = NULL;

   EXIT(rc);
}





