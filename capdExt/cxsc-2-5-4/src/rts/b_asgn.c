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

/* CVS $Id: b_asgn.c,v 1.22 2014/01/30 17:24:02 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lassign            Processor : C                   *
 *                                                                       *
 * Routine for Computing the Number of ulp's                             *
 * =========================================                             *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Used Functions :                                                      *
 *    Lerror          Error handling routine                             *
 *                    Intern arithmetic                                  *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhE, LhF, LhD                                                      *
 *    Lone              ( = 1 )                            (Constant)    *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Maxl            Length of INTERN Variables                         *
 *    Lcurrprec       Initial Length of INTERN Variables                 *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

/************************************************************************
 * Redefinition of Names of Intermediate Variables Used                 *
 ************************************************************************/

#define val       LhF           /* Function value                       */
#define err       LhE           /* Rounding error                       */
#define dummy     LhD           /* Variable for temporary use only      */

/************************************************************************
 * Constants for Routine Lassign                                        *
 ************************************************************************/

#define errsgn    ulp           /* Sign of error                        */
#define Currprec  Lcurrprec /* Maxl=Lcurrprec in contrary to not "Main" */
#define poolvars  1             /* Number of pool variables used        */

/************************************************************************
 * Algorithm                                                            *
 ************************************************************************/

#ifdef LINT_ARGS
int b_asgn(dynamic *res)
#else
int b_asgn(res)

dynamic *res;
#endif
#define  LRoutine    "Lassign"
{
   extern dynamic   Lone, val, err, dummy;
   extern a_btyp    Maxl;
   extern a_btyp    Lcurrprec;
   extern rounding  Lrnd;

          dynamic   *ub;
          int       sign, ulp, rc=0;

#ifdef Debug
   extern int       Ldebug;
#endif

#ifdef Debug
   Ldebug -= 1;                            /* Diminish Debug Level      */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lassign 102) Maxl = %d,   Lcurrprec = %d",Maxl,Lcurrprec);
      printf("\n  (Lassign 103) val = ");  Lprinti(&val);
      printf("\n  (Lassign 104) err = ");  Lprinti(&err);
      printf(",  err.z = %d",err.z);
      printf("\n  (Lassign 106)     rounding = %c",Lrnd);
      }
#endif

   if ((err.z == 1) || (Lrnd == exact)) {
      rc = copyii(&val,res);            /* Copy value                   */

      res->r = val.r;                   /* copy rounding bit explicitly */

      if (rc != 0) {
         Lerror(rc);                    /* Error message and handling   */
         }
      RETURN(rc);
      }

   errsgn = err.s;

   if (Lrnd != absolute) {                  /* Relative error           */
      Maxl = Minlerr;
      rc = 0;
      MUL_(&err,&val,&err,&dummy);          /* Compute absolute error   */
      if (dummy.z == 1)                     /*    err <- err*val        */
         NEXT_(&err,&err);
      }

   switch (Lrnd) {
      case signed  :  err.s = errsgn;
                      break;
      case rounded :
      case ibound  :  err.s = val.s;
                      break;
      case obound  :  err.s = ~val.s;
                      break;
      case lbound  :  err.s = 0;
                      break;
      case ubound  :  err.s = 1;
                      break;
      }

   Maxl = Lcurrprec;
   ub = Lvarget();

   if (Lrnd == rounded) {                   /* Sign of error not known  */
      SUB_(&val,&err,res);                  /* Inner bound = Result     */
      ADD_(&val,&err,ub);                   /* Outer bound              */
      if (ub->r != 0)
         NEXT_(ub,ub);
      }
    else {                                  /* Sign of error known      */
      if (err.s == val.s) {
         COPY_(&val,res);                   /* Inner bound = Result     */
         ADD_(&val,&err,ub);                /* Outer bound              */
         if (ub->r != 0)
            NEXT_(ub,ub);
         }
       else {
         ADD_(&val,&err,res);               /* Inner bound = Result     */
         COPY_(&val,ub);                    /* Outer bound              */
         if (ub->r != 0)
            NEXT_(ub,ub);
         }
      }

#if INT_HPREC
   if (b_case)
      {
      if (b_case & LESS_ABS_ARG)
         {
         if (GT_ABS_(ub,b_farg))
            {
            sign = ub->s;
            COPY_(b_farg,ub);
            ub->s = sign;
            }
         }
      else if (b_case & GREATER_ABS_ARG)
         {
         if (LT_ABS_(res,b_farg))
            {
            sign = res->s;
            COPY_(b_farg,res);
            res->s = sign;
            }
         }

      if (b_case & LESS_ABS_ONE)
         {
         if (GT_ABS_(ub,&Lone))
            {
            sign = ub->s;
            COPY_(&Lone,ub);
            ub->s = sign;
            }
         }
      else if (b_case & GREATER_ABS_ONE)
         {
         if (LT_ABS_(res,&Lone))
            {
            sign = res->s;
            COPY_(&Lone,res);
            res->s = sign;
            }
         }
      b_case = 0;
      }
#endif

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lassign 168) res = ");  Lprinti(res);
      printf("\n  (Lassign 169) ub  = ");  Lprinti(ub);
      }
#endif

   Maxl = Minlerr;
   SUB_(ub,res,&dummy);                     /* Get number of units      */
   if (dummy.r != 0) NEXT_(&dummy,&dummy);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lassign 179) ub->e = %d,   ub->l = %d",ub->e,ub->l);
      }
#endif

   dummy.e -= (ub->e - ub->l + 1);
   dummy.s  = 0;

   Lvardrop(poolvars);

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lassign 190) diff = ");  Lprinti(&dummy);
      }
#endif

   for(ulp=0; GE_(&dummy,&Lone) && ulp<4 ; ulp++)
      SUB_(&dummy,&Lone,&dummy);
   if (dummy.z == 0) ulp++;

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lassign 200) ulp = %d",ulp);
      }
#endif

   if (rc != 0) ERREXIT(ASSGN,rc,0);       /* Error message and handling */

   if (ulp > 2) {
      res->r = 3;                          /* Mark error bound invalid   */
      ERREXIT(ERRBD,UNITS,0);              /* Error message and handling */
      }
    else
      res->r = ulp;                        /* Set number of units (ulp's)*/

   RETURN(rc);
}





