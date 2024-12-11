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

/* CVS $Id: b_pool.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/*************************************************************************
 *                                                                       *
 * Descriptive Name : Lvarpool           Processor : C                   *
 *                                                                       *
 * Pool for Variables of type dynamic                                    *
 * ==================================                                    *
 *                                                                       *
 * Include files :  b_lari.h  - definitions for INTERN standard functions*
 *                                                                       *
 * Function value : int       - 0                                        *
 *                                                                       *
 * Used Number Base :  2**32                                             *
 *                                                                       *
 * Used Functions :                                                      *
 *                     Intern arithmetic                                 *
 *                                                                       *
 * Used Global INTERN Variables and Constants :                          *
 *    LhV[10]                                                            *
 *                                                                       *
 * Used UNSIGNED Global Variables :                                      *
 *    Ldebug          Flag for printing additional information           *
 *                                                                       *
 *************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

static dynamic   LhV[numvar];

#ifdef Debug
extern int       Ldebug;
#endif

/****************************************************************
 * Constants for Routine Lvarpool                               *
 ****************************************************************/

static int next = 0;

/****************************************************************
 * Algorithm for Fetching a new Variable                        *
 ****************************************************************/

#ifdef LINT_ARGS
dynamic *b_get_(void)
#else
dynamic *b_get_()
#endif
#define  LRoutine    "Lvarget"
{
#ifdef Debug
   Ldebug -= 1;                            /* Diminish Debug Level */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lvarget 78) assigning : %d\n",next);
      }
#endif

   if (next >= numvar-1) {
      Lerror(ALLOC);            /* Error Message                */
      }
   else {
      next++;
      EXIT(&(LhV[next])); /* Return address of next free variable */
      }
   EXIT(NULL);               /* Return NULL if no variable available */
}
#undef  LRoutine


/****************************************************************
 * Algorithm for Returning Variables                            *
 ****************************************************************/

#ifdef LINT_ARGS
int b_drop(int  n)
#else
int b_drop(n)
int  n;
#endif
#define  LRoutine    "Lvardrop"
{
#ifdef Debug
   Ldebug -= 1;                            /* Diminish Debug Level */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

   next = max(next-n,0);        /*   avoiding negative indices      */

#ifdef Debug
   if (Ldebug >= 0) {
      printf("\n  (Lvardrop 118) next = %d\n",next);
      }
#endif

   EXIT(0);
}





