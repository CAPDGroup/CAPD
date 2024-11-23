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

/* CVS $Id: b_gini.c,v 1.22 2014/01/30 17:24:04 cxsc Exp $ */

/************************************************************************
 *                                                                      *
 *   Global Constants and Variables common to all Standardfunctions     *
 *   ==============================================================     *
 *                                                                      *
 *   Initialization of constants and variables common to all standard   *
 *   function routines.                                                 *
 *                                                                      *
 *   The routine is called only once by the first standard function     *
 *   routine using one of the constants or variables defined in this    *
 *   routine.                                                           *
 *                                                                      *
 *   Future calls are avoided by setting the global flag variable       *
 *   Lgiflag.                                                           *
 *                                                                      *
 ************************************************************************/

#ifdef AIX
#include "/u/p88c/runtime/base/b_lari.h"
#else
#include "b_lari.h"
#endif

extern a_btyp  Maxl;

static a_btyp  mone[3] = { 1, 0, 0 };
static a_btyp  mmid[3] = { 4, 0, 0 };
static a_btyp  mmsq[3] = { Minsqrm, 0, 0 };
        
dynamic  LhF, LhE, LhD;
dynamic  LhV[numvar];
dynamic  Llnbase;
dynamic *baseptr = &Llnbase;
dynamic  Lone     = { 0, 0, 0, 0,       0, Minl, &mone[0] };
dynamic  Leps     = { 0, 0, 0, 0,       0, Minl, &mone[0] };
dynamic  Lmindbl  = { 0, 0, 0, 0,       0, Minl, &mmid[0] };
dynamic  LFminsq  = { 0, 0, 0, 0, Minsqre, Minl, &mmsq[0] };
dynamic  LMinreal = { 0, 0, 0, 0,-EFUFfac, Minl, &mone[0] };
        
a_real  FUFfac, FIUFfac, Fln2, Flnb;
        
a_btyp  Lcurrprec;
a_btyp  Lgiflag  = false;
a_btyp  LhI;
/* Uebergangsweise : */
a_btyp  LFulp    = 0;
/* */
a_btyp  EUFfac   = EFUFfac;
a_btyp  Lintern  = false;
a_btyp  Lversion = false;

/****************************************************************
 * Initial Values for PI, 4/PI and 2/PI                         *
 ****************************************************************/

#ifdef AIX
#include  "/u/p88c/runtime/base/b_lpi_.h"
#else
#include  "b_lpi_.h"
#endif

dynamic  LPiov4   = { 0, 0, 1, 0, -1, Lpi_init, &mPiov4[0] };
dynamic  L4ovPi   = { 0, 0, 1, 0,  0, Lpi_init, &m4ovPi[0] };

int   Ldebug    = 100;

#if INT_HPREC
int   b_case    = 0;
dynamic *b_farg = NULL;
#endif
        
rounding  Lrnd;

#ifdef LINT_ARGS
int b_gini(void)
#else
int b_gini()
#endif
#define  LRoutine    "Lginit"
{
   int       i;
   extern a_real *r_fln2,*r_flnb;

#ifdef Debug
   Ldebug -= 1;                 /* Diminish Debug Level      */
   if (Ldebug >= 0) printf("\n Entering Routine %s",LRoutine);
#endif

   Lgiflag = true;
        
   (void)b_bini(&LhF);
   (void)b_bini(&LhE);
   (void)b_bini(&LhD);

   for (i=0; i<numvar; i++) (void)b_bini(&LhV[i]);

   (void)b_bini(&Llnbase);

   R_ASSIGN(Fln2,*r_fln2);
   R_ASSIGN(Flnb,*r_flnb);

   (void)b_bcdi(Flnb,&baseptr,(a_intg)0);

   ((a_btyp *)&FUFfac)[B_HPART] = MFUFfac;
   ((a_btyp *)&FUFfac)[B_LPART] = ZERO;
   ((a_btyp *)&FIUFfac)[B_HPART] = MFIUFfac;
   ((a_btyp *)&FIUFfac)[B_LPART] = ZERO;
        
/* FMinsqr.u[0] = Minsqrd;  FMinsqr.u[1] = 0; */

   EXIT(0);
}
        
        





