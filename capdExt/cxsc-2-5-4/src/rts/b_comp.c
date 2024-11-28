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

/* CVS $Id: b_comp.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_comp.c                              */
/*                                                              */
/*      Entries         : void b_comp(x,expo,mant,vz)           */
/*                        a_real *x;                            */
/*                        a_intg expo;                          */
/*                        a_btyp *mant;                         */
/*                        a_bool vz;                            */
/*                                                              */
/*      Arguments       : x     - IEEE double format real       */
/*                        expo  - base 2 exponent               */
/*                        mant  - mantissa bits (a_btyp[2])     */
/*                        vz    - sign of number                */
/*                                                              */
/*      Description     : Composition of an IEEE double value   */
/*                        x from exponent, mantissa and sign    */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local void b_comp(a_real *x,a_intg expo,a_btyp *mant,a_bool  vz)
#else
local void b_comp(x,expo,mant,vz)

a_real     *x;
a_intg   expo;
a_btyp *mant;
a_bool  vz;
#endif
        {
        a_btyp *u;
#if C_P_3
#else
        size_t i;
#endif

        u = (a_btyp *)x;

        /* copy mantissa                                        */
#if C_P_3
        u[B_HPART] = mant[0]; u[B_LPART] = mant[1];
#else
        for (i=0;i<D_U_RATIO;i++)
            u[B_HPART-((B_HPART-B_LPART)/(D_U_RATIO-1))*i] = mant[i];
#endif

        /* denormalized number                                  */
        if (expo==EXPO_MIN && (mant[0] & HIDDEN_BIT)==ZERO)
           expo = -CHARAC;

        /* normalized number                                    */
        else
           u[B_HPART] &= NOT_HIDDEN_BIT;

        u[B_HPART] |= ((a_btyp)(expo+CHARAC))<<EXPO_SHIFT;

        /* generate sign                                        */
        if (vz)
           u[B_HPART] |= MSB;

        return;
        }





