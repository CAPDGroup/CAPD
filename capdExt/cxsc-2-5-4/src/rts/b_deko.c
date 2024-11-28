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

/* CVS $Id: b_deko.c,v 1.21 2014/01/30 17:24:04 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_deko.c                              */
/*                                                              */
/*      Entries         : a_bool b_deko(x,expo,mant,vz)         */
/*                        a_real x;                             */
/*                        a_btyp *mant;                         */
/*                        a_intg *expo;                         */
/*                        a_bool *vz;                           */
/*                                                              */
/*      Description     : Decomposition of an IEEE double value */
/*                        x into exponent, mantissa and sign    */
/*                                                              */
/*      Function value  : FALSE if x!=0                         */
/*                        TRUE if x==0                          */
/*                                                              */
/*                   statement sequence changed                 */
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
local a_bool b_deko(a_real x,a_intg *expo,a_btyp *mant,a_bool *vz)
#else
local a_bool b_deko(x,expo,mant,vz)

a_real x;
a_intg *expo;
a_btyp *mant;
a_bool *vz;
#endif
        {
        register a_btyp *u;

        u = (a_btyp *)&x;

        /* low part of mantissa */
        mant[1] = u[B_LPART];

        /* sign of real value   */
        *vz = (u[B_HPART] & MSB) ? TRUE : FALSE;

        /* high part of mantissa */
        mant[0] = (u[B_HPART] & MANT_HIGH) | HIDDEN_BIT;

        /* exponent of number   */
        if ((*expo = ((u[B_HPART] & EXPO_MASK)>>EXPO_SHIFT)-CHARAC)
            ==-CHARAC)
           {

           /* value is zero     */
           if (!((mant[0] &= NOT_HIDDEN_BIT) | mant[1]))
              return(TRUE);

           /* denormalized number  */
           *expo = -CHARAC+1;
           }

        return(FALSE);
        }





