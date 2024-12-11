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

/* CVS $Id: r_clss.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_clss.c                              */
/*                                                              */
/*      Entries         : a_intg r_clss(x)                      */
/*                                                              */
/*      Arguments       : real x;                               */
/*                                                              */
/*      Function value  : E_CLS0 - signaling NaN                */
/*                        E_CLS1 - quiet NaN                    */
/*                        E_CLS2 - -infinity                    */
/*                        E_CLS3 - negative normalized nonzero  */
/*                        E_CLS4 - negative denormalized        */
/*                        E_CLS5 - -0                           */
/*                        E_CLS6 - +0                           */
/*                        E_CLS7 - positive denormalized        */
/*                        E_CLS8 - positive normalized nonzero  */
/*                        E_CLS9 - +infinity                    */
/*                                                              */
/*      Description     : determine the class of an IEEE        */
/*                        coded real value.                     */
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
local a_intg r_clss(a_real x)
#else
local a_intg r_clss(x)

a_real x;
#endif
        {
        a_btyp mant[D_U_RATIO];
        a_intg expo,res;
        a_bool vz;

        E_TPUSH("r_clss")

        if (b_deko(x,&expo,mant,&vz))
           res = (vz) ? E_CLS5 : E_CLS6;
        else if (expo==EXPO_MAX+1)
           {
           if (MANT_INFINITY(mant))
              res = (vz) ? E_CLS2 : E_CLS9;
           else
              res = (SIGNALING(mant[0])) ? E_CLS0 : E_CLS1;
           }
        else if (expo==EXPO_MIN-1)
           res = (vz) ? E_CLS4 : E_CLS7;
        else
           res = (vz) ? E_CLS3 : E_CLS8;

        E_TPOPP("r_clss")
        return(res);
        }





