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

/* CVS $Id: r_comp.c,v 1.21 2014/01/30 17:24:11 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_comp.c                              */
/*                                                              */
/*      Entries         : a_real r_comp(m,e)                    */
/*                        a_real m;                             */
/*                        a_intg e;                             */
/*                                                              */
/*      Arguments       : m = mantissa                          */
/*                        e = exponent                          */
/*                                                              */
/*      Description     : Composition of an IEEE value.         */
/*                                                              */
/*      Note            : mantissa is 0.0 or in the range       */
/*                           0.5<=|mantissa|<1.0                */
/*                        sNaN causes INV_OP exception.         */
/*                        qNaN is returned unaltered.           */
/*                        infinity is returned unaltered.       */
/*                                                              */
/*                   argument of e_trap must be address         */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_minf;
extern a_real *r_pinf;
#endif

#ifdef LINT_ARGS
local a_real r_comp(a_real m,a_intg e)
#else
local a_real r_comp(m,e)

a_real   m;
a_intg e;
#endif
        {
        a_bool vz;
        a_intg expo;
        a_btyp mant[D_U_RATIO];

        E_TPUSH("r_comp")

        /* mantissa is zero */
        if (b_deko(m,&expo,mant,&vz))
           e = -CHARAC;

        /* mantissa out of range if expo not equal to -1*/
        else if (expo!=-1)
           {
           if (expo>EXPO_MAX)
              {
              if (MANT_INFINITY(mant))
                 e_trap(NO_ERROR+E_EMSG+E_ECNT,6,E_TMSG,13,
                        E_TDBL+E_TEXT(18),&m,E_TINT+E_TEXT(19),&e);
              else if (SIGNALING(mant[0]))
                 e_trap(INV_OP+E_IEEE,6,E_TMSG,5,E_TDBL+E_TEXT(18),&m,
                        E_TINT+E_TEXT(19),&e);
              else
                 e_trap(NO_ERROR+E_EMSG+E_ECNT,6,E_TMSG,14,
                        E_TDBL+E_TEXT(18),&m,E_TINT+E_TEXT(19),&e);
              }
           else
              e_trap(INV_ARG,6,E_TMSG,47,E_TDBL+E_TEXT(18),&m,
                     E_TINT+E_TEXT(19),&e);

           E_TPOPP("r_comp")
           return(m);
           }

        else
           {
           /* adjust exponent e */
           e--;

           /* exponent out of range */
           if (e>EXPO_MAX)
              {
              e_trap(OVERFLOW,6,E_TMSG,48,E_TDBL+E_TEXT(18),&m,
                     E_TINT+E_TEXT(19),&e);
              E_TPOPP("r_comp")
              return((vz) ? *r_minf : *r_pinf);
              }

           /* loss of accuracy due to small exponent */
           else if (e<EXPO_MIN-MANTL+1)
              {
              e_trap(INEXACT,6,E_TMSG,50,E_TDBL+E_TEXT(18),&m,
                     E_TINT+E_TEXT(19),&e);
              e = -CHARAC;
              }

           /* generate denormalized number                         */
           else if (e<EXPO_MIN)
              {

              /* expo used to indicate that message was displayed once */
              /* and for message displaying by e_trap() */
              expo = e;
              do {
                 if ((mant[D_U_RATIO-1] & LSB) && expo!=0)
                    {
                    e_trap(INEXACT,6,E_TMSG,49,E_TDBL+E_TEXT(18),&m,
                           E_TINT+E_TEXT(19),&expo);
                    expo = 0;
                    }
                 b_shr1(mant,D_U_RATIO);
                 }
              while (++e<EXPO_MIN);
              }
           }

        b_comp(&m,e,mant,vz);

        E_TPOPP("r_comp")
        return(m);
        }





