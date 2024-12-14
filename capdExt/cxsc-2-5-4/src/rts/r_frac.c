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

/* CVS $Id: r_frac.c,v 1.21 2014/01/30 17:24:12 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : r_frac.c                              */
/*                                                              */
/*      Entries         : a_real r_frac(a)                      */
/*                        a_real a;                             */
/*                                                              */
/*      Arguments       : a = IEEE value                        */
/*                                                              */
/*      Description     : Fraction part of IEEE number.         */
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
local a_real r_frac(a_real a)
#else
local a_real r_frac(a)

a_real a;
#endif
        {
        a_btyp expo;
        a_btyp *ptra;

        E_TPUSH("r_frac")

        /* special value */
        if ((expo = (ptra = ((a_btyp *)&a))[B_HPART] & EXPO_MASK)==EXPO_MASK)
           {
           if (((ptra[B_HPART] & EXPO_MASK) | ptra[B_LPART])!=0)
              {
              if (SIGNALING(ptra[B_HPART]))
                 {
                 e_trap(INV_OP+E_IEEE,4,E_TMSG,5,E_TDBL+E_TEXT(7),ptra);
                 ptra[B_HPART] = SET_QUIET(EXPO_MASK);
                 ptra[B_LPART] = INV_OP;
                 }
              }
           }

        /* number has fraction part */
        else if (expo>=0x3ff00000)
           {
           a_btyp a0,a1;

           /* number does not posses fraction part since >=2**52 */
           if (expo>=0x43300000) 
              {
              a0 = ptra[B_HPART] & MSB;
              expo = a1 = 0;
              }

           /* determine fraction */
           else
              {
              int s;
   
              a0 = ptra[B_HPART] & MANT_HIGH;
              a1 = ptra[B_LPART];
   
              /* clear all integer bits by shifting */
              s = ((int)(expo>>EXPO_SHIFT))-CHARAC;
              if (s>=B_LENGTH) 
                 {
                 a0 = (a1<<(s-B_LENGTH)) & MANT_HIGH;
                 a1 = 0;
                 }
              else if (s!=0)
                 {
                 a0 = ((a0<<s) & MANT_HIGH) | (a1>>(B_LENGTH-s));
                 a1 <<= s;
                 }
   
              /* determine exponent if non-zero fraction part */
              if ((a0|a1)!=0)
                 {
                 expo = 0x3ff00000;
                 while ((a0 & HIDDEN_BIT)==0)
                    {
                    expo -= HIDDEN_BIT;
                    a0 = (a0<<1) | (a1>>(B_LENGTH-1));
                    a1 <<= 1;
                    }
                 }
              else
                 {
                 expo = 0;
                 }
              }
   
           /* construct number */
           ptra[B_HPART] = (ptra[B_HPART] & MSB) | expo | (a0 & MANT_HIGH);
           ptra[B_LPART] = a1;
           }

        E_TPOPP("r_frac")
        return(a);
        }





