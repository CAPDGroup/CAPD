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

/* CVS $Id: b_bcdi.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bcdi.c          Processor : C                   */
/*                                                                      */
/* Conversion of an IEEE double variable to intern representation.      */
/* The conversion is done without rounding.                             */
/*                                                                      */
/* Function value : int     0 - normalized number converted             */
/*                      ALLOC - allocation error                        */
/*                      DENOR - denormalized number converted           */
/*                      PINFI - +infinity     (no conversion!!!)        */
/*                      MINFI - -infinity     (no conversion!!!)        */
/*                      NANDE - NAN           (no conversion!!!)        */
/*                                                                      */
/************************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
#endif

#ifdef LINT_ARGS
local int b_bcdi(a_real d,multiprecision *i,a_intg rnd)
#else
local int b_bcdi(d,i,rnd)

a_real d;                 /* value of double variable             */
multiprecision *i;       /* pointer to intern variable           */
a_intg rnd;
#endif
        {
        a_btyp mant[D_U_RATIO+1];
        a_bool vz;
        a_intg expo;
        a_intg c;
        a_intg k,l;
/*----------------------------------------------------------------------*/
        (*i)->r = 0;    /* operation is exact                           */
/*----------------------------------------------------------------------*/
        if (b_deko(d,&expo,mant,&vz))
           {
           (*i)->z = 1;
           if ((*i)->l)
              {
              (*i)->l = 0;
#ifdef HEAP_CHECK
b_freh((a_char *)&(*i)->m,(a_char *)(*i)->m,(a_char *)"b_bcdi");
#endif
              B_FREE((*i)->m)
              }
           return(0);
           }            /* double value is zero                         */
/*----------------------------------------------------------------------*/
        (*i)->z = 0;    /* number value is non-zero                     */
/*----------------------------------------------------------------------*/
        (*i)->s = vz;   /* determine sign of double value               */
/*----------------------------------------------------------------------*/
        if (expo>EXPO_MAX)
           {
           if (MANT_INFINITY(mant)) return( (vz) ? MINFI : PINFI );
           return(NANDE);
           }            /* double variable holds a NAN or infinity      */
/*----------------------------------------------------------------------*/
        (*i)->e = B_ASHR(expo,LOG_B_LENGTH);
                           /* exponent of (*i)->m[0]                   */
        c = ((*i)->e<<LOG_B_LENGTH)-expo+(B_LENGTH-ZERO_BITS);
                           /* shift value                               */
/*----------------------------------------------------------------------*/
        if (c<0)
           b_shlu(mant,(l = D_U_RATIO),-c);
        else if (c>0)
           {
           mant[D_U_RATIO] = ZERO;
           b_shru(mant,(l = D_U_RATIO+1),c);
           }
        else
           l = D_U_RATIO;
                                /* shift mantissa */
/*----------------------------------------------------------------------*/
        for (k=0;mant[k]==ZERO;k++) (*i)->e--;
/*----------------------------------------------------------------------*/
        while (mant[l-1]==ZERO) l--;
/*----------------------------------------------------------------------*/
        b_badj(l-k,*i);   /* adjust intern variable                      */
/*----------------------------------------------------------------------*/
        for (l=0;l<(*i)->l;l++) (*i)->m[l] = mant[l+k];
                         /* copy mantissa digits                        */
/*----------------------------------------------------------------------*/
        return(0);
        }





