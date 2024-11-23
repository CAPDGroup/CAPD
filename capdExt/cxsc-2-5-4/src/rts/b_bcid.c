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

/* CVS $Id: b_bcid.c,v 1.21 2014/01/30 17:24:03 cxsc Exp $ */

/************************************************************************/
/*                                                                      */
/* Descriptive Name : b_bcid.c          Processor : C                   */
/*                                                                      */
/* Conversion of intern representation to double variable               */
/*                                                                      */
/* Function value : int     0 - value is exactly representable          */
/*                      NALLO - data not allocated                      */
/*                      ROUND - value is rounded towards zero           */
/*                      UFLOW - underflow (result is zero)              */
/*                      OFLOW - overflow  (result is NAN)               */
/*                                                                      */
/* Rounding : towards zero (rnd argument ignored)                       */
/*                                                                      */
/************************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_zero;
#endif

#ifdef LINT_ARGS
local int b_bcid(multiprecision i,a_real *d,a_intg rnd)
#else
local int b_bcid(i,d,rnd)

multiprecision i;       /* pointer to intern variable           */
a_real *d;              /* pointer to double variable           */
a_intg rnd;
#endif
        {
        a_btyp mant[D_U_RATIO];
                      /* mantissa of IEEE value                 */
        a_intg c;     /* shift value determined from mantissa of*/
                      /* intern variable                        */
        a_intg k;     /* loop variable                          */
        a_intg expo;  /* exponent of IEEE value                 */
        a_bool vz;    /* sign of IEEE value                     */
        a_btyp uh;    /* unsigned value for double construction */
/*--------------------------------------------------------------*/
        R_ASSIGN(*d,*r_zero);
                        /* set double variable to zero          */
/*--------------------------------------------------------------*/
        if (i->z)       /* intern number is zero                */
           {
           return(0);
           }
/*--------------------------------------------------------------*/
        vz = i->s;
/*--------------------------------------------------------------*/
        if ((i->e*B_LENGTH>EXPO_MAX+1) ||
            ((i->e*B_LENGTH==EXPO_MAX+1) && (i->m[0]!=LSB)))
           {
           expo = EXPO_MAX;
           mant[0] = HIDDEN_BIT | MANT_HIGH;
           for (k=1;k<D_U_RATIO;k++) mant[k] = MAX_BASETYPE;
           b_comp(d,expo,mant,vz);
           return(OFLOW);
           }            /* overflow                                     */
/*----------------------------------------------------------------------*/
                        /* get most significant part of intern mantissa */
        uh = i->m[0];
        for (c=1-ZERO_BITS;!(uh & MSB);c++) uh <<= 1 ;
                        /* at least one bit is set !!!                  */
                        /* determine shift value                        */
/*----------------------------------------------------------------------*/
        if ((B_LENGTH-ZERO_BITS)-c+i->e*B_LENGTH<-CHARAC)
           {
           return(UFLOW);
           }            /* underflow                                    */
/*----------------------------------------------------------------------*/
        expo = i->e*B_LENGTH+(B_LENGTH-ZERO_BITS)-c;
                        /* determine exponent of double variable        */
/*----------------------------------------------------------------------*/
        for (k=0;k<D_U_RATIO;k++) mant[k] = ZERO;
                        /* initialize mantissa with zero                */
/*----------------------------------------------------------------------*/
        if (c==0)
           {
           for (k=0;k<D_U_RATIO;k++) mant[k] = i->m[k];
           b_comp(d,expo,mant,vz);
           if (i->l>1)
              {
              for (k=D_U_RATIO;k<i->l;k++) if (i->m[k]) return(ROUND);
              }
           return(0);
           }            /* no shift required                            */
/*----------------------------------------------------------------------*/
        if (c<0)
           {
           c = -c;
           mant[0] = i->m[0]>>c;
           if (i->l>1)
              {
              for (k=1;k<D_U_RATIO;k++)
                 mant[k] = (i->m[k-1] << (B_LENGTH-c)) |
                                         (i->m[k] >> c);
              b_comp(d,expo,mant,vz);
              if (i->m[D_U_RATIO-1] << (B_LENGTH-c)) return(ROUND);
              for (k=D_U_RATIO;k<i->l;k++) if (i->m[k]) return(ROUND);
              }
           else
              {
              mant[1] = i->m[0] << (B_LENGTH-c);
              b_comp(d,expo,mant,vz);
              }
           return(0);
           }            /* shift right                                  */
/*----------------------------------------------------------------------*/
        if (i->l>1)
           {
           for (k=0;k<D_U_RATIO-1;k++)
              mant[k] = (i->m[k] << c) |
                        (i->m[k+1] >> (B_LENGTH-c));
           mant[D_U_RATIO-1] = i->m[D_U_RATIO-1] << c;
           if (i->l>2)
              mant[D_U_RATIO-1] |= (i->m[D_U_RATIO] >> (B_LENGTH-c));
           b_comp(d,expo,mant,vz);
           if (i->m[D_U_RATIO-1] << c) return(ROUND);
           for (k=D_U_RATIO;k<i->l;k++) if (i->m[k]) return(ROUND);
           }
        else
           {
           mant[0] = i->m[0]<<c;
           b_comp(d,expo,mant,vz);
           }
/*----------------------------------------------------------------------*/
        return(0);      /* shift left                                   */
        }





