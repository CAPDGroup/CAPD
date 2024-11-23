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

/* CVS $Id: b_ltor.c,v 1.22 2014/02/27 15:29:42 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_ltor.c                              */
/*                                                              */
/*      Entries         : a_btyp b_ltor(i,r,rnd)                */
/*                        multiprecision i;                     */
/*                        a_real *r;                            */
/*                        a_intg rnd;                           */
/*                                                              */
/*      Arguments       : i    = long value                     */
/*                        r    = real value                     */
/*                        rnd  = rounding                       */
/*                               <0 : downwards                 */
/*                               =0 : nearest                   */
/*                               >0 : upwards                   */
/*                                                              */
/*      Function value  : error code                            */
/*                                                              */
/*      Description     : convert long to real value.           */
/*                                                              */
/*      Note            : long value mantissa must not exceed   */
/*                        three digits.                         */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern a_real *r_sero;
extern a_real *r_zero;
extern a_real *r_eps_;
extern a_real *r_meps;
extern a_btyp b_maxl;
#endif

#ifdef LINT_ARGS
local a_btyp b_ltor(multiprecision i,a_real *r,a_intg rnd)
#else
local a_btyp b_ltor(i,r,rnd)

multiprecision i;
a_real *r;
a_intg rnd;
#endif
        {
        /* a_btyp rc,rc1,mant[BSIZE]; */
        a_btyp rc,rc1,mant[BSIZE+1];
        a_intg expo,k;
        a_bool sign;

        E_TPUSH("b_ltor")

        /* number is zero                                       */
        if (i->z)
           {
           if (rnd<0)
              {
              if (i->r && i->s)
                 R_ASSIGN(*r,*r_meps);
              else
                 R_ASSIGN(*r,*r_zero);
              }
           else if (rnd>0)
              {
              if (i->r && i->s==0)
                 R_ASSIGN(*r,*r_eps_);
              else if (i->s)
                 R_ASSIGN(*r,*r_sero);
              else
                 R_ASSIGN(*r,*r_zero);
              }
           else
              R_ASSIGN(*r,*r_zero);

           E_TPOPP("b_ltor")
           return(NO_ERROR);
           }

        /* sign and exponent of real value */
        sign = i->s;
        expo = (i->e<<LOG_B_LENGTH)+EXPO_SHIFT;

        /* mantissa digits of real value */
        for (k=0;k<i->l && k<BSIZE-1;k++) mant[k] = i->m[k];
        if (i->l>=BSIZE)
           {

           /* rounding outward required */
           if ((rnd>0 && sign==0) || (rnd<0 && sign))
              {

              /* at least one trailing bit assumed */
              if (i->l<b_maxl) mant[BSIZE-1] = MSB;

              else
                 {

                 /* rounding bits may cause carry */
                 for (k = BSIZE;k<b_maxl-1;k++)
                    if (mant[k]!=MAX_BASETYPE) break;

                 /* no carry possible */
                 if (k<b_maxl-1) mant[BSIZE-1] = MSB;
                 else
                    {

                    /* no carry occurred */
                    if (mant[k]+i->r>=i->r) mant[BSIZE-1] = MSB;

                    /* carry occurred */
                    else
                       {
                       /* add carry to mantissa */
                       if (b_bcad((a_intg)(BSIZE-1),mant))
                          {
                          mant[0] = LSB;
                          expo += B_LENGTH;
                          }

                       /* non-zero remainder mantiassa */
                       if (k+1<i->l || mant[k]+i->r>ZERO)
                          mant[BSIZE-1] = MSB;
                       }
                    }
                 }
              }
           }
        else /* i->l<BSIZE */
           {

           /* clear all remaining mantissa positions in mant */
           while (k<BSIZE) mant[k++] = ZERO;

           /* rounding outward required */
           if ((rnd>0 && sign==0) || (rnd<0 && sign))
              {

              /* consider all rounding bits */
              if (i->r)
                 {
                 if ((mant[b_maxl-1] += i->r)<i->r)
                    {
                    if (b_bcad((a_intg)(b_maxl-1),mant))
                       {
                       mant[0] = LSB;
                       expo += B_LENGTH;
                       }
                    }
                 }
              }
           }

        /* normalization of mantissa                            */
        if (SHFT_MASK & *mant)
           {
           do {
              b_shr1(mant,BSIZE);
              expo++;
              }
           while (SHFT_MASK & *mant);
           }
        else if ((HIDDEN_BIT & *mant)==ZERO)
           {
           k = 1;
           while ((HIDDEN_BIT & (*mant<<k))==ZERO) k++;
           b_shlu(mant,BSIZE,k);
           expo -= k;
           }

        /* adjust denormalized number and check inexact data    */
        rc = b_adj(mant,&expo);

        /* round mantissa                                       */
        if (rnd<0) rc1 = b_rndd(mant,&expo,sign);
        else if (rnd>0) rc1 = b_rndu(mant,&expo,sign);
        else rc1 = b_rndn(mant,&expo);

        /* compose result                                       */
        b_comp(r,expo,mant,sign);

        if (rc==ZERO) rc = rc1;

        E_TPOPP("b_ltor")
        return(rc);
        }





