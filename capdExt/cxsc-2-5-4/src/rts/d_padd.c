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

/* CVS $Id: d_padd.c,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : d_padd.c                              */
/*                                                              */
/*      Entries         : void d_padd(c,a,b)                    */
/*                        d_otpr *c;                            */
/*                        a_real a,b;                           */
/*                                                              */
/*      Arguments       : c = dotprecision variable             */
/*                        a = real value                        */
/*                        b = real value                        */
/*                                                              */
/*      Description     : Add product to dotprecision variable  */
/*                        c = c+a*b                             */
/*                                                              */
/*                   by introducing variable length             */
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
local void d_padd(d_otpr *c,a_real a,a_real b)
#else
local void d_padd(c,a,b)

d_otpr *c;
a_real a;
a_real b;
#endif
        {
        register a_btyp *rega;
        a_btyp *d;
        a_btyp manta[D_U_RATIO];
        a_btyp mantb[D_U_RATIO];
        a_btyp lang[BSIZE];
        a_intg expoa,expob,expoc;
        a_bool vza,vzb,vzc;
        a_intg i,k,msb,lsb;
        a_bool zeroa,zerob;

        E_TPUSH("d_padd")

        rega = *c;

        /* decomposition                                        */
        zeroa = b_deko(a,&expoa,manta,&vza);
        zerob = b_deko(b,&expob,mantb,&vzb);

        /* sign of product                                      */
        vzc = vza^vzb;

        /* NaN or infinity                                      */
        if (expoa>EXPO_MAX || expob>EXPO_MAX)
           {
           rega[A_STATUS] |= A_PZERO+A_MZERO;

           /* a = infinity or NaN                               */
           if (expoa==EXPO_MAX+1) {

              /* a = infinity                                   */
              if (MANT_INFINITY(manta)) {

                 /* b = zero                                    */
                 if (zerob) {
                    e_trap(INV_OP+E_IEEE,6,E_TMSG,10,
                           E_TDBL,&a,E_TDBL,&b);
                    rega[A_STATUS] |= A_QUIETNAN;
                    rega[A_LOWNAN] = INV_OP;
                    E_TPOPP("d_padd")
                    return;
                    }

                 /* b = infinity or NaN                         */
                 if (expob>EXPO_MAX) {

                    /* a = b = infinity                         */
                    if (MANT_INFINITY(mantb)) {
                       if ((rega[A_STATUS] & A_QUIETNAN)==ZERO)
                          {
                          if (vzc && (rega[A_STATUS] & A_PINFINITY))
                             {
                             e_trap(INV_OP+E_IEEE,8,E_TMSG,9,
                                    E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                             rega[A_STATUS] |= A_QUIETNAN;
                             rega[A_LOWNAN] = INV_OP;
                             }
                          else if (NOT(vzc) &&
                                   (rega[A_STATUS] & A_MINFINITY))
                             {
                             e_trap(INV_OP+E_IEEE,8,E_TMSG,9,
                                    E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                             rega[A_STATUS] |= A_QUIETNAN;
                             rega[A_LOWNAN] = INV_OP;
                             }
                          else
                             {
                             e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,6,E_TMSG,13,
                                    E_TDBL,&a,E_TDBL,&b);
                             if (vzc) rega[A_STATUS] |= A_MINFINITY;
                             else rega[A_STATUS] |= A_PINFINITY;
                             }
                          }
                       E_TPOPP("d_padd")
                       return;
                       }

                    /* a = infinity b = NaN                     */
                    if (SIGNALING(mantb[0]))
                       {
                       e_trap(INV_OP+E_IEEE,6,E_TMSG,5,
                              E_TDBL,&a,E_TDBL,&b);
                       rega[A_STATUS] |= A_QUIETNAN;
                       rega[A_LOWNAN] = INV_OP;
                       }
                    else if ((rega[A_STATUS] & A_QUIETNAN)==ZERO)
                       {
                       e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,6,E_TMSG,14,
                              E_TDBL,&a,E_TDBL,&b);
                       rega[A_STATUS] |= A_QUIETNAN;
                       rega[A_LOWNAN] = ((a_btyp *)&b)[B_LPART];
                       }
                    E_TPOPP("d_padd")
                    return;
                    }

                 /* multiplication with infinity => a */
                 if ((rega[A_STATUS] & A_QUIETNAN)==ZERO)
                    {
                    if (vzc && (rega[A_STATUS] & A_PINFINITY))
                       {
                       e_trap(INV_OP+E_IEEE,8,E_TMSG,9,
                              E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                       rega[A_STATUS] |= A_QUIETNAN;
                       rega[A_LOWNAN] = INV_OP;
                       }
                    else if (NOT(vzc) && (rega[A_STATUS] & A_MINFINITY))
                       {
                       e_trap(INV_OP+E_IEEE,6,E_TMSG,9,
                              E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                       rega[A_STATUS] |= A_QUIETNAN;
                       rega[A_LOWNAN] = INV_OP;
                       }
                    else
                       {
                       e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,6,E_TMSG,13,
                              E_TDBL,&a,E_TDBL,&b);
                       if (vzc) rega[A_STATUS] |= A_MINFINITY;
                       else rega[A_STATUS] |= A_PINFINITY;
                       }
                    }
                 E_TPOPP("d_padd")
                 return;
                 }

              /* a = NaN                                        */
              else
                 {
                 if (SIGNALING(manta[0]) ||
                     (expob>EXPO_MAX && SIGNALING(mantb[0])))
                    {
                    e_trap(INV_OP+E_IEEE,6,E_TMSG,5,
                           E_TDBL,&a,E_TDBL,&b);
                    rega[A_STATUS] |= A_QUIETNAN;
                    rega[A_LOWNAN] = INV_OP;
                    }
                 else if ((rega[A_STATUS] & A_QUIETNAN)==ZERO)
                    {
                    e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,6,E_TMSG,14,
                           E_TDBL,&a,E_TDBL,&b);
                    rega[A_STATUS] |= A_QUIETNAN;
                    rega[A_LOWNAN] = ((a_btyp *)&a)[B_LPART];
                    }
                 E_TPOPP("d_padd")
                 return;
                 }
              }

           /* b = infinity or NaN                               */
           else {

              /* b = infinity                                   */
              if (MANT_INFINITY(mantb)) {

                 /* a = zero                                    */
                 if (zeroa)
                    {
                    e_trap(INV_OP+E_IEEE,8,E_TMSG,10,
                           E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                    rega[A_STATUS] |= A_QUIETNAN;
                    rega[A_LOWNAN] = INV_OP;
                    }
                 else if ((rega[A_STATUS] & A_QUIETNAN)==ZERO)
                    {
                    if (vzc && (rega[A_STATUS] & A_PINFINITY))
                       {
                       e_trap(INV_OP+E_IEEE,8,E_TMSG,9,
                              E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                       rega[A_STATUS] |= A_QUIETNAN;
                       rega[A_LOWNAN] = INV_OP;
                       }
                    else if (NOT(vzc) && (rega[A_STATUS] & A_MINFINITY))
                       {
                       e_trap(INV_OP+E_IEEE,8,E_TMSG,9,
                              E_TDTP,c,E_TDBL,&a,E_TDBL,&b);
                       rega[A_STATUS] |= A_QUIETNAN;
                       rega[A_LOWNAN] = INV_OP;
                       }
                    else
                       {
                       e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,6,E_TMSG,13,
                              E_TDBL,&a,E_TDBL,&b);
                       if (vzc) rega[A_STATUS] |= A_MINFINITY;
                       else rega[A_STATUS] |= A_PINFINITY;
                       }
                    }
                 E_TPOPP("d_padd")
                 return;
                 }

              /* b = NaN                                        */
              if (SIGNALING(mantb[0]))
                 {
                 e_trap(INV_OP+E_IEEE,6,E_TMSG,5,
                        E_TDBL,&a,E_TDBL,&b);
                 rega[A_STATUS] |= A_QUIETNAN;
                 rega[A_LOWNAN] = INV_OP;
                 }
              else if ((rega[A_STATUS] & A_QUIETNAN)==ZERO)
                 {
                 e_trap(NO_ERROR+E_IEEE+E_EMSG+E_ECNT,6,E_TMSG,14,
                        E_TDBL,&a,E_TDBL,&b);
                 rega[A_STATUS] |= A_QUIETNAN;
                 rega[A_LOWNAN] = ((a_btyp *)&b)[B_LPART];
                 }
              E_TPOPP("d_padd")
              return;
              }
           }

        /* product is zero */
        if (zeroa || zerob)
           {
           rega[A_STATUS] |= (vzc) ? A_PZERO : A_MZERO;
           E_TPOPP("d_padd")
           return;
           }

        rega[A_STATUS] |= A_PZERO+A_MZERO;

        /* determine exact product      */
        lang[0] = ZERO;
        b_prod(manta,mantb,lang+1);

        /* exponent of least significant bit of product */
        expoc = expoa+expob-2*(MANTL-1);

        /* allign mantissa to dotprecision      */
        lsb = A_D_P-B_ASHR(expoc,LOG_B_LENGTH);
        msb = lsb-(BSIZE-1);
        b_shlu(lang,BSIZE,expoc & (B_LENGTH-1));

        /* accu is empty        */
        if (rega[A_BEGIN]==ZERO)
           {
           rega[A_END] = lsb;
           rega[A_BEGIN] = msb;
           for (i=0;i<BSIZE;i++,msb++)
              rega[msb] = lang[i];
           rega[A_SIGN] = vzc;
           }

        else
           {
           if (rega[A_END]<lsb)
              {
              /* !!! assuming cleared  positions !!! */
              rega[A_END] = lsb;
              }

           /* add value to dotprecision variable of same sign   */
           if (vzc==rega[A_SIGN])
              {

              /* add product mantissa to accu mantissa */
              if (b_addm(BSIZE,&(rega[msb]),lang))
                 {

                 /* msb is updated here */
                 do
                    msb--;
                 while (++(rega[msb])==ZERO);
                 }
              }

           /* add value to dotprecision with different sign     */
           else
              {

              /* subtract product mantissa from accu mantissa */
              if (b_subm(BSIZE,&(rega[msb]),lang))
                 {
                 if ((k = msb-1)<rega[A_BEGIN])
                    {
                    i = rega[A_END];
                    while (++k<=i) rega[k] = ~rega[k];
                    while (++(rega[i])==ZERO) i--;
                    rega[A_SIGN] = 1-rega[A_SIGN];
                    }
                 else
                    {
                    while ((rega[k])--==ZERO) k--;
                    }
                 }
              }

           if (rega[A_BEGIN]>msb)
              {
              /* !!! assuming cleared position !!! */
              rega[A_BEGIN] = msb;
              }
           }

        /* eliminate leading zeros      */
        d = &(rega[rega[A_BEGIN]]);
        while (*d++==ZERO)
           {
           if (++(rega[A_BEGIN])>rega[A_END])
              {
              rega[A_BEGIN] = rega[A_END] = ZERO;
              break;
              }
           }

        /* eliminate trailing zeros     */
        if (rega[A_BEGIN])
           {
           d = &(rega[rega[A_END]]);
           while (*d--==ZERO) rega[A_END]--;
           }

        E_TPOPP("d_padd")
        return;
        }





