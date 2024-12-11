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

/* CVS $Id: mulbody.h,v 1.22 2014/01/30 17:24:10 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : mulbody.h                             */
/*                                                              */
/*      Description     : Multiplication of IEEE numbers.       */
/*                                                              */
/*      Notes           : included by r_muln,r_muld,r_mulu      */
/*                                                              */
/****************************************************************/

        if ( expoa>EXPO_MAX || expob>EXPO_MAX )
                {

                /* a = infinity or NaN                          */
                if (expoa==EXPO_MAX+1) {

                   /* a = infinity                              */
                   if (MANT_INFINITY(manta)) {

                      /* b = zero                               */
                      if (zerob) {
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,10,
                                                E_TDBL+E_TEXT(1),&a,
                                                E_TDBL+E_TEXT(2),&b,
                                                E_TRES+E_TDBL,&a);
                         E_TPOPP("mulbody")
                         return(a);
                         }

                      /* b = infinity or NaN                    */
                      if (expob==EXPO_MAX+1) {

                         /* a = b = infinity                    */
                         if (MANT_INFINITY(mantb)) {
                            if (vza!=vzc) b_comp(&a,expoa,manta,vzc);
                            E_TPOPP("mulbody")
                            return(a);
                            }

                         /* a = infinity b = NaN                */
                         if (SIGNALING(mantb[0]))
                            e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                                   E_TDBL+E_TEXT(1),&a,
                                                   E_TDBL+E_TEXT(2),&b,
                                                   E_TRES+E_TDBL,&b);
                         E_TPOPP("mulbody")
                         return(b);
                         }

                      /* multiplication with infinity => a      */
                      if (vza!=vzc) b_comp(&a,expoa,manta,vzc);
                      E_TPOPP("mulbody")
                      return(a);
                      }

                   /* a = NaN                                   */
                   else {
                      if (SIGNALING(manta[0]))
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                                E_TDBL+E_TEXT(1),&a,
                                                E_TDBL+E_TEXT(2),&b,
                                                E_TRES+E_TDBL,&a);

                      if (expob>EXPO_MAX && !MANT_INFINITY(mantb) && 
                          SIGNALING(mantb[0]))
                         {
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                                E_TDBL+E_TEXT(1),&a,
                                                E_TDBL+E_TEXT(2),&b,
                                                E_TRES+E_TDBL,&b);
                         E_TPOPP("mulbody")
                         return(b);
                         }

                      E_TPOPP("mulbody")
                      return(a);
                      }
                   }

                /* b = infinity or NaN                          */
                else {

                   /* b = infinity                             */
                   if (MANT_INFINITY(mantb)) {

                      /* a = zero                                  */
                      if (zeroa) {
                         e_trap(INV_OP+E_IEEE,8,E_TMSG,10,
                                              E_TDBL+E_TEXT(1),&a,
                                              E_TDBL+E_TEXT(2),&b,
                                              E_TRES+E_TDBL,&b);
                         }
                      else if (vzb!=vzc) b_comp(&b,expob,mantb,NOT(vzb));

                      E_TPOPP("mulbody")
                      return(b);
                      }

                    /* b = NaN                                  */
                    if (SIGNALING(mantb[0]))
                       e_trap(INV_OP+E_IEEE,8,E_TMSG,5,
                                              E_TDBL+E_TEXT(1),&a,
                                              E_TDBL+E_TEXT(2),&b,
                                              E_TRES+E_TDBL,&b);
                    E_TPOPP("mulbody")
                    return(b);
                    }
                }

        /* product is zero                                      */
        if (zeroa || zerob) {
            E_TPOPP("mulbody")
            return(*r_zero);
            }

        /* add exponents                                        */
        expoc = expoa+expob;

#if C_P_3
        /* determine exact product                              */
        lang[BSIZE-1] = ZERO;
        b_prod(manta,mantb,lang);
#else
        /* initialize exact mantissa array lang                 */
        B_CLEAR(lang,BSIZE)

        /* determine exact product                              */
        for ( i=2*D_U_RATIO-1; i>=D_U_RATIO; i-- ) {
            j = i-D_U_RATIO;
            while ( j<D_U_RATIO ) {
                b_muad(manta[j],mantb[i-j-1],lang+i);
                j++;
                }
            }
        for ( i=D_U_RATIO-1; i>0; i-- ) {
            j = i-1;
            while( j>=0 ) {
                b_muad(manta[j],mantb[i-j-1],lang+i);
                j--;
                }
            }
#endif

        /* adjust mantissa                                      */
        b_shlu(lang,BSIZE,ZERO_BITS);

        /* normalization                                        */
        if ( SHFT_MASK & lang[0]) {
            b_shr1(lang,BSIZE);
            expoc++;
            }
        else {
            while ( (HIDDEN_BIT & lang[0])==ZERO ) {
                b_shl1(lang,BSIZE);
                expoc--;
                }
            }

        /* adjust denormalized number                           */
        rc = b_adj(lang,&expoc);





