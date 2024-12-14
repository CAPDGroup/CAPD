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

/* CVS $Id: body.h,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : body.h                                */
/*                                                              */
/*      Description     : Addition or subtraction of a_btyp     */
/*                        mantissas                             */
/*                                                              */
/*      Notes           : included by addbody.h,subbody.h       */
/*                                                              */
/****************************************************************/

        /* order operands in order to get |a|>=|b|              */
        if ( expob > expoa ) {
                i = expoa; expoa = expob; expob = i;
                i = (a_intg) vza;   vza = vzb;  vzb = (a_bool) i;
#if C_P_3
                hu = manta[0];
                manta[0] = mantb[0];
                mantb[0] = hu;
                hu = manta[1];
                manta[1] = mantb[1];
                mantb[1] = hu;
#else
                for ( i=0 ; i<D_U_RATIO ; i++ )
                        {
                        hu = manta[i];
                        manta[i] = mantb[i];
                        mantb[i] = hu;
                        }
#endif
                }
        else if ( expoa == expob ) {
                for ( i=0; i<D_U_RATIO; i++ )
                    if ( manta[i]<mantb[i] ) {
                        i = (a_intg) vza; vza = vzb; vzb = (a_bool) i;
#if C_P_3
                        hu = manta[0];
                        manta[0] = mantb[0];
                        mantb[0] = hu;
                        hu = manta[1];
                        manta[1] = mantb[1];
                        mantb[1] = hu;
#else
                        for (i=0;i<D_U_RATIO;i++ ) {
                                hu = manta[i];
                                manta[i] = mantb[i];
                                mantb[i] = hu;
                                }
#endif
                        break;
                        }
                    else if ( manta[i]>mantb[i] ) break;
                }

        /* initialize unused mantissa positions                 */
        for ( i=D_U_RATIO; i<BSIZE; i++ ) manta[i] = ZERO;

        /* 1. case                                              */
        /* difference of exponents >= length of mantissa+2      */
        if ( expoa-expob>=MANTL+2 ) {

                /* propagate a borrow to mantissa               */
                if (vza!=vzb) b_subc(manta+D_U_RATIO);

                /* set lsb of of result mantissa                */
                else manta[D_U_RATIO] |= LSB;
                }

        /* 2. case                                              */
        /* difference of exponents < length of mantissa+2       */
        else {

                /* initialize unused mantissa positions         */
                for ( i=D_U_RATIO; i<BSIZE; i++ ) mantb[i] = ZERO;

                /* adjust mantissa of mantb to manta            */
                b_shru( mantb, BSIZE, expoa-expob );

                /* add a_btyp mantissas                         */
                if ( vza==vzb ) {
                        (void)b_addm(BSIZE,manta,mantb);
                        }

                /* subtract a_btyp mantissas                    */
                else {
                        (void)b_subm(BSIZE,manta,mantb);

                        /* return zero if no bits are set            */
                        if (b_test(BSIZE,manta))
                            {
                            b_comp(&a,-CHARAC,manta,vzz);
                            E_TPOPP("body")
                            return(a);
                            }
                        }
                }

        /* normalization of mantissa                            */
        if ( SHFT_MASK & *manta ) {
                b_shr1( manta, BSIZE );
                expoa++;
                }
        else while ( (HIDDEN_BIT & *manta)==ZERO ) {
                b_shl1( manta, BSIZE );
                expoa--;
                }

        /* adjust denormalized number and check inexact data    */
        rc = b_adj(manta,&expoa);





