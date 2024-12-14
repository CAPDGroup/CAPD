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

/* CVS $Id: divbody.h,v 1.21 2014/01/30 17:24:06 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : divbody.h                             */
/*                                                              */
/*      Description     : Division of a_btyp mantissas          */
/*                                                              */
/****************************************************************/

        /* denormalized numbers                                 */
        while ( !(HIDDEN_BIT & manta[0]) ) {
#if C_P_3
                manta[0] = (manta[0]<<1) |
                           (manta[1]>>(B_LENGTH-1));
                manta[1] <<= 1;
#else
                b_shl1( manta, D_U_RATIO );
#endif
                expoa--;
                }
        while ( !(HIDDEN_BIT & mantb[0]) ) {
#if C_P_3
                mantb[0] = (mantb[0]<<1) |
                           (mantb[1]>>(B_LENGTH-1));
                mantb[1] <<= 1;
#else
                b_shl1( mantb, D_U_RATIO );
#endif
                expob--;
                }

        /* exponent                                             */
        expoc = expoa-expob;

        /* division of mantissas                                */
        b_mdiv( manta, mantb, mantc, &expoc );

        /* adjust denormalized number                           */
        rc = b_adj(mantc,&expoc);





