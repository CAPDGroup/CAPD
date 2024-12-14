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

/* CVS $Id: b_tcom.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tcom.c                              */
/*                                                              */
/*      Entries         : void b_tcom(a,expo,mant,vz)           */
/*                        tenbyte *a;                           */
/*                        a_intg *expo;                         */
/*                        a_btyp mant[BSIZE];                   */
/*                        a_bool  *vz;                          */
/*                                                              */
/*      Arguments       : a = tenbyte value                     */
/*                        expo = exponent                       */
/*                        mant = mantissa                       */
/*                        vz = sign                             */
/*                                                              */
/*      Description     : Composition of a tenbyte value.       */
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
local void b_tcom(tenbyte *a,a_intg expo,a_btyp *mant,a_bool  vz)
#else
local void b_tcom(a,expo,mant,vz)

tenbyte *a;
a_btyp mant[BSIZE];
a_intg expo;
a_bool  vz;
#endif
        {
        E_TPUSH("b_tcom")

        /* mantissa of real value                               */
#if O_P_4
        a->c[0] = (unsigned char) (mant[1] & tMASK); mant[1] >>= tSHIFT;
        a->c[1] = (unsigned char) (mant[1] & tMASK); mant[1] >>= tSHIFT;
        a->c[2] = (unsigned char) (mant[1] & tMASK); mant[1] >>= tSHIFT;
        a->c[3] = (unsigned char) (mant[1] & tMASK);
        a->c[4] = (unsigned char) (mant[0] & tMASK); mant[0] >>= tSHIFT;
        a->c[5] = (unsigned char) (mant[0] & tMASK); mant[0] >>= tSHIFT;
        a->c[6] = (unsigned char) (mant[0] & tMASK); mant[0] >>= tSHIFT;
        a->c[7] = (unsigned char) (mant[0] & tMASK);
#else
        a->c[9] = (unsigned char) (mant[1] & tMASK) ; mant[1] >>= tSHIFT;
        a->c[8] = (unsigned char) (mant[1] & tMASK); mant[1] >>= tSHIFT;
        a->c[7] = (unsigned char) (mant[1] & tMASK); mant[1] >>= tSHIFT;
        a->c[6] = (unsigned char) (mant[1] & tMASK);
        a->c[5] = (unsigned char) (mant[0] & tMASK); mant[0] >>= tSHIFT;
        a->c[4] = (unsigned char) (mant[0] & tMASK); mant[0] >>= tSHIFT;
        a->c[3] = (unsigned char) (mant[0] & tMASK); mant[0] >>= tSHIFT;
        a->c[2] = (unsigned char) (mant[0] & tMASK);
#endif

        expo += tCHARAC;

        /* exponent of real value                               */
#if O_P_4
        a->c[8] = (unsigned char) (expo & tMASK);  expo >>= tSHIFT;
        a->c[9] = (unsigned char) (expo & tMASK);
#else
        a->c[1] = (unsigned char) (expo & tMASK);  expo >>= tSHIFT;
        a->c[0] = (unsigned char) (expo & tMASK);
#endif

        /* sign of real value                                   */
#if O_P_4
        if (vz) a->c[9] |= tMSB;
#else
        if (vz) a->c[0] |= tMSB;
#endif

        E_TPOPP("b_tcom")
        return;
        }





