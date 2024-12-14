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

/* CVS $Id: b_tdek.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tdek.c                              */
/*                                                              */
/*      Entries         : void b_tdek(a,expo,mant,vz)           */
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
/*      Description     : Decomposition of a tenbyte value.     */
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
local a_bool b_tdek(tenbyte *a,a_intg *expo,a_btyp *mant,a_bool  *vz)
#else
local a_bool b_tdek(a,expo,mant,vz)

tenbyte *a;
a_intg *expo;
a_btyp mant[BSIZE];
a_bool  *vz;
#endif
        {
        E_TPUSH("b_tdek")

        /* sign of real value                                   */
#if O_P_4
        *vz = (a->c[9] & tMSB) ? TRUE : FALSE;
#else
        *vz = (a->c[0] & tMSB) ? TRUE : FALSE;
#endif

        /* exponent of real value                               */
#if O_P_4
        *expo = (((a_intg)a->c[9] & ~tMSB)<<tSHIFT)+
                (a_intg)a->c[8]-tCHARAC;
#else
        *expo = (((a_intg)a->c[0] & ~tMSB)<<tSHIFT)+
                (a_intg)a->c[1]-tCHARAC;
#endif

        /* mantissa of real value                               */
#if O_P_4
        mant[0] = ((((a_btyp)a->c[7]<<tSHIFT)+
                     (a_btyp)a->c[6])<<tSHIFT)
                    +(a_btyp)a->c[5];
        mant[1] = (((((a->c[4]<<tSHIFT)+
                     (a_btyp)a->c[3])<<tSHIFT)
                    +(a_btyp)a->c[2])<<tSHIFT)+(a_btyp)a->c[1];
        mant[2] = (a_btyp)a->c[0]<<(3*tSHIFT);
        mant[3] = ZERO;
        mant[4] = ZERO;
#else
        mant[0] = ((((a_btyp)a->c[2]<<tSHIFT)
                    +(a_btyp)a->c[3])<<tSHIFT)
                    +(a_btyp)a->c[4];
        mant[1] = ((((((a_btyp)a->c[5]<<tSHIFT)
                    +(a_btyp)a->c[6])<<tSHIFT)
                    +(a_btyp)a->c[7])<<tSHIFT)+(a_btyp)a->c[8];
        mant[2] = (a_btyp)a->c[9]<<(3*tSHIFT);
        mant[3] = ZERO;
        mant[4] = ZERO;
#endif

        if (!(mant[0] | mant[1] | mant[2]))
           {

           /* real value is zero                                */
           if (*expo==-tCHARAC)
              {
              E_TPOPP("b_tdek")
              return(TRUE);
              }

           /* value is infinity or NaN                          */
           e_trap(INV_OP+E_IEEE,2,E_TMSG,70);
           }

        /* denormalized number                                  */
        else if (*expo==-tCHARAC)
           *expo = -tCHARAC+1;

        E_TPOPP("b_tdek")
        return(FALSE);
        }





