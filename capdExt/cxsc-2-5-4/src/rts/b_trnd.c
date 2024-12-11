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

/* CVS $Id: b_trnd.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_trnd.c                              */
/*                                                              */
/*      Entries         : void b_trnd(mant,expo,vz)             */
/*                        a_btyp *mant;                         */
/*                        a_intg *expo;                         */
/*                        a_bool  vz;                           */
/*                                                              */
/*      Arguments       : mant = mantissa to be rounded         */
/*                        expo = exponent                       */
/*                        vz = sign                             */
/*                                                              */
/*      Description     : Round mantissa according to           */
/*                        flag b_rflg.                          */
/*                                                              */
/****************************************************************/

#ifndef ALL_IN_ONE
#ifdef AIX
#include "/u/p88c/runtime/o_defs.h"
#else
#include "o_defs.h"
#endif
#define local
extern int b_rflg;
#endif

#ifdef LINT_ARGS
local void b_trnd(a_btyp *mant,a_intg *expo,a_bool  vz)
#else
local void b_trnd(mant,expo,vz)

a_btyp *mant;
a_intg *expo;
a_bool  vz;
#endif
        {
        a_intg i;

        E_TPUSH("b_trnd")

        if (b_rflg==NEAREST)
           {

           /* msb of mant[2] is set                     */
           if (mant[2] & MSB)
              {
              e_sieo();

              /* test trailing mantissa to be zero      */
              if ((mant[2] & NOT_MSB)==ZERO)
                 {
                 for (i=3;i<BSIZE;i++) if (mant[i]) break;

                 /* return if lsb of mant[1] is zero    */
                 if (i>=BSIZE)
                 if (!(mant[1] & LSB))
                    {
                    E_TPOPP("b_trnd")
                    return;
                    }
                 }
              }
           else
              {
              for (i=2;i<BSIZE;i++)
                 if (mant[i])
                    {
                    e_sieo();
                    break;
                    }

              E_TPOPP("b_trnd")
              return;
              }
           }
        else if (b_rflg==CHOP)
           {
           for (i=2;i<BSIZE;i++) if (mant[i]) break;
           if (i==BSIZE)
              {
              E_TPOPP("b_trnd")
              return;
              }
           e_sieo();

           E_TPOPP("b_trnd")
           return;
           }
        else
           {
           for (i=2;i<BSIZE;i++) if (mant[i]) break;
           if (i==BSIZE)
              {
              E_TPOPP("b_trnd")
              return;
              }
           e_sieo();

           /* add correction if trailing mantissa is non-zero   */
           if (!((b_rflg==DOWN && vz==TRUE) || (b_rflg==UP && vz==FALSE)))
              {
              E_TPOPP("b_trnd")
              return;
              }
           }

        /* propagate carry                              */
        mant[1]++;
        if (mant[1]==ZERO)
           {
           mant[0]++;

           /* normalization (only msb is set in mantissa)       */
           if (mant[0]==ZERO)
              {
              mant[0] = MSB;
              (*expo)++;

              /* overflow occurred                      */
              if (*expo>tEXPO_MAX)
                 {
                 if (e_of_e()) *expo -= tEXPO_ADJUST;
                 }
              }
           }

        E_TPOPP("b_trnd")
        return;
        }





