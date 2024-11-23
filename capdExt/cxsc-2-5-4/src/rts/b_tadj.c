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

/* CVS $Id: b_tadj.c,v 1.21 2014/01/30 17:24:05 cxsc Exp $ */

/****************************************************************/
/*                                                              */
/*      Filename        : b_tadj.c                              */
/*                                                              */
/*      Entries         : void b_tadj(mant,expo)                */
/*                        a_btyp mant[BSIZE];                   */
/*                        a_intg *expo;                         */
/*                                                              */
/*      Arguments       : mant = mantossa of tenbyte value      */
/*                        expo = exponent of tenbyte value      */
/*                                                              */
/*      Description     : Adjust denormalized tenbyte,          */
/*                        overflow and underflow.               */
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
local void b_tadj(a_btyp *mant,a_intg *expo)
#else
local void b_tadj(mant,expo)

a_btyp mant[BSIZE];
a_intg *expo;
#endif
        {
        a_intg i;
        a_btyp rc;

        E_TPUSH("b_tadj")

        if (*expo>tEXPO_MAX) {

           /* if overflow trap enabled                          */
           if (e_of_e()) {
              *expo -= tEXPO_ADJUST;

              E_TPOPP("b_tadj")
              return;
              }
           else    {

              /* set all sgnificant mantissa bits               */
              mant[0] = MAX_BASETYPE;
              mant[1] = MAX_BASETYPE;
              mant[2] = MSB;

              *expo = tEXPO_MAX;

              e_sieo();

              E_TPOPP("b_tadj")
              return;
              }
           }
        else if (*expo<tEXPO_MIN) {

           /* if underflow trap enabled                         */
           if (e_uf_e()) {
              *expo += tEXPO_ADJUST;

              E_TPOPP("b_tadj")
              return;
              }
           else    {
              rc = NO_ERROR;
              for (i=2; i<BSIZE; i++ )
              if (mant[i]) {
                    e_sufo();
                    rc = UNDERFLOW;
                    break;
                    }

              /* denormalized number will result or zero        */
              if (*expo<tEXPO_MIN-tMANTL)
                 b_shru( mant, BSIZE, tMANTL+1 );
              else
                 b_shru( mant, BSIZE, tEXPO_MIN-*expo );

              /* exponent set to denormalized exponent          */
              *expo = tEXPO_MIN-1;

              if (rc==NO_ERROR)
              for ( i=2; i<BSIZE; i++ )
              if (mant[i]) {
                       e_sufo();
                       break;
                       }
              }
           }

        E_TPOPP("b_tadj")
        return;
        }





